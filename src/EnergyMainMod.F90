module EnergyMainMod

!!! Main energy module including all energy relevant processes
!!! soil/snow thermal property -> radiation -> ground/vegtation heat flux -> snow/soil temperature solver -> soil/snow phase change
!
! --------------------------------------------------------------------------------------------------
! NoahMP uses different approaches to deal with subgrid features of radiation transfer and turbulent
! transfer. It uses 'tile' approach to compute turbulent fluxes, while it uses two-stream approx.
! to compute radiation transfer. Tile approach, assemblying vegetation canopies together,
! may expose too much ground surfaces (either covered by snow or grass) to solar radiation. The
! modified two-stream assumes vegetation covers fully the gridcell but with gaps between tree crowns.
! --------------------------------------------------------------------------------------------------
! turbulence transfer : 'tile' approach to compute energy fluxes in vegetated fraction and
!                         bare fraction separately and then sum them up weighted by fraction
!                     --------------------------------------
!                    / O  O  O  O  O  O  O  O  /          / 
!                   /  |  |  |  |  |  |  |  | /          /
!                  / O  O  O  O  O  O  O  O  /          /
!                 /  |  |  |tile1|  |  |  | /  tile2   /
!                / O  O  O  O  O  O  O  O  /  bare    /
!               /  |  |  | vegetated |  | /          /
!              / O  O  O  O  O  O  O  O  /          /
!             /  |  |  |  |  |  |  |  | /          /
!            --------------------------------------
! --------------------------------------------------------------------------------------------------
! radiation transfer : modified two-stream (Yang and Friedl, 2003, JGR; Niu ang Yang, 2004, JGR)
!                     --------------------------------------  two-stream treats leaves as
!                    /   O   O   O   O   O   O   O   O    /  cloud over the entire grid-cell,
!                   /    |   |   |   |   |   |   |   |   / while the modified two-stream 
!                  /   O   O   O   O   O   O   O   O    / aggregates cloudy leaves into  
!                 /    |   |   |   |   |   |   |   |   / tree crowns with gaps (as shown in
!                /   O   O   O   O   O   O   O   O    / the left figure). We assume these
!               /    |   |   |   |   |   |   |   |   / tree crowns are evenly distributed
!              /   O   O   O   O   O   O   O   O    / within the gridcell with 100% veg
!             /    |   |   |   |   |   |   |   |   / fraction, but with gaps. The 'tile'
!            -------------------------------------- approach overlaps too much shadows.
! --------------------------------------------------------------------------------------------------

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowCoverGroundNiu07Mod,        only : SnowCoverGroundNiu07
  use SnowCoverGroundAR25Mod,         only : SnowCoverGroundAR25
  use GroundRoughnessPropertyMod,     only : GroundRoughnessProperty
  use GroundThermalPropertyMod,       only : GroundThermalProperty
  use SurfaceAlbedoMod,               only : SurfaceAlbedo
  use SurfaceRadiationMod,            only : SurfaceRadiation
  use SurfaceEmissivityMod,           only : SurfaceEmissivity
  use SoilWaterTranspirationMod,      only : SoilWaterTranspiration
  use ResistanceGroundEvaporationMod, only : ResistanceGroundEvaporation
  use PsychrometricVariableMod,       only : PsychrometricVariable
  use SurfaceEnergyFluxVegetatedMod,  only : SurfaceEnergyFluxVegetated
  use SurfaceEnergyFluxBareGroundMod, only : SurfaceEnergyFluxBareGround
  use SoilSnowTemperatureMainMod,     only : SoilSnowTemperatureMain
  use SoilSnowWaterPhaseChangeMod,    only : SoilSnowWaterPhaseChange

  implicit none

contains

  subroutine EnergyMain(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ENERGY
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J          ! grid indices
    logical, allocatable             :: FlagVegSfc(:,:) ! flag: true if vegetated surface (2D array for subroutine call)

! --------------------------------------------------------------------

    ! Allocate local 2D FlagVegSfc array
    allocate(FlagVegSfc(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                        noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    !$acc data create(FlagVegSfc)

    ! Compute FlagVegSfc array and VegAreaIndEff before calling GroundRoughnessProperty
    !$acc parallel loop collapse(2) gang vector present(noahmp, FlagVegSfc)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                                             &
              LeafAreaIndEff          => noahmp%energy%state%LeafAreaIndEff(I,J)          ,& ! in,    leaf area index, after burying by snow
              StemAreaIndEff          => noahmp%energy%state%StemAreaIndEff(I,J)          ,& ! in,    stem area index, after burying by snow
              TemperatureAir2mVeg     => noahmp%energy%state%TemperatureAir2mVeg(I,J)     ,& ! out,   2 m height air temperature [K], vegetated
              VegAreaIndEff           => noahmp%energy%state%VegAreaIndEff(I,J)           ,& ! out,   one-sided leaf+stem area index [m2/m2]
              WindStressEwVeg         => noahmp%energy%state%WindStressEwVeg(I,J)         ,& ! out,   wind stress: east-west [N/m2] above canopy
              WindStressNsVeg         => noahmp%energy%state%WindStressNsVeg(I,J)         ,& ! out,   wind stress: north-south [N/m2] above canopy
              SpecHumidity2mVeg       => noahmp%energy%state%SpecHumidity2mVeg(I,J)       ,& ! out,   water vapor mixing ratio at 2m vegetated
              ExchCoeffShAbvCan       => noahmp%energy%state%ExchCoeffShAbvCan(I,J)       ,& ! out,   exchange coeff [m/s] for heat, above ZeroPlaneDisp, vegetated
              ExchCoeffShLeaf         => noahmp%energy%state%ExchCoeffShLeaf(I,J)         ,& ! out,   leaf sensible heat exchange coeff [m/s], leaf to canopy air
              ExchCoeffShUndCan       => noahmp%energy%state%ExchCoeffShUndCan(I,J)       ,& ! out,   under canopy sensible heat exchange coefficient [m/s]
              ExchCoeffSh2mVeg        => noahmp%energy%state%ExchCoeffSh2mVeg(I,J)        ,& ! out,   2m sensible heat exchange coefficient [m/s] vegetated
              HeatPrecipAdvSfc        => noahmp%energy%flux%HeatPrecipAdvSfc(I,J)         ,& ! out,   precipitation advected heat - total [W/m2]
              RadLwNetCanopy          => noahmp%energy%flux%RadLwNetCanopy(I,J)           ,& ! out,   canopy net longwave radiation [W/m2] (+ to atm)
              RadLwNetVegGrd          => noahmp%energy%flux%RadLwNetVegGrd(I,J)           ,& ! out,   ground net longwave radiation [W/m2] (+ to atm)
              HeatSensibleCanopy      => noahmp%energy%flux%HeatSensibleCanopy(I,J)       ,& ! out,   canopy sensible heat flux [W/m2] (+ to atm)
              HeatSensibleVegGrd      => noahmp%energy%flux%HeatSensibleVegGrd(I,J)       ,& ! out,   vegetated ground sensible heat flux [W/m2] (+ to atm)
              HeatLatentVegGrd        => noahmp%energy%flux%HeatLatentVegGrd(I,J)         ,& ! out,   ground evaporation heat flux [W/m2] (+ to atm)
              HeatLatentCanEvap       => noahmp%energy%flux%HeatLatentCanEvap(I,J)        ,& ! out,   canopy evaporation heat flux [W/m2] (+ to atm)
              HeatLatentCanTransp     => noahmp%energy%flux%HeatLatentCanTransp(I,J)      ,& ! out,   canopy transpiration heat flux [W/m2] (+ to atm)
              HeatGroundVegGrd        => noahmp%energy%flux%HeatGroundVegGrd(I,J)         ,& ! out,   vegetated ground heat [W/m2] (+ to soil/snow)
              HeatCanStorageChg       => noahmp%energy%flux%HeatCanStorageChg(I,J)        ,& ! out,   canopy heat storage change [W/m2]
              PhotosynLeafSunlit      => noahmp%biochem%flux%PhotosynLeafSunlit(I,J)      ,& ! out,   sunlit leaf photosynthesis [umol co2 /m2 /s]
              PhotosynLeafShade       => noahmp%biochem%flux%PhotosynLeafShade(I,J)        & ! out,   shaded leaf photosynthesis [umol co2 /m2 /s]
             )
! ----------------------------------------------------------------------

         ! initialization
         WindStressEwVeg     = 0.0
         WindStressNsVeg     = 0.0
         RadLwNetCanopy      = 0.0
         HeatSensibleCanopy  = 0.0
         RadLwNetVegGrd      = 0.0
         HeatSensibleVegGrd  = 0.0
         HeatLatentVegGrd    = 0.0
         HeatLatentCanEvap   = 0.0
         HeatLatentCanTransp = 0.0
         HeatGroundVegGrd    = 0.0
         PhotosynLeafSunlit  = 0.0
         PhotosynLeafShade   = 0.0
         TemperatureAir2mVeg = 0.0
         SpecHumidity2mVeg   = 0.0
         ExchCoeffShAbvCan   = 0.0
         ExchCoeffShLeaf     = 0.0
         ExchCoeffShUndCan   = 0.0
         ExchCoeffSh2mVeg    = 0.0
         HeatPrecipAdvSfc    = 0.0
         HeatCanStorageChg   = 0.0

         ! vegetated or non-vegetated
         VegAreaIndEff = LeafAreaIndEff + StemAreaIndEff
         FlagVegSfc(I,J)    = .false.
         if ( VegAreaIndEff > 0.0 ) FlagVegSfc(I,J) = .true.

         end associate
      end do
    end do
    !$acc end parallel loop

    ! ground snow cover fraction
    if ( noahmp%config%nmlist%OptSnowCoverGround == 1 ) call SnowCoverGroundNiu07(noahmp)
    if ( noahmp%config%nmlist%OptSnowCoverGround == 2 ) call SnowCoverGroundAR25(noahmp)

    ! ground and surface roughness length and reference height
    call GroundRoughnessProperty(noahmp, FlagVegSfc)

    ! Thermal properties of soil, snow, lake, and frozen soil
    call GroundThermalProperty(noahmp)

    ! Surface shortwave albedo: ground and canopy radiative transfer
    call SurfaceAlbedo(noahmp)

    ! Surface shortwave radiation: absorbed & reflected by the ground and canopy
    call SurfaceRadiation(noahmp)

    ! longwave emissivity for vegetation, ground, total net surface
    call SurfaceEmissivity(noahmp)

    ! soil water transpiration factor controlling stomatal resistance and evapotranspiration
    call SoilWaterTranspiration(noahmp)

    ! soil surface resistance for ground evaporation/sublimation
    call ResistanceGroundEvaporation(noahmp)

    ! set psychrometric variable/constant
    call PsychrometricVariable(noahmp)

    ! temperatures and energy fluxes of canopy, below-canopy ground, and bare ground
    !$acc parallel loop collapse(2) gang vector present(noahmp, FlagVegSfc)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         associate(                                                  &
              VegFrac                 => noahmp%energy%state%VegFrac(I,J)                 ,& ! in,    greeness vegetation fraction
              TemperatureGrd         => noahmp%energy%state%TemperatureGrd(I,J),         & ! inout, ground temperature [K]
              TemperatureGrdVeg    => noahmp%energy%state%TemperatureGrdVeg(I,J),     & ! inout, vegetation ground temperature [K]
              ExchCoeffMomSfc       => noahmp%energy%state%ExchCoeffMomSfc(I,J),       & ! inout, exchange coefficient [m/s] for momentum, surface, grid mean
              ExchCoeffMomAbvCan      => noahmp%energy%state%ExchCoeffMomAbvCan(I,J),      & ! out,   exchange coefficient [m/s] for momentum, above canopy
              ExchCoeffShSfc        => noahmp%energy%state%ExchCoeffShSfc(I,J)          ,& ! inout, exchange coefficient [m/s] for heat, surface, grid mean
              ExchCoeffShAbvCan       => noahmp%energy%state%ExchCoeffShAbvCan(I,J)         & ! out,   exchange coefficient [m/s] for heat, above canopy
             )
       if ( (FlagVegSfc(I,J) .eqv. .true.) .and. (VegFrac > 0) ) then ! vegetated portion of the grid
         TemperatureGrdVeg  = TemperatureGrd
         ExchCoeffMomAbvCan = ExchCoeffMomSfc
         ExchCoeffShAbvCan  = ExchCoeffShSfc
       endif

    end associate
      end do
    end do

    call SurfaceEnergyFluxVegetated(noahmp)

    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         associate(                                                  &
              TemperatureGrd         => noahmp%energy%state%TemperatureGrd(I,J),         & ! inout, ground temperature [K]
              TemperatureGrdBare    => noahmp%energy%state%TemperatureGrdBare(I,J),     & ! inout, bare ground temperature [K]
              ExchCoeffMomSfc       => noahmp%energy%state%ExchCoeffMomSfc(I,J),       & ! inout, exchange coefficient [m/s] for momentum, surface, grid mean
              ExchCoeffMomBare      => noahmp%energy%state%ExchCoeffMomBare(I,J),      & ! out,   exchange coefficient [m/s] for momentum, bare ground
              ExchCoeffShSfc        => noahmp%energy%state%ExchCoeffShSfc(I,J),         & ! inout, exchange coefficient [m/s] for heat, surface, grid mean
              ExchCoeffShBare       => noahmp%energy%state%ExchCoeffShBare(I,J)         & ! out,   exchange coefficient [m/s] for heat, bare ground
             )
    ! temperatures and energy fluxes of bare ground
    TemperatureGrdBare = TemperatureGrd
    ExchCoeffMomBare   = ExchCoeffMomSfc
    ExchCoeffShBare    = ExchCoeffShSfc

    end associate
      end do
    end do

    call SurfaceEnergyFluxBareGround(noahmp)

    ! Grid-level computations requiring 2D parallel loop
    !$acc parallel loop collapse(2) gang vector present(noahmp, FlagVegSfc)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                                             &
              NumSoilTimeStep         => noahmp%config%domain%NumSoilTimeStep        ,& ! in,    number of time step for calculating soil processes
              PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight(I,J)         ,& ! in,    air pressure [Pa] at reference height
              RadLwDownRefHeight      => noahmp%forcing%RadLwDownRefHeight(I,J)           ,& ! in,    downward longwave radiation [W/m2] at reference height
              RadSwDownRefHeight      => noahmp%forcing%RadSwDownRefHeight(I,J)           ,& ! in,    downward shortwave radiation [W/m2] at reference height
              OptSnowSoilTempTime     => noahmp%config%nmlist%OptSnowSoilTempTime         ,& ! in,    options for snow/soil temperature time scheme
              FlagCropland            => noahmp%config%domain%FlagCropland(I,J)           ,& ! in,    flag to identify croplands
              FlagSoilProcess         => noahmp%config%domain%FlagSoilProcess             ,& ! in,    flag to determine if calculating soil processes
              IrriFracThreshold       => noahmp%water%param%IrriFracThreshold             ,& ! in,    irrigation fraction parameter
              IrrigationFracGrid      => noahmp%water%state%IrrigationFracGrid(I,J)       ,& ! in,    total input irrigation fraction
              LeafAreaIndEff          => noahmp%energy%state%LeafAreaIndEff(I,J)          ,& ! in,    leaf area index, after burying by snow
              StemAreaIndEff          => noahmp%energy%state%StemAreaIndEff(I,J)          ,& ! in,    stem area index, after burying by snow
              VegFrac                 => noahmp%energy%state%VegFrac(I,J)                 ,& ! in,    greeness vegetation fraction
              HeatLatentIrriEvap      => noahmp%energy%flux%HeatLatentIrriEvap(I,J)       ,& ! in,    latent heating due to sprinkler evaporation [W/m2]
              HeatPrecipAdvCanopy     => noahmp%energy%flux%HeatPrecipAdvCanopy(I,J)      ,& ! in,    precipitation advected heat - vegetation net [W/m2]
              HeatPrecipAdvVegGrd     => noahmp%energy%flux%HeatPrecipAdvVegGrd(I,J)      ,& ! in,    precipitation advected heat - under canopy net [W/m2]
              HeatPrecipAdvBareGrd    => noahmp%energy%flux%HeatPrecipAdvBareGrd(I,J)     ,& ! in,    precipitation advected heat - bare ground net [W/m2]
              TemperatureSfc          => noahmp%energy%state%TemperatureSfc(I,J)          ,& ! inout, surface temperature [K]
              TemperatureGrd          => noahmp%energy%state%TemperatureGrd(I,J)          ,& ! inout, ground temperature [K]
              TemperatureCanopy       => noahmp%energy%state%TemperatureCanopy(I,J)       ,& ! inout, vegetation temperature [K]
              SpecHumiditySfc         => noahmp%energy%state%SpecHumiditySfc(I,J)         ,& ! inout, specific humidity [kg/kg] at bare/veg/urban surface
              SpecHumiditySfcMean     => noahmp%energy%state%SpecHumiditySfcMean(I,J)     ,& ! inout, specific humidity [kg/kg] at surface grid mean
              PressureVaporCanAir     => noahmp%energy%state%PressureVaporCanAir(I,J)     ,& ! inout, canopy air vapor pressure [Pa]
              ExchCoeffMomSfc         => noahmp%energy%state%ExchCoeffMomSfc(I,J)         ,& ! inout, exchange coefficient [m/s] for momentum, surface, grid mean
              ExchCoeffShSfc          => noahmp%energy%state%ExchCoeffShSfc(I,J)          ,& ! inout, exchange coefficient [m/s] for heat, surface, grid mean
              SnowDepth               => noahmp%water%state%SnowDepth(I,J)                ,& ! inout, snow depth [m]
              RoughLenMomSfcToAtm     => noahmp%energy%state%RoughLenMomSfcToAtm(I,J)     ,& ! out,   roughness length, momentum, surface, sent to coupled model
              WindStressEwSfc         => noahmp%energy%state%WindStressEwSfc(I,J)         ,& ! out,   wind stress: east-west [N/m2] grid mean
              WindStressNsSfc         => noahmp%energy%state%WindStressNsSfc(I,J)         ,& ! out,   wind stress: north-south [N/m2] grid mean
              TemperatureRadSfc       => noahmp%energy%state%TemperatureRadSfc(I,J)       ,& ! out,   surface radiative temperature [K]
              TemperatureAir2m        => noahmp%energy%state%TemperatureAir2m(I,J)        ,& ! out,   grid mean 2-m air temperature [K]
              ResistanceStomataSunlit => noahmp%energy%state%ResistanceStomataSunlit(I,J) ,& ! out,   sunlit leaf stomatal resistance [s/m]
              ResistanceStomataShade  => noahmp%energy%state%ResistanceStomataShade(I,J)  ,& ! out,   shaded leaf stomatal resistance [s/m]
              TemperatureAir2mVeg     => noahmp%energy%state%TemperatureAir2mVeg(I,J)     ,& ! out,   2 m height air temperature [K], vegetated
              TemperatureAir2mBare    => noahmp%energy%state%TemperatureAir2mBare(I,J)    ,& ! out,   2 m height air temperature [K] bare ground
              LeafAreaIndSunlit       => noahmp%energy%state%LeafAreaIndSunlit(I,J)       ,& ! out,   sunlit leaf area index, one-sided [m2/m2]
              LeafAreaIndShade        => noahmp%energy%state%LeafAreaIndShade(I,J)        ,& ! out,   shaded leaf area index, one-sided [m2/m2]
              EmissivitySfc           => noahmp%energy%state%EmissivitySfc(I,J)           ,& ! out,   surface emissivity
              VegAreaIndEff           => noahmp%energy%state%VegAreaIndEff(I,J)           ,& ! out,   one-sided leaf+stem area index [m2/m2]
              RoughLenMomSfc          => noahmp%energy%state%RoughLenMomSfc(I,J)          ,& ! out,   roughness length [m], momentum, surface
              RoughLenMomGrd          => noahmp%energy%state%RoughLenMomGrd(I,J)          ,& ! out,   roughness length [m], momentum, ground
              WindStressEwVeg         => noahmp%energy%state%WindStressEwVeg(I,J)         ,& ! out,   wind stress: east-west [N/m2] above canopy
              WindStressNsVeg         => noahmp%energy%state%WindStressNsVeg(I,J)         ,& ! out,   wind stress: north-south [N/m2] above canopy
              WindStressEwBare        => noahmp%energy%state%WindStressEwBare(I,J)        ,& ! out,   wind stress: east-west [N/m2] bare ground
              WindStressNsBare        => noahmp%energy%state%WindStressNsBare(I,J)        ,& ! out,   wind stress: north-south [N/m2] bare ground
              SpecHumidity2mVeg       => noahmp%energy%state%SpecHumidity2mVeg(I,J)       ,& ! out,   water vapor mixing ratio at 2m vegetated
              SpecHumidity2mBare      => noahmp%energy%state%SpecHumidity2mBare(I,J)      ,& ! out,   bare ground 2-m water vapor mixing ratio
              SpecHumidity2m          => noahmp%energy%state%SpecHumidity2m(I,J)          ,& ! out,   grid mean 2-m water vapor mixing ratio
              TemperatureGrdVeg       => noahmp%energy%state%TemperatureGrdVeg(I,J)       ,& ! out,   vegetated ground (below-canopy) temperature [K]
              TemperatureGrdBare      => noahmp%energy%state%TemperatureGrdBare(I,J)      ,& ! out,   bare ground temperature [K]
              ExchCoeffMomAbvCan      => noahmp%energy%state%ExchCoeffMomAbvCan(I,J)      ,& ! out,   exchange coeff [m/s] for momentum, above ZeroPlaneDisp, vegetated
              ExchCoeffMomBare        => noahmp%energy%state%ExchCoeffMomBare(I,J)        ,& ! out,   exchange coeff [m/s] for momentum, above ZeroPlaneDisp, bare ground
              ExchCoeffShAbvCan       => noahmp%energy%state%ExchCoeffShAbvCan(I,J)       ,& ! out,   exchange coeff [m/s] for heat, above ZeroPlaneDisp, vegetated
              ExchCoeffShBare         => noahmp%energy%state%ExchCoeffShBare(I,J)         ,& ! out,   exchange coeff [m/s] for heat, above ZeroPlaneDisp, bare ground
              ExchCoeffShLeaf         => noahmp%energy%state%ExchCoeffShLeaf(I,J)         ,& ! out,   leaf sensible heat exchange coeff [m/s], leaf to canopy air
              ExchCoeffShUndCan       => noahmp%energy%state%ExchCoeffShUndCan(I,J)       ,& ! out,   under canopy sensible heat exchange coefficient [m/s]
              ExchCoeffSh2mVeg        => noahmp%energy%state%ExchCoeffSh2mVeg(I,J)        ,& ! out,   2m sensible heat exchange coefficient [m/s] vegetated
              AlbedoSfc               => noahmp%energy%state%AlbedoSfc(I,J)               ,& ! out,   total shortwave surface albedo
              RadSwReflSfc            => noahmp%energy%flux%RadSwReflSfc(I,J)             ,& ! out,   total reflected solar radiation [W/m2]
              RadLwNetSfc             => noahmp%energy%flux%RadLwNetSfc(I,J)              ,& ! out,   total net longwave rad [W/m2] (+ to atm)
              HeatSensibleSfc         => noahmp%energy%flux%HeatSensibleSfc(I,J)          ,& ! out,   total sensible heat [W/m2] (+ to atm)
              HeatLatentGrd           => noahmp%energy%flux%HeatLatentGrd(I,J)            ,& ! out,   total ground latent heat [W/m2] (+ to atm)
              HeatLatentCanopy        => noahmp%energy%flux%HeatLatentCanopy(I,J)         ,& ! out,   canopy latent heat flux [W/m2] (+ to atm)
              HeatLatentTransp        => noahmp%energy%flux%HeatLatentTransp(I,J)         ,& ! out,   latent heat flux from transpiration [W/m2] (+ to atm)
              RadPhotoActAbsCan       => noahmp%energy%flux%RadPhotoActAbsCan(I,J)        ,& ! out,   total photosyn. active energy [W/m2) absorbed by canopy
              RadPhotoActAbsSunlit    => noahmp%energy%flux%RadPhotoActAbsSunlit(I,J)     ,& ! out,   average absorbed par for sunlit leaves [W/m2]
              RadPhotoActAbsShade     => noahmp%energy%flux%RadPhotoActAbsShade(I,J)      ,& ! out,   average absorbed par for shaded leaves [W/m2]
              HeatGroundTot           => noahmp%energy%flux%HeatGroundTot(I,J)            ,& ! out,   total ground heat flux [W/m2] (+ to soil/snow)
              HeatPrecipAdvSfc        => noahmp%energy%flux%HeatPrecipAdvSfc(I,J)         ,& ! out,   precipitation advected heat - total [W/m2]
              RadLwEmitSfc            => noahmp%energy%flux%RadLwEmitSfc(I,J)             ,& ! out,   emitted outgoing IR [W/m2]
              RadLwNetCanopy          => noahmp%energy%flux%RadLwNetCanopy(I,J)           ,& ! out,   canopy net longwave radiation [W/m2] (+ to atm)
              RadLwNetVegGrd          => noahmp%energy%flux%RadLwNetVegGrd(I,J)           ,& ! out,   ground net longwave radiation [W/m2] (+ to atm)
              RadLwNetBareGrd         => noahmp%energy%flux%RadLwNetBareGrd(I,J)          ,& ! out,   net longwave rad [W/m2] bare ground (+ to atm)
              HeatSensibleCanopy      => noahmp%energy%flux%HeatSensibleCanopy(I,J)       ,& ! out,   canopy sensible heat flux [W/m2] (+ to atm)
              HeatSensibleVegGrd      => noahmp%energy%flux%HeatSensibleVegGrd(I,J)       ,& ! out,   vegetated ground sensible heat flux [W/m2] (+ to atm)
              HeatSensibleBareGrd     => noahmp%energy%flux%HeatSensibleBareGrd(I,J)      ,& ! out,   sensible heat flux [W/m2] bare ground (+ to atm)
              HeatLatentVegGrd        => noahmp%energy%flux%HeatLatentVegGrd(I,J)         ,& ! out,   ground evaporation heat flux [W/m2] (+ to atm)
              HeatLatentBareGrd       => noahmp%energy%flux%HeatLatentBareGrd(I,J)        ,& ! out,   latent heat flux [W/m2] bare ground (+ to atm)
              HeatLatentCanEvap       => noahmp%energy%flux%HeatLatentCanEvap(I,J)        ,& ! out,   canopy evaporation heat flux [W/m2] (+ to atm)
              HeatLatentCanTransp     => noahmp%energy%flux%HeatLatentCanTransp(I,J)      ,& ! out,   canopy transpiration heat flux [W/m2] (+ to atm)
              HeatGroundVegGrd        => noahmp%energy%flux%HeatGroundVegGrd(I,J)         ,& ! out,   vegetated ground heat [W/m2] (+ to soil/snow)
              HeatGroundBareGrd       => noahmp%energy%flux%HeatGroundBareGrd(I,J)        ,& ! out,   bare ground heat flux [W/m2] (+ to soil/snow)
              HeatCanStorageChg       => noahmp%energy%flux%HeatCanStorageChg(I,J)        ,& ! out,   canopy heat storage change [W/m2]
              HeatFromSoilBot         => noahmp%energy%flux%HeatFromSoilBot(I,J)          ,& ! out,   energy influx from soil bottom [J/m2] during soil timestep
              HeatGroundTotMean       => noahmp%energy%flux%HeatGroundTotMean(I,J)        ,& ! out,   mean ground heat flux during soil timestep [W/m2]
              HeatGroundTotAcc        => noahmp%energy%flux%HeatGroundTotAcc(I,J)         ,& ! inout, accumulated total ground heat flux per soil timestep [W/m2 * dt_soil/dt_main]
              PhotosynTotal           => noahmp%biochem%flux%PhotosynTotal(I,J)           ,& ! out,   total leaf photosynthesis [umol co2 /m2 /s]
              PhotosynLeafSunlit      => noahmp%biochem%flux%PhotosynLeafSunlit(I,J)      ,& ! out,   sunlit leaf photosynthesis [umol co2 /m2 /s]
              PhotosynLeafShade       => noahmp%biochem%flux%PhotosynLeafShade(I,J)        & ! out,   shaded leaf photosynthesis [umol co2 /m2 /s]
             )
! ----------------------------------------------------------------------

    ! compute grid mean quantities by weighting vegetated and bare portions
    ! Energy balance at vege canopy: 
    ! RadSwAbsVeg = (RadLwNetCanopy + HeatSensibleCanopy + HeatLatentCanEvap + HeatLatentCanTransp) * VegFrac at VegFrac 
    ! Energy balance at vege ground: 
    ! RadSwAbsGrd * VegFrac = (RadLwNetVegGrd + HeatSensibleVegGrd + HeatLatentVegGrd + HeatGroundVegGrd) * VegFrac at VegFrac
    ! Energy balance at bare ground: 
    ! RadSwAbsGrd * (1-VegFrac) = (RadLwNetBareGrd + HeatSensibleBareGrd + HeatLatentBareGrd + HeatGroundBareGrd) * (1-VegFrac) at 1-VegFrac
    if ( (FlagVegSfc(I,J) .eqv. .true.) .and. (VegFrac > 0) ) then
       WindStressEwSfc     = VegFrac * WindStressEwVeg     + (1.0 - VegFrac) * WindStressEwBare
       WindStressNsSfc     = VegFrac * WindStressNsVeg     + (1.0 - VegFrac) * WindStressNsBare
       RadLwNetSfc         = VegFrac * RadLwNetVegGrd      + (1.0 - VegFrac) * RadLwNetBareGrd     + RadLwNetCanopy
       HeatSensibleSfc     = VegFrac * HeatSensibleVegGrd  + (1.0 - VegFrac) * HeatSensibleBareGrd + HeatSensibleCanopy
       HeatLatentGrd       = VegFrac * HeatLatentVegGrd    + (1.0 - VegFrac) * HeatLatentBareGrd
       HeatGroundTot       = VegFrac * HeatGroundVegGrd    + (1.0 - VegFrac) * HeatGroundBareGrd
       HeatLatentCanopy    = HeatLatentCanEvap
       HeatLatentTransp    = HeatLatentCanTransp
       HeatPrecipAdvSfc    = VegFrac * HeatPrecipAdvVegGrd + (1.0 - VegFrac) * HeatPrecipAdvBareGrd + HeatPrecipAdvCanopy
       TemperatureGrd      = VegFrac * TemperatureGrdVeg   + (1.0 - VegFrac) * TemperatureGrdBare
       TemperatureAir2m    = VegFrac * TemperatureAir2mVeg + (1.0 - VegFrac) * TemperatureAir2mBare
       TemperatureSfc      = VegFrac * TemperatureCanopy   + (1.0 - VegFrac) * TemperatureGrdBare
       ExchCoeffMomSfc     = VegFrac * ExchCoeffMomAbvCan  + (1.0 - VegFrac) * ExchCoeffMomBare     ! better way to average?
       ExchCoeffShSfc      = VegFrac * ExchCoeffShAbvCan   + (1.0 - VegFrac) * ExchCoeffShBare
       SpecHumidity2m      = VegFrac * SpecHumidity2mVeg   + (1.0 - VegFrac) * SpecHumidity2mBare 
       SpecHumiditySfcMean = VegFrac * (PressureVaporCanAir * 0.622 / &
                             (PressureAirRefHeight - 0.378*PressureVaporCanAir)) + (1.0 - VegFrac) * SpecHumiditySfc
       RoughLenMomSfcToAtm = RoughLenMomSfc
    else
       WindStressEwSfc         = WindStressEwBare
       WindStressNsSfc         = WindStressNsBare
       RadLwNetSfc             = RadLwNetBareGrd
       HeatSensibleSfc         = HeatSensibleBareGrd
       HeatLatentGrd           = HeatLatentBareGrd
       HeatGroundTot           = HeatGroundBareGrd
       TemperatureGrd          = TemperatureGrdBare
       TemperatureAir2m        = TemperatureAir2mBare
       HeatLatentCanopy        = 0.0
       HeatLatentTransp        = 0.0
       HeatPrecipAdvSfc        = HeatPrecipAdvBareGrd
       TemperatureSfc          = TemperatureGrd
       ExchCoeffMomSfc         = ExchCoeffMomBare
       ExchCoeffShSfc          = ExchCoeffShBare
       SpecHumiditySfcMean     = SpecHumiditySfc
       SpecHumidity2m          = SpecHumidity2mBare
       ResistanceStomataSunlit = 0.0
       ResistanceStomataShade  = 0.0
       TemperatureGrdVeg       = TemperatureGrdBare
       ExchCoeffShAbvCan       = ExchCoeffShBare
       RoughLenMomSfcToAtm     = RoughLenMomGrd
    endif

    ! emitted longwave radiation and physical check
    RadLwEmitSfc = RadLwDownRefHeight + RadLwNetSfc
#ifndef _OPENACC
    if ( RadLwEmitSfc <= 0.0 ) then
       write(*,*) "emitted longwave <0; skin T may be wrong due to inconsistent"
       write(*,*) "input of VegFracGreen with LeafAreaIndex"
       write(*,*) "VegFrac = ", VegFrac, "VegAreaIndEff = ", VegAreaIndEff, &
                  "TemperatureCanopy = ", TemperatureCanopy, "TemperatureGrd = ", TemperatureGrd
       write(*,*) "RadLwDownRefHeight = ", RadLwDownRefHeight, "RadLwNetSfc = ", RadLwNetSfc, "SnowDepth = ", SnowDepth
       stop "Error: Longwave radiation budget problem in NoahMP LSM"
    endif
#endif

    ! radiative temperature: subtract from the emitted IR the
    ! reflected portion of the incoming longwave radiation, so just
    ! considering the IR originating/emitted in the canopy/ground system.
    ! Old TemperatureRadSfc calculation not taking into account Emissivity:
    ! TemperatureRadSfc = (RadLwEmitSfc/ConstStefanBoltzmann)**0.25
    TemperatureRadSfc = ((RadLwEmitSfc - (1.0-EmissivitySfc)*RadLwDownRefHeight) / (EmissivitySfc*ConstStefanBoltzmann))**0.25

    ! other photosynthesis related quantities for biochem process
    RadPhotoActAbsCan = RadPhotoActAbsSunlit * LeafAreaIndSunlit + RadPhotoActAbsShade * LeafAreaIndShade
    PhotosynTotal     = PhotosynLeafSunlit   * LeafAreaIndSunlit + PhotosynLeafShade   * LeafAreaIndShade

    ! compute snow and soil layer temperature at soil timestep
    HeatFromSoilBot = 0.0
    HeatGroundTotAcc = HeatGroundTotAcc + HeatGroundTot
    if ( FlagSoilProcess .eqv. .true. ) HeatGroundTotMean = HeatGroundTotAcc / NumSoilTimeStep

    end associate

      end do
    end do

    if ( noahmp%config%domain%FlagSoilProcess .eqv. .true. ) then
       call SoilSnowTemperatureMain(noahmp)
    endif ! FlagSoilProcess

    !$acc parallel loop collapse(2) gang vector present(noahmpm, FlagVegSfc)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         associate(                                                  &
              OptSnowSoilTempTime     => noahmp%config%nmlist%OptSnowSoilTempTime         ,& ! in,    options for snow/soil temperature time scheme
              VegFrac                 => noahmp%energy%state%VegFrac(I,J)                 ,& ! in,    greeness vegetation fraction
              TemperatureGrdVeg       => noahmp%energy%state%TemperatureGrdVeg(I,J)       ,& ! out,   vegetated ground (below-canopy) temperature [K]
              TemperatureGrdBare      => noahmp%energy%state%TemperatureGrdBare(I,J)      ,& ! out,   bare ground temperature [K]
              TemperatureSfc          => noahmp%energy%state%TemperatureSfc(I,J)          ,& ! inout, surface temperature [K]
              TemperatureGrd          => noahmp%energy%state%TemperatureGrd(I,J)          ,& ! inout, ground temperature [K]
              TemperatureCanopy       => noahmp%energy%state%TemperatureCanopy(I,J)       ,& ! inout, vegetation temperature [K]
              SnowDepth               => noahmp%water%state%SnowDepth(I,J)                ,& ! inout, snow depth [m]
              RadSwReflSfc        => noahmp%energy%flux%RadSwReflSfc(I,J)           & ! out,   total reflected solar radiation [W/m2]
             )

    ! adjusting suface temperature based on snow condition
    if ( OptSnowSoilTempTime == 2 ) then
       if ( (SnowDepth > 0.05) .and. (TemperatureGrd > ConstFreezePoint) ) then
          TemperatureGrdVeg  = ConstFreezePoint
          TemperatureGrdBare = ConstFreezePoint
          if ( (FlagVegSfc(I,J) .eqv. .true.) .and. (VegFrac > 0) ) then
             TemperatureGrd  = VegFrac * TemperatureGrdVeg + (1.0 - VegFrac) * TemperatureGrdBare
             TemperatureSfc  = VegFrac * TemperatureCanopy + (1.0 - VegFrac) * TemperatureGrdBare
          else
             TemperatureGrd  = TemperatureGrdBare
             TemperatureSfc  = TemperatureGrdBare
          endif
       endif
    endif

    end associate

      end do
    end do
    !$acc end parallel loop

    ! Phase change and Energy released or consumed by snow & frozen soil
    call SoilSnowWaterPhaseChange(noahmp)

    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         associate(                                                  &
              FlagCropland        => noahmp%config%domain%FlagCropland(I,J)       ,& ! in,    flag to identify croplands
              IrrigationFracGrid  => noahmp%water%state%IrrigationFracGrid(I,J)   ,& ! in,    total input irrigation fraction
              IrriFracThreshold   => noahmp%water%param%IrriFracThreshold(I,J)         ,& ! in,    irrigation fraction parameter
              HeatLatentIrriEvap  => noahmp%energy%flux%HeatLatentIrriEvap(I,J)   ,& ! in,    latent heating due to sprinkler evaporation [W/m2]
              HeatSensibleSfc     => noahmp%energy%flux%HeatSensibleSfc(I,J)      ,& ! inout, total sensible heat [W/m2] (+ to atm)
              RadSwDownRefHeight  => noahmp%forcing%RadSwDownRefHeight(I,J)       ,& ! in,    downward shortwave radiation [W/m2] at reference height
              AlbedoSfc           => noahmp%energy%state%AlbedoSfc(I,J)           ,& ! out,   total shortwave surface albedo
              RadSwReflSfc        => noahmp%energy%flux%RadSwReflSfc(I,J)          & ! out,   total reflected solar radiation [W/m2]
             )
    ! update sensible heat flux due to sprinkler irrigation evaporation
    if ( (FlagCropland .eqv. .true.) .and. (IrrigationFracGrid >= IrriFracThreshold) ) &
       HeatSensibleSfc = HeatSensibleSfc - HeatLatentIrriEvap

    ! update total surface albedo
    if ( RadSwDownRefHeight > 0.0 ) then
       AlbedoSfc = RadSwReflSfc / RadSwDownRefHeight
    else
       AlbedoSfc = undefined_real
    endif

    end associate

      end do
    end do
    !$acc end parallel loop
    !$acc end data

    ! Deallocate local array
    deallocate(FlagVegSfc)

  end subroutine EnergyMain

end module EnergyMainMod
