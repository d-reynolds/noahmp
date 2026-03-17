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
  use GlacierPhaseChangeMod,             only : GlacierPhaseChange
  use GlacierTemperatureMainMod,          only : GlacierTemperatureMain
  use SnowCoverGlacierMod,               only : SnowCoverGlacier
  use GroundRoughnessPropertyGlacierMod, only : GroundRoughnessPropertyGlacier
  use GroundThermalPropertyGlacierMod,   only : GroundThermalPropertyGlacier
  use GroundAlbedoGlacierMod,            only : GroundAlbedoGlacier

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
    associate(                                                                      &
              LeafAreaIndEff          => noahmp%energy%state%LeafAreaIndEff ,& ! in,    leaf area index, after burying by snow
              StemAreaIndEff          => noahmp%energy%state%StemAreaIndEff ,& ! in,    stem area index, after burying by snow
              TemperatureAir2mVeg     => noahmp%energy%state%TemperatureAir2mVeg ,& ! out,   2 m height air temperature [K], vegetated
              VegAreaIndEff           => noahmp%energy%state%VegAreaIndEff ,& ! out,   one-sided leaf+stem area index [m2/m2]
              WindStressEwVeg         => noahmp%energy%state%WindStressEwVeg ,& ! out,   wind stress: east-west [N/m2] above canopy
              WindStressNsVeg         => noahmp%energy%state%WindStressNsVeg ,& ! out,   wind stress: north-south [N/m2] above canopy
              SpecHumidity2mVeg       => noahmp%energy%state%SpecHumidity2mVeg ,& ! out,   water vapor mixing ratio at 2m vegetated
              ExchCoeffShAbvCan       => noahmp%energy%state%ExchCoeffShAbvCan ,& ! out,   exchange coeff [m/s] for heat, above ZeroPlaneDisp, vegetated
              ExchCoeffShLeaf         => noahmp%energy%state%ExchCoeffShLeaf ,& ! out,   leaf sensible heat exchange coeff [m/s], leaf to canopy air
              ExchCoeffShUndCan       => noahmp%energy%state%ExchCoeffShUndCan ,& ! out,   under canopy sensible heat exchange coefficient [m/s]
              ExchCoeffSh2mVeg        => noahmp%energy%state%ExchCoeffSh2mVeg ,& ! out,   2m sensible heat exchange coefficient [m/s] vegetated
              HeatPrecipAdvSfc        => noahmp%energy%flux%HeatPrecipAdvSfc ,& ! out,   precipitation advected heat - total [W/m2]
              RadLwNetCanopy          => noahmp%energy%flux%RadLwNetCanopy ,& ! out,   canopy net longwave radiation [W/m2] (+ to atm)
              RadLwNetVegGrd          => noahmp%energy%flux%RadLwNetVegGrd ,& ! out,   ground net longwave radiation [W/m2] (+ to atm)
              HeatSensibleCanopy      => noahmp%energy%flux%HeatSensibleCanopy ,& ! out,   canopy sensible heat flux [W/m2] (+ to atm)
              HeatSensibleVegGrd      => noahmp%energy%flux%HeatSensibleVegGrd ,& ! out,   vegetated ground sensible heat flux [W/m2] (+ to atm)
              HeatLatentVegGrd        => noahmp%energy%flux%HeatLatentVegGrd ,& ! out,   ground evaporation heat flux [W/m2] (+ to atm)
              HeatLatentCanEvap       => noahmp%energy%flux%HeatLatentCanEvap ,& ! out,   canopy evaporation heat flux [W/m2] (+ to atm)
              HeatLatentCanTransp     => noahmp%energy%flux%HeatLatentCanTransp ,& ! out,   canopy transpiration heat flux [W/m2] (+ to atm)
              HeatGroundVegGrd        => noahmp%energy%flux%HeatGroundVegGrd ,& ! out,   vegetated ground heat [W/m2] (+ to soil/snow)
              HeatCanStorageChg       => noahmp%energy%flux%HeatCanStorageChg ,& ! out,   canopy heat storage change [W/m2]
              PhotosynLeafSunlit      => noahmp%biochem%flux%PhotosynLeafSunlit ,& ! out,   sunlit leaf photosynthesis [umol co2 /m2 /s]
              PhotosynLeafShade       => noahmp%biochem%flux%PhotosynLeafShade ,& ! out,   shaded leaf photosynthesis [umol co2 /m2 /s]
              VegFrac                 => noahmp%energy%state%VegFrac ,& ! in,    greeness vegetation fraction
              TemperatureGrd         => noahmp%energy%state%TemperatureGrd ,& ! inout, ground temperature [K]
              TemperatureGrdVeg    => noahmp%energy%state%TemperatureGrdVeg ,& ! inout, vegetation ground temperature [K]
              ExchCoeffMomSfc       => noahmp%energy%state%ExchCoeffMomSfc ,& ! inout, exchange coefficient [m/s] for momentum, surface, grid mean
              ExchCoeffMomAbvCan      => noahmp%energy%state%ExchCoeffMomAbvCan ,& ! out,   exchange coefficient [m/s] for momentum, above canopy
              ExchCoeffShSfc        => noahmp%energy%state%ExchCoeffShSfc ,& ! inout, exchange coefficient [m/s] for heat, surface, grid mean
              TemperatureGrdBare    => noahmp%energy%state%TemperatureGrdBare ,& ! inout, bare ground temperature [K]
              ExchCoeffMomBare      => noahmp%energy%state%ExchCoeffMomBare ,& ! out,   exchange coefficient [m/s] for momentum, bare ground
              ExchCoeffShBare       => noahmp%energy%state%ExchCoeffShBare ,& ! out,   exchange coefficient [m/s] for heat, bare ground
              NumSoilTimeStep         => noahmp%config%domain%NumSoilTimeStep ,& ! in,    number of time step for calculating soil processes
              PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight ,& ! in,    air pressure [Pa] at reference height
              RadLwDownRefHeight      => noahmp%forcing%RadLwDownRefHeight ,& ! in,    downward longwave radiation [W/m2] at reference height
              RadSwDownRefHeight      => noahmp%forcing%RadSwDownRefHeight ,& ! in,    downward shortwave radiation [W/m2] at reference height
              OptSnowSoilTempTime     => noahmp%config%nmlist%OptSnowSoilTempTime ,& ! in,    options for snow/soil temperature time scheme
              FlagCropland            => noahmp%config%domain%FlagCropland ,& ! in,    flag to identify croplands
              FlagSoilProcess         => noahmp%config%domain%FlagSoilProcess ,& ! in,    flag to determine if calculating soil processes
              IrriFracThreshold       => noahmp%water%param%IrriFracThreshold ,& ! in,    irrigation fraction parameter
              IrrigationFracGrid      => noahmp%water%state%IrrigationFracGrid ,& ! in,    total input irrigation fraction
              HeatLatentIrriEvap      => noahmp%energy%flux%HeatLatentIrriEvap ,& ! in,    latent heating due to sprinkler evaporation [W/m2]
              HeatPrecipAdvCanopy     => noahmp%energy%flux%HeatPrecipAdvCanopy ,& ! in,    precipitation advected heat - vegetation net [W/m2]
              HeatPrecipAdvVegGrd     => noahmp%energy%flux%HeatPrecipAdvVegGrd ,& ! in,    precipitation advected heat - under canopy net [W/m2]
              HeatPrecipAdvBareGrd    => noahmp%energy%flux%HeatPrecipAdvBareGrd ,& ! in,    precipitation advected heat - bare ground net [W/m2]
              TemperatureSfc          => noahmp%energy%state%TemperatureSfc ,& ! inout, surface temperature [K]
              TemperatureCanopy       => noahmp%energy%state%TemperatureCanopy ,& ! inout, vegetation temperature [K]
              SpecHumiditySfc         => noahmp%energy%state%SpecHumiditySfc ,& ! inout, specific humidity [kg/kg] at bare/veg/urban surface
              SpecHumiditySfcMean     => noahmp%energy%state%SpecHumiditySfcMean ,& ! inout, specific humidity [kg/kg] at surface grid mean
              PressureVaporCanAir     => noahmp%energy%state%PressureVaporCanAir ,& ! inout, canopy air vapor pressure [Pa]
              SnowDepth               => noahmp%water%state%SnowDepth ,& ! inout, snow depth [m]
              RoughLenMomSfcToAtm     => noahmp%energy%state%RoughLenMomSfcToAtm ,& ! out,   roughness length, momentum, surface, sent to coupled model
              WindStressEwSfc         => noahmp%energy%state%WindStressEwSfc ,& ! out,   wind stress: east-west [N/m2] grid mean
              WindStressNsSfc         => noahmp%energy%state%WindStressNsSfc ,& ! out,   wind stress: north-south [N/m2] grid mean
              TemperatureRadSfc       => noahmp%energy%state%TemperatureRadSfc ,& ! out,   surface radiative temperature [K]
              TemperatureAir2m        => noahmp%energy%state%TemperatureAir2m ,& ! out,   grid mean 2-m air temperature [K]
              ResistanceStomataSunlit => noahmp%energy%state%ResistanceStomataSunlit ,& ! out,   sunlit leaf stomatal resistance [s/m]
              ResistanceStomataShade  => noahmp%energy%state%ResistanceStomataShade ,& ! out,   shaded leaf stomatal resistance [s/m]
              TemperatureAir2mBare    => noahmp%energy%state%TemperatureAir2mBare ,& ! out,   2 m height air temperature [K] bare ground
              LeafAreaIndSunlit       => noahmp%energy%state%LeafAreaIndSunlit ,& ! out,   sunlit leaf area index, one-sided [m2/m2]
              LeafAreaIndShade        => noahmp%energy%state%LeafAreaIndShade ,& ! out,   shaded leaf area index, one-sided [m2/m2]
              EmissivitySfc           => noahmp%energy%state%EmissivitySfc ,& ! out,   surface emissivity
              RoughLenMomSfc          => noahmp%energy%state%RoughLenMomSfc ,& ! out,   roughness length [m], momentum, surface
              RoughLenMomGrd          => noahmp%energy%state%RoughLenMomGrd ,& ! out,   roughness length [m], momentum, ground
              WindStressEwBare        => noahmp%energy%state%WindStressEwBare ,& ! out,   wind stress: east-west [N/m2] bare ground
              WindStressNsBare        => noahmp%energy%state%WindStressNsBare ,& ! out,   wind stress: north-south [N/m2] bare ground
              SpecHumidity2mBare      => noahmp%energy%state%SpecHumidity2mBare ,& ! out,   bare ground 2-m water vapor mixing ratio
              SpecHumidity2m          => noahmp%energy%state%SpecHumidity2m ,& ! out,   grid mean 2-m water vapor mixing ratio
              AlbedoSfc               => noahmp%energy%state%AlbedoSfc ,& ! out,   total shortwave surface albedo
              RadSwReflSfc            => noahmp%energy%flux%RadSwReflSfc ,& ! out,   total reflected solar radiation [W/m2]
              RadLwNetSfc             => noahmp%energy%flux%RadLwNetSfc ,& ! out,   total net longwave rad [W/m2] (+ to atm)
              HeatSensibleSfc         => noahmp%energy%flux%HeatSensibleSfc ,& ! out,   total sensible heat [W/m2] (+ to atm)
              HeatLatentGrd           => noahmp%energy%flux%HeatLatentGrd ,& ! out,   total ground latent heat [W/m2] (+ to atm)
              HeatLatentCanopy        => noahmp%energy%flux%HeatLatentCanopy ,& ! out,   canopy latent heat flux [W/m2] (+ to atm)
              HeatLatentTransp        => noahmp%energy%flux%HeatLatentTransp ,& ! out,   latent heat flux from transpiration [W/m2] (+ to atm)
              RadPhotoActAbsCan       => noahmp%energy%flux%RadPhotoActAbsCan ,& ! out,   total photosyn. active energy [W/m2) absorbed by canopy
              IndicatorIceSfc         => noahmp%config%domain%IndicatorIceSfc ,& ! in,    flag to identify ice surface
              RadPhotoActAbsSunlit    => noahmp%energy%flux%RadPhotoActAbsSunlit ,& ! out,   average absorbed par for sunlit leaves [W/m2]
              RadPhotoActAbsShade     => noahmp%energy%flux%RadPhotoActAbsShade ,& ! out,   average absorbed par for shaded leaves [W/m2]
              HeatGroundTot           => noahmp%energy%flux%HeatGroundTot ,& ! out,   total ground heat flux [W/m2] (+ to soil/snow)
              RadLwEmitSfc            => noahmp%energy%flux%RadLwEmitSfc ,& ! out,   emitted outgoing IR [W/m2]
              RadLwNetBareGrd         => noahmp%energy%flux%RadLwNetBareGrd ,& ! out,   net longwave rad [W/m2] bare ground (+ to atm)
              HeatSensibleBareGrd     => noahmp%energy%flux%HeatSensibleBareGrd ,& ! out,   sensible heat flux [W/m2] bare ground (+ to atm)
              HeatLatentBareGrd       => noahmp%energy%flux%HeatLatentBareGrd ,& ! out,   latent heat flux [W/m2] bare ground (+ to atm)
              HeatGroundBareGrd       => noahmp%energy%flux%HeatGroundBareGrd ,& ! out,   bare ground heat flux [W/m2] (+ to soil/snow)
              HeatFromSoilBot         => noahmp%energy%flux%HeatFromSoilBot ,& ! out,   energy influx from soil bottom [J/m2] during soil timestep
              HeatGroundTotMean       => noahmp%energy%flux%HeatGroundTotMean ,& ! out,   mean ground heat flux during soil timestep [W/m2]
              HeatGroundTotAcc        => noahmp%energy%flux%HeatGroundTotAcc ,& ! inout, accumulated total ground heat flux per soil timestep [W/m2 * dt_soil/dt_main]
              PhotosynTotal           => noahmp%biochem%flux%PhotosynTotal  & ! out,   total leaf photosynthesis [umol co2 /m2 /s]
             )

    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


         ! initialization
         WindStressEwVeg(I,J)     = 0.0
         WindStressNsVeg(I,J)     = 0.0
         RadLwNetCanopy(I,J)      = 0.0
         HeatSensibleCanopy(I,J)  = 0.0
         RadLwNetVegGrd(I,J)      = 0.0
         HeatSensibleVegGrd(I,J)  = 0.0
         HeatLatentVegGrd(I,J)    = 0.0
         HeatLatentCanEvap(I,J)   = 0.0
         HeatLatentCanTransp(I,J) = 0.0
         HeatGroundVegGrd(I,J)    = 0.0
         PhotosynLeafSunlit(I,J)  = 0.0
         PhotosynLeafShade(I,J)   = 0.0
         TemperatureAir2mVeg(I,J) = 0.0
         SpecHumidity2mVeg(I,J)   = 0.0
         ExchCoeffShAbvCan(I,J)   = 0.0
         ExchCoeffShLeaf(I,J)     = 0.0
         ExchCoeffShUndCan(I,J)   = 0.0
         ExchCoeffSh2mVeg(I,J)    = 0.0
         HeatPrecipAdvSfc(I,J)    = 0.0
         HeatCanStorageChg(I,J)   = 0.0

         ! vegetated or non-vegetated
         VegAreaIndEff(I,J) = LeafAreaIndEff(I,J) + StemAreaIndEff(I,J)
         FlagVegSfc(I,J)    = .false.
         if ( VegAreaIndEff(I,J) > 0.0 ) FlagVegSfc(I,J) = .true.

      end do
    end do
    !$acc end parallel loop

    ! ground snow cover fraction
    if ( noahmp%config%nmlist%OptSnowCoverGround == 1 ) call SnowCoverGroundNiu07(noahmp)
    if ( noahmp%config%nmlist%OptSnowCoverGround == 2 ) call SnowCoverGroundAR25(noahmp)

    ! glaicer snow cover fraction
    call SnowCoverGlacier(noahmp)

    ! ground and surface roughness length and reference height
    call GroundRoughnessProperty(noahmp, FlagVegSfc)
    ! ground and surface roughness length and reference height
    call GroundRoughnessPropertyGlacier(noahmp)

    ! Thermal properties of soil, snow, lake, and frozen soil
    call GroundThermalProperty(noahmp)
    ! Thermal properties of snow and glacier ice
    call GroundThermalPropertyGlacier(noahmp)

    ! Surface shortwave albedo: ground and canopy radiative transfer
    call SurfaceAlbedo(noahmp)
    ! Glacier surface shortwave abeldo
    call GroundAlbedoGlacier(noahmp)

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
    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

       if ( (FlagVegSfc(I,J) .eqv. .true.) .and. (VegFrac(I,J) > 0) ) then ! vegetated portion of the grid
         TemperatureGrdVeg(I,J)  = TemperatureGrd(I,J)
         ExchCoeffMomAbvCan(I,J) = ExchCoeffMomSfc(I,J)
         ExchCoeffShAbvCan(I,J)  = ExchCoeffShSfc(I,J)
       endif

      end do
    end do

    call SurfaceEnergyFluxVegetated(noahmp)

    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! temperatures and energy fluxes of bare ground
    TemperatureGrdBare(I,J) = TemperatureGrd(I,J)
    ExchCoeffMomBare(I,J)   = ExchCoeffMomSfc(I,J)
    ExchCoeffShBare(I,J)    = ExchCoeffShSfc(I,J)

      end do
    end do

    call SurfaceEnergyFluxBareGround(noahmp)

    ! Grid-level computations requiring 2D parallel loop
    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! compute grid mean quantities by weighting vegetated and bare portions
    ! Energy balance at vege canopy: 
    ! RadSwAbsVeg = (RadLwNetCanopy + HeatSensibleCanopy + HeatLatentCanEvap + HeatLatentCanTransp) * VegFrac at VegFrac 
    ! Energy balance at vege ground: 
    ! RadSwAbsGrd * VegFrac = (RadLwNetVegGrd + HeatSensibleVegGrd + HeatLatentVegGrd + HeatGroundVegGrd) * VegFrac at VegFrac
    ! Energy balance at bare ground: 
    ! RadSwAbsGrd * (1-VegFrac) = (RadLwNetBareGrd + HeatSensibleBareGrd + HeatLatentBareGrd + HeatGroundBareGrd) * (1-VegFrac) at 1-VegFrac
    if ( (FlagVegSfc(I,J) .eqv. .true.) .and. (VegFrac(I,J) > 0) ) then
       WindStressEwSfc(I,J)     = VegFrac(I,J) * WindStressEwVeg(I,J)     + (1.0 - VegFrac(I,J)) * WindStressEwBare(I,J)
       WindStressNsSfc(I,J)     = VegFrac(I,J) * WindStressNsVeg(I,J)     + (1.0 - VegFrac(I,J)) * WindStressNsBare(I,J)
       RadLwNetSfc(I,J)         = VegFrac(I,J) * RadLwNetVegGrd(I,J)      + (1.0 - VegFrac(I,J)) * RadLwNetBareGrd(I,J) + RadLwNetCanopy(I,J)
       HeatSensibleSfc(I,J)     = VegFrac(I,J) * HeatSensibleVegGrd(I,J)  + (1.0 - VegFrac(I,J)) * HeatSensibleBareGrd(I,J) + HeatSensibleCanopy(I,J)
       HeatLatentGrd(I,J)       = VegFrac(I,J) * HeatLatentVegGrd(I,J)    + (1.0 - VegFrac(I,J)) * HeatLatentBareGrd(I,J)
       HeatGroundTot(I,J)  = VegFrac(I,J) * HeatGroundVegGrd(I,J)    + (1.0 - VegFrac(I,J)) * HeatGroundBareGrd(I,J)
       HeatLatentCanopy(I,J)    = HeatLatentCanEvap(I,J)
       HeatLatentTransp(I,J)    = HeatLatentCanTransp(I,J)
       HeatPrecipAdvSfc(I,J)    = VegFrac(I,J) * HeatPrecipAdvVegGrd(I,J) + (1.0 - VegFrac(I,J)) * HeatPrecipAdvBareGrd(I,J) + HeatPrecipAdvCanopy(I,J)
       TemperatureGrd(I,J)      = VegFrac(I,J) * TemperatureGrdVeg(I,J)   + (1.0 - VegFrac(I,J)) * TemperatureGrdBare(I,J)
       TemperatureAir2m(I,J)    = VegFrac(I,J) * TemperatureAir2mVeg(I,J) + (1.0 - VegFrac(I,J)) * TemperatureAir2mBare(I,J)
       TemperatureSfc(I,J)      = VegFrac(I,J) * TemperatureCanopy(I,J)   + (1.0 - VegFrac(I,J)) * TemperatureGrdBare(I,J)
       ExchCoeffMomSfc(I,J)     = VegFrac(I,J) * ExchCoeffMomAbvCan(I,J)  + (1.0 - VegFrac(I,J)) * ExchCoeffMomBare(I,J)     ! better way to average?
       ExchCoeffShSfc(I,J)      = VegFrac(I,J) * ExchCoeffShAbvCan(I,J)   + (1.0 - VegFrac(I,J)) * ExchCoeffShBare(I,J)
       SpecHumidity2m(I,J)      = VegFrac(I,J) * SpecHumidity2mVeg(I,J)   + (1.0 - VegFrac(I,J)) * SpecHumidity2mBare(I,J) 
       SpecHumiditySfcMean(I,J) = VegFrac(I,J) * (PressureVaporCanAir(I,J) * 0.622 / &
                             (PressureAirRefHeight(I,J) - 0.378*PressureVaporCanAir(I,J))) + (1.0 - VegFrac(I,J)) * SpecHumiditySfc(I,J)
       RoughLenMomSfcToAtm(I,J) = RoughLenMomSfc(I,J)
    else
       WindStressEwSfc(I,J)         = WindStressEwBare(I,J)
       WindStressNsSfc(I,J)         = WindStressNsBare(I,J)
       RadLwNetSfc(I,J)             = RadLwNetBareGrd(I,J)
       HeatSensibleSfc(I,J)         = HeatSensibleBareGrd(I,J)
       HeatLatentGrd(I,J)           = HeatLatentBareGrd(I,J)
       HeatGroundTot(I,J)      = HeatGroundBareGrd(I,J)
       TemperatureGrd(I,J)          = TemperatureGrdBare(I,J)
       TemperatureAir2m(I,J)        = TemperatureAir2mBare(I,J)
       HeatLatentCanopy(I,J)        = 0.0
       HeatLatentTransp(I,J)        = 0.0
       HeatPrecipAdvSfc(I,J)        = HeatPrecipAdvBareGrd(I,J)
       TemperatureSfc(I,J)          = TemperatureGrd(I,J)
       ExchCoeffMomSfc(I,J)         = ExchCoeffMomBare(I,J)
       ExchCoeffShSfc(I,J)          = ExchCoeffShBare(I,J)
       SpecHumiditySfcMean(I,J)     = SpecHumiditySfc(I,J)
       SpecHumidity2m(I,J)          = SpecHumidity2mBare(I,J)
       ResistanceStomataSunlit(I,J) = 0.0
       ResistanceStomataShade(I,J)  = 0.0
       TemperatureGrdVeg(I,J)       = TemperatureGrdBare(I,J)
       ExchCoeffShAbvCan(I,J)       = ExchCoeffShBare(I,J)
       RoughLenMomSfcToAtm(I,J)     = RoughLenMomGrd(I,J)
    endif

    ! emitted longwave radiation and physical check
    RadLwEmitSfc(I,J) = RadLwDownRefHeight(I,J) + RadLwNetSfc(I,J)
#ifndef _OPENACC
    if ( RadLwEmitSfc(I,J) <= 0.0 ) then
       write(*,*) "emitted longwave <0; skin T may be wrong due to inconsistent"
       write(*,*) "input of VegFracGreen with LeafAreaIndex"
       write(*,*) "VegFrac(I,J) = ", VegFrac(I,J), "VegAreaIndEff(I,J) = ", VegAreaIndEff(I,J), &
                  "TemperatureCanopy(I,J) = ", TemperatureCanopy(I,J), "TemperatureGrd(I,J) = ", TemperatureGrd(I,J)
       write(*,*) "RadLwDownRefHeight(I,J) = ", RadLwDownRefHeight(I,J), "RadLwNetSfc(I,J) = ", RadLwNetSfc(I,J), "SnowDepth(I,J) = ", SnowDepth(I,J)
       stop "Error: Longwave radiation budget problem in NoahMP LSM"
    endif
#endif

    ! radiative temperature: subtract from the emitted IR the
    ! reflected portion of the incoming longwave radiation, so just
    ! considering the IR originating/emitted in the canopy/ground system.
    ! Old TemperatureRadSfc calculation not taking into account Emissivity:
    ! TemperatureRadSfc = (RadLwEmitSfc/ConstStefanBoltzmann)**0.25
    TemperatureRadSfc(I,J) = ((RadLwEmitSfc(I,J) - (1.0-EmissivitySfc(I,J))*RadLwDownRefHeight(I,J)) / (EmissivitySfc(I,J)*ConstStefanBoltzmann))**0.25

    ! other photosynthesis related quantities for biochem process
    RadPhotoActAbsCan(I,J) = RadPhotoActAbsSunlit(I,J) * LeafAreaIndSunlit(I,J) + RadPhotoActAbsShade(I,J) * LeafAreaIndShade(I,J)
    PhotosynTotal(I,J) = PhotosynLeafSunlit(I,J)  * LeafAreaIndSunlit(I,J) + PhotosynLeafShade(I,J)   * LeafAreaIndShade(I,J)

    ! compute snow and soil layer temperature at soil timestep
    HeatFromSoilBot(I,J) = 0.0
    HeatGroundTotAcc(I,J) = HeatGroundTotAcc(I,J) + HeatGroundTot(I,J)
    if ( FlagSoilProcess .eqv. .true. ) HeatGroundTotMean(I,J) = HeatGroundTotAcc(I,J) / NumSoilTimeStep


      end do
    end do

    if ( noahmp%config%domain%FlagSoilProcess .eqv. .true. ) then
       call SoilSnowTemperatureMain(noahmp)
    endif ! FlagSoilProcess

    ! compute snow and glacier ice temperature (every main timestep)
    call GlacierTemperatureMain(noahmp)

    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! adjusting suface temperature based on snow condition
    if ( OptSnowSoilTempTime == 2 ) then
       if ( (SnowDepth(I,J) > 0.05) .and. (TemperatureGrd(I,J) > ConstFreezePoint) ) then
          TemperatureGrdVeg(I,J)  = ConstFreezePoint
          TemperatureGrdBare(I,J) = ConstFreezePoint
          if ( (FlagVegSfc(I,J) .eqv. .true.) .and. (VegFrac(I,J) > 0) ) then
             TemperatureGrd(I,J)  = VegFrac(I,J) * TemperatureGrdVeg(I,J) + (1.0 - VegFrac(I,J)) * TemperatureGrdBare(I,J)
             TemperatureSfc(I,J)  = VegFrac(I,J) * TemperatureCanopy(I,J) + (1.0 - VegFrac(I,J)) * TemperatureGrdBare(I,J)
          else
             if (IndicatorIceSfc(I,J) == -1) TemperatureGrdBare(I,J) = ConstFreezePoint
             TemperatureGrd(I,J)  = TemperatureGrdBare(I,J)
             TemperatureSfc(I,J)  = TemperatureGrdBare(I,J)
          endif
       endif
    endif


      end do
    end do
    !$acc end parallel loop

    ! Phase change and Energy released or consumed by snow & frozen soil
    call SoilSnowWaterPhaseChange(noahmp)
    ! Phase change and Energy released or consumed by snow & glacier ice
    call GlacierPhaseChange(noahmp)

    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! update sensible heat flux due to sprinkler irrigation evaporation
    if ( (FlagCropland(I,J) .eqv. .true.) .and. (IrrigationFracGrid(I,J) >= IrriFracThreshold(I,J)) ) &
       HeatSensibleSfc(I,J) = HeatSensibleSfc(I,J) - HeatLatentIrriEvap(I,J)

    ! update total surface albedo
    if ( RadSwDownRefHeight(I,J) > 0.0 ) then
       AlbedoSfc(I,J) = RadSwReflSfc(I,J) / RadSwDownRefHeight(I,J)
    else
       AlbedoSfc(I,J) = undefined_real
    endif


      end do
    end do
    !$acc end parallel loop

    end associate

    !$acc end data

    ! Deallocate local array
    deallocate(FlagVegSfc)



  end subroutine EnergyMain

end module EnergyMainMod
