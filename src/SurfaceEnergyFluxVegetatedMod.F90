module SurfaceEnergyFluxVegetatedMod

!!! Compute surface energy fluxes and budget for vegetated surface
!!! Use newton-raphson iteration to solve for vegetation and ground temperatures
!!! Surface energy balance:
!!! Canopy level: -RadSwAbsVeg - HeatPrecipAdvCanopy + RadLwNetCanopy + HeatSensibleCanopy + HeatLatentCanEvap + HeatLatentCanTransp + HeatCanStorageChg = 0
!!! Ground level: -RadSwAbsGrd - HeatPrecipAdvVegGrd + RadLwNetVegGrd + HeatSensibleVegGrd + HeatLatentVegGrd + HeatGroundVegGrd = 0

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use VaporPressureSaturationMod,          only : VaporPressureSaturation
  use ResistanceAboveCanopyMostMod,        only : ResistanceAboveCanopyMOST
  use ResistanceAboveCanopyChen97Mod,      only : ResistanceAboveCanopyChen97
  use ResistanceLeafToGroundMod,           only : ResistanceLeafToGround
  use ResistanceCanopyStomataBallBerryMod, only : ResistanceCanopyStomataBallBerry
  use ResistanceCanopyStomataJarvisMod,    only : ResistanceCanopyStomataJarvis

  implicit none

contains

  subroutine SurfaceEnergyFluxVegetated(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: VEGE_FLUX
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type)     , intent(inout) :: noahmp

! local variable
    integer                               :: I, J                    ! grid indices
    integer                               :: IndIter                 ! iteration index
    integer                               :: LastIter                ! Last iteration
    integer, allocatable, dimension(:,:)  :: MoStabParaSgn           ! number of times MoStabParaAbvCan changes sign
    integer                               :: IndexShade              ! index for sunlit/shaded (0=sunlit;1=shaded)
    integer, parameter                    :: NumIterC = 20           ! number of iterations for surface temperature (5~20)
    integer, parameter                    :: NumIterG = 5            ! number of iterations for ground temperature (3~5)
    real(kind=kind_noahmp)                :: ExchCoeffShAbvCanTmp    ! sensible heat conductance, canopy air to reference height air [m/s]
    real(kind=kind_noahmp)                :: TemperatureCanChg       ! change in tv, last iteration [K]
    real(kind=kind_noahmp)                :: TemperatureGrdChg       ! change in tg, last iteration [K]
    real(kind=kind_noahmp)                :: LwCoeffAir              ! coefficients for longwave emission as function of ts**4
    real(kind=kind_noahmp)                :: LwCoeffCan              ! coefficients for longwave emission as function of ts**4
    real(kind=kind_noahmp)                :: ShCoeff                 ! coefficients for sensible heat as function of ts
    real(kind=kind_noahmp)                :: LhCoeff                 ! coefficients for latent heat as function of ts
    real(kind=kind_noahmp)                :: GrdHeatCoeff            ! coefficients for ground heat as function of ts
    real(kind=kind_noahmp)                :: TranspHeatCoeff         ! coefficients for transpiration heat as function of ts
    real(kind=kind_noahmp)                :: TempShGhTmp             ! partial temperature by sensible and ground heat
    real(kind=kind_noahmp)                :: ExchCoeffShFrac         ! exchange coefficient fraction for sensible heat 
    real(kind=kind_noahmp)                :: VapPresLhTot            ! vapor pressure related to total latent heat
    real(kind=kind_noahmp)                :: ExchCoeffEtFrac         ! exchange coefficient fraction for evapotranspiration heat
    real(kind=kind_noahmp)                :: VapPresSatWatTmp        ! saturated vapor pressure for water
    real(kind=kind_noahmp)                :: VapPresSatIceTmp        ! saturated vapor pressure for ice
    real(kind=kind_noahmp)                :: VapPresSatWatTmpD       ! saturated vapor pressure gradient with ground temp. [Pa/K] for water
    real(kind=kind_noahmp)                :: VapPresSatIceTmpD       ! saturated vapor pressure gradient with ground temp. [Pa/K] for ice
    real(kind=kind_noahmp)                :: FluxTotCoeff            ! temporary total coefficients for all energy flux
    real(kind=kind_noahmp)                :: EnergyResTmp            ! temporary energy residual
    real(kind=kind_noahmp)                :: ExchCoeffShLeafTmp      ! sensible heat conductance, leaf surface to canopy air [m/s]
    real(kind=kind_noahmp)                :: ExchCoeffTot            ! sum of conductances [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ShCanTmp                ! temporary sensible heat flux [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ShGrdTmp                ! temporary sensible heat flux [W/m2]
    real(kind=kind_noahmp)                :: MoistureFluxSfc         ! moisture flux
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VegAreaIndTmp           ! total leaf area index + stem area index,effective
    real(kind=kind_noahmp)                :: LeafAreaIndSunEff       ! sunlit leaf area index, one-sided [m2/m2],effective
    real(kind=kind_noahmp)                :: LeafAreaIndShdEff       ! shaded leaf area index, one-sided [m2/m2],effective
    real(kind=kind_noahmp)                :: TempTmp                 ! temporary temperature
    real(kind=kind_noahmp)                :: TempUnitConv            ! Kelvin to degree Celsius with limit -50 to +50
    real(kind=kind_noahmp)                :: HeatCapacCan            ! canopy heat capacity [J/m2/K]
! local statement function
    TempUnitConv(TempTmp) = min(50.0, max(-50.0, (TempTmp - ConstFreezePoint)))

    allocate(VegAreaIndTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                             noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MoStabParaSgn(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                              noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(ShCanTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                      noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(ShGrdTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                      noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    !$acc data create(ShGrdTmp, ShCanTmp, MoStabParaSgn, VegAreaIndTmp)
! --------------------------------------------------------------------

    ! Initialization parallel region
    !$acc parallel loop collapse(2) gang vector present(noahmp, ShGrdTmp, ShCanTmp, VegAreaIndTmp, MoStabParaSgn) &
    !$acc private(LastIter, TemperatureCanChg, TemperatureGrdChg) &
    !$acc private(VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD, TempTmp) &
    !$acc private(LwCoeffAir, LwCoeffCan, MoistureFluxSfc)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
       if ( .not. ((noahmp%energy%state%VegAreaIndEff(I,J) > 0.0 ) .and. (noahmp%energy%state%VegFrac(I,J) > 0)) ) cycle ! skip non-vegetated surface

        associate(                                                                        &
                  SnowDepth               => noahmp%water%state%SnowDepth(I,J)           ,& ! in,    snow depth [m]
                  HeightCanopyTop         => noahmp%energy%param%HeightCanopyTop(I,J)    ,& ! in,    top of canopy [m]
                  RadLwDownRefHeight      => noahmp%forcing%RadLwDownRefHeight(I,J)      ,& ! in,    downward longwave radiation [W/m2] at reference height
                  VegAreaIndEff           => noahmp%energy%state%VegAreaIndEff(I,J)      ,& ! in,    one-sided leaf+stem area index [m2/m2]
                  LeafAreaIndSunlit       => noahmp%energy%state%LeafAreaIndSunlit(I,J)  ,& ! in,    sunlit leaf area index, one-sided [m2/m2]
                  LeafAreaIndShade        => noahmp%energy%state%LeafAreaIndShade(I,J)   ,& ! in,    shaded leaf area index, one-sided [m2/m2]
                  ZeroPlaneDispSfc        => noahmp%energy%state%ZeroPlaneDispSfc(I,J)   ,& ! in,    zero plane displacement [m]
                  RoughLenMomSfc          => noahmp%energy%state%RoughLenMomSfc(I,J)     ,& ! in,    roughness length [m], momentum, surface
                  EmissivityVeg           => noahmp%energy%state%EmissivityVeg(I,J)      ,& ! in,    vegetation emissivity
                  EmissivityGrd           => noahmp%energy%state%EmissivityGrd(I,J)      ,& ! in,    ground emissivity
                  RefHeightAboveGrd       => noahmp%energy%state%RefHeightAboveGrd(I,J)  ,& ! in,    surface reference height [m]
                  WindSpdRefHeight        => noahmp%energy%state%WindSpdRefHeight(I,J)   ,& ! in,    wind speed [m/s] at reference height
                  TemperatureGrdVeg       => noahmp%energy%state%TemperatureGrdVeg(I,J)  ,& ! inout, vegetated ground (below-canopy) temperature [K]
                  FrictionVelVeg          => noahmp%energy%state%FrictionVelVeg(I,J)     ,& ! out,   friction velocity [m/s], vegetated
                  MoStabParaAbvCan        => noahmp%energy%state%MoStabParaAbvCan(I,J)   ,& ! out,   Monin-Obukhov stability (z/L), above ZeroPlaneDispSfc, vegetated
                  MoStabCorrShVeg2m       => noahmp%energy%state%MoStabCorrShVeg2m(I,J)  ,& ! out,   M-O sen heat stability correction, 2m, vegetated
                  VapPresSatGrdVeg        => noahmp%energy%state%VapPresSatGrdVeg(I,J)   ,& ! out,   saturation vapor pressure at TemperatureGrd [Pa]
                  CanopyHeight            => noahmp%energy%state%CanopyHeight(I,J)       ,& ! out,   canopy height [m]
                  WindSpdCanopyTop        => noahmp%energy%state%WindSpdCanopyTop(I,J)    & ! out,   wind speed at top of canopy [m/s]
                 )

        ! initialization (including variables that do not depend on stability iteration)
        LastIter          = 0
        FrictionVelVeg    = 0.1
        TemperatureCanChg = 0.0
        TemperatureGrdChg = 0.0
        MoStabParaAbvCan  = 0.0
        MoStabParaSgn(I,J)     = 0
        MoStabCorrShVeg2m = 0.0
        ShGrdTmp(I,J)     = 0.0
        ShCanTmp(I,J)     = 0.0
        MoistureFluxSfc   = 0.0
        ! limit LeafAreaIndex
        VegAreaIndTmp(I,J)     = min(6.0, VegAreaIndEff)

        ! saturation vapor pressure at ground temperature
        TempTmp = TempUnitConv(TemperatureGrdVeg)
        call VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD)
        if ( TempTmp > 0.0 ) then
           VapPresSatGrdVeg = VapPresSatWatTmp
        else
           VapPresSatGrdVeg = VapPresSatIceTmp
        endif

        ! canopy height
        CanopyHeight = HeightCanopyTop
        ! wind speed at canopy height
       !WindSpdCanopyTop = WindSpdRefHeight * log(CanopyHeight/RoughLenMomSfc) / log(RefHeightAboveGrd/RoughLenMomSfc)
        WindSpdCanopyTop = WindSpdRefHeight * log((CanopyHeight - ZeroPlaneDispSfc + RoughLenMomSfc)/RoughLenMomSfc) / &
                           log(RefHeightAboveGrd/RoughLenMomSfc)                                           ! MB: add ZeroPlaneDispSfc v3.7
#ifndef _OPENACC
        if ( (CanopyHeight-ZeroPlaneDispSfc) <= 0.0 ) then
           print*, "CRITICAL PROBLEM: CanopyHeight <= ZeroPlaneDispSfc"
           print*, "CanopyHeight = "         , CanopyHeight
           print*, "ZeroPlaneDispSfc = "     , ZeroPlaneDispSfc
           print*, "SnowDepth = "            , SnowDepth
           stop "Error: ZeroPlaneDisp problem in NoahMP LSM"
        endif
#endif

        end associate
      end do
    end do
    !$acc end parallel loop

    ! begin stability iteration for canopy temperature and flux
    loop1: do IndIter = 1, NumIterC

       ! Roughness length calculation
       !$acc parallel loop collapse(2) gang vector present(noahmp, ShGrdTmp, ShCanTmp, VegAreaIndTmp, MoStabParaSgn)
       do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
            
         if ( .not. ((noahmp%energy%state%VegAreaIndEff(I,J) > 0.0 ) .and. (noahmp%energy%state%VegFrac(I,J) > 0)) ) cycle ! skip non-vegetated surface

           associate(                                                                        &
                     RoughLenMomSfc          => noahmp%energy%state%RoughLenMomSfc(I,J)     ,& ! in,    roughness length [m], momentum, surface
                     RoughLenMomGrd          => noahmp%energy%state%RoughLenMomGrd(I,J)     ,& ! in,    roughness length [m], momentum, ground
                     RoughLenShCanopy        => noahmp%energy%state%RoughLenShCanopy(I,J)  ,& ! out,   roughness length [m], sensible heat, vegetated
                     RoughLenShVegGrd        => noahmp%energy%state%RoughLenShVegGrd(I,J)   & ! out,   roughness length [m], sensible heat ground, below canopy
                    )

           ! ground and surface roughness length
           if ( IndIter == 1 ) then
              RoughLenShCanopy = RoughLenMomSfc
              RoughLenShVegGrd = RoughLenMomGrd
           else
              RoughLenShCanopy = RoughLenMomSfc  !* exp(-ZilitinkevichCoeff*0.4*258.2*sqrt(FrictionVelVeg*RoughLenMomSfc))
              RoughLenShVegGrd = RoughLenMomGrd  !* exp(-ZilitinkevichCoeff*0.4*258.2*sqrt(FrictionVelVeg*RoughLenMomGrd))
           endif

           end associate
         end do
       end do
       !$acc end parallel loop

       ! aerodyn resistances between RefHeightAboveGrd and d+z0v
       if ( noahmp%config%nmlist%OptSurfaceDrag == 1 ) call ResistanceAboveCanopyMOST(noahmp, IndIter, ShCanTmp, MoStabParaSgn)
       if ( noahmp%config%nmlist%OptSurfaceDrag == 2 ) call ResistanceAboveCanopyChen97(noahmp, IndIter)

       ! aerodyn resistance between z0g and d+z0v, and leaf boundary layer resistance
       call ResistanceLeafToGround(noahmp, IndIter, VegAreaIndTmp, ShGrdTmp)

       ! stomatal resistance (only on first iteration)
       if ( IndIter == 1 ) then
          if ( noahmp%config%nmlist%OptStomataResistance == 1 ) then  ! Ball-Berry
             IndexShade = 0 ! sunlit case
             call ResistanceCanopyStomataBallBerry(noahmp, IndexShade)
             IndexShade = 1 ! shaded case
             call ResistanceCanopyStomataBallBerry(noahmp, IndexShade)
          endif
          if ( noahmp%config%nmlist%OptStomataResistance == 2 ) then  ! Jarvis
             IndexShade = 0 ! sunlit case
             call ResistanceCanopyStomataJarvis(noahmp, IndexShade)
             IndexShade = 1 ! shaded case
             call ResistanceCanopyStomataJarvis(noahmp, IndexShade)
          endif
       endif

       ! Canopy flux calculations
       !$acc parallel loop collapse(2) gang vector present(noahmp, VegAreaIndTmp) &
       !$acc private(LastIter, TemperatureCanChg, TemperatureGrdChg) &
       !$acc private(LeafAreaIndSunEff, LeafAreaIndShdEff) &
       !$acc private(VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD, TempTmp) &
       !$acc private(LwCoeffAir, LwCoeffCan, ShCoeff, LhCoeff, GrdHeatCoeff, TranspHeatCoeff) &
       !$acc private(ExchCoeffShAbvCanTmp, ExchCoeffShLeafTmp, ExchCoeffTot, TempShGhTmp) &
       !$acc private(ExchCoeffShFrac, VapPresLhTot, ExchCoeffEtFrac, FluxTotCoeff, EnergyResTmp) &
       !$acc private(MoistureFluxSfc, HeatCapacCan)
       do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if ( .not. ((noahmp%energy%state%VegAreaIndEff(I,J) > 0.0 ) .and. (noahmp%energy%state%VegFrac(I,J) > 0)) ) cycle ! skip non-vegetated surface

           associate(                                                                        &
                     MainTimeStep            => noahmp%config%domain%MainTimeStep           ,& ! in,    main noahmp timestep [s]
                     RadLwDownRefHeight      => noahmp%forcing%RadLwDownRefHeight(I,J)      ,& ! in,    downward longwave radiation [W/m2] at reference height
                     TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight(I,J) ,& ! in,    air temperature [K] at reference height
                     PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight(I,J)    ,& ! in,    air pressure [Pa] at reference height
                     CanopyWetFrac           => noahmp%water%state%CanopyWetFrac(I,J)       ,& ! in,    wetted or snowed fraction of the canopy
                     CanopyLiqWater          => noahmp%water%state%CanopyLiqWater(I,J)      ,& ! in,    canopy intercepted liquid water [mm]
                     CanopyIce               => noahmp%water%state%CanopyIce(I,J)           ,& ! in,    canopy intercepted ice [mm]
                     HeatCapacCanFac         => noahmp%energy%param%HeatCapacCanFac(I,J)    ,& ! in,    canopy biomass heat capacity parameter [m]
                     RadSwAbsVeg             => noahmp%energy%flux%RadSwAbsVeg(I,J)         ,& ! in,    solar radiation absorbed by vegetation [W/m2]
                     HeatPrecipAdvCanopy     => noahmp%energy%flux%HeatPrecipAdvCanopy(I,J) ,& ! in,    precipitation advected heat - vegetation net [W/m2]
                     VegFrac                 => noahmp%energy%state%VegFrac(I,J)            ,& ! in,    greeness vegetation fraction
                     PressureVaporRefHeight  => noahmp%energy%state%PressureVaporRefHeight(I,J),& ! in,    vapor pressure air [Pa] at reference height
                     DensityAirRefHeight     => noahmp%energy%state%DensityAirRefHeight(I,J),& ! in,    density air [kg/m3]
                     VegAreaIndEff           => noahmp%energy%state%VegAreaIndEff(I,J)      ,& ! in,    one-sided leaf+stem area index [m2/m2]
                     LeafAreaIndSunlit       => noahmp%energy%state%LeafAreaIndSunlit(I,J)  ,& ! in,    sunlit leaf area index, one-sided [m2/m2]
                     LeafAreaIndShade        => noahmp%energy%state%LeafAreaIndShade(I,J)   ,& ! in,    shaded leaf area index, one-sided [m2/m2]
                     EmissivityVeg           => noahmp%energy%state%EmissivityVeg(I,J)      ,& ! in,    vegetation emissivity
                     EmissivityGrd           => noahmp%energy%state%EmissivityGrd(I,J)      ,& ! in,    ground emissivity
                     ResistanceGrdEvap       => noahmp%energy%state%ResistanceGrdEvap(I,J)  ,& ! in,    ground surface resistance [s/m] to evaporation
                     PsychConstCanopy        => noahmp%energy%state%PsychConstCanopy(I,J)   ,& ! in,    psychrometric constant [Pa/K], canopy
                     LatHeatVapCanopy        => noahmp%energy%state%LatHeatVapCanopy(I,J)   ,& ! in,    latent heat of vaporization/subli [J/kg], canopy
                     ResistanceStomataSunlit => noahmp%energy%state%ResistanceStomataSunlit(I,J),& ! in,  sunlit leaf stomatal resistance [s/m]
                     ResistanceStomataShade  => noahmp%energy%state%ResistanceStomataShade(I,J),& ! in,   shaded leaf stomatal resistance [s/m]
                     ResistanceLeafBoundary  => noahmp%energy%state%ResistanceLeafBoundary(I,J),& ! in,   bulk leaf boundary layer resistance [s/m]
                     ResistanceShAbvCan      => noahmp%energy%state%ResistanceShAbvCan(I,J) ,& ! in,    aerodynamic resistance for sensible heat [s/m], above canopy
                     ResistanceLhAbvCan      => noahmp%energy%state%ResistanceLhAbvCan(I,J) ,& ! in,    aerodynamic resistance for water vapor [s/m], above canopy
                     ResistanceShUndCan      => noahmp%energy%state%ResistanceShUndCan(I,J) ,& ! in,    ground aerodynamic resistance for sensible heat [s/m]
                     ResistanceLhUndCan      => noahmp%energy%state%ResistanceLhUndCan(I,J) ,& ! in,    ground aerodynamic resistance for water vapor [s/m]
                     TemperatureGrdVeg       => noahmp%energy%state%TemperatureGrdVeg(I,J)  ,& ! in,    vegetated ground (below-canopy) temperature [K]
                     VapPresSatGrdVeg        => noahmp%energy%state%VapPresSatGrdVeg(I,J)   ,& ! in,    saturation vapor pressure at TemperatureGrd [Pa]
                     SpecHumiditySfc         => noahmp%energy%state%SpecHumiditySfc(I,J)    ,& ! inout, specific humidity at vegetated surface
                     PressureVaporCanAir     => noahmp%energy%state%PressureVaporCanAir(I,J),& ! inout, canopy air vapor pressure [Pa]
                     TemperatureCanopyAir    => noahmp%energy%state%TemperatureCanopyAir(I,J),& ! inout, canopy air temperature [K]
                     TemperatureCanopy       => noahmp%energy%state%TemperatureCanopy(I,J)  ,& ! inout, vegetation temperature [K]
                     ExchCoeffLhAbvCan       => noahmp%energy%state%ExchCoeffLhAbvCan(I,J)  ,& ! out,   latent heat conductance, canopy air to reference height [m/s]
                     ExchCoeffLhTransp       => noahmp%energy%state%ExchCoeffLhTransp(I,J)  ,& ! out,   transpiration conductance, leaf to canopy air [m/s]
                     ExchCoeffLhEvap         => noahmp%energy%state%ExchCoeffLhEvap(I,J)    ,& ! out,   evaporation conductance, leaf to canopy air [m/s]
                     ExchCoeffLhUndCan       => noahmp%energy%state%ExchCoeffLhUndCan(I,J)  ,& ! out,   latent heat conductance, ground to canopy air [m/s]
                     VapPresSatCanopy        => noahmp%energy%state%VapPresSatCanopy(I,J)   ,& ! out,   saturation vapor pressure at TemperatureCanopy [Pa]
                     VapPresSatCanTempD      => noahmp%energy%state%VapPresSatCanTempD(I,J) ,& ! out,   d(VapPresSatCanopy)/dt at TemperatureCanopy [Pa/K]
                     RadLwNetCanopy          => noahmp%energy%flux%RadLwNetCanopy(I,J)      ,& ! out,   canopy net longwave radiation [W/m2] (+ to atm)
                     HeatSensibleCanopy      => noahmp%energy%flux%HeatSensibleCanopy(I,J)  ,& ! out,   canopy sensible heat flux [W/m2] (+ to atm)
                     HeatLatentCanEvap       => noahmp%energy%flux%HeatLatentCanEvap(I,J)   ,& ! out,   canopy evaporation heat flux [W/m2] (+ to atm)
                     HeatLatentCanTransp     => noahmp%energy%flux%HeatLatentCanTransp(I,J) ,& ! out,   canopy transpiration heat flux [W/m2] (+ to atm)
                     HeatCanStorageChg       => noahmp%energy%flux%HeatCanStorageChg(I,J)    & ! out,   canopy heat storage change [W/m2]
                    )

           ! limit LeafAreaIndex
           LeafAreaIndSunEff = min(6.0, LeafAreaIndSunlit)
           LeafAreaIndShdEff = min(6.0, LeafAreaIndShade)

           ! prepare for longwave rad.
           LwCoeffAir = -EmissivityVeg * (1.0 + (1.0-EmissivityVeg)*(1.0-EmissivityGrd)) * RadLwDownRefHeight - &
                         EmissivityVeg * EmissivityGrd * ConstStefanBoltzmann * TemperatureGrdVeg**4
           LwCoeffCan = (2.0 - EmissivityVeg * (1.0-EmissivityGrd)) * EmissivityVeg * ConstStefanBoltzmann

           ! ES and d(ES)/dt evaluated at TemperatureCanopy
           TempTmp = TempUnitConv(TemperatureCanopy)
           call VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD)
           if ( TempTmp > 0.0 ) then
              VapPresSatCanopy   = VapPresSatWatTmp
              VapPresSatCanTempD = VapPresSatWatTmpD
           else
              VapPresSatCanopy   = VapPresSatIceTmp
              VapPresSatCanTempD = VapPresSatIceTmpD
           endif

           ! sensible heat conductance and coeff above veg.
           ExchCoeffShAbvCanTmp = 1.0 / ResistanceShAbvCan
           ExchCoeffShLeafTmp   = 2.0 * VegAreaIndTmp(I,J) / ResistanceLeafBoundary
           GrdHeatCoeff         = 1.0 / ResistanceShUndCan
           ExchCoeffTot         = ExchCoeffShAbvCanTmp + ExchCoeffShLeafTmp + GrdHeatCoeff
           TempShGhTmp          = (TemperatureAirRefHeight*ExchCoeffShAbvCanTmp + TemperatureGrdVeg*GrdHeatCoeff) / ExchCoeffTot
           ExchCoeffShFrac      = ExchCoeffShLeafTmp / ExchCoeffTot
           ShCoeff              = (1.0 - ExchCoeffShFrac) * DensityAirRefHeight * ConstHeatCapacAir * ExchCoeffShLeafTmp

           ! latent heat conductance and coeff above veg.
           ExchCoeffLhAbvCan = 1.0 / ResistanceLhAbvCan
           ExchCoeffLhEvap   = CanopyWetFrac * VegAreaIndTmp(I,J) / ResistanceLeafBoundary
           ExchCoeffLhTransp = (1.0 - CanopyWetFrac) * (LeafAreaIndSunEff/(ResistanceLeafBoundary+ResistanceStomataSunlit) + &
                                                        LeafAreaIndShdEff/(ResistanceLeafBoundary+ResistanceStomataShade))
           ExchCoeffLhUndCan = 1.0 / (ResistanceLhUndCan + ResistanceGrdEvap)
           ExchCoeffTot      = ExchCoeffLhAbvCan + ExchCoeffLhEvap + ExchCoeffLhTransp + ExchCoeffLhUndCan
           VapPresLhTot      = (PressureVaporRefHeight*ExchCoeffLhAbvCan + VapPresSatGrdVeg*ExchCoeffLhUndCan ) / ExchCoeffTot
           ExchCoeffEtFrac   = (ExchCoeffLhEvap + ExchCoeffLhTransp) / ExchCoeffTot
           LhCoeff           = (1.0 - ExchCoeffEtFrac) * ExchCoeffLhEvap * DensityAirRefHeight * &
                               ConstHeatCapacAir / PsychConstCanopy
           TranspHeatCoeff   = (1.0 - ExchCoeffEtFrac) * ExchCoeffLhTransp * DensityAirRefHeight * &
                               ConstHeatCapacAir / PsychConstCanopy

           ! evaluate surface fluxes with current temperature and solve for temperature change
           TemperatureCanopyAir = TempShGhTmp + ExchCoeffShFrac * TemperatureCanopy
           PressureVaporCanAir  = VapPresLhTot + ExchCoeffEtFrac * VapPresSatCanopy
           RadLwNetCanopy       = VegFrac * (LwCoeffAir + LwCoeffCan * TemperatureCanopy**4)
           HeatSensibleCanopy   = VegFrac * DensityAirRefHeight * ConstHeatCapacAir * &
                                  ExchCoeffShLeafTmp * (TemperatureCanopy - TemperatureCanopyAir)
           HeatLatentCanEvap    = VegFrac * DensityAirRefHeight * ConstHeatCapacAir * ExchCoeffLhEvap * &
                                  (VapPresSatCanopy - PressureVaporCanAir) / PsychConstCanopy
           HeatLatentCanTransp  = VegFrac * DensityAirRefHeight * ConstHeatCapacAir * ExchCoeffLhTransp * &
                                  (VapPresSatCanopy - PressureVaporCanAir) / PsychConstCanopy
           if ( TemperatureCanopy > ConstFreezePoint ) then
              HeatLatentCanEvap = min(CanopyLiqWater*LatHeatVapCanopy/MainTimeStep, HeatLatentCanEvap)
           else
              HeatLatentCanEvap = min(CanopyIce*LatHeatVapCanopy/MainTimeStep, HeatLatentCanEvap)
           endif
           ! canopy heat capacity
           HeatCapacCan         = HeatCapacCanFac*VegAreaIndTmp(I,J)*ConstHeatCapacWater + CanopyLiqWater*ConstHeatCapacWater/ConstDensityWater + &
                                  CanopyIce*ConstHeatCapacIce/ConstDensityIce
           ! compute vegetation temperature change
           EnergyResTmp         = RadSwAbsVeg - RadLwNetCanopy - HeatSensibleCanopy - &
                                  HeatLatentCanEvap - HeatLatentCanTransp + HeatPrecipAdvCanopy
           FluxTotCoeff         = VegFrac * (4.0*LwCoeffCan*TemperatureCanopy**3 + ShCoeff + &
                                            (LhCoeff+TranspHeatCoeff)*VapPresSatCanTempD + HeatCapacCan/MainTimeStep)
           TemperatureCanChg    = EnergyResTmp / FluxTotCoeff
           ! update fluxes with temperature change
           RadLwNetCanopy       = RadLwNetCanopy      + VegFrac * 4.0 * LwCoeffCan * TemperatureCanopy**3 * TemperatureCanChg
           HeatSensibleCanopy   = HeatSensibleCanopy  + VegFrac * ShCoeff * TemperatureCanChg
           HeatLatentCanEvap    = HeatLatentCanEvap   + VegFrac * LhCoeff * VapPresSatCanTempD * TemperatureCanChg
           HeatLatentCanTransp  = HeatLatentCanTransp + VegFrac * TranspHeatCoeff * VapPresSatCanTempD * TemperatureCanChg
           HeatCanStorageChg    = VegFrac * HeatCapacCan / MainTimeStep * TemperatureCanChg
           ! update vegetation temperature
           TemperatureCanopy    = TemperatureCanopy + TemperatureCanChg
          !TemperatureCanopyAir = TempShGhTmp + ExchCoeffShFrac * TemperatureCanopy                        ! canopy air T; update here for consistency

            ! for computing M-O length in the next iteration
            ShCanTmp(I,J) = DensityAirRefHeight * ConstHeatCapacAir * (TemperatureCanopyAir-TemperatureAirRefHeight) / ResistanceShAbvCan
            ShGrdTmp(I,J) = DensityAirRefHeight * ConstHeatCapacAir * (TemperatureGrdVeg-TemperatureCanopyAir) / ResistanceShUndCan

           ! consistent specific humidity from canopy air vapor pressure
           SpecHumiditySfc = (0.622 * PressureVaporCanAir) / (PressureAirRefHeight - 0.378 * PressureVaporCanAir)
         !   if ( LastIter == 1 ) then
         !      exit loop1
         !   endif
         !   if ( (IndIter >= 5) .and. (abs(TemperatureCanChg) <= 0.01) .and. (LastIter == 0) ) then
         !      LastIter = 1
         !   endif

           end associate
         end do
       end do
       !$acc end parallel loop

    enddo loop1  ! end stability iteration

    ! Ground temperature iteration (loop2) and final calculations
    !$acc parallel loop collapse(2) gang vector present(noahmp, VegAreaIndTmp) &
    !$acc private(LastIter, TemperatureCanChg, TemperatureGrdChg, IndIter) &
    !$acc private(VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD, TempTmp) &
    !$acc private(LwCoeffAir, LwCoeffCan, ShCoeff, LhCoeff, GrdHeatCoeff) &
    !$acc private(FluxTotCoeff, EnergyResTmp, ExchCoeffShAbvCanTmp, ExchCoeffShLeafTmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
        if ( .not. ((noahmp%energy%state%VegAreaIndEff(I,J) > 0.0 ) .and. (noahmp%energy%state%VegFrac(I,J) > 0)) ) cycle ! skip non-vegetated surface

        associate(                                                                        &
                  NumSnowLayerNeg         => noahmp%config%domain%NumSnowLayerNeg(I,J)   ,& ! in,    actual number of snow layers (negative)
                  ThicknessSnowSoilLayer  => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
                  OptSurfaceDrag          => noahmp%config%nmlist%OptSurfaceDrag         ,& ! in,    options for surface layer drag/exchange coefficient
                  OptSnowSoilTempTime     => noahmp%config%nmlist%OptSnowSoilTempTime    ,& ! in,    options for snow/soil temperature time scheme (only layer 1)
                  WindEastwardRefHeight   => noahmp%forcing%WindEastwardRefHeight(I,J)   ,& ! in,    wind speed [m/s] in eastward direction at reference height
                  WindNorthwardRefHeight  => noahmp%forcing%WindNorthwardRefHeight(I,J)  ,& ! in,    wind speed [m/s] in northward direction at reference height
                  RadLwDownRefHeight      => noahmp%forcing%RadLwDownRefHeight(I,J)      ,& ! in,    downward longwave radiation [W/m2] at reference height
                  TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight(I,J) ,& ! in,    air temperature [K] at reference height
                  PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight(I,J)    ,& ! in,    air pressure [Pa] at reference height
                  SnowDepth               => noahmp%water%state%SnowDepth(I,J)           ,& ! in,    snow depth [m]
                  SnowCoverFrac           => noahmp%water%state%SnowCoverFrac(I,J)       ,& ! in,    snow cover fraction
                  RadSwAbsGrd             => noahmp%energy%flux%RadSwAbsGrd(I,J)         ,& ! in,    solar radiation absorbed by ground [W/m2]
                  HeatPrecipAdvVegGrd     => noahmp%energy%flux%HeatPrecipAdvVegGrd(I,J) ,& ! in,    precipitation advected heat - under canopy net [W/m2]
                  VegFrac                 => noahmp%energy%state%VegFrac(I,J)            ,& ! in,    greeness vegetation fraction
                  WindSpdRefHeight        => noahmp%energy%state%WindSpdRefHeight(I,J)   ,& ! in,    wind speed [m/s] at reference height
                  DensityAirRefHeight     => noahmp%energy%state%DensityAirRefHeight(I,J),& ! in,    density air [kg/m3]
                  VegAreaIndEff           => noahmp%energy%state%VegAreaIndEff(I,J)      ,& ! in,    one-sided leaf+stem area index [m2/m2]
                  EmissivityVeg           => noahmp%energy%state%EmissivityVeg(I,J)      ,& ! in,    vegetation emissivity
                  EmissivityGrd           => noahmp%energy%state%EmissivityGrd(I,J)      ,& ! in,    ground emissivity
                  TemperatureSoilSnow     => noahmp%energy%state%TemperatureSoilSnow     ,& ! in,    snow and soil layer temperature [K]
                  ThermConductSoilSnow    => noahmp%energy%state%ThermConductSoilSnow    ,& ! in,    thermal conductivity [W/m/K] for all soil & snow
                  ResistanceGrdEvap       => noahmp%energy%state%ResistanceGrdEvap(I,J)  ,& ! in,    ground surface resistance [s/m] to evaporation
                  LatHeatVapCanopy        => noahmp%energy%state%LatHeatVapCanopy(I,J)   ,& ! in,    latent heat of vaporization/subli [J/kg], canopy
                  PsychConstGrd           => noahmp%energy%state%PsychConstGrd(I,J)      ,& ! in,    psychrometric constant [Pa/K], ground
                  RelHumidityGrd          => noahmp%energy%state%RelHumidityGrd(I,J)     ,& ! in,    raltive humidity in surface soil/snow air space
                  RoughLenShCanopy        => noahmp%energy%state%RoughLenShCanopy(I,J)   ,& ! in,    roughness length [m], sensible heat, vegetated
                  ResistanceLeafBoundary  => noahmp%energy%state%ResistanceLeafBoundary(I,J),& ! in,   bulk leaf boundary layer resistance [s/m]
                  ResistanceShAbvCan      => noahmp%energy%state%ResistanceShAbvCan(I,J) ,& ! in,    aerodynamic resistance for sensible heat [s/m], above canopy
                  ResistanceShUndCan      => noahmp%energy%state%ResistanceShUndCan(I,J) ,& ! in,    ground aerodynamic resistance for sensible heat [s/m]
                  ResistanceLhUndCan      => noahmp%energy%state%ResistanceLhUndCan(I,J) ,& ! in,    ground aerodynamic resistance for water vapor [s/m]
                  FrictionVelVeg          => noahmp%energy%state%FrictionVelVeg(I,J)     ,& ! in,    friction velocity [m/s], vegetated
                  MoStabCorrShVeg2m       => noahmp%energy%state%MoStabCorrShVeg2m(I,J)  ,& ! in,    M-O sen heat stability correction, 2m, vegetated
                  SpecHumiditySfc         => noahmp%energy%state%SpecHumiditySfc(I,J)    ,& ! inout, specific humidity at vegetated surface
                  PressureVaporCanAir     => noahmp%energy%state%PressureVaporCanAir(I,J),& ! inout, canopy air vapor pressure [Pa]
                  TemperatureCanopyAir    => noahmp%energy%state%TemperatureCanopyAir(I,J),& ! inout, canopy air temperature [K]
                  TemperatureCanopy       => noahmp%energy%state%TemperatureCanopy(I,J)  ,& ! inout, vegetation temperature [K]
                  TemperatureGrdVeg       => noahmp%energy%state%TemperatureGrdVeg(I,J)  ,& ! inout, vegetated ground (below-canopy) temperature [K]
                  ExchCoeffMomAbvCan      => noahmp%energy%state%ExchCoeffMomAbvCan(I,J) ,& ! inout, momentum exchange coeff [m/s], above ZeroPlaneDisp, vegetated
                  ExchCoeffShAbvCan       => noahmp%energy%state%ExchCoeffShAbvCan(I,J)  ,& ! inout, heat exchange coeff [m/s], above ZeroPlaneDisp, vegetated
                  VapPresSatGrdVeg        => noahmp%energy%state%VapPresSatGrdVeg(I,J)   ,& ! inout, saturation vapor pressure at TemperatureGrd [Pa]
                  VapPresSatGrdVegTempD   => noahmp%energy%state%VapPresSatGrdVegTempD(I,J),& ! out,   d(VapPresSatGrdVeg)/dt at TemperatureGrd [Pa/K]
                  WindStressEwVeg         => noahmp%energy%state%WindStressEwVeg(I,J)    ,& ! out,   wind stress: east-west [N/m2] above canopy
                  WindStressNsVeg         => noahmp%energy%state%WindStressNsVeg(I,J)    ,& ! out,   wind stress: north-south [N/m2] above canopy
                  TemperatureAir2mVeg     => noahmp%energy%state%TemperatureAir2mVeg(I,J),& ! out,   2 m height air temperature [K], vegetated
                  ExchCoeffShLeaf         => noahmp%energy%state%ExchCoeffShLeaf(I,J)    ,& ! out,   sensible heat exchange coeff [m/s],leaf surface to canopy air
                  ExchCoeffShUndCan       => noahmp%energy%state%ExchCoeffShUndCan(I,J)  ,& ! out,   under canopy sensible heat exchange coefficient [m/s]
                  ExchCoeffSh2mVeg        => noahmp%energy%state%ExchCoeffSh2mVeg(I,J)   ,& ! out,   2m sensible heat exchange coefficient [m/s]
                  SpecHumidity2mVeg       => noahmp%energy%state%SpecHumidity2mVeg(I,J)  ,& ! out,   specific humidity [kg/kg] at 2m vegetated
                  RadLwNetVegGrd          => noahmp%energy%flux%RadLwNetVegGrd(I,J)      ,& ! out,   ground net longwave radiation [W/m2] (+ to atm)
                  HeatSensibleCanopy      => noahmp%energy%flux%HeatSensibleCanopy(I,J)  ,& ! in,    canopy sensible heat flux [W/m2] (+ to atm)
                  HeatLatentCanEvap       => noahmp%energy%flux%HeatLatentCanEvap(I,J)   ,& ! in,    canopy evaporation heat flux [W/m2] (+ to atm)
                  HeatLatentCanTransp     => noahmp%energy%flux%HeatLatentCanTransp(I,J) ,& ! in,    canopy transpiration heat flux [W/m2] (+ to atm)
                  HeatSensibleVegGrd      => noahmp%energy%flux%HeatSensibleVegGrd(I,J)  ,& ! out,   vegetated ground sensible heat flux [W/m2] (+ to atm)
                  HeatLatentVegGrd        => noahmp%energy%flux%HeatLatentVegGrd(I,J)    ,& ! out,   ground evaporation heat flux [W/m2] (+ to atm)
                  HeatGroundVegGrd        => noahmp%energy%flux%HeatGroundVegGrd(I,J)     & ! out,   vegetated ground heat [W/m2] (+ to soil/snow)
                 )

        ! under-canopy fluxes and ground temperature
        LwCoeffAir   = -EmissivityGrd * (1.0 - EmissivityVeg) * RadLwDownRefHeight - &
                        EmissivityGrd * EmissivityVeg * ConstStefanBoltzmann * TemperatureCanopy**4
        LwCoeffCan   = EmissivityGrd * ConstStefanBoltzmann
        ShCoeff      = DensityAirRefHeight * ConstHeatCapacAir / ResistanceShUndCan
        LhCoeff      = DensityAirRefHeight * ConstHeatCapacAir / (PsychConstGrd * (ResistanceLhUndCan+ResistanceGrdEvap))  ! Barlage: change to ground v3.6
        GrdHeatCoeff = 2.0 * ThermConductSoilSnow(I,NumSnowLayerNeg+1,J) / ThicknessSnowSoilLayer(I,NumSnowLayerNeg+1,J)

        ! begin stability iteration
        !$acc loop seq
        loop2: do IndIter = 1, NumIterG
           TempTmp = TempUnitConv(TemperatureGrdVeg)
           call VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD)
           if ( TempTmp > 0.0 ) then
              VapPresSatGrdVeg      = VapPresSatWatTmp
              VapPresSatGrdVegTempD = VapPresSatWatTmpD
           else
              VapPresSatGrdVeg      = VapPresSatIceTmp
              VapPresSatGrdVegTempD = VapPresSatIceTmpD
           endif
           RadLwNetVegGrd     = LwCoeffCan * TemperatureGrdVeg**4 + LwCoeffAir
           HeatSensibleVegGrd = ShCoeff * (TemperatureGrdVeg - TemperatureCanopyAir)
           HeatLatentVegGrd   = LhCoeff * (VapPresSatGrdVeg*RelHumidityGrd - PressureVaporCanAir)
           HeatGroundVegGrd   = GrdHeatCoeff * (TemperatureGrdVeg - TemperatureSoilSnow(I,NumSnowLayerNeg+1,J))
           EnergyResTmp       = RadSwAbsGrd - RadLwNetVegGrd - HeatSensibleVegGrd - &
                                HeatLatentVegGrd - HeatGroundVegGrd + HeatPrecipAdvVegGrd
           FluxTotCoeff       = 4.0 * LwCoeffCan * TemperatureGrdVeg**3 + ShCoeff + LhCoeff*VapPresSatGrdVegTempD + GrdHeatCoeff
           TemperatureGrdChg  = EnergyResTmp / FluxTotCoeff
           RadLwNetVegGrd     = RadLwNetVegGrd + 4.0 * LwCoeffCan * TemperatureGrdVeg**3 * TemperatureGrdChg
           HeatSensibleVegGrd = HeatSensibleVegGrd + ShCoeff * TemperatureGrdChg
           HeatLatentVegGrd   = HeatLatentVegGrd + LhCoeff * VapPresSatGrdVegTempD * TemperatureGrdChg
           HeatGroundVegGrd   = HeatGroundVegGrd + GrdHeatCoeff * TemperatureGrdChg
           TemperatureGrdVeg  = TemperatureGrdVeg + TemperatureGrdChg
        enddo loop2
        !TemperatureCanopyAir = (ExchCoeffShAbvCanTmp*TemperatureAirRefHeight + ExchCoeffShLeafTmp*TemperatureCanopy + &
        !                        GrdHeatCoeff*TemperatureGrdVeg)/(ExchCoeffShAbvCanTmp + ExchCoeffShLeafTmp + GrdHeatCoeff)

        ! if snow on ground and TemperatureGrdVeg > freezing point: reset TemperatureGrdVeg = freezing point
        if ( (OptSnowSoilTempTime == 1) .or. (OptSnowSoilTempTime == 3) ) then
           if ( (SnowDepth > 0.05) .and. (TemperatureGrdVeg > ConstFreezePoint) ) then
              if ( OptSnowSoilTempTime == 1 ) &
                 TemperatureGrdVeg = ConstFreezePoint
              if ( OptSnowSoilTempTime == 3 ) &
                 TemperatureGrdVeg = (1.0 - SnowCoverFrac) * TemperatureGrdVeg + SnowCoverFrac * ConstFreezePoint

              RadLwNetVegGrd     = LwCoeffCan * TemperatureGrdVeg**4 - EmissivityGrd * (1.0-EmissivityVeg) * RadLwDownRefHeight - &
                                   EmissivityGrd * EmissivityVeg * ConstStefanBoltzmann * TemperatureCanopy**4
              HeatSensibleVegGrd = ShCoeff * (TemperatureGrdVeg - TemperatureCanopyAir)
              HeatLatentVegGrd   = LhCoeff * (VapPresSatGrdVeg*RelHumidityGrd - PressureVaporCanAir)
              HeatGroundVegGrd   = RadSwAbsGrd + HeatPrecipAdvVegGrd - (RadLwNetVegGrd + HeatSensibleVegGrd + HeatLatentVegGrd)
           endif
        endif

        ! wind stresses
        WindStressEwVeg = -DensityAirRefHeight * ExchCoeffMomAbvCan * WindSpdRefHeight * WindEastwardRefHeight
        WindStressNsVeg = -DensityAirRefHeight * ExchCoeffMomAbvCan * WindSpdRefHeight * WindNorthwardRefHeight

        ! consistent vegetation air temperature and vapor pressure 
        ! since TemperatureGrdVeg is not consistent with the TemperatureCanopyAir/PressureVaporCanAir calculation.
        !TemperatureCanopyAir = TemperatureAirRefHeight + (HeatSensibleVegGrd + HeatSensibleCanopy) / &
        !                       (DensityAirRefHeight*ConstHeatCapacAir*ExchCoeffShAbvCanTmp) 
        !TemperatureCanopyAir = TemperatureAirRefHeight + (HeatSensibleVegGrd * VegFrac + HeatSensibleCanopy) / &
        !                       (DensityAirRefHeight*ConstHeatCapacAir*ExchCoeffShAbvCanTmp)                     ! ground flux need fveg
        !PressureVaporCanAir  = PressureVaporRefHeight + (HeatLatentCanEvap+VegFrac*(HeatLatentCanTransp+HeatLatentVegGrd)) / &
        !                       (DensityAirRefHeight*ExchCoeffLhAbvCan*ConstHeatCapacAir/PsychConstGrd)
        !MoistureFluxSfc      = (SpecHumiditySfc - SpecHumidityRefHeight) * DensityAirRefHeight * ExchCoeffLhAbvCan !*ConstHeatCapacAir/PsychConstGrd
    
        ! 2m temperature over vegetation ( corrected for low LH exchange coeff values )
        if ( (OptSurfaceDrag == 1) .or. (OptSurfaceDrag == 2) ) then
           !ExchCoeffSh2mVeg = FrictionVelVeg * 1.0 / ConstVonKarman * log((2.0+RoughLenShCanopy)/RoughLenShCanopy)
           !ExchCoeffSh2mVeg = FrictionVelVeg * ConstVonKarman / log((2.0+RoughLenShCanopy)/RoughLenShCanopy)
           ExchCoeffSh2mVeg = FrictionVelVeg * ConstVonKarman / (log((2.0+RoughLenShCanopy)/RoughLenShCanopy) - MoStabCorrShVeg2m)
           if ( ExchCoeffSh2mVeg < 1.0e-5 ) then
              TemperatureAir2mVeg = TemperatureCanopyAir
              !SpecHumidity2mVeg   = (PressureVaporCanAir*0.622/(PressureAirRefHeight - 0.378*PressureVaporCanAir))
              SpecHumidity2mVeg   = SpecHumiditySfc
           else
              TemperatureAir2mVeg = TemperatureCanopyAir - (HeatSensibleVegGrd + HeatSensibleCanopy/VegFrac) / &
                                    (DensityAirRefHeight * ConstHeatCapacAir) * 1.0 / ExchCoeffSh2mVeg
              !SpecHumidity2mVeg   = (PressureVaporCanAir*0.622/(PressureAirRefHeight - 0.378*PressureVaporCanAir)) - &
              !                      MoistureFluxSfc/(DensityAirRefHeight*FrictionVelVeg)* 1.0/ConstVonKarman * &
              !                      log((2.0+RoughLenShCanopy)/RoughLenShCanopy)
              SpecHumidity2mVeg   = SpecHumiditySfc - ((HeatLatentCanEvap+HeatLatentCanTransp)/VegFrac + HeatLatentVegGrd) / &
                                                      (LatHeatVapCanopy * DensityAirRefHeight) * 1.0 / ExchCoeffSh2mVeg
           endif
        endif

        ! update ExchCoeffSh for output
        ExchCoeffShAbvCanTmp = 1.0 / ResistanceShAbvCan
        ExchCoeffShLeafTmp   = 2.0 * VegAreaIndTmp(I,J) / ResistanceLeafBoundary
        ExchCoeffShAbvCan    = ExchCoeffShAbvCanTmp
        ExchCoeffShLeaf      = ExchCoeffShLeafTmp
        ExchCoeffShUndCan    = 1.0 / ResistanceShUndCan

        end associate
      end do
    end do
    !$acc end parallel loop

    !$acc end data

  end subroutine SurfaceEnergyFluxVegetated

end module SurfaceEnergyFluxVegetatedMod
