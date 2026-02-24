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
        associate(                                                                      &
                  SnowDepth               => noahmp%water%state%SnowDepth ,& ! in,    snow depth [m]
                  HeightCanopyTop         => noahmp%energy%param%HeightCanopyTop ,& ! in,    top of canopy [m]
                  RadLwDownRefHeight      => noahmp%forcing%RadLwDownRefHeight ,& ! in,    downward longwave radiation [W/m2] at reference height
                  VegAreaIndEff           => noahmp%energy%state%VegAreaIndEff ,& ! in,    one-sided leaf+stem area index [m2/m2]
                  LeafAreaIndSunlit       => noahmp%energy%state%LeafAreaIndSunlit ,& ! in,    sunlit leaf area index, one-sided [m2/m2]
                  LeafAreaIndShade        => noahmp%energy%state%LeafAreaIndShade ,& ! in,    shaded leaf area index, one-sided [m2/m2]
                  ZeroPlaneDispSfc        => noahmp%energy%state%ZeroPlaneDispSfc ,& ! in,    zero plane displacement [m]
                  RoughLenMomSfc          => noahmp%energy%state%RoughLenMomSfc ,& ! in,    roughness length [m], momentum, surface
                  EmissivityVeg           => noahmp%energy%state%EmissivityVeg ,& ! in,    vegetation emissivity
                  EmissivityGrd           => noahmp%energy%state%EmissivityGrd ,& ! in,    ground emissivity
                  RefHeightAboveGrd       => noahmp%energy%state%RefHeightAboveGrd ,& ! in,    surface reference height [m]
                  WindSpdRefHeight        => noahmp%energy%state%WindSpdRefHeight ,& ! in,    wind speed [m/s] at reference height
                  TemperatureGrdVeg       => noahmp%energy%state%TemperatureGrdVeg ,& ! inout, vegetated ground (below-canopy) temperature [K]
                  FrictionVelVeg          => noahmp%energy%state%FrictionVelVeg ,& ! out,   friction velocity [m/s], vegetated
                  MoStabParaAbvCan        => noahmp%energy%state%MoStabParaAbvCan ,& ! out,   Monin-Obukhov stability (z/L), above ZeroPlaneDispSfc, vegetated
                  MoStabCorrShVeg2m       => noahmp%energy%state%MoStabCorrShVeg2m ,& ! out,   M-O sen heat stability correction, 2m, vegetated
                  VapPresSatGrdVeg        => noahmp%energy%state%VapPresSatGrdVeg ,& ! out,   saturation vapor pressure at TemperatureGrd [Pa]
                  CanopyHeight            => noahmp%energy%state%CanopyHeight ,& ! out,   canopy height [m]
                  WindSpdCanopyTop        => noahmp%energy%state%WindSpdCanopyTop ,& ! out,   wind speed at top of canopy [m/s]
                     RoughLenMomGrd          => noahmp%energy%state%RoughLenMomGrd ,& ! in,    roughness length [m], momentum, ground
                     RoughLenShCanopy        => noahmp%energy%state%RoughLenShCanopy ,& ! out,   roughness length [m], sensible heat, vegetated
                     RoughLenShVegGrd        => noahmp%energy%state%RoughLenShVegGrd ,& ! out,   roughness length [m], sensible heat ground, below canopy
                     MainTimeStep            => noahmp%config%domain%MainTimeStep ,& ! in,    main noahmp timestep [s]
                     TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight ,& ! in,    air temperature [K] at reference height
                     PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight ,& ! in,    air pressure [Pa] at reference height
                     CanopyWetFrac           => noahmp%water%state%CanopyWetFrac ,& ! in,    wetted or snowed fraction of the canopy
                     CanopyLiqWater          => noahmp%water%state%CanopyLiqWater ,& ! in,    canopy intercepted liquid water [mm]
                     CanopyIce               => noahmp%water%state%CanopyIce ,& ! in,    canopy intercepted ice [mm]
                     HeatCapacCanFac         => noahmp%energy%param%HeatCapacCanFac ,& ! in,    canopy biomass heat capacity parameter [m]
                     RadSwAbsVeg             => noahmp%energy%flux%RadSwAbsVeg ,& ! in,    solar radiation absorbed by vegetation [W/m2]
                     HeatPrecipAdvCanopy     => noahmp%energy%flux%HeatPrecipAdvCanopy ,& ! in,    precipitation advected heat - vegetation net [W/m2]
                     VegFrac                 => noahmp%energy%state%VegFrac ,& ! in,    greeness vegetation fraction
                     PressureVaporRefHeight  => noahmp%energy%state%PressureVaporRefHeight ,& ! in,    vapor pressure air [Pa] at reference height
                     DensityAirRefHeight     => noahmp%energy%state%DensityAirRefHeight ,& ! in,    density air [kg/m3]
                     ResistanceGrdEvap       => noahmp%energy%state%ResistanceGrdEvap ,& ! in,    ground surface resistance [s/m] to evaporation
                     PsychConstCanopy        => noahmp%energy%state%PsychConstCanopy ,& ! in,    psychrometric constant [Pa/K], canopy
                     LatHeatVapCanopy        => noahmp%energy%state%LatHeatVapCanopy ,& ! in,    latent heat of vaporization/subli [J/kg], canopy
                     ResistanceStomataSunlit => noahmp%energy%state%ResistanceStomataSunlit ,& ! in,  sunlit leaf stomatal resistance [s/m]
                     ResistanceStomataShade  => noahmp%energy%state%ResistanceStomataShade ,& ! in,   shaded leaf stomatal resistance [s/m]
                     ResistanceLeafBoundary  => noahmp%energy%state%ResistanceLeafBoundary ,& ! in,   bulk leaf boundary layer resistance [s/m]
                     ResistanceShAbvCan      => noahmp%energy%state%ResistanceShAbvCan ,& ! in,    aerodynamic resistance for sensible heat [s/m], above canopy
                     ResistanceLhAbvCan      => noahmp%energy%state%ResistanceLhAbvCan ,& ! in,    aerodynamic resistance for water vapor [s/m], above canopy
                     ResistanceShUndCan      => noahmp%energy%state%ResistanceShUndCan ,& ! in,    ground aerodynamic resistance for sensible heat [s/m]
                     ResistanceLhUndCan      => noahmp%energy%state%ResistanceLhUndCan ,& ! in,    ground aerodynamic resistance for water vapor [s/m]
                     SpecHumiditySfc         => noahmp%energy%state%SpecHumiditySfc ,& ! inout, specific humidity at vegetated surface
                     PressureVaporCanAir     => noahmp%energy%state%PressureVaporCanAir ,& ! inout, canopy air vapor pressure [Pa]
                     TemperatureCanopyAir    => noahmp%energy%state%TemperatureCanopyAir ,& ! inout, canopy air temperature [K]
                     TemperatureCanopy       => noahmp%energy%state%TemperatureCanopy ,& ! inout, vegetation temperature [K]
                     ExchCoeffLhAbvCan       => noahmp%energy%state%ExchCoeffLhAbvCan ,& ! out,   latent heat conductance, canopy air to reference height [m/s]
                     ExchCoeffLhTransp       => noahmp%energy%state%ExchCoeffLhTransp ,& ! out,   transpiration conductance, leaf to canopy air [m/s]
                     ExchCoeffLhEvap         => noahmp%energy%state%ExchCoeffLhEvap ,& ! out,   evaporation conductance, leaf to canopy air [m/s]
                     ExchCoeffLhUndCan       => noahmp%energy%state%ExchCoeffLhUndCan ,& ! out,   latent heat conductance, ground to canopy air [m/s]
                     VapPresSatCanopy        => noahmp%energy%state%VapPresSatCanopy ,& ! out,   saturation vapor pressure at TemperatureCanopy [Pa]
                     VapPresSatCanTempD      => noahmp%energy%state%VapPresSatCanTempD ,& ! out,   d(VapPresSatCanopy)/dt at TemperatureCanopy [Pa/K]
                     RadLwNetCanopy          => noahmp%energy%flux%RadLwNetCanopy ,& ! out,   canopy net longwave radiation [W/m2] (+ to atm)
                     HeatSensibleCanopy      => noahmp%energy%flux%HeatSensibleCanopy ,& ! out,   canopy sensible heat flux [W/m2] (+ to atm)
                     HeatLatentCanEvap       => noahmp%energy%flux%HeatLatentCanEvap ,& ! out,   canopy evaporation heat flux [W/m2] (+ to atm)
                     HeatLatentCanTransp     => noahmp%energy%flux%HeatLatentCanTransp ,& ! out,   canopy transpiration heat flux [W/m2] (+ to atm)
                     HeatCanStorageChg       => noahmp%energy%flux%HeatCanStorageChg ,& ! out,   canopy heat storage change [W/m2]
                  NumSnowLayerNeg         => noahmp%config%domain%NumSnowLayerNeg ,& ! in,    actual number of snow layers (negative)
                  ThicknessSnowSoilLayer  => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
                  OptSurfaceDrag          => noahmp%config%nmlist%OptSurfaceDrag ,& ! in,    options for surface layer drag/exchange coefficient
                  OptSnowSoilTempTime     => noahmp%config%nmlist%OptSnowSoilTempTime ,& ! in,    options for snow/soil temperature time scheme (only layer 1)
                  WindEastwardRefHeight   => noahmp%forcing%WindEastwardRefHeight ,& ! in,    wind speed [m/s] in eastward direction at reference height
                  WindNorthwardRefHeight  => noahmp%forcing%WindNorthwardRefHeight ,& ! in,    wind speed [m/s] in northward direction at reference height
                  SnowCoverFrac           => noahmp%water%state%SnowCoverFrac ,& ! in,    snow cover fraction
                  RadSwAbsGrd             => noahmp%energy%flux%RadSwAbsGrd ,& ! in,    solar radiation absorbed by ground [W/m2]
                  HeatPrecipAdvVegGrd     => noahmp%energy%flux%HeatPrecipAdvVegGrd ,& ! in,    precipitation advected heat - under canopy net [W/m2]
                  TemperatureSoilSnow     => noahmp%energy%state%TemperatureSoilSnow ,& ! in,    snow and soil layer temperature [K]
                  ThermConductSoilSnow    => noahmp%energy%state%ThermConductSoilSnow ,& ! in,    thermal conductivity [W/m/K] for all soil & snow
                  PsychConstGrd           => noahmp%energy%state%PsychConstGrd ,& ! in,    psychrometric constant [Pa/K], ground
                  RelHumidityGrd          => noahmp%energy%state%RelHumidityGrd ,& ! in,    raltive humidity in surface soil/snow air space
                  ExchCoeffMomAbvCan      => noahmp%energy%state%ExchCoeffMomAbvCan ,& ! inout, momentum exchange coeff [m/s], above ZeroPlaneDisp, vegetated
                  ExchCoeffShAbvCan       => noahmp%energy%state%ExchCoeffShAbvCan ,& ! inout, heat exchange coeff [m/s], above ZeroPlaneDisp, vegetated
                  VapPresSatGrdVegTempD   => noahmp%energy%state%VapPresSatGrdVegTempD ,& ! out,   d(VapPresSatGrdVeg)/dt at TemperatureGrd [Pa/K]
                  WindStressEwVeg         => noahmp%energy%state%WindStressEwVeg ,& ! out,   wind stress: east-west [N/m2] above canopy
                  WindStressNsVeg         => noahmp%energy%state%WindStressNsVeg ,& ! out,   wind stress: north-south [N/m2] above canopy
                  TemperatureAir2mVeg     => noahmp%energy%state%TemperatureAir2mVeg ,& ! out,   2 m height air temperature [K], vegetated
                  ExchCoeffShLeaf         => noahmp%energy%state%ExchCoeffShLeaf ,& ! out,   sensible heat exchange coeff [m/s],leaf surface to canopy air
                  ExchCoeffShUndCan       => noahmp%energy%state%ExchCoeffShUndCan ,& ! out,   under canopy sensible heat exchange coefficient [m/s]
                  ExchCoeffSh2mVeg        => noahmp%energy%state%ExchCoeffSh2mVeg ,& ! out,   2m sensible heat exchange coefficient [m/s]
                  SpecHumidity2mVeg       => noahmp%energy%state%SpecHumidity2mVeg ,& ! out,   specific humidity [kg/kg] at 2m vegetated
                  RadLwNetVegGrd          => noahmp%energy%flux%RadLwNetVegGrd ,& ! out,   ground net longwave radiation [W/m2] (+ to atm)
                  HeatSensibleVegGrd      => noahmp%energy%flux%HeatSensibleVegGrd ,& ! out,   vegetated ground sensible heat flux [W/m2] (+ to atm)
                  HeatLatentVegGrd        => noahmp%energy%flux%HeatLatentVegGrd ,& ! out,   ground evaporation heat flux [W/m2] (+ to atm)
                  HeatGroundVegGrd        => noahmp%energy%flux%HeatGroundVegGrd  & ! out,   vegetated ground heat [W/m2] (+ to soil/snow)
                 )

    !$acc parallel loop collapse(2) gang vector default(present) &
    !$acc private(LastIter, TemperatureCanChg, TemperatureGrdChg) &
    !$acc private(VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD, TempTmp) &
    !$acc private(LwCoeffAir, LwCoeffCan, MoistureFluxSfc)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
       if ( .not. ((noahmp%energy%state%VegAreaIndEff(I,J) > 0.0 ) .and. (noahmp%energy%state%VegFrac(I,J) > 0)) ) cycle ! skip non-vegetated surface


        ! initialization (including variables that do not depend on stability iteration)
        LastIter          = 0
        FrictionVelVeg(I,J)    = 0.1
        TemperatureCanChg = 0.0
        TemperatureGrdChg = 0.0
        MoStabParaAbvCan(I,J)  = 0.0
        MoStabParaSgn(I,J)     = 0
        MoStabCorrShVeg2m(I,J) = 0.0
        ShGrdTmp(I,J)     = 0.0
        ShCanTmp(I,J)     = 0.0
        MoistureFluxSfc   = 0.0
        ! limit LeafAreaIndex
        VegAreaIndTmp(I,J)     = min(6.0, VegAreaIndEff(I,J))

        ! saturation vapor pressure at ground temperature
        TempTmp = TempUnitConv(TemperatureGrdVeg(I,J))
        call VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD)
        if ( TempTmp > 0.0 ) then
           VapPresSatGrdVeg(I,J) = VapPresSatWatTmp
        else
           VapPresSatGrdVeg(I,J) = VapPresSatIceTmp
        endif

        ! canopy height
        CanopyHeight(I,J) = HeightCanopyTop(I,J)
        ! wind speed at canopy height
       !WindSpdCanopyTop = WindSpdRefHeight * log(CanopyHeight/RoughLenMomSfc) / log(RefHeightAboveGrd/RoughLenMomSfc)
        WindSpdCanopyTop(I,J) = WindSpdRefHeight(I,J) * log((CanopyHeight(I,J) - ZeroPlaneDispSfc(I,J) + RoughLenMomSfc(I,J))/RoughLenMomSfc(I,J)) / &
                           log(RefHeightAboveGrd(I,J)/RoughLenMomSfc(I,J))                                           ! MB: add ZeroPlaneDispSfc(I,J) v3.7
#ifndef _OPENACC
        if ( (CanopyHeight(I,J)-ZeroPlaneDispSfc(I,J)) <= 0.0 ) then
           print*, "CRITICAL PROBLEM: CanopyHeight(I,J) <= ZeroPlaneDispSfc(I,J)"
           print*, "CanopyHeight(I,J) = "         , CanopyHeight(I,J)
           print*, "ZeroPlaneDispSfc(I,J) = "     , ZeroPlaneDispSfc(I,J)
           print*, "SnowDepth(I,J) = "            , SnowDepth(I,J)
           stop "Error: ZeroPlaneDisp problem in NoahMP LSM"
        endif
#endif

      end do
    end do
    !$acc end parallel loop

   !$acc parallel default(present) firstprivate(IndIter) &
   !$acc private(LastIter, TemperatureCanChg, TemperatureGrdChg) &
   !$acc private(LeafAreaIndSunEff, LeafAreaIndShdEff) &
   !$acc private(VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD, TempTmp) &
   !$acc private(LwCoeffAir, LwCoeffCan, ShCoeff, LhCoeff, GrdHeatCoeff, TranspHeatCoeff) &
   !$acc private(ExchCoeffShAbvCanTmp, ExchCoeffShLeafTmp, ExchCoeffTot, TempShGhTmp) &
   !$acc private(ExchCoeffShFrac, VapPresLhTot, ExchCoeffEtFrac, FluxTotCoeff, EnergyResTmp) &
   !$acc private(MoistureFluxSfc, HeatCapacCan)

    ! begin stability iteration for canopy temperature and flux
    loop1: do IndIter = 1, NumIterC

       ! Roughness length calculation
       !$acc loop gang vector collapse(2) 
       do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
            
         if ( .not. ((noahmp%energy%state%VegAreaIndEff(I,J) > 0.0 ) .and. (noahmp%energy%state%VegFrac(I,J) > 0)) ) cycle ! skip non-vegetated surface


           ! ground and surface roughness length
           if ( IndIter == 1 ) then
              RoughLenShCanopy(I,J) = RoughLenMomSfc(I,J)
              RoughLenShVegGrd(I,J) = RoughLenMomGrd(I,J)
           else
              RoughLenShCanopy(I,J) = RoughLenMomSfc(I,J)  !* exp(-ZilitinkevichCoeff*0.4*258.2*sqrt(FrictionVelVeg(I,J)*RoughLenMomSfc(I,J)))
              RoughLenShVegGrd(I,J) = RoughLenMomGrd(I,J)  !* exp(-ZilitinkevichCoeff*0.4*258.2*sqrt(FrictionVelVeg(I,J)*RoughLenMomGrd(I,J)))
           endif

         end do
       end do

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
       !$acc loop gang vector collapse(2) 
       do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if ( .not. ((noahmp%energy%state%VegAreaIndEff(I,J) > 0.0 ) .and. (noahmp%energy%state%VegFrac(I,J) > 0)) ) cycle ! skip non-vegetated surface


           ! limit LeafAreaIndex
           LeafAreaIndSunEff = min(6.0, LeafAreaIndSunlit(I,J))
           LeafAreaIndShdEff = min(6.0, LeafAreaIndShade(I,J))

           ! prepare for longwave rad.
           LwCoeffAir = -EmissivityVeg(I,J) * (1.0 + (1.0-EmissivityVeg(I,J))*(1.0-EmissivityGrd(I,J))) * RadLwDownRefHeight(I,J) - &
                         EmissivityVeg(I,J) * EmissivityGrd(I,J) * ConstStefanBoltzmann * TemperatureGrdVeg(I,J)**4
           LwCoeffCan = (2.0 - EmissivityVeg(I,J) * (1.0-EmissivityGrd(I,J))) * EmissivityVeg(I,J) * ConstStefanBoltzmann

           ! ES and d(ES)/dt evaluated at TemperatureCanopy
           TempTmp = TempUnitConv(TemperatureCanopy(I,J))
           call VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD)
           if ( TempTmp > 0.0 ) then
              VapPresSatCanopy(I,J)   = VapPresSatWatTmp
              VapPresSatCanTempD(I,J) = VapPresSatWatTmpD
           else
              VapPresSatCanopy(I,J)   = VapPresSatIceTmp
              VapPresSatCanTempD(I,J) = VapPresSatIceTmpD
           endif

           ! sensible heat conductance and coeff above veg.
           ExchCoeffShAbvCanTmp = 1.0 / ResistanceShAbvCan(I,J)
           ExchCoeffShLeafTmp   = 2.0 * VegAreaIndTmp(I,J) / ResistanceLeafBoundary(I,J)
           GrdHeatCoeff         = 1.0 / ResistanceShUndCan(I,J)
           ExchCoeffTot         = ExchCoeffShAbvCanTmp + ExchCoeffShLeafTmp + GrdHeatCoeff
           TempShGhTmp          = (TemperatureAirRefHeight(I,J)*ExchCoeffShAbvCanTmp + TemperatureGrdVeg(I,J)*GrdHeatCoeff) / ExchCoeffTot
           ExchCoeffShFrac      = ExchCoeffShLeafTmp / ExchCoeffTot
           ShCoeff              = (1.0 - ExchCoeffShFrac) * DensityAirRefHeight(I,J) * ConstHeatCapacAir * ExchCoeffShLeafTmp

           ! latent heat conductance and coeff above veg.
           ExchCoeffLhAbvCan(I,J) = 1.0 / ResistanceLhAbvCan(I,J)
           ExchCoeffLhEvap(I,J)   = CanopyWetFrac(I,J) * VegAreaIndTmp(I,J) / ResistanceLeafBoundary(I,J)
           ExchCoeffLhTransp(I,J) = (1.0 - CanopyWetFrac(I,J)) * (LeafAreaIndSunEff/(ResistanceLeafBoundary(I,J)+ResistanceStomataSunlit(I,J)) + &
                                                        LeafAreaIndShdEff/(ResistanceLeafBoundary(I,J)+ResistanceStomataShade(I,J)))
           ExchCoeffLhUndCan(I,J) = 1.0 / (ResistanceLhUndCan(I,J) + ResistanceGrdEvap(I,J))
           ExchCoeffTot      = ExchCoeffLhAbvCan(I,J) + ExchCoeffLhEvap(I,J) + ExchCoeffLhTransp(I,J) + ExchCoeffLhUndCan(I,J)
           VapPresLhTot      = (PressureVaporRefHeight(I,J)*ExchCoeffLhAbvCan(I,J) + VapPresSatGrdVeg(I,J)*ExchCoeffLhUndCan(I,J) ) / ExchCoeffTot
           ExchCoeffEtFrac   = (ExchCoeffLhEvap(I,J) + ExchCoeffLhTransp(I,J)) / ExchCoeffTot
           LhCoeff           = (1.0 - ExchCoeffEtFrac) * ExchCoeffLhEvap(I,J) * DensityAirRefHeight(I,J) * &
                               ConstHeatCapacAir / PsychConstCanopy(I,J)
           TranspHeatCoeff   = (1.0 - ExchCoeffEtFrac) * ExchCoeffLhTransp(I,J) * DensityAirRefHeight(I,J) * &
                               ConstHeatCapacAir / PsychConstCanopy(I,J)

           ! evaluate surface fluxes with current temperature and solve for temperature change
           TemperatureCanopyAir(I,J) = TempShGhTmp + ExchCoeffShFrac * TemperatureCanopy(I,J)
           PressureVaporCanAir(I,J)  = VapPresLhTot + ExchCoeffEtFrac * VapPresSatCanopy(I,J)
           RadLwNetCanopy(I,J)       = VegFrac(I,J) * (LwCoeffAir + LwCoeffCan * TemperatureCanopy(I,J)**4)
           HeatSensibleCanopy(I,J)   = VegFrac(I,J) * DensityAirRefHeight(I,J) * ConstHeatCapacAir * &
                                  ExchCoeffShLeafTmp * (TemperatureCanopy(I,J) - TemperatureCanopyAir(I,J))
           HeatLatentCanEvap(I,J)    = VegFrac(I,J) * DensityAirRefHeight(I,J) * ConstHeatCapacAir * ExchCoeffLhEvap(I,J) * &
                                  (VapPresSatCanopy(I,J) - PressureVaporCanAir(I,J)) / PsychConstCanopy(I,J)
           HeatLatentCanTransp(I,J)  = VegFrac(I,J) * DensityAirRefHeight(I,J) * ConstHeatCapacAir * ExchCoeffLhTransp(I,J) * &
                                  (VapPresSatCanopy(I,J) - PressureVaporCanAir(I,J)) / PsychConstCanopy(I,J)
           if ( TemperatureCanopy(I,J) > ConstFreezePoint ) then
              HeatLatentCanEvap(I,J) = min(CanopyLiqWater(I,J)*LatHeatVapCanopy(I,J)/MainTimeStep, HeatLatentCanEvap(I,J))
           else
              HeatLatentCanEvap(I,J) = min(CanopyIce(I,J)*LatHeatVapCanopy(I,J)/MainTimeStep, HeatLatentCanEvap(I,J))
           endif
           ! canopy heat capacity
           HeatCapacCan         = HeatCapacCanFac(I,J)*VegAreaIndTmp(I,J)*ConstHeatCapacWater + CanopyLiqWater(I,J)*ConstHeatCapacWater/ConstDensityWater + &
                                  CanopyIce(I,J)*ConstHeatCapacIce/ConstDensityIce
           ! compute vegetation temperature change
           EnergyResTmp         = RadSwAbsVeg(I,J) - RadLwNetCanopy(I,J) - HeatSensibleCanopy(I,J) - &
                                  HeatLatentCanEvap(I,J) - HeatLatentCanTransp(I,J) + HeatPrecipAdvCanopy(I,J)
           FluxTotCoeff         = VegFrac(I,J) * (4.0*LwCoeffCan*TemperatureCanopy(I,J)**3 + ShCoeff + &
                                            (LhCoeff+TranspHeatCoeff)*VapPresSatCanTempD(I,J) + HeatCapacCan/MainTimeStep)
           TemperatureCanChg    = EnergyResTmp / FluxTotCoeff
           ! update fluxes with temperature change
           RadLwNetCanopy(I,J)       = RadLwNetCanopy(I,J)      + VegFrac(I,J) * 4.0 * LwCoeffCan * TemperatureCanopy(I,J)**3 * TemperatureCanChg
           HeatSensibleCanopy(I,J)   = HeatSensibleCanopy(I,J)  + VegFrac(I,J) * ShCoeff * TemperatureCanChg
           HeatLatentCanEvap(I,J)    = HeatLatentCanEvap(I,J)   + VegFrac(I,J) * LhCoeff * VapPresSatCanTempD(I,J) * TemperatureCanChg
           HeatLatentCanTransp(I,J)  = HeatLatentCanTransp(I,J) + VegFrac(I,J) * TranspHeatCoeff * VapPresSatCanTempD(I,J) * TemperatureCanChg
           HeatCanStorageChg(I,J)    = VegFrac(I,J) * HeatCapacCan / MainTimeStep * TemperatureCanChg
           ! update vegetation temperature
           TemperatureCanopy(I,J)    = TemperatureCanopy(I,J) + TemperatureCanChg
          !TemperatureCanopyAir = TempShGhTmp + ExchCoeffShFrac * TemperatureCanopy                        ! canopy air T; update here for consistency

            ! for computing M-O length in the next iteration
            ShCanTmp(I,J) = DensityAirRefHeight(I,J) * ConstHeatCapacAir * (TemperatureCanopyAir(I,J)-TemperatureAirRefHeight(I,J)) / ResistanceShAbvCan(I,J)
            ShGrdTmp(I,J) = DensityAirRefHeight(I,J) * ConstHeatCapacAir * (TemperatureGrdVeg(I,J)-TemperatureCanopyAir(I,J)) / ResistanceShUndCan(I,J)

           ! consistent specific humidity from canopy air vapor pressure
           SpecHumiditySfc(I,J) = (0.622 * PressureVaporCanAir(I,J)) / (PressureAirRefHeight(I,J) - 0.378 * PressureVaporCanAir(I,J))
         !   if ( LastIter == 1 ) then
         !      exit loop1
         !   endif
         !   if ( (IndIter >= 5) .and. (abs(TemperatureCanChg) <= 0.01) .and. (LastIter == 0) ) then
         !      LastIter = 1
         !   endif

         end do
       end do

    enddo loop1  ! end stability iteration

    ! Ground temperature iteration (loop2) and final calculations
    !$acc loop gang vector collapse(2)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
        if ( .not. ((noahmp%energy%state%VegAreaIndEff(I,J) > 0.0 ) .and. (noahmp%energy%state%VegFrac(I,J) > 0)) ) cycle ! skip non-vegetated surface


        ! under-canopy fluxes and ground temperature
        LwCoeffAir   = -EmissivityGrd(I,J) * (1.0 - EmissivityVeg(I,J)) * RadLwDownRefHeight(I,J) - &
                        EmissivityGrd(I,J) * EmissivityVeg(I,J) * ConstStefanBoltzmann * TemperatureCanopy(I,J)**4
        LwCoeffCan   = EmissivityGrd(I,J) * ConstStefanBoltzmann
        ShCoeff      = DensityAirRefHeight(I,J) * ConstHeatCapacAir / ResistanceShUndCan(I,J)
        LhCoeff      = DensityAirRefHeight(I,J) * ConstHeatCapacAir / (PsychConstGrd(I,J) * (ResistanceLhUndCan(I,J)+ResistanceGrdEvap(I,J)))  ! Barlage: change to ground v3.6
        GrdHeatCoeff = 2.0 * ThermConductSoilSnow(I,NumSnowLayerNeg(I,J)+1,J) / ThicknessSnowSoilLayer(I,NumSnowLayerNeg(I,J)+1,J)

        ! begin stability iteration
        !$acc loop seq
        loop2: do IndIter = 1, NumIterG
           TempTmp = TempUnitConv(TemperatureGrdVeg(I,J))
           call VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD)
           if ( TempTmp > 0.0 ) then
              VapPresSatGrdVeg(I,J)      = VapPresSatWatTmp
              VapPresSatGrdVegTempD(I,J) = VapPresSatWatTmpD
           else
              VapPresSatGrdVeg(I,J)      = VapPresSatIceTmp
              VapPresSatGrdVegTempD(I,J) = VapPresSatIceTmpD
           endif
           RadLwNetVegGrd(I,J)     = LwCoeffCan * TemperatureGrdVeg(I,J)**4 + LwCoeffAir
           HeatSensibleVegGrd(I,J) = ShCoeff * (TemperatureGrdVeg(I,J) - TemperatureCanopyAir(I,J))
           HeatLatentVegGrd(I,J)   = LhCoeff * (VapPresSatGrdVeg(I,J)*RelHumidityGrd(I,J) - PressureVaporCanAir(I,J))
           HeatGroundVegGrd(I,J)   = GrdHeatCoeff * (TemperatureGrdVeg(I,J) - TemperatureSoilSnow(I,NumSnowLayerNeg(I,J)+1,J))
           EnergyResTmp       = RadSwAbsGrd(I,J) - RadLwNetVegGrd(I,J) - HeatSensibleVegGrd(I,J) - &
                                HeatLatentVegGrd(I,J) - HeatGroundVegGrd(I,J) + HeatPrecipAdvVegGrd(I,J)
           FluxTotCoeff       = 4.0 * LwCoeffCan * TemperatureGrdVeg(I,J)**3 + ShCoeff + LhCoeff*VapPresSatGrdVegTempD(I,J) + GrdHeatCoeff
           TemperatureGrdChg  = EnergyResTmp / FluxTotCoeff
           RadLwNetVegGrd(I,J)     = RadLwNetVegGrd(I,J) + 4.0 * LwCoeffCan * TemperatureGrdVeg(I,J)**3 * TemperatureGrdChg
           HeatSensibleVegGrd(I,J) = HeatSensibleVegGrd(I,J) + ShCoeff * TemperatureGrdChg
           HeatLatentVegGrd(I,J)   = HeatLatentVegGrd(I,J) + LhCoeff * VapPresSatGrdVegTempD(I,J) * TemperatureGrdChg
           HeatGroundVegGrd(I,J)   = HeatGroundVegGrd(I,J) + GrdHeatCoeff * TemperatureGrdChg
           TemperatureGrdVeg(I,J)  = TemperatureGrdVeg(I,J) + TemperatureGrdChg
        enddo loop2
        !TemperatureCanopyAir = (ExchCoeffShAbvCanTmp*TemperatureAirRefHeight + ExchCoeffShLeafTmp*TemperatureCanopy + &
        !                        GrdHeatCoeff*TemperatureGrdVeg)/(ExchCoeffShAbvCanTmp + ExchCoeffShLeafTmp + GrdHeatCoeff)

        ! if snow on ground and TemperatureGrdVeg > freezing point: reset TemperatureGrdVeg = freezing point
        if ( (OptSnowSoilTempTime == 1) .or. (OptSnowSoilTempTime == 3) ) then
           if ( (SnowDepth(I,J) > 0.05) .and. (TemperatureGrdVeg(I,J) > ConstFreezePoint) ) then
              if ( OptSnowSoilTempTime == 1 ) &
                 TemperatureGrdVeg(I,J) = ConstFreezePoint
              if ( OptSnowSoilTempTime == 3 ) &
                 TemperatureGrdVeg(I,J) = (1.0 - SnowCoverFrac(I,J)) * TemperatureGrdVeg(I,J) + SnowCoverFrac(I,J) * ConstFreezePoint

              RadLwNetVegGrd(I,J)     = LwCoeffCan * TemperatureGrdVeg(I,J)**4 - EmissivityGrd(I,J) * (1.0-EmissivityVeg(I,J)) * RadLwDownRefHeight(I,J) - &
                                   EmissivityGrd(I,J) * EmissivityVeg(I,J) * ConstStefanBoltzmann * TemperatureCanopy(I,J)**4
              HeatSensibleVegGrd(I,J) = ShCoeff * (TemperatureGrdVeg(I,J) - TemperatureCanopyAir(I,J))
              HeatLatentVegGrd(I,J)   = LhCoeff * (VapPresSatGrdVeg(I,J)*RelHumidityGrd(I,J) - PressureVaporCanAir(I,J))
              HeatGroundVegGrd(I,J)   = RadSwAbsGrd(I,J) + HeatPrecipAdvVegGrd(I,J) - (RadLwNetVegGrd(I,J) + HeatSensibleVegGrd(I,J) + HeatLatentVegGrd(I,J))
           endif
        endif

        ! wind stresses
        WindStressEwVeg(I,J) = -DensityAirRefHeight(I,J) * ExchCoeffMomAbvCan(I,J) * WindSpdRefHeight(I,J) * WindEastwardRefHeight(I,J)
        WindStressNsVeg(I,J) = -DensityAirRefHeight(I,J) * ExchCoeffMomAbvCan(I,J) * WindSpdRefHeight(I,J) * WindNorthwardRefHeight(I,J)

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
           ExchCoeffSh2mVeg(I,J) = FrictionVelVeg(I,J) * ConstVonKarman / (log((2.0+RoughLenShCanopy(I,J))/RoughLenShCanopy(I,J)) - MoStabCorrShVeg2m(I,J))
           if ( ExchCoeffSh2mVeg(I,J) < 1.0e-5 ) then
              TemperatureAir2mVeg(I,J) = TemperatureCanopyAir(I,J)
              !SpecHumidity2mVeg   = (PressureVaporCanAir*0.622/(PressureAirRefHeight - 0.378*PressureVaporCanAir))
              SpecHumidity2mVeg(I,J)   = SpecHumiditySfc(I,J)
           else
              TemperatureAir2mVeg(I,J) = TemperatureCanopyAir(I,J) - (HeatSensibleVegGrd(I,J) + HeatSensibleCanopy(I,J)/VegFrac(I,J)) / &
                                    (DensityAirRefHeight(I,J) * ConstHeatCapacAir) * 1.0 / ExchCoeffSh2mVeg(I,J)
              !SpecHumidity2mVeg   = (PressureVaporCanAir*0.622/(PressureAirRefHeight - 0.378*PressureVaporCanAir)) - &
              !                      MoistureFluxSfc/(DensityAirRefHeight*FrictionVelVeg)* 1.0/ConstVonKarman * &
              !                      log((2.0+RoughLenShCanopy)/RoughLenShCanopy)
              SpecHumidity2mVeg(I,J)   = SpecHumiditySfc(I,J) - ((HeatLatentCanEvap(I,J)+HeatLatentCanTransp(I,J))/VegFrac(I,J) + HeatLatentVegGrd(I,J)) / &
                                                      (LatHeatVapCanopy(I,J) * DensityAirRefHeight(I,J)) * 1.0 / ExchCoeffSh2mVeg(I,J)
           endif
        endif

        ! update ExchCoeffSh for output
        ExchCoeffShAbvCanTmp = 1.0 / ResistanceShAbvCan(I,J)
        ExchCoeffShLeafTmp   = 2.0 * VegAreaIndTmp(I,J) / ResistanceLeafBoundary(I,J)
        ExchCoeffShAbvCan(I,J)    = ExchCoeffShAbvCanTmp
        ExchCoeffShLeaf(I,J)      = ExchCoeffShLeafTmp
        ExchCoeffShUndCan(I,J)    = 1.0 / ResistanceShUndCan(I,J)

      end do
    end do
    !$acc end parallel



    end associate
    !$acc end data

  end subroutine SurfaceEnergyFluxVegetated

end module SurfaceEnergyFluxVegetatedMod
