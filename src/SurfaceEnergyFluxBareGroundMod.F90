module SurfaceEnergyFluxBareGroundMod

!!! Compute surface energy fluxes and budget for bare ground
!!! Use newton-raphson iteration to solve for ground temperatures
!!! Surface energy balance (bare soil):
!!! Ground level: -RadSwAbsGrd - HeatPrecipAdvBareGrd + RadLwNetBareGrd + HeatSensibleBareGrd + HeatLatentBareGrd + HeatGroundBareGrd = 0

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use VaporPressureSaturationMod,    only : VaporPressureSaturation
  use ResistanceBareGroundMostMod,   only : ResistanceBareGroundMOST
  use ResistanceBareGroundChen97Mod, only : ResistanceBareGroundChen97

  implicit none

contains

  subroutine SurfaceEnergyFluxBareGround(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: BARE_FLUX
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type)     , intent(inout) :: noahmp

! local variable
    integer                               :: I, J                 ! grid indices
    integer                               :: IndIter              ! iteration index
    integer, allocatable, dimension(:,:)  :: MoStabParaSgn        ! number of times MoStabParaBare changes sign
    integer, parameter                    :: NumIter = 5          ! number of iterations for surface temperature
    real(kind=kind_noahmp)                :: TemperatureGrdChg    ! change in ground temperature, last iteration [K]
    real(kind=kind_noahmp)                :: LwRadCoeff           ! coefficients for LW radiation as function of ts**4
    real(kind=kind_noahmp)                :: ShCoeff              ! coefficients for sensible heat as function of ts
    real(kind=kind_noahmp)                :: LhCoeff              ! coefficients for latent heat as function of ts
    real(kind=kind_noahmp)                :: GrdHeatCoeff         ! coefficients for ground heat as function of ts
    real(kind=kind_noahmp)                :: ExchCoeffShTmp       ! temporary sensible heat exchange coefficient [m/s]
    real(kind=kind_noahmp)                :: ExchCoeffMomTmp      ! temporary momentum heat exchange coefficient [m/s]
    real(kind=kind_noahmp)                :: MoistureFluxSfc      ! moisture flux
    real(kind=kind_noahmp)                :: VapPresSatWatTmp     ! saturated vapor pressure for water [Pa]
    real(kind=kind_noahmp)                :: VapPresSatIceTmp     ! saturated vapor pressure for ice [Pa]
    real(kind=kind_noahmp)                :: VapPresSatWatTmpD    ! saturated vapor pressure gradient with ground temp. [Pa/K] for water
    real(kind=kind_noahmp)                :: VapPresSatIceTmpD    ! saturated vapor pressure gradient with ground temp. [Pa/K] for ice
    real(kind=kind_noahmp)                :: FluxTotCoeff         ! temporary total coefficients for all energy flux
    real(kind=kind_noahmp)                :: EnergyResTmp         ! temporary energy residual
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatSensibleTmp      ! temporary sensible heat flux [W/m2]
    real(kind=kind_noahmp)                :: TempTmp              ! temporary temperature
    real(kind=kind_noahmp)                :: SoilIceTmp           ! temporary variable for soil ice content
    real(kind=kind_noahmp)                :: TempUnitConv         ! Kelvin to degree Celsius with limit -50 to +50
! local statement function
    TempUnitConv(TempTmp) = min(50.0, max(-50.0, (TempTmp-ConstFreezePoint)))

! --------------------------------------------------------------------
    allocate(MoStabParaSgn(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                            noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(HeatSensibleTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                             noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    associate(                                                                      &
              RoughLenMomGrd          => noahmp%energy%state%RoughLenMomGrd ,& ! in,    roughness length, momentum, ground [m]
              RoughLenShBareGrd       => noahmp%energy%state%RoughLenShBareGrd ,& ! out,   roughness length [m], sensible heat, bare ground
              FrictionVelBare         => noahmp%energy%state%FrictionVelBare ,& ! inout, friction velocity [m/s], bare ground
              MoStabParaBare          => noahmp%energy%state%MoStabParaBare ,& ! out,   Monin-Obukhov stability (z/L), above ZeroPlaneDisp, bare ground
              MoStabCorrShBare2m      => noahmp%energy%state%MoStabCorrShBare2m ,& ! out,   M-O sen heat stability correction, 2m, bare ground
              NumSnowLayerNeg         => noahmp%config%domain%NumSnowLayerNeg ,& ! in,    actual number of snow layers (negative)
              OptGlacierTreatment     => noahmp%config%nmlist%OptGlacierTreatment ,& ! in,    options for glacier treatment
              IndicatorIceSfc         => noahmp%config%domain%IndicatorIceSfc ,& ! in,    flag to identify ice surface
              SnowDepth               => noahmp%water%state%SnowDepth ,& ! in,    snow depth [m]
              ThicknessSnowSoilLayer  => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
              RadLwDownRefHeight      => noahmp%forcing%RadLwDownRefHeight ,& ! in,    downward longwave radiation [W/m2] at reference height
              TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight ,& ! in,    air temperature [K] at reference height
              PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight ,& ! in,    air pressure [Pa] at surface reference height
              SpecHumidityRefHeight   => noahmp%forcing%SpecHumidityRefHeight ,& ! in,    specific humidity [kg/kg] at reference height
              DensityAirRefHeight     => noahmp%energy%state%DensityAirRefHeight ,& ! in,    density air [kg/m3]
              RelHumidityGrd          => noahmp%energy%state%RelHumidityGrd ,& ! in,    raltive humidity in surface soil/snow air space
              EmissivityGrd           => noahmp%energy%state%EmissivityGrd ,& ! in,    ground emissivity
              TemperatureSoilSnow     => noahmp%energy%state%TemperatureSoilSnow ,& ! in,    snow and soil layer temperature [K]
              ThermConductSoilSnow    => noahmp%energy%state%ThermConductSoilSnow ,& ! in,    thermal conductivity [W/m/K] for all soil & snow
              ResistanceGrdEvap       => noahmp%energy%state%ResistanceGrdEvap ,& ! in,    ground surface resistance [s/m] to evaporation
              PressureVaporRefHeight  => noahmp%energy%state%PressureVaporRefHeight ,& ! in,    vapor pressure air [Pa] at reference height
              PsychConstGrd           => noahmp%energy%state%PsychConstGrd ,& ! in,    psychrometric constant [Pa/K], ground
              RadSwAbsGrd             => noahmp%energy%flux%RadSwAbsGrd ,& ! in,    solar radiation absorbed by ground [W/m2]
              HeatPrecipAdvBareGrd    => noahmp%energy%flux%HeatPrecipAdvBareGrd ,& ! in,    precipitation advected heat - bare ground net [W/m2]
              ResistanceLhBareGrd     => noahmp%energy%state%ResistanceLhBareGrd ,& ! in,    aerodynamic resistance for water vapor [s/m], bare ground
              ResistanceShBareGrd     => noahmp%energy%state%ResistanceShBareGrd ,& ! in,    aerodynamic resistance for sensible heat [s/m], bare ground
              ResistanceMomBareGrd    => noahmp%energy%state%ResistanceMomBareGrd ,& ! in,   aerodynamic resistance for momentum [s/m], bare ground
              SpecHumiditySfc         => noahmp%energy%state%SpecHumiditySfc ,& ! inout, specific humidity [kg/kg] at bare surface
              TemperatureGrdBare      => noahmp%energy%state%TemperatureGrdBare ,& ! inout, bare ground temperature (K)
              VapPresSatGrdBare       => noahmp%energy%state%VapPresSatGrdBare ,& ! out,   bare ground saturation vapor pressure [Pa]
              VapPresSatGrdBareTempD  => noahmp%energy%state%VapPresSatGrdBareTempD ,& ! out, bare ground d(VapPresSat)/dt [Pa/K]
              RadLwNetBareGrd         => noahmp%energy%flux%RadLwNetBareGrd ,& ! out,   net longwave rad [W/m2] bare ground (+ to atm)
              HeatSensibleBareGrd     => noahmp%energy%flux%HeatSensibleBareGrd ,& ! out,   sensible heat flux [W/m2] bare ground (+ to atm)
              HeatLatentBareGrd       => noahmp%energy%flux%HeatLatentBareGrd ,& ! out,   latent heat flux [W/m2] bare ground (+ to atm)
              HeatGroundBareGrd       => noahmp%energy%flux%HeatGroundBareGrd ,& ! out,   bare ground heat flux [W/m2] (+ to soil/snow)
              NumSoilLayer            => noahmp%config%domain%NumSoilLayer ,& ! in,    number of soil layers
              FlagUrban               => noahmp%config%domain%FlagUrban ,& ! in,    logical flag for urban grid
              OptSurfaceDrag          => noahmp%config%nmlist%OptSurfaceDrag ,& ! in,    options for surface layer drag/exchange coefficient
              OptSnowSoilTempTime     => noahmp%config%nmlist%OptSnowSoilTempTime ,& ! in,    options for snow/soil temperature time scheme (only layer 1)
              WindEastwardRefHeight   => noahmp%forcing%WindEastwardRefHeight ,& ! in,    wind speed [m/s] in eastward direction at reference height
              WindNorthwardRefHeight  => noahmp%forcing%WindNorthwardRefHeight ,& ! in,    wind speed [m/s] in northward direction at reference height
              SnowCoverFrac           => noahmp%water%state%SnowCoverFrac ,& ! in,    snow cover fraction
              SoilMoisture            => noahmp%water%state%SoilMoisture ,& ! in,    soil moisture for all soil layers [m3/m3]
              SoilLiqWater            => noahmp%water%state%SoilLiqWater ,& ! in,    soil liquid water content for all soil layers [m3/m3]
              WindSpdRefHeight        => noahmp%energy%state%WindSpdRefHeight ,& ! in,    wind speed [m/s] at reference height
              LatHeatVapGrd           => noahmp%energy%state%LatHeatVapGrd ,& ! in,    latent heat of vaporization/subli [J/kg], ground
              ExchCoeffMomBare        => noahmp%energy%state%ExchCoeffMomBare ,& ! inout, momentum exchange coeff [m/s], above ZeroPlaneDisp, bare ground
              ExchCoeffShBare         => noahmp%energy%state%ExchCoeffShBare ,& ! inout, heat exchange coeff [m/s], above ZeroPlaneDisp, bare ground
              WindStressEwBare        => noahmp%energy%state%WindStressEwBare ,& ! out,   wind stress: east-west [N/m2] bare ground
              WindStressNsBare        => noahmp%energy%state%WindStressNsBare ,& ! out,   wind stress: north-south [N/m2] bare ground
              TemperatureAir2mBare    => noahmp%energy%state%TemperatureAir2mBare ,& ! out,   2 m height air temperature [K] bare ground
              SpecHumidity2mBare      => noahmp%energy%state%SpecHumidity2mBare ,& ! out,   bare ground 2-m specific humidity [kg/kg]
              ExchCoeffSh2mBare       => noahmp%energy%state%ExchCoeffSh2mBare  & ! out,   bare ground 2-m sensible heat exchange coefficient [m/s]
                )

    !$acc data create(MoStabParaSgn, HeatSensibleTmp)
    ! begin stability iteration for ground temperature and flux
    loop3: do IndIter = 1, NumIter

       !$acc parallel loop collapse(2) gang vector default(present) private(TemperatureGrdChg, LwRadCoeff, ShCoeff, &
       !$acc LhCoeff, GrdHeatCoeff) private(ExchCoeffShTmp, ExchCoeffMomTmp, MoistureFluxSfc, VapPresSatWatTmp, &
       !$acc VapPresSatIceTmp) private(VapPresSatWatTmpD, VapPresSatIceTmpD, FluxTotCoeff, EnergyResTmp, TempTmp) &
       !$acc firstprivate(IndIter)
       do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE



           ! ground roughness length
           if ( IndIter == 1 ) then
              RoughLenShBareGrd(I,J) = RoughLenMomGrd(I,J)
              ! initialization (including variables that do not depend on stability iteration)
              TemperatureGrdChg  = 0.0
              MoStabParaBare(I,J)     = 0.0
              MoStabParaSgn(I,J)      = 0
              MoStabCorrShBare2m(I,J) = 0.0
              MoistureFluxSfc    = 0.0
              FrictionVelBare(I,J)    = 0.1
              HeatSensibleTmp(I,J)    = 0.0

           else
              RoughLenShBareGrd(I,J) = RoughLenMomGrd(I,J) !* exp(-ZilitinkevichCoeff*0.4*258.2*sqrt(FrictionVelBare(I,J)*RoughLenMomGrd(I,J)))
           endif

         end do
       end do
       !$acc end parallel loop

       ! aerodyn resistances between reference heigths and d+z0v
       if ( noahmp%config%nmlist%OptSurfaceDrag == 1 ) call ResistanceBareGroundMOST(noahmp, IndIter, HeatSensibleTmp, MoStabParaSgn)
       if ( noahmp%config%nmlist%OptSurfaceDrag == 2 ) call ResistanceBareGroundChen97(noahmp, IndIter)

       ! Second parallel region: flux calculations
       !$acc parallel loop collapse(2) gang vector default(present) &
       !$acc private(TemperatureGrdChg, LwRadCoeff, ShCoeff, LhCoeff, GrdHeatCoeff) &
       !$acc private(ExchCoeffShTmp, ExchCoeffMomTmp, MoistureFluxSfc, VapPresSatWatTmp, VapPresSatIceTmp) &
       !$acc private(VapPresSatWatTmpD, VapPresSatIceTmpD, FluxTotCoeff, EnergyResTmp, TempTmp)
       do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


           ! conductance variables for diagnostics
           ExchCoeffMomTmp = 1.0 / ResistanceMomBareGrd(I,J)
           ExchCoeffShTmp  = 1.0 / ResistanceShBareGrd(I,J)

           ! coefficients
           LwRadCoeff   = EmissivityGrd(I,J) * ConstStefanBoltzmann
           GrdHeatCoeff = 2.0*ThermConductSoilSnow(I,NumSnowLayerNeg(I,J)+1,J)/ThicknessSnowSoilLayer(I,NumSnowLayerNeg(I,J)+1,J)

           ! ES and d(ES)/dt evaluated at ground temperature
           TempTmp = min(50.0, max(-50.0, (TemperatureGrdBare(I,J)-ConstFreezePoint)))
           call VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD)
           if ( TempTmp > 0.0 ) then
              VapPresSatGrdBare(I,J)      = VapPresSatWatTmp
              VapPresSatGrdBareTempD(I,J) = VapPresSatWatTmpD
           else
              VapPresSatGrdBare(I,J)      = VapPresSatIceTmp
              VapPresSatGrdBareTempD(I,J) = VapPresSatIceTmpD
           endif

           ! ground fluxes and temperature change
           ShCoeff             = DensityAirRefHeight(I,J) * ConstHeatCapacAir / ResistanceShBareGrd(I,J)
           LhCoeff             = DensityAirRefHeight(I,J) * ConstHeatCapacAir / PsychConstGrd(I,J) / (ResistanceGrdEvap(I,J)+ResistanceLhBareGrd(I,J))
           if ( .not.((SnowDepth(I,J) > 0.0) .or. (OptGlacierTreatment == 1)) .and. (IndicatorIceSfc(I,J) == -1) ) then
              LhCoeff = 0.0   ! don't allow any sublimation of glacier in OptGlacierTreatment=2
           endif

           RadLwNetBareGrd(I,J)     = LwRadCoeff * TemperatureGrdBare(I,J)**4 - EmissivityGrd(I,J) * RadLwDownRefHeight(I,J)
           HeatSensibleBareGrd(I,J) = ShCoeff * (TemperatureGrdBare(I,J) - TemperatureAirRefHeight(I,J))
           HeatLatentBareGrd(I,J)   = LhCoeff * (VapPresSatGrdBare(I,J)*RelHumidityGrd(I,J) - PressureVaporRefHeight(I,J))
           HeatGroundBareGrd(I,J)   = GrdHeatCoeff * (TemperatureGrdBare(I,J) - TemperatureSoilSnow(I,NumSnowLayerNeg(I,J)+1,J))
           EnergyResTmp        = RadSwAbsGrd(I,J) - RadLwNetBareGrd(I,J) - HeatSensibleBareGrd(I,J) - HeatLatentBareGrd(I,J) - &
                                 HeatGroundBareGrd(I,J) + HeatPrecipAdvBareGrd(I,J)
           FluxTotCoeff        = 4.0*LwRadCoeff*TemperatureGrdBare(I,J)**3 + ShCoeff + LhCoeff*VapPresSatGrdBareTempD(I,J) + GrdHeatCoeff
           TemperatureGrdChg   = EnergyResTmp / FluxTotCoeff
           RadLwNetBareGrd(I,J)     = RadLwNetBareGrd(I,J) + 4.0 * LwRadCoeff * TemperatureGrdBare(I,J)**3 * TemperatureGrdChg
           HeatSensibleBareGrd(I,J) = HeatSensibleBareGrd(I,J) + ShCoeff * TemperatureGrdChg
           HeatLatentBareGrd(I,J)   = HeatLatentBareGrd(I,J) + LhCoeff * VapPresSatGrdBareTempD(I,J) * TemperatureGrdChg
           HeatGroundBareGrd(I,J)   = HeatGroundBareGrd(I,J) + GrdHeatCoeff * TemperatureGrdChg
           TemperatureGrdBare(I,J)  = TemperatureGrdBare(I,J) + TemperatureGrdChg

           ! for computing M-O length
           HeatSensibleTmp(I,J) = ShCoeff * (TemperatureGrdBare(I,J) - TemperatureAirRefHeight(I,J))

           ! update specific humidity
           TempTmp = min(50.0, max(-50.0, (TemperatureGrdBare(I,J)-ConstFreezePoint)))
           call VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD)
           if ( TempTmp > 0.0 ) then
              VapPresSatGrdBare(I,J) = VapPresSatWatTmp
           else
              VapPresSatGrdBare(I,J) = VapPresSatIceTmp
           endif
           SpecHumiditySfc(I,J)      = 0.622 * (VapPresSatGrdBare(I,J)*RelHumidityGrd(I,J)) / &
                                  (PressureAirRefHeight(I,J) - 0.378 * (VapPresSatGrdBare(I,J)*RelHumidityGrd(I,J)))
           MoistureFluxSfc      = (SpecHumiditySfc(I,J) - SpecHumidityRefHeight(I,J)) * LhCoeff * PsychConstGrd(I,J) / ConstHeatCapacAir

         end do
       end do
       !$acc end parallel loop

    enddo loop3 ! end stability iteration

    ! Post-iteration parallel region: snow check, wind stresses, 2m temperature
    !$acc parallel loop collapse(2) gang vector default(present) private(LwRadCoeff, ShCoeff, LhCoeff, ExchCoeffShTmp) &
    !$acc private(SoilIceTmp, VapPresSatIceTmp, VapPresSatIceTmpD, VapPresSatWatTmp, VapPresSatWatTmpD) &
    !$acc firstprivate(GrdHeatCoeff, IndIter, MoistureFluxSfc, TempTmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


        ! coefficients needed for snow check
        LwRadCoeff   = EmissivityGrd(I,J) * ConstStefanBoltzmann
        GrdHeatCoeff = 2.0*ThermConductSoilSnow(I,NumSnowLayerNeg(I,J)+1,J)/ThicknessSnowSoilLayer(I,NumSnowLayerNeg(I,J)+1,J)
        ShCoeff      = DensityAirRefHeight(I,J) * ConstHeatCapacAir / ResistanceShBareGrd(I,J)
        LhCoeff      = DensityAirRefHeight(I,J) * ConstHeatCapacAir / PsychConstGrd(I,J) / (ResistanceGrdEvap(I,J)+ResistanceLhBareGrd(I,J))
        if ( .not.((SnowDepth(I,J) > 0.0) .or. (OptGlacierTreatment == 1)) .and. (IndicatorIceSfc(I,J) == -1) ) then
           LhCoeff = 0.0   ! don't allow any sublimation of glacier in OptGlacierTreatment=2
        endif

        ! if snow on ground and TemperatureGrdBare > freezing point: reset TemperatureGrdBare = freezing point. reevaluate ground fluxes.
        if ( (OptSnowSoilTempTime == 1) .or. (OptSnowSoilTempTime == 3) ) then
         if (IndicatorIceSfc(I,J) == 0) then
           if ( (SnowDepth(I,J) > 0.05) .and. (TemperatureGrdBare(I,J) > ConstFreezePoint) ) then
              if ( OptSnowSoilTempTime == 1 ) &
                 TemperatureGrdBare(I,J) = ConstFreezePoint
              if ( OptSnowSoilTempTime == 3 ) &
                 TemperatureGrdBare(I,J) = (1.0-SnowCoverFrac(I,J)) * TemperatureGrdBare(I,J) + SnowCoverFrac(I,J) * ConstFreezePoint  ! MB: allow TemperatureGrd>0C during melt v3.7

              RadLwNetBareGrd(I,J)     = LwRadCoeff * TemperatureGrdBare(I,J)**4 - EmissivityGrd(I,J) * RadLwDownRefHeight(I,J)
              HeatSensibleBareGrd(I,J) = ShCoeff * (TemperatureGrdBare(I,J) - TemperatureAirRefHeight(I,J))
              HeatLatentBareGrd(I,J)   = LhCoeff * (VapPresSatGrdBare(I,J)*RelHumidityGrd(I,J) - PressureVaporRefHeight(I,J))
              HeatGroundBareGrd(I,J)   = RadSwAbsGrd(I,J) + HeatPrecipAdvBareGrd(I,J) - &
                                    (RadLwNetBareGrd(I,J) + HeatSensibleBareGrd(I,J) + HeatLatentBareGrd(I,J))
           endif
         else if (IndicatorIceSfc(I,J) == -1) then ! if a glacier
            SoilIceTmp = 0.0
            !$acc loop seq
            do IndIter = 1, NumSoilLayer
               SoilIceTmp = max(SoilIceTmp, SoilMoisture(I,IndIter,J) - SoilLiqWater(I,IndIter,J))
            enddo
            if ( (SoilIceTmp > 0.0 .or. SnowDepth(I,J) > 0.05) .and. &
                  (TemperatureGrdBare(I,J) > ConstFreezePoint) .and. (OptGlacierTreatment == 1) ) then
               TemperatureGrdBare(I,J)  = ConstFreezePoint
               TempTmp                 = min(50.0, max(-50.0, (TemperatureGrdBare(I,J)-ConstFreezePoint)))
               call VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD)
               VapPresSatGrdBare(I,J)   = VapPresSatIceTmp
               SpecHumiditySfc(I,J)     = 0.622 * (VapPresSatGrdBare(I,J)*RelHumidityGrd(I,J)) / &
                                    (PressureAirRefHeight(I,J) - 0.378 * (VapPresSatGrdBare(I,J)*RelHumidityGrd(I,J)))
               MoistureFluxSfc     = (SpecHumiditySfc(I,J) - SpecHumidityRefHeight(I,J)) * LhCoeff * PsychConstGrd(I,J) / ConstHeatCapacAir
               RadLwNetBareGrd(I,J)     = LwRadCoeff * TemperatureGrdBare(I,J)**4 - EmissivityGrd(I,J) * RadLwDownRefHeight(I,J)
               HeatSensibleBareGrd(I,J) = ShCoeff * (TemperatureGrdBare(I,J) - TemperatureAirRefHeight(I,J))
               HeatLatentBareGrd(I,J)   = LhCoeff * (VapPresSatGrdBare(I,J)*RelHumidityGrd(I,J) - PressureVaporRefHeight(I,J))
               HeatGroundBareGrd(I,J)   = RadSwAbsGrd(I,J) + HeatPrecipAdvBareGrd(I,J) - &
                                    (RadLwNetBareGrd(I,J) + HeatSensibleBareGrd(I,J) + HeatLatentBareGrd(I,J))
            endif
         endif
        endif

        ! wind stresses
        WindStressEwBare(I,J) = -DensityAirRefHeight(I,J) * ExchCoeffMomBare(I,J) * WindSpdRefHeight(I,J) * WindEastwardRefHeight(I,J)
        WindStressNsBare(I,J) = -DensityAirRefHeight(I,J) * ExchCoeffMomBare(I,J) * WindSpdRefHeight(I,J) * WindNorthwardRefHeight(I,J)

        ! 2m air temperature
        if ( (OptSurfaceDrag == 1) .or. (OptSurfaceDrag == 2) ) then
          !ExchCoeffSh2mBare = FrictionVelBare * ConstVonKarman / log((2.0+RoughLenShBareGrd)/RoughLenShBareGrd)
           ExchCoeffSh2mBare(I,J) = FrictionVelBare(I,J) * ConstVonKarman / &
                               (log((2.0+RoughLenShBareGrd(I,J))/RoughLenShBareGrd(I,J)) - MoStabCorrShBare2m(I,J))
           if ( ExchCoeffSh2mBare(I,J) < 1.0e-5 ) then
              TemperatureAir2mBare(I,J) = TemperatureGrdBare(I,J)
              SpecHumidity2mBare(I,J)   = SpecHumiditySfc(I,J)
           else
              TemperatureAir2mBare(I,J) = TemperatureGrdBare(I,J) - HeatSensibleBareGrd(I,J) / &
                                     (DensityAirRefHeight(I,J)*ConstHeatCapacAir) * 1.0 / ExchCoeffSh2mBare(I,J)
              SpecHumidity2mBare(I,J)   = SpecHumiditySfc(I,J) - HeatLatentBareGrd(I,J) /  &
                                     (LatHeatVapGrd(I,J)*DensityAirRefHeight(I,J)) * (1.0/ExchCoeffSh2mBare(I,J) + ResistanceGrdEvap(I,J))
           endif
           if ( FlagUrban(I,J) .eqv. .true. ) SpecHumidity2mBare(I,J) = SpecHumiditySfc(I,J)
        endif

        ! update ExchCoeffShBare
        ExchCoeffShTmp  = 1.0 / ResistanceShBareGrd(I,J)
        ExchCoeffShBare(I,J) = ExchCoeffShTmp

      end do
    end do
    !$acc end parallel loop
    !$acc end data

    end associate

  end subroutine SurfaceEnergyFluxBareGround

end module SurfaceEnergyFluxBareGroundMod
