module WaterMainMod

!!! Main water module including all water relevant processes
!!! canopy water -> snowpack water -> soil water -> ground water

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use CanopyHydrologyMod,     only : CanopyHydrology
  use SnowWaterMainMod,       only : SnowWaterMain
  use IrrigationFloodMod,     only : IrrigationFlood
  use IrrigationMicroMod,     only : IrrigationMicro
  use SoilWaterMainMod,       only : SoilWaterMain
  use WetlandWaterZhang22Mod, only : WetlandWaterZhang22

  implicit none

contains

  subroutine WaterMain(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: WATER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd      ! loop index
    integer                          :: I, J         ! grid indices
    real(kind=kind_noahmp)           :: WatReplaceSublim  ! replacement water due to sublimation of glacier
    real(kind=kind_noahmp)  :: SoilIceTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE,1:noahmp%config%domain%NumSoilLayer,noahmp%config%domain%JTS:noahmp%config%domain%JTE)       ! temporary glacier ice content [m3/m3]
    real(kind=kind_noahmp)  :: SoilLiqWaterTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE,1:noahmp%config%domain%NumSoilLayer,noahmp%config%domain%JTS:noahmp%config%domain%JTE)  ! temporary glacier liquid water content [m3/m3]

    associate(                                                                      &
              SoilTimeStep           => noahmp%config%domain%SoilTimeStep ,& ! in,    soil process timestep [s]
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer ,& ! in,    number of soil layers
              SurfaceType            => noahmp%config%domain%SurfaceType ,& ! in,    surface type 1-soil; 2-lake
              FlagCropland           => noahmp%config%domain%FlagCropland ,& ! in,    flag to identify croplands
              FlagWetland            => noahmp%config%domain%FlagWetland ,& ! in,    flag to identify wetlands
              FlagUrban              => noahmp%config%domain%FlagUrban ,& ! in,    urban point flag
              FlagSoilProcess        => noahmp%config%domain%FlagSoilProcess ,& ! in,    flag to calculate soil processes
              NumSoilTimeStep        => noahmp%config%domain%NumSoilTimeStep ,& ! in,    number of timesteps for soil process calculation
              OptWetlandModel        => noahmp%config%nmlist%OptWetlandModel ,& ! in,    options for wetland model
              VaporizeGrd            => noahmp%water%flux%VaporizeGrd ,& ! in,    ground vaporize rate total (evap+sublim) [mm/s]
              CondenseVapGrd         => noahmp%water%flux%CondenseVapGrd ,& ! in,    ground vapor condense rate total (dew+frost) [mm/s]
              RainfallGround         => noahmp%water%flux%RainfallGround ,& ! in,    ground surface rain rate [mm/s]
              SoilTranspFac          => noahmp%water%state%SoilTranspFac ,& ! in,    soil water transpiration factor (0 to 1)
              WaterStorageLakeMax    => noahmp%water%param%WaterStorageLakeMax ,& ! in,    maximum lake water storage [mm]
              NumSoilLayerRoot       => noahmp%water%param%NumSoilLayerRoot ,& ! in,    number of soil layers with root present
              FlagFrozenGround       => noahmp%energy%state%FlagFrozenGround ,& ! in,    frozen ground (logical) to define latent heat pathway
              LatHeatVapGrd          => noahmp%energy%state%LatHeatVapGrd ,& ! in,    latent heat of vaporization/subli [J/kg], ground
              DensityAirRefHeight    => noahmp%energy%state%DensityAirRefHeight ,& ! in,    density air [kg/m3]
              ExchCoeffShSfc         => noahmp%energy%state%ExchCoeffShSfc ,& ! in,    exchange coefficient [m/s] for heat, surface, grid mean
              SpecHumidityRefHeight  => noahmp%forcing%SpecHumidityRefHeight ,& ! in,    specific humidity [kg/kg] at reference height
              HeatLatentGrd          => noahmp%energy%flux%HeatLatentGrd ,& ! in,    total ground latent heat [W/m2] (+ to atm)
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg ,& ! inout, actual number of snow layers (negative)
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv ,& ! inout, snow water equivalent [mm]
              SnowWaterEquivPrev     => noahmp%water%state%SnowWaterEquivPrev ,& ! inout, snow water equivalent at last time step [mm]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater ,& ! inout, soil water content [m3/m3]
              SoilIce                => noahmp%water%state%SoilIce ,& ! inout, soil ice moisture [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture ,& ! inout, total soil moisture [m3/m3]
              WaterStorageLake       => noahmp%water%state%WaterStorageLake ,& ! inout, water storage in lake (can be negative) [mm]
              PondSfcThinSnwMelt     => noahmp%water%state%PondSfcThinSnwMelt ,& ! inout, surface ponding [mm] from snowmelt when thin snow has no layer
              WaterHeadSfc           => noahmp%water%state%WaterHeadSfc ,& ! inout, surface water head (mm)
              IrrigationAmtFlood     => noahmp%water%state%IrrigationAmtFlood ,& ! inout, flood irrigation water amount [m]
              IrrigationAmtMicro     => noahmp%water%state%IrrigationAmtMicro ,& ! inout, micro irrigation water amount [m]
              SoilSfcInflow          => noahmp%water%flux%SoilSfcInflow ,& ! inout, water input on soil surface [m/s]
              EvapSoilSfcLiq         => noahmp%water%flux%EvapSoilSfcLiq ,& ! inout, evaporation from soil surface [m/s]
              DewSoilSfcLiq          => noahmp%water%flux%DewSoilSfcLiq ,& ! inout, soil surface dew rate [mm/s]
              FrostSnowSfcIce        => noahmp%water%flux%FrostSnowSfcIce ,& ! inout, snow surface frost rate[mm/s]
              SublimSnowSfcIce       => noahmp%water%flux%SublimSnowSfcIce ,& ! inout, snow surface sublimation rate[mm/s]
              TranspWatLossSoil      => noahmp%water%flux%TranspWatLossSoil ,& ! inout, transpiration water loss from soil layers [m/s]
              GlacierExcessFlow      => noahmp%water%flux%GlacierExcessFlow ,& ! inout, glacier excess flow [mm/s]
              GlacierExcessFlowAcc   => noahmp%water%flux%GlacierExcessFlowAcc ,& ! inout, accumulated glacier excess flow [mm]
              SoilSfcInflowAcc       => noahmp%water%flux%SoilSfcInflowAcc ,& ! inout, accumulated water flux into soil during soil timestep [m/s * dt_soil/dt_main]
              EvapSoilSfcLiqAcc      => noahmp%water%flux%EvapSoilSfcLiqAcc ,& ! inout, accumulated soil surface evaporation during soil timestep [m/s * dt_soil/dt_main]
              TranspWatLossSoilAcc   => noahmp%water%flux%TranspWatLossSoilAcc ,& ! inout, accumualted transpiration water loss during soil timestep [m/s * dt_soil/dt_main]
              SpecHumidity2mBare     => noahmp%energy%state%SpecHumidity2mBare ,& ! out,   bare ground 2-m specific humidity [kg/kg]
              SpecHumiditySfc        => noahmp%energy%state%SpecHumiditySfc ,& ! out,   specific humidity at surface [kg/kg]
              EvapGroundNet          => noahmp%water%flux%EvapGroundNet ,& ! out,   net ground (soil/snow) evaporation [mm/s]
              Transpiration          => noahmp%water%flux%Transpiration ,& ! out,   transpiration rate [mm/s]
              RunoffSurface          => noahmp%water%flux%RunoffSurface ,& ! out,   surface runoff [mm/dt_soil] per soil timestep
              RunoffSubsurface       => noahmp%water%flux%RunoffSubsurface ,& ! out,   subsurface runoff [mm/dt_soil] per soil timestep
              TileDrain              => noahmp%water%flux%TileDrain ,& ! out,   tile drainage per soil timestep [mm/dt_soil]
              SnowBotOutflow         => noahmp%water%flux%SnowBotOutflow ,& ! out,   total water (snowmelt+rain through pack) out of snow bottom [mm/s]
              WaterToAtmosTotal      => noahmp%water%flux%WaterToAtmosTotal ,& ! out,   total water vapor flux to atmosphere [mm/s]
              SoilSfcInflowMean      => noahmp%water%flux%SoilSfcInflowMean ,& ! out,   mean water flux into soil during soil timestep [m/s]
              TranspWatLossSoilMean  => noahmp%water%flux%TranspWatLossSoilMean ,& ! out,   mean transpiration water loss during soil timestep [m/s]
              PondSfcThinSnwComb     => noahmp%water%state%PondSfcThinSnwComb ,& ! out,   surface ponding [mm] from liquid in thin snow layer combination
              PondSfcThinSnwTrans    => noahmp%water%state%PondSfcThinSnwTrans ,& ! out,   surface ponding [mm] from thin snow liquid during transition from multilayer to no layer
            MainTimeStep           => noahmp%config%domain%MainTimeStep ,& ! in,    noahmp main time step [s]
            IndicatorIceSfc        => noahmp%config%domain%IndicatorIceSfc ,& ! in,    indicator for ice surface (1-ice surface; 0-non ice surface)
            OptGlacierTreatment     => noahmp%config%nmlist%OptGlacierTreatment ,& ! in,    option for glacier treatment (0-no glacier; 1-glacier with excess flow; 2-glacier without excess flow)
            EvapSoilSfcLiqMean     => noahmp%water%flux%EvapSoilSfcLiqMean ,& ! out,   mean soil surface evaporation during soil timestep [m/s]
              EvapCanopyNet        => noahmp%water%flux%EvapCanopyNet  & ! in,    evaporation of intercepted water [mm/s]
             )

      !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
      do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! initialize
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       TranspWatLossSoil(I,LoopInd,J)   = 0.0
       ! prepare for water process
       SoilIce(I,LoopInd,J)         = max(0.0, SoilMoisture(I,LoopInd,J)-SoilLiqWater(I,LoopInd,J))
       SoilIceTmp(I,LoopInd,J)      = SoilIce(I,LoopInd,J)
       SoilLiqWaterTmp(I,LoopInd,J) = SoilLiqWater(I,LoopInd,J)
    enddo
    GlacierExcessFlow(I,J)  = 0.0
    SoilSfcInflow(I,J)      = 0.0
    RunoffSurface(I,J)      = 0.0
    RunoffSubsurface(I,J)   = 0.0
    TileDrain(I,J)          = 0.0

    SnowWaterEquivPrev(I,J) = SnowWaterEquiv(I,J)
    ! compute soil/snow surface evap/dew rate based on energy flux
    VaporizeGrd(I,J)        = max(HeatLatentGrd(I,J)/LatHeatVapGrd(I,J), 0.0)       ! positive part of ground latent heat; Barlage change to ground v3.6
    CondenseVapGrd(I,J)     = abs(min(HeatLatentGrd(I,J)/LatHeatVapGrd(I,J), 0.0))  ! negative part of ground latent heat
    EvapGroundNet(I,J)      = VaporizeGrd(I,J) - CondenseVapGrd(I,J)


   enddo
   enddo

    ! canopy-intercepted snowfall/rainfall, drips, and throughfall
    call CanopyHydrology(noahmp)

   !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
   do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
    ! ground sublimation and evaporation
    SublimSnowSfcIce(I,J)    = 0.0
    if ( SnowWaterEquiv(I,J) > 0.0 ) then
       SublimSnowSfcIce(I,J) = min(VaporizeGrd(I,J), SnowWaterEquiv(I,J)/MainTimeStep)
    endif
    if ( IndicatorIceSfc(I,J) == -1 ) then
      SublimSnowSfcIce(I,J)   = VaporizeGrd(I,J)
    endif
    EvapSoilSfcLiq(I,J)      = VaporizeGrd(I,J) - SublimSnowSfcIce(I,J)

    ! ground frost and dew
    FrostSnowSfcIce(I,J)     = 0.0
    if ( (SnowWaterEquiv(I,J) > 0.0) .or. (IndicatorIceSfc(I,J) == -1) ) then
       FrostSnowSfcIce(I,J)  = CondenseVapGrd(I,J)
    endif
    DewSoilSfcLiq(I,J)       = CondenseVapGrd(I,J) - FrostSnowSfcIce(I,J)

   enddo
enddo


    ! snowpack water processs
    call SnowWaterMain(noahmp)

   !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd) private(WatReplaceSublim)
   do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    if (IndicatorIceSfc(I,J) == 0) then
      ! accumulate glacier excessive flow [mm]
      GlacierExcessFlowAcc(I,J) = GlacierExcessFlowAcc(I,J) + GlacierExcessFlow(I,J) * MainTimeStep

      ! treat frozen ground/soil
      if ( FlagFrozenGround(I,J) .eqv. .true. ) then
         SoilIce(I,1,J)     = SoilIce(I,1,J) + (DewSoilSfcLiq(I,J)-EvapSoilSfcLiq(I,J)) * MainTimeStep / &
                                       (ThicknessSnowSoilLayer(I,1,J)*1000.0)
         DewSoilSfcLiq(I,J)  = 0.0
         EvapSoilSfcLiq(I,J) = 0.0
         if ( SoilIce(I,1,J) < 0.0 ) then
            SoilLiqWater(I,1,J) = SoilLiqWater(I,1,J) + SoilIce(I,1,J)
            SoilIce(I,1,J)      = 0.0
         endif
         SoilMoisture(I,1,J) = SoilLiqWater(I,1,J) + SoilIce(I,1,J)
      endif
      EvapSoilSfcLiq(I,J) = EvapSoilSfcLiq(I,J) * 0.001 ! mm/s -> m/s

      ! transpiration mm/s -> m/s
      !$acc loop seq
      do LoopInd = 1, NumSoilLayerRoot(I,J)
         TranspWatLossSoil(I,LoopInd,J) = Transpiration(I,J) * SoilTranspFac(I,LoopInd,J) * 0.001
      enddo

      ! total surface input water to soil mm/s -> m/s
      SoilSfcInflow(I,J)    = (PondSfcThinSnwMelt(I,J) + PondSfcThinSnwComb(I,J) + PondSfcThinSnwTrans(I,J)) / &
                        MainTimeStep * 0.001  ! convert units (mm/s -> m/s)
      if ( NumSnowLayerNeg(I,J) == 0 ) then
         SoilSfcInflow(I,J) = SoilSfcInflow(I,J) + (SnowBotOutflow(I,J) + DewSoilSfcLiq(I,J) + RainfallGround(I,J)) * 0.001
      else
         SoilSfcInflow(I,J) = SoilSfcInflow(I,J) + (SnowBotOutflow(I,J) + DewSoilSfcLiq(I,J)) * 0.001
      endif

#ifdef WRF_HYDRO
      SoilSfcInflow(I,J)    = SoilSfcInflow(I,J) + WaterHeadSfc(I,J) / MainTimeStep * 0.001
#endif

      ! calculate soil process only at soil timestep
      SoilSfcInflowAcc(I,J)     = SoilSfcInflowAcc(I,J)     + SoilSfcInflow(I,J)
      EvapSoilSfcLiqAcc(I,J)    = EvapSoilSfcLiqAcc(I,J)    + EvapSoilSfcLiq(I,J)
      !$acc loop seq
      do LoopInd = 1, NumSoilLayer
         TranspWatLossSoilAcc(I,LoopInd,J) = TranspWatLossSoilAcc(I,LoopInd,J) + TranspWatLossSoil(I,LoopInd,J)
      enddo

   else if (IndicatorIceSfc(I,J) == -1) then
      ! total surface input water to glacier ice
      SoilSfcInflow(I,J) = (PondSfcThinSnwMelt(I,J) + PondSfcThinSnwComb(I,J) + PondSfcThinSnwTrans(I,J)) / MainTimeStep * 0.001  ! convert units (mm/s -> m/s)
      if ( NumSnowLayerNeg(I,J) == 0 ) then
         SoilSfcInflow(I,J) = SoilSfcInflow(I,J) + (SnowBotOutflow(I,J) + RainfallGround(I,J)) * 0.001
      else
         SoilSfcInflow(I,J) = SoilSfcInflow(I,J) + SnowBotOutflow(I,J) * 0.001
      endif
#ifdef WRF_HYDRO
      SoilSfcInflow(I,J) = SoilSfcInflow(I,J) + WaterHeadSfc(I,J) / MainTimeStep * 0.001
#endif

      ! surface runoff
      RunoffSurface(I,J) = SoilSfcInflow(I,J) * 1000.0   ! mm/s

      ! glacier ice water
      if ( OptGlacierTreatment == 1 ) then
         WatReplaceSublim = 0.0
         !$acc loop seq
         do LoopInd = 1, NumSoilLayer
            WatReplaceSublim = WatReplaceSublim + ThicknessSnowSoilLayer(I,LoopInd,J)*(SoilIce(I,LoopInd,J) - &
                              SoilIceTmp(I,LoopInd,J) + SoilLiqWater(I,LoopInd,J) - SoilLiqWaterTmp(I,LoopInd,J))
         enddo
         WatReplaceSublim = WatReplaceSublim * 1000.0 / MainTimeStep     ! convert to [mm/s]
         !$acc loop seq
         do LoopInd = 1, NumSoilLayer
            SoilIce(I,LoopInd,J) = min(1.0, SoilIceTmp(I,LoopInd,J))
         enddo
      elseif ( OptGlacierTreatment == 2 ) then
         WatReplaceSublim = 0.0
         !$acc loop seq
         do LoopInd = 1, NumSoilLayer
            SoilIce(I,LoopInd,J) = 1.0
         enddo
      endif

      !$acc loop seq
      do LoopInd = 1, NumSoilLayer
         SoilLiqWater(I,LoopInd,J) = 1.0 - SoilIce(I,LoopInd,J)
      enddo

      ! use RunoffSubsurface as a water balancer, GlacierExcessFlow is snow that disappears, WatReplaceSublim is
      ! water from below that replaces glacier loss
      if ( OptGlacierTreatment == 1 ) then
         RunoffSubsurface(I,J) = GlacierExcessFlow(I,J) + WatReplaceSublim
      elseif ( OptGlacierTreatment == 2 ) then
         RunoffSubsurface(I,J) = GlacierExcessFlow(I,J)
         VaporizeGrd(I,J)      = SublimSnowSfcIce(I,J)
         CondenseVapGrd(I,J)   = FrostSnowSfcIce(I,J)
      endif

      if ( OptGlacierTreatment == 2 ) then
         EvapGroundNet(I,J) = VaporizeGrd(I,J) - CondenseVapGrd(I,J)
         HeatLatentGrd(I,J) = EvapGroundNet(I,J) * LatHeatVapGrd(I,J)
      endif

#ifndef _OPENACC
      if ( maxval(SoilIce(I,:,J)) < 0.0001 ) then
         write(*,*) "GLACIER HAS MELTED AT: ", I, J, " ARE YOU SURE THIS SHOULD BE A GLACIER POINT?"
      endif
#endif
   endif ! IndicatorIceSfc(I,J) == -1
   enddo
enddo

    ! start soil water processes
    if ( noahmp%config%domain%FlagSoilProcess .eqv. .true. ) then

       ! irrigation: call flood irrigation and add to SoilSfcInflowAcc
       ! condition if ( (FlagCropland .eqv. .true.) .and. (IrrigationAmtFlood > 0.0) ) moved inside of function
       call IrrigationFlood(noahmp)

       ! irrigation: call micro irrigation assuming we implement drip in first layer
       ! of the Noah-MP. Change layer 1 moisture wrt to MI rate
       call IrrigationMicro(noahmp)
    endif

   !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
   do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points
    ! start soil water processes
    if ( FlagSoilProcess .eqv. .true. ) then
       ! compute mean water flux during soil timestep
       SoilSfcInflowMean(I,J)     = SoilSfcInflowAcc(I,J) / NumSoilTimeStep
       EvapSoilSfcLiqMean(I,J)    = EvapSoilSfcLiqAcc(I,J) / NumSoilTimeStep

       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
         TranspWatLossSoilMean(I,LoopInd,J) = TranspWatLossSoilAcc(I,LoopInd,J) / NumSoilTimeStep
       enddo

       ! lake/soil water balances
       if ( SurfaceType(I,J) == 2 ) then   ! lake
          RunoffSurface(I,J) = 0.0
          if ( WaterStorageLake(I,J) >= WaterStorageLakeMax(I,J) ) RunoffSurface(I,J) = SoilSfcInflowMean(I,J)*1000.0*SoilTimeStep             ! mm per soil timestep
          WaterStorageLake(I,J) = WaterStorageLake(I,J) + (SoilSfcInflowMean(I,J)-EvapSoilSfcLiqMean(I,J))*1000.0*SoilTimeStep - RunoffSurface(I,J) ! mm per soil timestep
       endif


    endif ! FlagSoilProcess soil timestep
   
   enddo
   enddo

    ! soil water processes (including Top model groundwater and shallow water MMF groundwater)
    if (noahmp%config%domain%FlagSoilProcess .eqv. .True.) call SoilWaterMain(noahmp)

   !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
   do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points

    !copied from above in code, moved here for compactness of GPU port
    if (FlagSoilProcess .eqv. .true.) then
      ! merge excess glacier snow flow to subsurface runoff
      RunoffSubsurface(I,J) = RunoffSubsurface(I,J) + GlacierExcessFlowAcc(I,J)  ! mm per soil timestep
    endif

    ! update surface water vapor flux ! urban - jref
    WaterToAtmosTotal(I,J) = Transpiration(I,J) + EvapCanopyNet(I,J) + EvapGroundNet(I,J)
    if ( (FlagUrban(I,J) .eqv. .true.) ) then
       SpecHumiditySfc(I,J)    = WaterToAtmosTotal(I,J) / (DensityAirRefHeight(I,J)*ExchCoeffShSfc(I,J)) + SpecHumidityRefHeight(I,J)
       SpecHumidity2mBare(I,J) = SpecHumiditySfc(I,J)
    endif



   enddo
enddo

    ! call surface wetland scheme (due to subgrid wetland treatment, currently no flag control)
    !if ( (FlagWetland .eqv. .true.) .and. (OptWetlandModel > 0) ) then
    if ( noahmp%config%nmlist%OptWetlandModel > 0 ) then
       call WetlandWaterZhang22(noahmp,noahmp%config%domain%MainTimeStep)
    endif 


    end associate

  end subroutine WaterMain

end module WaterMainMod
