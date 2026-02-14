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

      !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
      do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
! --------------------------------------------------------------------
    associate(                                                                       &
              SoilTimeStep           => noahmp%config%domain%SoilTimeStep           ,& ! in,    soil process timestep [s]
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              SurfaceType            => noahmp%config%domain%SurfaceType(I,J)            ,& ! in,    surface type 1-soil; 2-lake 
              FlagCropland           => noahmp%config%domain%FlagCropland(I,J)           ,& ! in,    flag to identify croplands
              FlagWetland            => noahmp%config%domain%FlagWetland(I,J)            ,& ! in,    flag to identify wetlands
              FlagUrban              => noahmp%config%domain%FlagUrban(I,J)              ,& ! in,    urban point flag
              FlagSoilProcess        => noahmp%config%domain%FlagSoilProcess        ,& ! in,    flag to calculate soil processes
              NumSoilTimeStep        => noahmp%config%domain%NumSoilTimeStep        ,& ! in,    number of timesteps for soil process calculation
              OptWetlandModel        => noahmp%config%nmlist%OptWetlandModel        ,& ! in,    options for wetland model
              VaporizeGrd            => noahmp%water%flux%VaporizeGrd(I,J)               ,& ! in,    ground vaporize rate total (evap+sublim) [mm/s]
              CondenseVapGrd         => noahmp%water%flux%CondenseVapGrd(I,J)            ,& ! in,    ground vapor condense rate total (dew+frost) [mm/s]
              RainfallGround         => noahmp%water%flux%RainfallGround(I,J)            ,& ! in,    ground surface rain rate [mm/s]
              SoilTranspFac          => noahmp%water%state%SoilTranspFac            ,& ! in,    soil water transpiration factor (0 to 1)
              WaterStorageLakeMax    => noahmp%water%param%WaterStorageLakeMax(I,J)      ,& ! in,    maximum lake water storage [mm]
              NumSoilLayerRoot       => noahmp%water%param%NumSoilLayerRoot(I,J)         ,& ! in,    number of soil layers with root present
              FlagFrozenGround       => noahmp%energy%state%FlagFrozenGround(I,J)        ,& ! in,    frozen ground (logical) to define latent heat pathway
              LatHeatVapGrd          => noahmp%energy%state%LatHeatVapGrd(I,J)           ,& ! in,    latent heat of vaporization/subli [J/kg], ground
              DensityAirRefHeight    => noahmp%energy%state%DensityAirRefHeight(I,J)     ,& ! in,    density air [kg/m3]
              ExchCoeffShSfc         => noahmp%energy%state%ExchCoeffShSfc(I,J)          ,& ! in,    exchange coefficient [m/s] for heat, surface, grid mean
              SpecHumidityRefHeight  => noahmp%forcing%SpecHumidityRefHeight(I,J)        ,& ! in,    specific humidity [kg/kg] at reference height
              HeatLatentGrd          => noahmp%energy%flux%HeatLatentGrd(I,J)            ,& ! in,    total ground latent heat [W/m2] (+ to atm)
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg(I,J)        ,& ! inout, actual number of snow layers (negative)
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv(I,J)           ,& ! inout, snow water equivalent [mm]
              SnowWaterEquivPrev     => noahmp%water%state%SnowWaterEquivPrev(I,J)       ,& ! inout, snow water equivalent at last time step [mm]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! inout, soil water content [m3/m3]
              SoilIce                => noahmp%water%state%SoilIce                  ,& ! inout, soil ice moisture [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! inout, total soil moisture [m3/m3]
              WaterStorageLake       => noahmp%water%state%WaterStorageLake(I,J)         ,& ! inout, water storage in lake (can be negative) [mm]
              PondSfcThinSnwMelt     => noahmp%water%state%PondSfcThinSnwMelt(I,J)       ,& ! inout, surface ponding [mm] from snowmelt when thin snow has no layer
              WaterHeadSfc           => noahmp%water%state%WaterHeadSfc(I,J)             ,& ! inout, surface water head (mm) 
              IrrigationAmtFlood     => noahmp%water%state%IrrigationAmtFlood(I,J)       ,& ! inout, flood irrigation water amount [m]
              IrrigationAmtMicro     => noahmp%water%state%IrrigationAmtMicro(I,J)       ,& ! inout, micro irrigation water amount [m]
              SoilSfcInflow          => noahmp%water%flux%SoilSfcInflow(I,J)             ,& ! inout, water input on soil surface [m/s]
              EvapSoilSfcLiq         => noahmp%water%flux%EvapSoilSfcLiq(I,J)            ,& ! inout, evaporation from soil surface [m/s]
              DewSoilSfcLiq          => noahmp%water%flux%DewSoilSfcLiq(I,J)             ,& ! inout, soil surface dew rate [mm/s]
              FrostSnowSfcIce        => noahmp%water%flux%FrostSnowSfcIce(I,J)           ,& ! inout, snow surface frost rate[mm/s]
              SublimSnowSfcIce       => noahmp%water%flux%SublimSnowSfcIce(I,J)          ,& ! inout, snow surface sublimation rate[mm/s]
              TranspWatLossSoil      => noahmp%water%flux%TranspWatLossSoil         ,& ! inout, transpiration water loss from soil layers [m/s]
              GlacierExcessFlow      => noahmp%water%flux%GlacierExcessFlow(I,J)         ,& ! inout, glacier excess flow [mm/s]
              GlacierExcessFlowAcc   => noahmp%water%flux%GlacierExcessFlowAcc(I,J)      ,& ! inout, accumulated glacier excess flow [mm]
              SoilSfcInflowAcc       => noahmp%water%flux%SoilSfcInflowAcc(I,J)          ,& ! inout, accumulated water flux into soil during soil timestep [m/s * dt_soil/dt_main]
              EvapSoilSfcLiqAcc      => noahmp%water%flux%EvapSoilSfcLiqAcc(I,J)         ,& ! inout, accumulated soil surface evaporation during soil timestep [m/s * dt_soil/dt_main]
              TranspWatLossSoilAcc   => noahmp%water%flux%TranspWatLossSoilAcc      ,& ! inout, accumualted transpiration water loss during soil timestep [m/s * dt_soil/dt_main]
              SpecHumidity2mBare     => noahmp%energy%state%SpecHumidity2mBare(I,J)      ,& ! out,   bare ground 2-m specific humidity [kg/kg]
              SpecHumiditySfc        => noahmp%energy%state%SpecHumiditySfc(I,J)         ,& ! out,   specific humidity at surface [kg/kg]
              EvapGroundNet          => noahmp%water%flux%EvapGroundNet(I,J)             ,& ! out,   net ground (soil/snow) evaporation [mm/s]
              Transpiration          => noahmp%water%flux%Transpiration(I,J)             ,& ! out,   transpiration rate [mm/s]
              RunoffSurface          => noahmp%water%flux%RunoffSurface(I,J)             ,& ! out,   surface runoff [mm/dt_soil] per soil timestep
              RunoffSubsurface       => noahmp%water%flux%RunoffSubsurface(I,J)          ,& ! out,   subsurface runoff [mm/dt_soil] per soil timestep
              TileDrain              => noahmp%water%flux%TileDrain(I,J)                 ,& ! out,   tile drainage per soil timestep [mm/dt_soil]
              SnowBotOutflow         => noahmp%water%flux%SnowBotOutflow(I,J)            ,& ! out,   total water (snowmelt+rain through pack) out of snow bottom [mm/s]
              WaterToAtmosTotal      => noahmp%water%flux%WaterToAtmosTotal(I,J)         ,& ! out,   total water vapor flux to atmosphere [mm/s]
              SoilSfcInflowMean      => noahmp%water%flux%SoilSfcInflowMean(I,J)         ,& ! out,   mean water flux into soil during soil timestep [m/s]
              TranspWatLossSoilMean  => noahmp%water%flux%TranspWatLossSoilMean     ,& ! out,   mean transpiration water loss during soil timestep [m/s]
              PondSfcThinSnwComb     => noahmp%water%state%PondSfcThinSnwComb(I,J)       ,& ! out,   surface ponding [mm] from liquid in thin snow layer combination
              PondSfcThinSnwTrans    => noahmp%water%state%PondSfcThinSnwTrans(I,J)       & ! out,   surface ponding [mm] from thin snow liquid during transition from multilayer to no layer
             )
! ----------------------------------------------------------------------

    ! initialize
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       TranspWatLossSoil(I,LoopInd,J)   = 0.0
       ! prepare for water process
       SoilIce(I,LoopInd,J)         = max(0.0, SoilMoisture(I,LoopInd,J)-SoilLiqWater(I,LoopInd,J))
    enddo
    SoilSfcInflow      = 0.0
    RunoffSurface      = 0.0
    RunoffSubsurface   = 0.0
    TileDrain          = 0.0

    SnowWaterEquivPrev = SnowWaterEquiv
    ! compute soil/snow surface evap/dew rate based on energy flux
    VaporizeGrd        = max(HeatLatentGrd/LatHeatVapGrd, 0.0)       ! positive part of ground latent heat; Barlage change to ground v3.6
    CondenseVapGrd     = abs(min(HeatLatentGrd/LatHeatVapGrd, 0.0))  ! negative part of ground latent heat
    EvapGroundNet      = VaporizeGrd - CondenseVapGrd

    end associate

   enddo
   enddo

    ! canopy-intercepted snowfall/rainfall, drips, and throughfall
    call CanopyHydrology(noahmp)

   !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
   do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
! --------------------------------------------------------------------
         associate(                                                                       &
            MainTimeStep           => noahmp%config%domain%MainTimeStep           ,& ! in,    noahmp main time step [s]
            SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv(I,J)           ,& ! inout, snow water equivalent [mm]
            EvapSoilSfcLiq         => noahmp%water%flux%EvapSoilSfcLiq(I,J)            ,& ! inout, evaporation from soil surface [m/s]
            DewSoilSfcLiq          => noahmp%water%flux%DewSoilSfcLiq(I,J)             ,& ! inout, soil surface dew rate [mm/s]
            FrostSnowSfcIce        => noahmp%water%flux%FrostSnowSfcIce(I,J)           ,& ! inout, snow surface frost rate[mm/s]
            VaporizeGrd            => noahmp%water%flux%VaporizeGrd(I,J)               ,& ! in,    ground vaporize rate total (evap+sublim) [mm/s]
            CondenseVapGrd         => noahmp%water%flux%CondenseVapGrd(I,J)            ,& ! in,    ground vapor condense rate total (dew+frost) [mm/s]
            SublimSnowSfcIce       => noahmp%water%flux%SublimSnowSfcIce(I,J)          & ! inout, snow surface sublimation rate[mm/s]
           )
    ! ground sublimation and evaporation
    SublimSnowSfcIce    = 0.0
    if ( SnowWaterEquiv > 0.0 ) then
       SublimSnowSfcIce = min(VaporizeGrd, SnowWaterEquiv/MainTimeStep)
    endif
    EvapSoilSfcLiq      = VaporizeGrd - SublimSnowSfcIce

    ! ground frost and dew
    FrostSnowSfcIce     = 0.0
    if ( SnowWaterEquiv > 0.0 ) then
       FrostSnowSfcIce  = CondenseVapGrd
    endif
    DewSoilSfcLiq       = CondenseVapGrd - FrostSnowSfcIce

    end associate
   enddo
enddo


    ! snowpack water processs
    call SnowWaterMain(noahmp)

   !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
   do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

! --------------------------------------------------------------------
         associate(                                                                       &
            MainTimeStep           => noahmp%config%domain%MainTimeStep           ,& ! in,    noahmp main time step [s]
            NumSoilTimeStep        => noahmp%config%domain%NumSoilTimeStep        ,& ! in,    number of timesteps for soil process calculation
            SoilTimeStep           => noahmp%config%domain%SoilTimeStep           ,& ! in,    soil process timestep [s]
            SurfaceType            => noahmp%config%domain%SurfaceType(I,J)            ,& ! in,    surface type 1-soil; 2-lake 
            FlagUrban            => noahmp%config%domain%FlagUrban(I,J)              ,& ! in,    urban point flag
            NumSoilLayer         => noahmp%config%domain%NumSoilLayer              ,& ! in,    number of soil layers
            FlagCropland         => noahmp%config%domain%FlagCropland(I,J)           ,& ! in,    flag to identify croplands
            FlagFrozenGround     => noahmp%energy%state%FlagFrozenGround(I,J)        ,& ! in,    frozen ground (logical) to define latent heat pathway
            FlagSoilProcess     => noahmp%config%domain%FlagSoilProcess        ,& ! in,    flag to calculate soil processes
            NumSoilLayerRoot     => noahmp%water%param%NumSoilLayerRoot(I,J)         ,& ! in,    number of soil layers with root present
            IrrigationAmtFlood   => noahmp%water%state%IrrigationAmtFlood(I,J)       ,& ! inout, flood irrigation water amount [m]
            IrrigationAmtMicro   => noahmp%water%state%IrrigationAmtMicro(I,J)       ,& ! inout, micro irrigation water amount [m]
            NumSnowLayerNeg      => noahmp%config%domain%NumSnowLayerNeg(I,J)        ,& ! inout, actual number of snow layers (negative)
            PondSfcThinSnwMelt   => noahmp%water%state%PondSfcThinSnwMelt(I,J)       ,& ! inout, surface ponding [mm] from snowmelt when thin snow has no layer
            PondSfcThinSnwComb   => noahmp%water%state%PondSfcThinSnwComb(I,J)       ,& ! out,   surface ponding [mm] from liquid in thin snow layer combination
            PondSfcThinSnwTrans  => noahmp%water%state%PondSfcThinSnwTrans(I,J)       ,& ! out,   surface ponding [mm] from thin snow liquid during transition from multilayer to no layer
            ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m]
            SoilLiqWater         => noahmp%water%state%SoilLiqWater               ,& ! inout, soil water content [m3/m3]
            SoilIce              => noahmp%water%state%SoilIce                  ,& ! inout, soil ice moisture [m3/m3]
            SoilMoisture         => noahmp%water%state%SoilMoisture             ,& ! inout, total soil moisture [m3/m3]
            DewSoilSfcLiq        => noahmp%water%flux%DewSoilSfcLiq(I,J)             ,& ! inout, soil surface dew rate [mm/s]
            EvapSoilSfcLiq       => noahmp%water%flux%EvapSoilSfcLiq(I,J)            ,& ! inout, evaporation from soil surface [m/s]
            Transpiration        => noahmp%water%flux%Transpiration(I,J)             ,& ! out,   transpiration rate [mm/s]
            TranspWatLossSoil    => noahmp%water%flux%TranspWatLossSoil         ,& ! inout, transpiration water loss from soil layers [m/s]
            WaterStorageLake     => noahmp%water%state%WaterStorageLake(I,J)         ,& ! inout, water storage in lake (can be negative) [mm]
            WaterStorageLakeMax  => noahmp%water%param%WaterStorageLakeMax(I,J)      ,& ! in,    maximum lake water storage [mm]
            SoilSfcInflowMean    => noahmp%water%flux%SoilSfcInflowMean(I,J)         ,& ! out,   mean water flux into soil during soil timestep [m/s]
            TranspWatLossSoilMean  => noahmp%water%flux%TranspWatLossSoilMean     ,& ! out,   mean transpiration water loss during soil timestep [m/s]
            SoilTranspFac        => noahmp%water%state%SoilTranspFac            ,& ! in,    soil water transpiration factor (0 to 1)
            GlacierExcessFlow    => noahmp%water%flux%GlacierExcessFlow(I,J)         ,& ! inout, glacier excess flow [mm/s]
            GlacierExcessFlowAcc => noahmp%water%flux%GlacierExcessFlowAcc(I,J)      ,& ! inout, accumulated glacier excess flow [mm]
            SoilSfcInflow        => noahmp%water%flux%SoilSfcInflow(I,J)             ,& ! inout, water input on soil surface [m/s]
            SoilSfcInflowAcc     => noahmp%water%flux%SoilSfcInflowAcc(I,J)          ,& ! inout, accumulated water flux into soil during soil timestep [m/s * dt_soil/dt_main]
            EvapSoilSfcLiqAcc    => noahmp%water%flux%EvapSoilSfcLiqAcc(I,J)         ,& ! inout, accumulated soil surface evaporation during soil timestep [m/s * dt_soil/dt_main]
            TranspWatLossSoilAcc => noahmp%water%flux%TranspWatLossSoilAcc           ,& ! inout, accumualted transpiration water loss during soil timestep [m/s * dt_soil/dt_main]
            RunoffSurface        => noahmp%water%flux%RunoffSurface(I,J)             ,& ! out,   surface runoff [mm/dt_soil] per soil timestep
            SnowBotOutflow       => noahmp%water%flux%SnowBotOutflow(I,J)            ,& ! out,   total water (snowmelt+rain through pack) out of snow bottom [mm/s]
            EvapSoilSfcLiqMean     => noahmp%water%flux%EvapSoilSfcLiqMean(I,J)        ,& ! out,   mean soil surface evaporation during soil timestep [m/s]
            RainfallGround       => noahmp%water%flux%RainfallGround(I,J)             & ! in,    ground surface rain rate [mm/s]
           )
    ! accumulate glacier excessive flow [mm]
    GlacierExcessFlowAcc = GlacierExcessFlowAcc + GlacierExcessFlow * MainTimeStep

    ! treat frozen ground/soil
    if ( FlagFrozenGround .eqv. .true. ) then
       SoilIce(I,1,J)     = SoilIce(I,1,J) + (DewSoilSfcLiq-EvapSoilSfcLiq) * MainTimeStep / &
                                     (ThicknessSnowSoilLayer(I,1,J)*1000.0)
       DewSoilSfcLiq  = 0.0
       EvapSoilSfcLiq = 0.0
       if ( SoilIce(I,1,J) < 0.0 ) then
          SoilLiqWater(I,1,J) = SoilLiqWater(I,1,J) + SoilIce(I,1,J)
          SoilIce(I,1,J)      = 0.0
       endif
       SoilMoisture(I,1,J) = SoilLiqWater(I,1,J) + SoilIce(I,1,J)
    endif
    EvapSoilSfcLiq = EvapSoilSfcLiq * 0.001 ! mm/s -> m/s

    ! transpiration mm/s -> m/s
    do LoopInd = 1, NumSoilLayerRoot
       TranspWatLossSoil(I,LoopInd,J) = Transpiration * SoilTranspFac(I,LoopInd,J) * 0.001
    enddo

    ! total surface input water to soil mm/s -> m/s
    SoilSfcInflow    = (PondSfcThinSnwMelt + PondSfcThinSnwComb + PondSfcThinSnwTrans) / &
                       MainTimeStep * 0.001  ! convert units (mm/s -> m/s)
    if ( NumSnowLayerNeg == 0 ) then
       SoilSfcInflow = SoilSfcInflow + (SnowBotOutflow + DewSoilSfcLiq + RainfallGround) * 0.001
    else
       SoilSfcInflow = SoilSfcInflow + (SnowBotOutflow + DewSoilSfcLiq) * 0.001
    endif

#ifdef WRF_HYDRO
    SoilSfcInflow    = SoilSfcInflow + WaterHeadSfc / MainTimeStep * 0.001
#endif

    ! calculate soil process only at soil timestep
    SoilSfcInflowAcc     = SoilSfcInflowAcc     + SoilSfcInflow
    EvapSoilSfcLiqAcc    = EvapSoilSfcLiqAcc    + EvapSoilSfcLiq
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
      TranspWatLossSoilAcc(I,LoopInd,J) = TranspWatLossSoilAcc(I,LoopInd,J) + TranspWatLossSoil(I,LoopInd,J)
    enddo
    end associate
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

   !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
   do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         associate(                                                                       &
            NumSoilTimeStep        => noahmp%config%domain%NumSoilTimeStep        ,& ! in,    number of timesteps for soil process calculation
            SoilTimeStep           => noahmp%config%domain%SoilTimeStep           ,& ! in,    soil process timestep [s]
            FlagSoilProcess     => noahmp%config%domain%FlagSoilProcess        ,& ! in,    flag to calculate soil processes
            SurfaceType            => noahmp%config%domain%SurfaceType(I,J)            ,& ! in,    surface type 1-soil; 2-lake 
            NumSoilLayer         => noahmp%config%domain%NumSoilLayer              ,& ! in,    number of soil layers
            WaterStorageLake     => noahmp%water%state%WaterStorageLake(I,J)         ,& ! inout, water storage in lake (can be negative) [mm]
            WaterStorageLakeMax  => noahmp%water%param%WaterStorageLakeMax(I,J)      ,& ! in,    maximum lake water storage [mm]
            SoilSfcInflowMean    => noahmp%water%flux%SoilSfcInflowMean(I,J)         ,& ! out,   mean water flux into soil during soil timestep [m/s]
            TranspWatLossSoilMean  => noahmp%water%flux%TranspWatLossSoilMean     ,& ! out,   mean transpiration water loss during soil timestep [m/s]
            SoilSfcInflowAcc     => noahmp%water%flux%SoilSfcInflowAcc(I,J)          ,& ! inout, accumulated water flux into soil during soil timestep [m/s * dt_soil/dt_main]
            EvapSoilSfcLiqAcc    => noahmp%water%flux%EvapSoilSfcLiqAcc(I,J)         ,& ! inout, accumulated soil surface evaporation during soil timestep [m/s * dt_soil/dt_main]
            TranspWatLossSoilAcc => noahmp%water%flux%TranspWatLossSoilAcc           ,& ! inout, accumualted transpiration water loss during soil timestep [m/s * dt_soil/dt_main]
            RunoffSurface        => noahmp%water%flux%RunoffSurface(I,J)             ,& ! out,   surface runoff [mm/dt_soil] per soil timestep
            EvapSoilSfcLiqMean     => noahmp%water%flux%EvapSoilSfcLiqMean(I,J)       & ! out,   mean soil surface evaporation during soil timestep [m/s]
           )
    ! start soil water processes
    if ( FlagSoilProcess .eqv. .true. ) then
       ! compute mean water flux during soil timestep
       SoilSfcInflowMean     = SoilSfcInflowAcc / NumSoilTimeStep
       EvapSoilSfcLiqMean    = EvapSoilSfcLiqAcc / NumSoilTimeStep

       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
         TranspWatLossSoilMean(I,LoopInd,J) = TranspWatLossSoilAcc(I,LoopInd,J) / NumSoilTimeStep
       enddo

       ! lake/soil water balances
       if ( SurfaceType == 2 ) then   ! lake
          RunoffSurface = 0.0
          if ( WaterStorageLake >= WaterStorageLakeMax ) RunoffSurface = SoilSfcInflowMean*1000.0*SoilTimeStep             ! mm per soil timestep
          WaterStorageLake = WaterStorageLake + (SoilSfcInflowMean-EvapSoilSfcLiqMean)*1000.0*SoilTimeStep - RunoffSurface ! mm per soil timestep
       endif


    endif ! FlagSoilProcess soil timestep
   
    end associate
   enddo
   enddo

    ! soil water processes (including Top model groundwater and shallow water MMF groundwater)
    if (noahmp%config%domain%FlagSoilProcess .eqv. .True.) call SoilWaterMain(noahmp)

   !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
   do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

! --------------------------------------------------------------------
    associate(                                                                       &
              FlagSoilProcess      => noahmp%config%domain%FlagSoilProcess        ,& ! in,    flag to calculate soil processes
              RunoffSubsurface     => noahmp%water%flux%RunoffSubsurface(I,J)          ,& ! inout,   subsurface runoff [mm/dt_soil] per soil timestep
              GlacierExcessFlowAcc  => noahmp%water%flux%GlacierExcessFlowAcc(I,J)     ,& ! inout,   accumulated glacier excess flow [mm]
              FlagUrban            => noahmp%config%domain%FlagUrban(I,J)              ,& ! in,    urban point flag
              Transpiration        => noahmp%water%flux%Transpiration(I,J)             ,& ! in,    transpiration rate [mm/s]
              EvapCanopyNet        => noahmp%water%flux%EvapCanopyNet(I,J)             ,& ! in,    evaporation of intercepted water [mm/s]
              EvapGroundNet      => noahmp%water%flux%EvapGroundNet(I,J)             ,& ! in,    net ground (soil/snow) evaporation [mm/s]
              DensityAirRefHeight  => noahmp%energy%state%DensityAirRefHeight(I,J)     ,& ! in,    density air [kg/m3]
              ExchCoeffShSfc       => noahmp%energy%state%ExchCoeffShSfc(I,J)          ,& ! in,    exchange coefficient [m/s] for heat, surface, grid mean
              SpecHumidityRefHeight=> noahmp%forcing%SpecHumidityRefHeight(I,J)        ,& ! in,    specific humidity [kg/kg] at reference height
              SpecHumiditySfc      => noahmp%energy%state%SpecHumiditySfc(I,J)         ,& ! out,   specific humidity at surface [kg/kg]
              WaterToAtmosTotal    => noahmp%water%flux%WaterToAtmosTotal(I,J)         ,& ! out,   total water vapor flux to atmosphere [mm/s]
              SpecHumidity2mBare   => noahmp%energy%state%SpecHumidity2mBare(I,J)      & ! out,   bare ground 2-m specific humidity [kg/kg]
             )
    !copied from above in code, moved here for compactness of GPU port
    if (FlagSoilProcess .eqv. .true.) then
      ! merge excess glacier snow flow to subsurface runoff
      RunoffSubsurface = RunoffSubsurface + GlacierExcessFlowAcc  ! mm per soil timestep
    endif

    ! update surface water vapor flux ! urban - jref
    WaterToAtmosTotal = Transpiration + EvapCanopyNet + EvapGroundNet
    if ( (FlagUrban .eqv. .true.) ) then
       SpecHumiditySfc    = WaterToAtmosTotal / (DensityAirRefHeight*ExchCoeffShSfc) + SpecHumidityRefHeight
       SpecHumidity2mBare = SpecHumiditySfc
    endif


    end associate

   enddo
enddo

    ! call surface wetland scheme (due to subgrid wetland treatment, currently no flag control)
    !if ( (FlagWetland .eqv. .true.) .and. (OptWetlandModel > 0) ) then
    if ( noahmp%config%nmlist%OptWetlandModel > 0 ) then
       call WetlandWaterZhang22(noahmp,noahmp%config%domain%MainTimeStep)
    endif 

  end subroutine WaterMain

end module WaterMainMod
