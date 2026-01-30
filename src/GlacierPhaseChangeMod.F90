module GlacierPhaseChangeMod

!!! Compute the phase change (melting/freezing) of snow and glacier ice

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine GlacierPhaseChange(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: PHASECHANGE_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)     , intent(inout) :: noahmp

! local variable
    integer                               :: I, J                       ! grid indices
    integer                               :: LoopInd1, LoopInd2         ! loop index
    real(kind=kind_noahmp)                :: SnowWaterPrev              ! old/previous snow water equivalent [kg/m2]
    real(kind=kind_noahmp)                :: SnowWaterRatio             ! ratio of previous vs updated snow water equivalent
    real(kind=kind_noahmp)                :: HeatLhTotPhsChg            ! total latent heat of phase change
    logical                               :: FlagAnyAboveFreeze         ! flag for any layer above freezing
    logical                               :: FlagAnyBelowFreeze         ! flag for any layer below freezing
    logical                               :: FlagAnyIce                 ! flag for any layer with ice
    logical                               :: FlagAnyLiq                 ! flag for any layer with liquid
    real(kind=kind_noahmp)                :: EnergyRes      (-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)      ! energy residual [W/m2]
    real(kind=kind_noahmp)                :: GlacierPhaseChg(-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)      ! melting or freezing glacier water [kg/m2]
    real(kind=kind_noahmp)                :: MassWatTotInit (-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)      ! initial total water (ice + liq) mass
    real(kind=kind_noahmp)                :: MassWatIceInit (-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)      ! initial ice content
    real(kind=kind_noahmp)                :: MassWatLiqInit (-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)      ! initial liquid content
    real(kind=kind_noahmp)                :: MassWatIceTmp  (-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)      ! soil/snow ice mass [mm]
    real(kind=kind_noahmp)                :: MassWatLiqTmp  (-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)      ! soil/snow liquid water mass [mm]
    real(kind=kind_noahmp)                :: EnergyResLeft  (-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)      ! energy residual or loss after melting/freezing

! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp) &
    !$acc private(LoopInd1,LoopInd2,SnowWaterPrev,SnowWaterRatio,HeatLhTotPhsChg) &
    !$acc private(FlagAnyAboveFreeze,FlagAnyBelowFreeze,FlagAnyIce,FlagAnyLiq) &
    !$acc private(EnergyRes,GlacierPhaseChg,MassWatTotInit,MassWatIceInit,MassWatLiqInit) &
    !$acc private(MassWatIceTmp,MassWatLiqTmp,EnergyResLeft)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                                            &
              OptGlacierTreatment    => noahmp%config%nmlist%OptGlacierTreatment          ,& ! in,    options for glacier treatment
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer                 ,& ! in,    number of soil layers
              NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax              ,& ! in,    maximum number of snow layers
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg(I,J)         ,& ! in,    actual number of snow layers (negative)
              MainTimeStep           => noahmp%config%domain%MainTimeStep                 ,& ! in,    main noahmp timestep [s]
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer       ,& ! in,    thickness of snow/soil layers [m]
              PhaseChgFacSoilSnow    => noahmp%energy%state%PhaseChgFacSoilSnow           ,& ! in,    energy factor for soil & snow phase change
              TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow           ,& ! inout, snow and soil layer temperature [K]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater                   ,& ! inout, soil water content [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture                   ,& ! inout, total soil moisture [m3/m3]
              SnowIce                => noahmp%water%state%SnowIce                        ,& ! inout, snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater                   ,& ! inout, snow layer liquid water [mm]
              SnowDepth              => noahmp%water%state%SnowDepth(I,J)                 ,& ! inout, snow depth [m]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv(I,J)            ,& ! inout, snow water equivalent [mm]
              IndexPhaseChange       => noahmp%water%state%IndexPhaseChange              ,& ! out,   phase change index [0-none;1-melt;2-refreeze]
              MeltGroundSnow         => noahmp%water%flux%MeltGroundSnow(I,J)             ,& ! out,   ground snowmelt rate [mm/s]
              PondSfcThinSnwMelt     => noahmp%water%state%PondSfcThinSnwMelt(I,J)         & ! out,   surface ponding [mm] from snowmelt when thin snow has no layer
             )
! ----------------------------------------------------------------------

    !--- Initialization
    !$acc loop seq
    do LoopInd1 = -NumSnowLayerMax+1, NumSoilLayer
       EnergyRes(LoopInd1)       = 0.0
       GlacierPhaseChg(LoopInd1) = 0.0
       MassWatTotInit(LoopInd1)  = 0.0
       MassWatIceInit(LoopInd1)  = 0.0
       MassWatLiqInit(LoopInd1)  = 0.0
       MassWatIceTmp(LoopInd1)   = 0.0
       MassWatLiqTmp(LoopInd1)   = 0.0
       EnergyResLeft(LoopInd1)   = 0.0
    enddo
    MeltGroundSnow     = 0.0
    PondSfcThinSnwMelt = 0.0
    HeatLhTotPhsChg    = 0.0

    !--- treat snowpack over glacier ice first

    ! snow layer water mass
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg+1, 0
       MassWatIceTmp(LoopInd1) = SnowIce(I,LoopInd1,J)
       MassWatLiqTmp(LoopInd1) = SnowLiqWater(I,LoopInd1,J)
    enddo

    ! other required variables
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg+1, 0
       IndexPhaseChange(I,LoopInd1,J) = 0
       EnergyRes       (LoopInd1) = 0.0
       GlacierPhaseChg (LoopInd1) = 0.0
       EnergyResLeft   (LoopInd1) = 0.0
       MassWatIceInit  (LoopInd1) = MassWatIceTmp(LoopInd1)
       MassWatLiqInit  (LoopInd1) = MassWatLiqTmp(LoopInd1)
       MassWatTotInit  (LoopInd1) = MassWatIceTmp(LoopInd1) + MassWatLiqTmp(LoopInd1)
    enddo

    ! determine melting or freezing state
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg+1, 0
       if ( (MassWatIceTmp(LoopInd1) > 0.0) .and. (TemperatureSoilSnow(I,LoopInd1,J) >= ConstFreezePoint) ) then
          IndexPhaseChange(I,LoopInd1,J) = 1  ! melting
       endif
       if ( (MassWatLiqTmp(LoopInd1) > 0.0) .and. (TemperatureSoilSnow(I,LoopInd1,J) < ConstFreezePoint) ) then
          IndexPhaseChange(I,LoopInd1,J) = 2  ! freezing
       endif
    enddo

    ! Calculate the energy surplus and loss for melting and freezing
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg+1, 0
       if ( IndexPhaseChange(I,LoopInd1,J) > 0 ) then
          EnergyRes(LoopInd1)           = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
          TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint
       endif
       if ( (IndexPhaseChange(I,LoopInd1,J) == 1) .and. (EnergyRes(LoopInd1) < 0.0) ) then
          EnergyRes(LoopInd1)           = 0.0
          IndexPhaseChange(I,LoopInd1,J)    = 0
       endif
       if ( (IndexPhaseChange(I,LoopInd1,J) == 2) .and. (EnergyRes(LoopInd1) > 0.0) ) then
          EnergyRes(LoopInd1)           = 0.0
          IndexPhaseChange(I,LoopInd1,J)    = 0
       endif
       GlacierPhaseChg(LoopInd1) = EnergyRes(LoopInd1) * MainTimeStep / ConstLatHeatFusion
    enddo

    ! The rate of melting for snow without a layer, needs more work.
    if ( OptGlacierTreatment == 2 ) then
       if ( (NumSnowLayerNeg == 0) .and. (SnowWaterEquiv > 0.0) .and. (TemperatureSoilSnow(I,1,J) > ConstFreezePoint) ) then
          EnergyRes(1)           = (TemperatureSoilSnow(I,1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,1,J)             ! available heat
          TemperatureSoilSnow(I,1,J) = ConstFreezePoint                                                                     ! set T to freezing
          GlacierPhaseChg(1)     = EnergyRes(1) * MainTimeStep / ConstLatHeatFusion                                         ! total snow melt possible
          SnowWaterPrev          = SnowWaterEquiv
          SnowWaterEquiv         = max(0.0, SnowWaterPrev-GlacierPhaseChg(1))                                               ! snow remaining
          SnowWaterRatio         = SnowWaterEquiv / SnowWaterPrev                                                           ! fraction melted
          SnowDepth              = max(0.0, SnowWaterRatio*SnowDepth)                                                       ! new snow height
          SnowDepth              = min(max(SnowDepth,SnowWaterEquiv/500.0), SnowWaterEquiv/50.0)                            ! limit to a reasonable snow density
          EnergyResLeft(1)       = EnergyRes(1) - ConstLatHeatFusion * (SnowWaterPrev - SnowWaterEquiv) / MainTimeStep      ! excess heat
          if ( EnergyResLeft(1) > 0.0 ) then
             GlacierPhaseChg(1)         = EnergyResLeft(1) * MainTimeStep / ConstLatHeatFusion
             TemperatureSoilSnow(I,1,J) = TemperatureSoilSnow(I,1,J) + PhaseChgFacSoilSnow(I,1,J) * EnergyResLeft(1)        ! re-heat ice
          else
             GlacierPhaseChg(1) = 0.0
             EnergyRes(1)       = 0.0
          endif
          MeltGroundSnow     = max(0.0, SnowWaterPrev-SnowWaterEquiv) / MainTimeStep                                        ! melted snow rate
          HeatLhTotPhsChg    = ConstLatHeatFusion * MeltGroundSnow                                                          ! melted snow energy
          PondSfcThinSnwMelt = SnowWaterPrev - SnowWaterEquiv                                                               ! melt water
       endif
    endif ! OptGlacierTreatment==2

    ! The rate of melting and freezing for multi-layer snow
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg+1, 0
       if ( (IndexPhaseChange(I,LoopInd1,J) > 0) .and. (abs(EnergyRes(LoopInd1)) > 0.0) ) then
          EnergyResLeft(LoopInd1)    = 0.0
          if ( GlacierPhaseChg(LoopInd1) > 0.0 ) then
             MassWatIceTmp(LoopInd1) = max(0.0, MassWatIceInit(LoopInd1)-GlacierPhaseChg(LoopInd1))
             EnergyResLeft(LoopInd1) = EnergyRes(LoopInd1) - ConstLatHeatFusion * &
                                       (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / MainTimeStep
          elseif ( GlacierPhaseChg(LoopInd1) < 0.0 ) then
             MassWatIceTmp(LoopInd1) = min(MassWatTotInit(LoopInd1), MassWatIceInit(LoopInd1)-GlacierPhaseChg(LoopInd1))
             EnergyResLeft(LoopInd1) = EnergyRes(LoopInd1) - ConstLatHeatFusion * &
                                       (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / MainTimeStep
          endif
          MassWatLiqTmp(LoopInd1)    = max(0.0, MassWatTotInit(LoopInd1)-MassWatIceTmp(LoopInd1))                           ! update liquid water mass

          ! update snow temperature and energy surplus/loss
          if ( abs(EnergyResLeft(LoopInd1)) > 0.0 ) then
             TemperatureSoilSnow(I,LoopInd1,J) = TemperatureSoilSnow(I,LoopInd1,J) + &
                                                 PhaseChgFacSoilSnow(I,LoopInd1,J) * EnergyResLeft(LoopInd1)
             if ( (MassWatLiqTmp(LoopInd1)*MassWatIceTmp(LoopInd1)) > 0.0 ) &
                TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint
          endif
          HeatLhTotPhsChg = HeatLhTotPhsChg + &
                            ConstLatHeatFusion * (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / MainTimeStep

          ! snow melting rate
          MeltGroundSnow  = MeltGroundSnow + max(0.0, (MassWatIceInit(LoopInd1)-MassWatIceTmp(LoopInd1))) / MainTimeStep
       endif
    enddo

    !---- glacier ice layer treatment

    if ( OptGlacierTreatment == 1 ) then

       ! ice layer water mass
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          MassWatLiqTmp(LoopInd1) = SoilLiqWater(I,LoopInd1,J) * ThicknessSnowSoilLayer(I,LoopInd1,J) * 1000.0
          MassWatIceTmp(LoopInd1) = (SoilMoisture(I,LoopInd1,J) - SoilLiqWater(I,LoopInd1,J)) * ThicknessSnowSoilLayer(I,LoopInd1,J) * 1000.0
       enddo

       ! other required variables
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          IndexPhaseChange(I,LoopInd1,J) = 0
          EnergyRes(LoopInd1)        = 0.0
          GlacierPhaseChg(LoopInd1)  = 0.0
          EnergyResLeft(LoopInd1)    = 0.0
          MassWatIceInit(LoopInd1)   = MassWatIceTmp(LoopInd1)
          MassWatLiqInit(LoopInd1)   = MassWatLiqTmp(LoopInd1)
          MassWatTotInit(LoopInd1)   = MassWatIceTmp(LoopInd1) + MassWatLiqTmp(LoopInd1)
       enddo

       ! determine melting or freezing state
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( (MassWatIceTmp(LoopInd1) > 0.0) .and. (TemperatureSoilSnow(I,LoopInd1,J) >= ConstFreezePoint) ) then
             IndexPhaseChange(I,LoopInd1,J) = 1  ! melting
          endif
          if ( (MassWatLiqTmp(LoopInd1) > 0.0) .and. (TemperatureSoilSnow(I,LoopInd1,J) < ConstFreezePoint) ) then
             IndexPhaseChange(I,LoopInd1,J) = 2  ! freezing
          endif
          ! If snow exists, but its thickness is not enough to create a layer
          if ( (NumSnowLayerNeg == 0) .and. (SnowWaterEquiv > 0.0) .and. (LoopInd1 == 1) ) then
             if ( TemperatureSoilSnow(I,LoopInd1,J) >= ConstFreezePoint ) then
                IndexPhaseChange(I,LoopInd1,J) = 1
             endif
          endif
       enddo

       ! Calculate the energy surplus and loss for melting and freezing
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( IndexPhaseChange(I,LoopInd1,J) > 0 ) then
             EnergyRes(LoopInd1)           = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
             TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint
          endif
          if ( (IndexPhaseChange(I,LoopInd1,J) == 1) .and. (EnergyRes(LoopInd1) < 0.0) ) then
             EnergyRes(LoopInd1)        = 0.0
             IndexPhaseChange(I,LoopInd1,J) = 0
          endif
          if ( (IndexPhaseChange(I,LoopInd1,J) == 2) .and. (EnergyRes(LoopInd1) > 0.0) ) then
             EnergyRes(LoopInd1)        = 0.0
             IndexPhaseChange(I,LoopInd1,J) = 0
          endif
          GlacierPhaseChg(LoopInd1) = EnergyRes(LoopInd1) * MainTimeStep / ConstLatHeatFusion
       enddo

       ! The rate of melting for snow without a layer, needs more work.
       if ( (NumSnowLayerNeg == 0) .and. (SnowWaterEquiv > 0.0) .and. (GlacierPhaseChg(1) > 0.0) ) then
          SnowWaterPrev = SnowWaterEquiv
          SnowWaterEquiv     = max(0.0, SnowWaterPrev-GlacierPhaseChg(1))
          SnowWaterRatio   = SnowWaterEquiv / SnowWaterPrev
          SnowDepth          = max(0.0, SnowWaterRatio*SnowDepth)
          SnowDepth          = min(max(SnowDepth,SnowWaterEquiv/500.0), SnowWaterEquiv/50.0)  ! limit to a reasonable snow density
          EnergyResLeft(1)   = EnergyRes(1) - ConstLatHeatFusion * (SnowWaterPrev - SnowWaterEquiv) / MainTimeStep
          if ( EnergyResLeft(1) > 0.0 ) then
             GlacierPhaseChg(1)  = EnergyResLeft(1) * MainTimeStep / ConstLatHeatFusion
             EnergyRes(1)        = EnergyResLeft(1)
             IndexPhaseChange(I,1,J)    = 1
          else
             GlacierPhaseChg(1)  = 0.0
             EnergyRes(1)        = 0.0
             IndexPhaseChange(I,1,J)    = 0
          endif
          MeltGroundSnow         = max(0.0, (SnowWaterPrev-SnowWaterEquiv)) / MainTimeStep
          HeatLhTotPhsChg        = ConstLatHeatFusion * MeltGroundSnow
          PondSfcThinSnwMelt     = SnowWaterPrev - SnowWaterEquiv
       endif

       ! The rate of melting and freezing for glacier ice
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( (IndexPhaseChange(I,LoopInd1,J) > 0) .and. (abs(EnergyRes(LoopInd1)) > 0.0) ) then
             EnergyResLeft(LoopInd1) = 0.0
             if ( GlacierPhaseChg(LoopInd1) > 0.0 ) then
                MassWatIceTmp(LoopInd1) = max(0.0, MassWatIceInit(LoopInd1)-GlacierPhaseChg(LoopInd1))
                EnergyResLeft(LoopInd1) = EnergyRes(LoopInd1) - ConstLatHeatFusion * &
                                          (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / MainTimeStep
             elseif ( GlacierPhaseChg(LoopInd1) < 0.0 ) then
                MassWatIceTmp(LoopInd1) = min(MassWatTotInit(LoopInd1), MassWatIceInit(LoopInd1)-GlacierPhaseChg(LoopInd1))
                EnergyResLeft(LoopInd1) = EnergyRes(LoopInd1) - ConstLatHeatFusion * &
                                          (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / MainTimeStep
             endif
             MassWatLiqTmp(LoopInd1)    = max(0.0, MassWatTotInit(LoopInd1)-MassWatIceTmp(LoopInd1)) ! update liquid water mass

             ! update ice temperature and energy surplus/loss
             if ( abs(EnergyResLeft(LoopInd1)) > 0.0 ) then
                TemperatureSoilSnow(I,LoopInd1,J) = TemperatureSoilSnow(I,LoopInd1,J) + &
                                                    PhaseChgFacSoilSnow(I,LoopInd1,J) * EnergyResLeft(LoopInd1)
             endif
             HeatLhTotPhsChg = HeatLhTotPhsChg + &
                               ConstLatHeatFusion * (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / MainTimeStep
          endif
       enddo
       !$acc loop seq
       do LoopInd1 = -NumSnowLayerMax+1, NumSoilLayer
          EnergyResLeft(LoopInd1)   = 0.0
          GlacierPhaseChg(LoopInd1) = 0.0
       enddo

       !--- Deal with residuals in ice/soil

       ! Check if any layer above/below freezing (replace any() intrinsic)
       FlagAnyAboveFreeze = .false.
       FlagAnyBelowFreeze = .false.
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( TemperatureSoilSnow(I,LoopInd1,J) > ConstFreezePoint ) FlagAnyAboveFreeze = .true.
          if ( TemperatureSoilSnow(I,LoopInd1,J) < ConstFreezePoint ) FlagAnyBelowFreeze = .true.
       enddo

       ! first remove excess heat by reducing layer temperature
       if ( FlagAnyAboveFreeze .and. FlagAnyBelowFreeze ) then
          !$acc loop seq
          do LoopInd1 = 1, NumSoilLayer
             if ( TemperatureSoilSnow(I,LoopInd1,J) > ConstFreezePoint ) then
                EnergyResLeft(LoopInd1) = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
                !$acc loop seq
                do LoopInd2 = 1, NumSoilLayer
                   if ( (LoopInd1 /= LoopInd2) .and. (TemperatureSoilSnow(I,LoopInd2,J) < ConstFreezePoint) .and. &
                        (EnergyResLeft(LoopInd1) > 0.1) ) then
                      EnergyResLeft(LoopInd2) = (TemperatureSoilSnow(I,LoopInd2,J) - ConstFreezePoint) / &
                                                PhaseChgFacSoilSnow(I,LoopInd2,J)
                      if ( abs(EnergyResLeft(LoopInd2)) > EnergyResLeft(LoopInd1) ) then ! LAYER ABSORBS ALL
                         EnergyResLeft(LoopInd2)       = EnergyResLeft(LoopInd2) + EnergyResLeft(LoopInd1)
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint + &
                                                             EnergyResLeft(LoopInd2) * PhaseChgFacSoilSnow(I,LoopInd2,J)
                         EnergyResLeft(LoopInd1)       = 0.0
                      else
                         EnergyResLeft(LoopInd1)       = EnergyResLeft(LoopInd1) + EnergyResLeft(LoopInd2)
                         EnergyResLeft(LoopInd2)       = 0.0
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                      endif
                   endif
                enddo
                TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint + EnergyResLeft(LoopInd1) * PhaseChgFacSoilSnow(I,LoopInd1,J)
             endif
          enddo
       endif

       ! Re-check flags after temperature adjustment
       FlagAnyAboveFreeze = .false.
       FlagAnyBelowFreeze = .false.
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( TemperatureSoilSnow(I,LoopInd1,J) > ConstFreezePoint ) FlagAnyAboveFreeze = .true.
          if ( TemperatureSoilSnow(I,LoopInd1,J) < ConstFreezePoint ) FlagAnyBelowFreeze = .true.
       enddo

       ! now remove excess cold by increasing temperture (may not be necessary with above loop)
       if ( FlagAnyAboveFreeze .and. FlagAnyBelowFreeze ) then
          !$acc loop seq
          do LoopInd1 = 1, NumSoilLayer
             if ( TemperatureSoilSnow(I,LoopInd1,J) < ConstFreezePoint ) then
                EnergyResLeft(LoopInd1) = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
                !$acc loop seq
                do LoopInd2 = 1, NumSoilLayer
                   if ( (LoopInd1 /= LoopInd2) .and. (TemperatureSoilSnow(I,LoopInd2,J) > ConstFreezePoint) .and. &
                        (EnergyResLeft(LoopInd1) < -0.1) ) then
                      EnergyResLeft(LoopInd2) = (TemperatureSoilSnow(I,LoopInd2,J) - ConstFreezePoint) / &
                                                PhaseChgFacSoilSnow(I,LoopInd2,J)
                      if ( EnergyResLeft(LoopInd2) > abs(EnergyResLeft(LoopInd1)) ) then  ! LAYER ABSORBS ALL
                         EnergyResLeft(LoopInd2)       = EnergyResLeft(LoopInd2) + EnergyResLeft(LoopInd1)
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint + &
                                                             EnergyResLeft(LoopInd2) * PhaseChgFacSoilSnow(I,LoopInd2,J)
                         EnergyResLeft(LoopInd1)       = 0.0
                      else
                         EnergyResLeft(LoopInd1)       = EnergyResLeft(LoopInd1) + EnergyResLeft(LoopInd2)
                         EnergyResLeft(LoopInd2)       = 0.0
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                      endif
                   endif
                enddo
                TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint + EnergyResLeft(LoopInd1) * PhaseChgFacSoilSnow(I,LoopInd1,J)
             endif
          enddo
       endif

       ! Check for ice and temperature flags
       FlagAnyAboveFreeze = .false.
       FlagAnyIce = .false.
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( TemperatureSoilSnow(I,LoopInd1,J) > ConstFreezePoint ) FlagAnyAboveFreeze = .true.
          if ( MassWatIceTmp(LoopInd1) > 0.0 ) FlagAnyIce = .true.
       enddo

       ! now remove excess heat by melting ice
       if ( FlagAnyAboveFreeze .and. FlagAnyIce ) then
          !$acc loop seq
          do LoopInd1 = 1, NumSoilLayer
             if ( TemperatureSoilSnow(I,LoopInd1,J) > ConstFreezePoint ) then
                EnergyResLeft(LoopInd1)   = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
                GlacierPhaseChg(LoopInd1) = EnergyResLeft(LoopInd1) * MainTimeStep / ConstLatHeatFusion
                !$acc loop seq
                do LoopInd2 = 1, NumSoilLayer
                   if ( (LoopInd1 /= LoopInd2) .and. (MassWatIceTmp(LoopInd2) > 0.0) .and. &
                        (GlacierPhaseChg(LoopInd1) > 0.1) ) then
                      if ( MassWatIceTmp(LoopInd2) > GlacierPhaseChg(LoopInd1) ) then  ! LAYER ABSORBS ALL
                         MassWatIceTmp(LoopInd2)       = MassWatIceTmp(LoopInd2) - GlacierPhaseChg(LoopInd1)
                         HeatLhTotPhsChg               = HeatLhTotPhsChg + &
                                                         ConstLatHeatFusion * GlacierPhaseChg(LoopInd1)/MainTimeStep
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                         GlacierPhaseChg(LoopInd1)     = 0.0
                      else
                         GlacierPhaseChg(LoopInd1)     = GlacierPhaseChg(LoopInd1) - MassWatIceTmp(LoopInd2)
                         HeatLhTotPhsChg               = HeatLhTotPhsChg + &
                                                         ConstLatHeatFusion * MassWatIceTmp(LoopInd2) / MainTimeStep
                         MassWatIceTmp(LoopInd2)       = 0.0
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                      endif
                      MassWatLiqTmp(LoopInd2) = max(0.0, MassWatTotInit(LoopInd2)-MassWatIceTmp(LoopInd2))
                   endif
                enddo
                EnergyResLeft(LoopInd1)       = GlacierPhaseChg(LoopInd1) * ConstLatHeatFusion / MainTimeStep
                TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint + EnergyResLeft(LoopInd1) * PhaseChgFacSoilSnow(I,LoopInd1,J)
             endif
          enddo
       endif

       ! Check for liquid and temperature flags
       FlagAnyBelowFreeze = .false.
       FlagAnyLiq = .false.
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( TemperatureSoilSnow(I,LoopInd1,J) < ConstFreezePoint ) FlagAnyBelowFreeze = .true.
          if ( MassWatLiqTmp(LoopInd1) > 0.0 ) FlagAnyLiq = .true.
       enddo

       ! snow remove excess cold by refreezing liquid (may not be necessary with above loop)
       if ( FlagAnyBelowFreeze .and. FlagAnyLiq ) then
          !$acc loop seq
          do LoopInd1 = 1, NumSoilLayer
             if ( TemperatureSoilSnow(I,LoopInd1,J) < ConstFreezePoint ) then
                EnergyResLeft(LoopInd1)   = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
                GlacierPhaseChg(LoopInd1) = EnergyResLeft(LoopInd1) * MainTimeStep / ConstLatHeatFusion
                !$acc loop seq
                do LoopInd2 = 1, NumSoilLayer
                   if ( (LoopInd1 /= LoopInd2) .and. (MassWatLiqTmp(LoopInd2) > 0.0) .and. &
                        (GlacierPhaseChg(LoopInd1) < -0.1) ) then
                      if ( MassWatLiqTmp(LoopInd2) > abs(GlacierPhaseChg(LoopInd1)) ) then  ! LAYER ABSORBS ALL
                         MassWatIceTmp(LoopInd2)       = MassWatIceTmp(LoopInd2) - GlacierPhaseChg(LoopInd1)
                         HeatLhTotPhsChg               = HeatLhTotPhsChg + &
                                                         ConstLatHeatFusion * GlacierPhaseChg(LoopInd1) / MainTimeStep
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                         GlacierPhaseChg(LoopInd1)     = 0.0
                      else
                         GlacierPhaseChg(LoopInd1)     = GlacierPhaseChg(LoopInd1) + MassWatLiqTmp(LoopInd2)
                         HeatLhTotPhsChg               = HeatLhTotPhsChg - &
                                                         ConstLatHeatFusion * MassWatLiqTmp(LoopInd2) / MainTimeStep
                         MassWatIceTmp(LoopInd2)       = MassWatTotInit(LoopInd2)
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                      endif
                      MassWatLiqTmp(LoopInd2) = max(0.0, MassWatTotInit(LoopInd2)-MassWatIceTmp(LoopInd2))
                   endif
                enddo
                EnergyResLeft(LoopInd1)           = GlacierPhaseChg(LoopInd1) * ConstLatHeatFusion / MainTimeStep
                TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint + EnergyResLeft(LoopInd1) * PhaseChgFacSoilSnow(I,LoopInd1,J)
             endif
          enddo
       endif

    endif ! OptGlacierTreatment==1

    !--- update snow and soil ice and liquid content
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg+1, 0     ! snow
       SnowLiqWater(I,LoopInd1,J) = MassWatLiqTmp(LoopInd1)
       SnowIce(I,LoopInd1,J)      = MassWatIceTmp(LoopInd1)
    enddo
    !$acc loop seq
    do LoopInd1 = 1, NumSoilLayer       ! glacier ice
       if ( OptGlacierTreatment == 1 ) then
          SoilLiqWater(I,LoopInd1,J) = MassWatLiqTmp(LoopInd1) / (1000.0 * ThicknessSnowSoilLayer(I,LoopInd1,J))
          SoilLiqWater(I,LoopInd1,J) = max(0.0, min(1.0,SoilLiqWater(I,LoopInd1,J)))
       elseif ( OptGlacierTreatment == 2 ) then
          SoilLiqWater(I,LoopInd1,J) = 0.0             ! ice, assume all frozen forever
       endif
       SoilMoisture(I,LoopInd1,J) = 1.0
    enddo

    end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine GlacierPhaseChange

end module GlacierPhaseChangeMod
