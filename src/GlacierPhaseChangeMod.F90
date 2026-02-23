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
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: EnergyRes          ! energy residual [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: GlacierPhaseChg    ! melting or freezing glacier water [kg/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassWatTotInit     ! initial total water (ice + liq) mass
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassWatIceInit     ! initial ice content
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassWatLiqInit     ! initial liquid content
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassWatIceTmp      ! soil/snow ice mass [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassWatLiqTmp      ! soil/snow liquid water mass [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: EnergyResLeft      ! energy residual or loss after melting/freezing

! --------------------------------------------------------------------
    associate(                                                                            &
              OptGlacierTreatment    => noahmp%config%nmlist%OptGlacierTreatment          ,& ! in,    options for glacier treatment
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer                 ,& ! in,    number of soil layers
              NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax              ,& ! in,    maximum number of snow layers
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg         ,& ! in,    actual number of snow layers (negative)
              MainTimeStep           => noahmp%config%domain%MainTimeStep                 ,& ! in,    main noahmp timestep [s]
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer       ,& ! in,    thickness of snow/soil layers [m]
              PhaseChgFacSoilSnow    => noahmp%energy%state%PhaseChgFacSoilSnow           ,& ! in,    energy factor for soil & snow phase change
              TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow           ,& ! inout, snow and soil layer temperature [K]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater                   ,& ! inout, soil water content [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture                   ,& ! inout, total soil moisture [m3/m3]
              SnowIce                => noahmp%water%state%SnowIce                        ,& ! inout, snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater                   ,& ! inout, snow layer liquid water [mm]
              SnowDepth              => noahmp%water%state%SnowDepth                 ,& ! inout, snow depth [m]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv            ,& ! inout, snow water equivalent [mm]
              IndexPhaseChange       => noahmp%water%state%IndexPhaseChange              ,& ! out,   phase change index [0-none;1-melt;2-refreeze]
              MeltGroundSnow         => noahmp%water%flux%MeltGroundSnow             ,& ! out,   ground snowmelt rate [mm/s]
              PondSfcThinSnwMelt     => noahmp%water%state%PondSfcThinSnwMelt         & ! out,   surface ponding [mm] from snowmelt when thin snow has no layer
             )

    allocate(EnergyRes      (noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(GlacierPhaseChg(noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassWatTotInit (noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassWatIceInit (noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassWatLiqInit (noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassWatIceTmp  (noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassWatLiqTmp  (noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(EnergyResLeft  (noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    !$acc data create(EnergyRes, GlacierPhaseChg, MassWatTotInit, MassWatIceInit, MassWatLiqInit, MassWatIceTmp, MassWatLiqTmp, EnergyResLeft)

    !$acc parallel loop collapse(2) gang vector default(present) &
    !$acc private(LoopInd1,LoopInd2,SnowWaterPrev,SnowWaterRatio,HeatLhTotPhsChg) &
    !$acc private(FlagAnyAboveFreeze,FlagAnyBelowFreeze,FlagAnyIce,FlagAnyLiq)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         if (noahmp%config%domain%IndicatorIceSfc(I,J) /= -1) cycle
         

    !--- Initialization
    !$acc loop seq
    do LoopInd1 = -NumSnowLayerMax+1, NumSoilLayer
       EnergyRes(I,LoopInd1,J)       = 0.0
       GlacierPhaseChg(I,LoopInd1,J) = 0.0
       MassWatTotInit(I,LoopInd1,J)  = 0.0
       MassWatIceInit(I,LoopInd1,J)  = 0.0
       MassWatLiqInit(I,LoopInd1,J)  = 0.0
       MassWatIceTmp(I,LoopInd1,J)   = 0.0
       MassWatLiqTmp(I,LoopInd1,J)   = 0.0
       EnergyResLeft(I,LoopInd1,J)   = 0.0
    enddo
    MeltGroundSnow(I,J)     = 0.0
    PondSfcThinSnwMelt(I,J) = 0.0
    HeatLhTotPhsChg    = 0.0

    !--- treat snowpack over glacier ice first

    ! snow layer water mass
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg(I,J)+1, 0
       MassWatIceTmp(I,LoopInd1,J) = SnowIce(I,LoopInd1,J)
       MassWatLiqTmp(I,LoopInd1,J) = SnowLiqWater(I,LoopInd1,J)
    enddo

    ! other required variables
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg(I,J)+1, 0
       IndexPhaseChange(I,LoopInd1,J) = 0
       EnergyRes       (I,LoopInd1,J) = 0.0
       GlacierPhaseChg (I,LoopInd1,J) = 0.0
       EnergyResLeft   (I,LoopInd1,J) = 0.0
       MassWatIceInit  (I,LoopInd1,J) = MassWatIceTmp(I,LoopInd1,J)
       MassWatLiqInit  (I,LoopInd1,J) = MassWatLiqTmp(I,LoopInd1,J)
       MassWatTotInit  (I,LoopInd1,J) = MassWatIceTmp(I,LoopInd1,J) + MassWatLiqTmp(I,LoopInd1,J)
    enddo

    ! determine melting or freezing state
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg(I,J)+1, 0
       if ( (MassWatIceTmp(I,LoopInd1,J) > 0.0) .and. (TemperatureSoilSnow(I,LoopInd1,J) >= ConstFreezePoint) ) then
          IndexPhaseChange(I,LoopInd1,J) = 1  ! melting
       endif
       if ( (MassWatLiqTmp(I,LoopInd1,J) > 0.0) .and. (TemperatureSoilSnow(I,LoopInd1,J) < ConstFreezePoint) ) then
          IndexPhaseChange(I,LoopInd1,J) = 2  ! freezing
       endif
    enddo

    ! Calculate the energy surplus and loss for melting and freezing
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg(I,J)+1, 0
       if ( IndexPhaseChange(I,LoopInd1,J) > 0 ) then
          EnergyRes(I,LoopInd1,J)           = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
          TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint
       endif
       if ( (IndexPhaseChange(I,LoopInd1,J) == 1) .and. (EnergyRes(I,LoopInd1,J) < 0.0) ) then
          EnergyRes(I,LoopInd1,J)           = 0.0
          IndexPhaseChange(I,LoopInd1,J)    = 0
       endif
       if ( (IndexPhaseChange(I,LoopInd1,J) == 2) .and. (EnergyRes(I,LoopInd1,J) > 0.0) ) then
          EnergyRes(I,LoopInd1,J)           = 0.0
          IndexPhaseChange(I,LoopInd1,J)    = 0
       endif
       GlacierPhaseChg(I,LoopInd1,J) = EnergyRes(I,LoopInd1,J) * MainTimeStep / ConstLatHeatFusion
    enddo

    ! The rate of melting for snow without a layer, needs more work.
    if ( OptGlacierTreatment == 2 ) then
       if ( (NumSnowLayerNeg(I,J) == 0) .and. (SnowWaterEquiv(I,J) > 0.0) .and. (TemperatureSoilSnow(I,1,J) > ConstFreezePoint) ) then
          EnergyRes(I,1,J)           = (TemperatureSoilSnow(I,1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,1,J)             ! available heat
          TemperatureSoilSnow(I,1,J) = ConstFreezePoint                                                                     ! set T to freezing
          GlacierPhaseChg(I,1,J)     = EnergyRes(I,1,J) * MainTimeStep / ConstLatHeatFusion                                         ! total snow melt possible
          SnowWaterPrev          = SnowWaterEquiv(I,J)
          SnowWaterEquiv(I,J)         = max(0.0, SnowWaterPrev-GlacierPhaseChg(I,1,J))                                               ! snow remaining
          SnowWaterRatio         = SnowWaterEquiv(I,J) / SnowWaterPrev                                                           ! fraction melted
          SnowDepth(I,J)              = max(0.0, SnowWaterRatio*SnowDepth(I,J))                                                       ! new snow height
          SnowDepth(I,J)              = min(max(SnowDepth(I,J),SnowWaterEquiv(I,J)/500.0), SnowWaterEquiv(I,J)/50.0)                            ! limit to a reasonable snow density
          EnergyResLeft(I,1,J)       = EnergyRes(I,1,J) - ConstLatHeatFusion * (SnowWaterPrev - SnowWaterEquiv(I,J)) / MainTimeStep      ! excess heat
          if ( EnergyResLeft(I,1,J) > 0.0 ) then
             GlacierPhaseChg(I,1,J)         = EnergyResLeft(I,1,J) * MainTimeStep / ConstLatHeatFusion
             TemperatureSoilSnow(I,1,J) = TemperatureSoilSnow(I,1,J) + PhaseChgFacSoilSnow(I,1,J) * EnergyResLeft(I,1,J)        ! re-heat ice
          else
             GlacierPhaseChg(I,1,J) = 0.0
             EnergyRes(I,1,J)       = 0.0
          endif
          MeltGroundSnow(I,J)     = max(0.0, SnowWaterPrev-SnowWaterEquiv(I,J)) / MainTimeStep                                        ! melted snow rate
          HeatLhTotPhsChg    = ConstLatHeatFusion * MeltGroundSnow(I,J)                                                          ! melted snow energy
          PondSfcThinSnwMelt(I,J) = SnowWaterPrev - SnowWaterEquiv(I,J)                                                               ! melt water
       endif
    endif ! OptGlacierTreatment==2

    ! The rate of melting and freezing for multi-layer snow
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg(I,J)+1, 0
       if ( (IndexPhaseChange(I,LoopInd1,J) > 0) .and. (abs(EnergyRes(I,LoopInd1,J)) > 0.0) ) then
          EnergyResLeft(I,LoopInd1,J)    = 0.0
          if ( GlacierPhaseChg(I,LoopInd1,J) > 0.0 ) then
             MassWatIceTmp(I,LoopInd1,J) = max(0.0, MassWatIceInit(I,LoopInd1,J)-GlacierPhaseChg(I,LoopInd1,J))
             EnergyResLeft(I,LoopInd1,J) = EnergyRes(I,LoopInd1,J) - ConstLatHeatFusion * &
                                       (MassWatIceInit(I,LoopInd1,J) - MassWatIceTmp(I,LoopInd1,J)) / MainTimeStep
          elseif ( GlacierPhaseChg(I,LoopInd1,J) < 0.0 ) then
             MassWatIceTmp(I,LoopInd1,J) = min(MassWatTotInit(I,LoopInd1,J), MassWatIceInit(I,LoopInd1,J)-GlacierPhaseChg(I,LoopInd1,J))
             EnergyResLeft(I,LoopInd1,J) = EnergyRes(I,LoopInd1,J) - ConstLatHeatFusion * &
                                       (MassWatIceInit(I,LoopInd1,J) - MassWatIceTmp(I,LoopInd1,J)) / MainTimeStep
          endif
          MassWatLiqTmp(I,LoopInd1,J)    = max(0.0, MassWatTotInit(I,LoopInd1,J)-MassWatIceTmp(I,LoopInd1,J))                           ! update liquid water mass

          ! update snow temperature and energy surplus/loss
          if ( abs(EnergyResLeft(I,LoopInd1,J)) > 0.0 ) then
             TemperatureSoilSnow(I,LoopInd1,J) = TemperatureSoilSnow(I,LoopInd1,J) + &
                                                 PhaseChgFacSoilSnow(I,LoopInd1,J) * EnergyResLeft(I,LoopInd1,J)
             if ( (MassWatLiqTmp(I,LoopInd1,J)*MassWatIceTmp(I,LoopInd1,J)) > 0.0 ) &
                TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint
          endif
          HeatLhTotPhsChg = HeatLhTotPhsChg + &
                            ConstLatHeatFusion * (MassWatIceInit(I,LoopInd1,J) - MassWatIceTmp(I,LoopInd1,J)) / MainTimeStep

          ! snow melting rate
          MeltGroundSnow(I,J)  = MeltGroundSnow(I,J) + max(0.0, (MassWatIceInit(I,LoopInd1,J)-MassWatIceTmp(I,LoopInd1,J))) / MainTimeStep
       endif
    enddo

    !---- glacier ice layer treatment

    if ( OptGlacierTreatment == 1 ) then

       ! ice layer water mass
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          MassWatLiqTmp(I,LoopInd1,J) = SoilLiqWater(I,LoopInd1,J) * ThicknessSnowSoilLayer(I,LoopInd1,J) * 1000.0
          MassWatIceTmp(I,LoopInd1,J) = (SoilMoisture(I,LoopInd1,J) - SoilLiqWater(I,LoopInd1,J)) * ThicknessSnowSoilLayer(I,LoopInd1,J) * 1000.0
       enddo

       ! other required variables
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          IndexPhaseChange(I,LoopInd1,J) = 0
          EnergyRes(I,LoopInd1,J)        = 0.0
          GlacierPhaseChg(I,LoopInd1,J)  = 0.0
          EnergyResLeft(I,LoopInd1,J)    = 0.0
          MassWatIceInit(I,LoopInd1,J)   = MassWatIceTmp(I,LoopInd1,J)
          MassWatLiqInit(I,LoopInd1,J)   = MassWatLiqTmp(I,LoopInd1,J)
          MassWatTotInit(I,LoopInd1,J)   = MassWatIceTmp(I,LoopInd1,J) + MassWatLiqTmp(I,LoopInd1,J)
       enddo

       ! determine melting or freezing state
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( (MassWatIceTmp(I,LoopInd1,J) > 0.0) .and. (TemperatureSoilSnow(I,LoopInd1,J) >= ConstFreezePoint) ) then
             IndexPhaseChange(I,LoopInd1,J) = 1  ! melting
          endif
          if ( (MassWatLiqTmp(I,LoopInd1,J) > 0.0) .and. (TemperatureSoilSnow(I,LoopInd1,J) < ConstFreezePoint) ) then
             IndexPhaseChange(I,LoopInd1,J) = 2  ! freezing
          endif
          ! If snow exists, but its thickness is not enough to create a layer
          if ( (NumSnowLayerNeg(I,J) == 0) .and. (SnowWaterEquiv(I,J) > 0.0) .and. (LoopInd1 == 1) ) then
             if ( TemperatureSoilSnow(I,LoopInd1,J) >= ConstFreezePoint ) then
                IndexPhaseChange(I,LoopInd1,J) = 1
             endif
          endif
       enddo

       ! Calculate the energy surplus and loss for melting and freezing
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( IndexPhaseChange(I,LoopInd1,J) > 0 ) then
             EnergyRes(I,LoopInd1,J)           = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
             TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint
          endif
          if ( (IndexPhaseChange(I,LoopInd1,J) == 1) .and. (EnergyRes(I,LoopInd1,J) < 0.0) ) then
             EnergyRes(I,LoopInd1,J)        = 0.0
             IndexPhaseChange(I,LoopInd1,J) = 0
          endif
          if ( (IndexPhaseChange(I,LoopInd1,J) == 2) .and. (EnergyRes(I,LoopInd1,J) > 0.0) ) then
             EnergyRes(I,LoopInd1,J)        = 0.0
             IndexPhaseChange(I,LoopInd1,J) = 0
          endif
          GlacierPhaseChg(I,LoopInd1,J) = EnergyRes(I,LoopInd1,J) * MainTimeStep / ConstLatHeatFusion
       enddo

       ! The rate of melting for snow without a layer, needs more work.
       if ( (NumSnowLayerNeg(I,J) == 0) .and. (SnowWaterEquiv(I,J) > 0.0) .and. (GlacierPhaseChg(I,1,J) > 0.0) ) then
          SnowWaterPrev = SnowWaterEquiv(I,J)
          SnowWaterEquiv(I,J)     = max(0.0, SnowWaterPrev-GlacierPhaseChg(I,1,J))
          SnowWaterRatio   = SnowWaterEquiv(I,J) / SnowWaterPrev
          SnowDepth(I,J)          = max(0.0, SnowWaterRatio*SnowDepth(I,J))
          SnowDepth(I,J)          = min(max(SnowDepth(I,J),SnowWaterEquiv(I,J)/500.0), SnowWaterEquiv(I,J)/50.0)  ! limit to a reasonable snow density
          EnergyResLeft(I,1,J)   = EnergyRes(I,1,J) - ConstLatHeatFusion * (SnowWaterPrev - SnowWaterEquiv(I,J)) / MainTimeStep
          if ( EnergyResLeft(I,1,J) > 0.0 ) then
             GlacierPhaseChg(I,1,J)  = EnergyResLeft(I,1,J) * MainTimeStep / ConstLatHeatFusion
             EnergyRes(I,1,J)        = EnergyResLeft(I,1,J)
             IndexPhaseChange(I,1,J)    = 1
          else
             GlacierPhaseChg(I,1,J)  = 0.0
             EnergyRes(I,1,J)        = 0.0
             IndexPhaseChange(I,1,J)    = 0
          endif
          MeltGroundSnow(I,J)         = max(0.0, (SnowWaterPrev-SnowWaterEquiv(I,J))) / MainTimeStep
          HeatLhTotPhsChg        = ConstLatHeatFusion * MeltGroundSnow(I,J)
          PondSfcThinSnwMelt(I,J)     = SnowWaterPrev - SnowWaterEquiv(I,J)
       endif

       ! The rate of melting and freezing for glacier ice
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( (IndexPhaseChange(I,LoopInd1,J) > 0) .and. (abs(EnergyRes(I,LoopInd1,J)) > 0.0) ) then
             EnergyResLeft(I,LoopInd1,J) = 0.0
             if ( GlacierPhaseChg(I,LoopInd1,J) > 0.0 ) then
                MassWatIceTmp(I,LoopInd1,J) = max(0.0, MassWatIceInit(I,LoopInd1,J)-GlacierPhaseChg(I,LoopInd1,J))
                EnergyResLeft(I,LoopInd1,J) = EnergyRes(I,LoopInd1,J) - ConstLatHeatFusion * &
                                          (MassWatIceInit(I,LoopInd1,J) - MassWatIceTmp(I,LoopInd1,J)) / MainTimeStep
             elseif ( GlacierPhaseChg(I,LoopInd1,J) < 0.0 ) then
                MassWatIceTmp(I,LoopInd1,J) = min(MassWatTotInit(I,LoopInd1,J), MassWatIceInit(I,LoopInd1,J)-GlacierPhaseChg(I,LoopInd1,J))
                EnergyResLeft(I,LoopInd1,J) = EnergyRes(I,LoopInd1,J) - ConstLatHeatFusion * &
                                          (MassWatIceInit(I,LoopInd1,J) - MassWatIceTmp(I,LoopInd1,J)) / MainTimeStep
             endif
             MassWatLiqTmp(I,LoopInd1,J)    = max(0.0, MassWatTotInit(I,LoopInd1,J)-MassWatIceTmp(I,LoopInd1,J)) ! update liquid water mass

             ! update ice temperature and energy surplus/loss
             if ( abs(EnergyResLeft(I,LoopInd1,J)) > 0.0 ) then
                TemperatureSoilSnow(I,LoopInd1,J) = TemperatureSoilSnow(I,LoopInd1,J) + &
                                                    PhaseChgFacSoilSnow(I,LoopInd1,J) * EnergyResLeft(I,LoopInd1,J)
             endif
             HeatLhTotPhsChg = HeatLhTotPhsChg + &
                               ConstLatHeatFusion * (MassWatIceInit(I,LoopInd1,J) - MassWatIceTmp(I,LoopInd1,J)) / MainTimeStep
          endif
       enddo
       !$acc loop seq
       do LoopInd1 = -NumSnowLayerMax+1, NumSoilLayer
          EnergyResLeft(I,LoopInd1,J)   = 0.0
          GlacierPhaseChg(I,LoopInd1,J) = 0.0
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
                EnergyResLeft(I,LoopInd1,J) = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
                !$acc loop seq
                do LoopInd2 = 1, NumSoilLayer
                   if ( (LoopInd1 /= LoopInd2) .and. (TemperatureSoilSnow(I,LoopInd2,J) < ConstFreezePoint) .and. &
                        (EnergyResLeft(I,LoopInd1,J) > 0.1) ) then
                      EnergyResLeft(I,LoopInd2,J) = (TemperatureSoilSnow(I,LoopInd2,J) - ConstFreezePoint) / &
                                                PhaseChgFacSoilSnow(I,LoopInd2,J)
                      if ( abs(EnergyResLeft(I,LoopInd2,J)) > EnergyResLeft(I,LoopInd1,J) ) then ! LAYER ABSORBS ALL
                         EnergyResLeft(I,LoopInd2,J)       = EnergyResLeft(I,LoopInd2,J) + EnergyResLeft(I,LoopInd1,J)
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint + &
                                                             EnergyResLeft(I,LoopInd2,J) * PhaseChgFacSoilSnow(I,LoopInd2,J)
                         EnergyResLeft(I,LoopInd1,J)       = 0.0
                      else
                         EnergyResLeft(I,LoopInd1,J)       = EnergyResLeft(I,LoopInd1,J) + EnergyResLeft(I,LoopInd2,J)
                         EnergyResLeft(I,LoopInd2,J)       = 0.0
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                      endif
                   endif
                enddo
                TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint + EnergyResLeft(I,LoopInd1,J) * PhaseChgFacSoilSnow(I,LoopInd1,J)
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
                EnergyResLeft(I,LoopInd1,J) = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
                !$acc loop seq
                do LoopInd2 = 1, NumSoilLayer
                   if ( (LoopInd1 /= LoopInd2) .and. (TemperatureSoilSnow(I,LoopInd2,J) > ConstFreezePoint) .and. &
                        (EnergyResLeft(I,LoopInd1,J) < -0.1) ) then
                      EnergyResLeft(I,LoopInd2,J) = (TemperatureSoilSnow(I,LoopInd2,J) - ConstFreezePoint) / &
                                                PhaseChgFacSoilSnow(I,LoopInd2,J)
                      if ( EnergyResLeft(I,LoopInd2,J) > abs(EnergyResLeft(I,LoopInd1,J)) ) then  ! LAYER ABSORBS ALL
                         EnergyResLeft(I,LoopInd2,J)       = EnergyResLeft(I,LoopInd2,J) + EnergyResLeft(I,LoopInd1,J)
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint + &
                                                             EnergyResLeft(I,LoopInd2,J) * PhaseChgFacSoilSnow(I,LoopInd2,J)
                         EnergyResLeft(I,LoopInd1,J)       = 0.0
                      else
                         EnergyResLeft(I,LoopInd1,J)       = EnergyResLeft(I,LoopInd1,J) + EnergyResLeft(I,LoopInd2,J)
                         EnergyResLeft(I,LoopInd2,J)       = 0.0
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                      endif
                   endif
                enddo
                TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint + EnergyResLeft(I,LoopInd1,J) * PhaseChgFacSoilSnow(I,LoopInd1,J)
             endif
          enddo
       endif

       ! Check for ice and temperature flags
       FlagAnyAboveFreeze = .false.
       FlagAnyIce = .false.
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( TemperatureSoilSnow(I,LoopInd1,J) > ConstFreezePoint ) FlagAnyAboveFreeze = .true.
          if ( MassWatIceTmp(I,LoopInd1,J) > 0.0 ) FlagAnyIce = .true.
       enddo

       ! now remove excess heat by melting ice
       if ( FlagAnyAboveFreeze .and. FlagAnyIce ) then
          !$acc loop seq
          do LoopInd1 = 1, NumSoilLayer
             if ( TemperatureSoilSnow(I,LoopInd1,J) > ConstFreezePoint ) then
                EnergyResLeft(I,LoopInd1,J)   = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
                GlacierPhaseChg(I,LoopInd1,J) = EnergyResLeft(I,LoopInd1,J) * MainTimeStep / ConstLatHeatFusion
                !$acc loop seq
                do LoopInd2 = 1, NumSoilLayer
                   if ( (LoopInd1 /= LoopInd2) .and. (MassWatIceTmp(I,LoopInd2,J) > 0.0) .and. &
                        (GlacierPhaseChg(I,LoopInd1,J) > 0.1) ) then
                      if ( MassWatIceTmp(I,LoopInd2,J) > GlacierPhaseChg(I,LoopInd1,J) ) then  ! LAYER ABSORBS ALL
                         MassWatIceTmp(I,LoopInd2,J)       = MassWatIceTmp(I,LoopInd2,J) - GlacierPhaseChg(I,LoopInd1,J)
                         HeatLhTotPhsChg               = HeatLhTotPhsChg + &
                                                         ConstLatHeatFusion * GlacierPhaseChg(I,LoopInd1,J)/MainTimeStep
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                         GlacierPhaseChg(I,LoopInd1,J)     = 0.0
                      else
                         GlacierPhaseChg(I,LoopInd1,J)     = GlacierPhaseChg(I,LoopInd1,J) - MassWatIceTmp(I,LoopInd2,J)
                         HeatLhTotPhsChg               = HeatLhTotPhsChg + &
                                                         ConstLatHeatFusion * MassWatIceTmp(I,LoopInd2,J) / MainTimeStep
                         MassWatIceTmp(I,LoopInd2,J)       = 0.0
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                      endif
                      MassWatLiqTmp(I,LoopInd2,J) = max(0.0, MassWatTotInit(I,LoopInd2,J)-MassWatIceTmp(I,LoopInd2,J))
                   endif
                enddo
                EnergyResLeft(I,LoopInd1,J)       = GlacierPhaseChg(I,LoopInd1,J) * ConstLatHeatFusion / MainTimeStep
                TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint + EnergyResLeft(I,LoopInd1,J) * PhaseChgFacSoilSnow(I,LoopInd1,J)
             endif
          enddo
       endif

       ! Check for liquid and temperature flags
       FlagAnyBelowFreeze = .false.
       FlagAnyLiq = .false.
       !$acc loop seq
       do LoopInd1 = 1, NumSoilLayer
          if ( TemperatureSoilSnow(I,LoopInd1,J) < ConstFreezePoint ) FlagAnyBelowFreeze = .true.
          if ( MassWatLiqTmp(I,LoopInd1,J) > 0.0 ) FlagAnyLiq = .true.
       enddo

       ! snow remove excess cold by refreezing liquid (may not be necessary with above loop)
       if ( FlagAnyBelowFreeze .and. FlagAnyLiq ) then
          !$acc loop seq
          do LoopInd1 = 1, NumSoilLayer
             if ( TemperatureSoilSnow(I,LoopInd1,J) < ConstFreezePoint ) then
                EnergyResLeft(I,LoopInd1,J)   = (TemperatureSoilSnow(I,LoopInd1,J) - ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd1,J)
                GlacierPhaseChg(I,LoopInd1,J) = EnergyResLeft(I,LoopInd1,J) * MainTimeStep / ConstLatHeatFusion
                !$acc loop seq
                do LoopInd2 = 1, NumSoilLayer
                   if ( (LoopInd1 /= LoopInd2) .and. (MassWatLiqTmp(I,LoopInd2,J) > 0.0) .and. &
                        (GlacierPhaseChg(I,LoopInd1,J) < -0.1) ) then
                      if ( MassWatLiqTmp(I,LoopInd2,J) > abs(GlacierPhaseChg(I,LoopInd1,J)) ) then  ! LAYER ABSORBS ALL
                         MassWatIceTmp(I,LoopInd2,J)       = MassWatIceTmp(I,LoopInd2,J) - GlacierPhaseChg(I,LoopInd1,J)
                         HeatLhTotPhsChg               = HeatLhTotPhsChg + &
                                                         ConstLatHeatFusion * GlacierPhaseChg(I,LoopInd1,J) / MainTimeStep
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                         GlacierPhaseChg(I,LoopInd1,J)     = 0.0
                      else
                         GlacierPhaseChg(I,LoopInd1,J)     = GlacierPhaseChg(I,LoopInd1,J) + MassWatLiqTmp(I,LoopInd2,J)
                         HeatLhTotPhsChg               = HeatLhTotPhsChg - &
                                                         ConstLatHeatFusion * MassWatLiqTmp(I,LoopInd2,J) / MainTimeStep
                         MassWatIceTmp(I,LoopInd2,J)       = MassWatTotInit(I,LoopInd2,J)
                         TemperatureSoilSnow(I,LoopInd2,J) = ConstFreezePoint
                      endif
                      MassWatLiqTmp(I,LoopInd2,J) = max(0.0, MassWatTotInit(I,LoopInd2,J)-MassWatIceTmp(I,LoopInd2,J))
                   endif
                enddo
                EnergyResLeft(I,LoopInd1,J)           = GlacierPhaseChg(I,LoopInd1,J) * ConstLatHeatFusion / MainTimeStep
                TemperatureSoilSnow(I,LoopInd1,J) = ConstFreezePoint + EnergyResLeft(I,LoopInd1,J) * PhaseChgFacSoilSnow(I,LoopInd1,J)
             endif
          enddo
       endif

    endif ! OptGlacierTreatment==1

    !--- update snow and soil ice and liquid content
    !$acc loop seq
    do LoopInd1 = NumSnowLayerNeg(I,J)+1, 0     ! snow
       SnowLiqWater(I,LoopInd1,J) = MassWatLiqTmp(I,LoopInd1,J)
       SnowIce(I,LoopInd1,J)      = MassWatIceTmp(I,LoopInd1,J)
    enddo
    !$acc loop seq
    do LoopInd1 = 1, NumSoilLayer       ! glacier ice
       if ( OptGlacierTreatment == 1 ) then
          SoilLiqWater(I,LoopInd1,J) = MassWatLiqTmp(I,LoopInd1,J) / (1000.0 * ThicknessSnowSoilLayer(I,LoopInd1,J))
          SoilLiqWater(I,LoopInd1,J) = max(0.0, min(1.0,SoilLiqWater(I,LoopInd1,J)))
       elseif ( OptGlacierTreatment == 2 ) then
          SoilLiqWater(I,LoopInd1,J) = 0.0             ! ice, assume all frozen forever
       endif
       SoilMoisture(I,LoopInd1,J) = 1.0
    enddo


      end do
    end do
    !$acc end parallel loop
    !$acc end data
    deallocate(EnergyRes, GlacierPhaseChg, MassWatTotInit, MassWatIceInit, MassWatLiqInit)
    deallocate(MassWatIceTmp, MassWatLiqTmp, EnergyResLeft)


    end associate

  end subroutine GlacierPhaseChange

end module GlacierPhaseChangeMod
