module SoilSnowWaterPhaseChangeMod

!!! Compute the phase change (melting/freezing) of snow water and soil water

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SoilWaterSupercoolKoren99Mod, only : SoilWaterSupercoolKoren99
  use SoilWaterSupercoolNiu06Mod,   only : SoilWaterSupercoolNiu06

  implicit none

contains

  subroutine SoilSnowWaterPhaseChange(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: PHASECHANGE
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)     , intent(inout) :: noahmp

! local variable
    integer                               :: LoopInd                        ! do loop index
    real(kind=kind_noahmp)                :: EnergyResLeft                  ! energy residual or loss after melting/freezing
    real(kind=kind_noahmp)                :: SnowWaterPrev                  ! old/previous snow water equivalent [kg/m2]
    real(kind=kind_noahmp)                :: SnowWaterRatio                 ! ratio of previous vs updated snow water equivalent 
    real(kind=kind_noahmp)                :: HeatLhTotPhsChg                ! total latent heat of phase change
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: EnergyRes          ! energy residual [w/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: WaterPhaseChg      ! melting or freezing water [kg/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassWatTotInit     ! initial total water (ice + liq) mass
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassWatIceInit     ! initial ice content
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassWatLiqInit     ! initial liquid content
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassWatIceTmp      ! soil/snow ice mass [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassWatLiqTmp      ! soil/snow liquid water mass [mm]
    integer                               :: I, J                           ! grid indices
    associate(                                                                       &
              OptSnowAlbedo          => noahmp%config%nmlist%OptSnowAlbedo          ,& ! in,    options for ground snow surface albedo
              OptSoilSupercoolWater  => noahmp%config%nmlist%OptSoilSupercoolWater  ,& ! in,    options for soil supercooled liquid water
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax        ,& ! in,    maximum number of snow layers
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg        ,& ! in,    actual number of snow layers (negative)
              MainTimeStep           => noahmp%config%domain%MainTimeStep           ,& ! in,    main noahmp timestep [s]
              SurfaceType            => noahmp%config%domain%SurfaceType            ,& ! in,    surface type 1-soil; 2-lake
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
              PhaseChgFacSoilSnow    => noahmp%energy%state%PhaseChgFacSoilSnow     ,& ! in,    energy factor for soil & snow phase change
              TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow     ,& ! inout, snow and soil layer temperature [K]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! inout, soil water content [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! inout, total soil moisture [m3/m3]
              SnowIce                => noahmp%water%state%SnowIce                  ,& ! inout, snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater             ,& ! inout, snow layer liquid water [mm]
              SnowDepth              => noahmp%water%state%SnowDepth                ,& ! inout, snow depth [m]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv           ,& ! inout, snow water equivalent [mm]
              IndexPhaseChange       => noahmp%water%state%IndexPhaseChange         ,& ! out,   phase change index [0-none;1-melt;2-refreeze]
              SoilSupercoolWater     => noahmp%water%state%SoilSupercoolWater       ,& ! out,   supercooled water in soil [kg/m2]
              PondSfcThinSnwMelt     => noahmp%water%state%PondSfcThinSnwMelt       ,& ! out,   surface ponding [mm] from melt when thin snow w/o layer
              MeltGroundSnow         => noahmp%water%flux%MeltGroundSnow            ,& ! out,   ground snowmelt rate [mm/s]
              SnowFreezeRate         => noahmp%water%flux%SnowFreezeRate             & ! out,   rate of snow freezing [mm/s]
             )

    allocate(EnergyRes(noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(WaterPhaseChg(noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassWatTotInit(noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassWatIceInit(noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassWatLiqInit(noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassWatIceTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassWatLiqTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, -NumSnowLayerMax+1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    !$acc data create(EnergyRes, WaterPhaseChg, MassWatTotInit, MassWatIceInit, MassWatLiqInit, MassWatIceTmp, MassWatLiqTmp)
    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd, EnergyResLeft, SnowWaterPrev, SnowWaterRatio, HeatLhTotPhsChg)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip land ice points

    !--- Initialization
    MeltGroundSnow(I,J)     = 0.0
    PondSfcThinSnwMelt(I,J) = 0.0
    HeatLhTotPhsChg    = 0.0

    ! supercooled water content
    !$acc loop seq
    do LoopInd = -NumSnowLayerMax+1, NumSoilLayer 
         EnergyRes(I,LoopInd,J)          = 0.0
         WaterPhaseChg(I,LoopInd,J)      = 0.0
         MassWatTotInit(I,LoopInd,J)     = 0.0
         MassWatIceInit(I,LoopInd,J)     = 0.0
         MassWatLiqInit(I,LoopInd,J)     = 0.0
         MassWatIceTmp(I,LoopInd,J)      = 0.0
         MassWatLiqTmp(I,LoopInd,J)      = 0.0
         SoilSupercoolWater(I,LoopInd,J) = 0.0
         if (OptSnowAlbedo == 3) then
            SnowFreezeRate(I,LoopInd,J) = 0.0
         endif
    enddo

    ! snow layer water mass
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, 0
       MassWatIceTmp(I,LoopInd,J) = SnowIce(I,LoopInd,J)
       MassWatLiqTmp(I,LoopInd,J) = SnowLiqWater(I,LoopInd,J)
    enddo

    ! soil layer water mass
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       MassWatLiqTmp(I,LoopInd,J) = SoilLiqWater(I,LoopInd,J) * ThicknessSnowSoilLayer(I,LoopInd,J) * 1000.0
       MassWatIceTmp(I,LoopInd,J) = (SoilMoisture(I,LoopInd,J) - SoilLiqWater(I,LoopInd,J)) * ThicknessSnowSoilLayer(I,LoopInd,J) * 1000.0
    enddo

    ! other required variables
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       IndexPhaseChange(I,LoopInd,J) = 0
       EnergyRes(I,LoopInd,J)        = 0.0
       WaterPhaseChg(I,LoopInd,J)    = 0.0
       MassWatIceInit(I,LoopInd,J)   = MassWatIceTmp(I,LoopInd,J)
       MassWatLiqInit(I,LoopInd,J)   = MassWatLiqTmp(I,LoopInd,J)
       MassWatTotInit(I,LoopInd,J)   = MassWatIceTmp(I,LoopInd,J) + MassWatLiqTmp(I,LoopInd,J)
    enddo

    !--- compute soil supercool water content
    if ( SurfaceType(I,J) == 1 ) then ! land points
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          if ( OptSoilSupercoolWater == 1 ) then
             if ( TemperatureSoilSnow(I,LoopInd,J) < ConstFreezePoint ) then
                call SoilWaterSupercoolNiu06(noahmp, LoopInd, SoilSupercoolWater(I,LoopInd,J),TemperatureSoilSnow(I,LoopInd,J), I, J)
                SoilSupercoolWater(I,LoopInd,J) = SoilSupercoolWater(I,LoopInd,J) * ThicknessSnowSoilLayer(I,LoopInd,J) * 1000.0
             endif
          endif
          if ( OptSoilSupercoolWater == 2 ) then
             if ( TemperatureSoilSnow(I,LoopInd,J) < ConstFreezePoint ) then
                call SoilWaterSupercoolKoren99(noahmp, LoopInd, SoilSupercoolWater(I,LoopInd,J), &
                                               TemperatureSoilSnow(I,LoopInd,J), SoilMoisture(I,LoopInd,J), SoilLiqWater(I,LoopInd,J), I, J)
                SoilSupercoolWater(I,LoopInd,J) = SoilSupercoolWater(I,LoopInd,J) * ThicknessSnowSoilLayer(I,LoopInd,J) * 1000.0
             endif
          endif
       enddo
    endif

    !--- determine melting or freezing state
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       if ( (MassWatIceTmp(I,LoopInd,J) > 0.0) .and. (TemperatureSoilSnow(I,LoopInd,J) >= ConstFreezePoint) ) then
          IndexPhaseChange(I,LoopInd,J) = 1  ! melting
       endif
       if ( (MassWatLiqTmp(I,LoopInd,J) > SoilSupercoolWater(I,LoopInd,J)) .and. &
            (TemperatureSoilSnow(I,LoopInd,J) < ConstFreezePoint) ) then
          IndexPhaseChange(I,LoopInd,J) = 2  ! freezing
       endif
       ! If snow exists, but its thickness is not enough to create a layer
       if ( (NumSnowLayerNeg(I,J) == 0) .and. (SnowWaterEquiv(I,J) > 0.0) .and. (LoopInd == 1) ) then
          if ( TemperatureSoilSnow(I,LoopInd,J) >= ConstFreezePoint ) then
             IndexPhaseChange(I,LoopInd,J) = 1
          endif
       endif
    enddo

    !--- Calculate the energy surplus and loss for melting and freezing
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       if ( IndexPhaseChange(I,LoopInd,J) > 0 ) then
          EnergyRes(I,LoopInd,J)           = (TemperatureSoilSnow(I,LoopInd,J)-ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd,J)
          TemperatureSoilSnow(I,LoopInd,J) = ConstFreezePoint
       endif
       if ( (IndexPhaseChange(I,LoopInd,J) == 1) .and. (EnergyRes(I,LoopInd,J) < 0.0) ) then
          EnergyRes(I,LoopInd,J)        = 0.0
          IndexPhaseChange(I,LoopInd,J) = 0
       endif
       if ( (IndexPhaseChange(I,LoopInd,J) == 2) .and. (EnergyRes(I,LoopInd,J) > 0.0) ) then
          EnergyRes(I,LoopInd,J)        = 0.0
          IndexPhaseChange(I,LoopInd,J) = 0
       endif
       WaterPhaseChg(I,LoopInd,J) = EnergyRes(I,LoopInd,J) * MainTimeStep / ConstLatHeatFusion
    enddo

    !--- The rate of melting for snow without a layer, needs more work.
    if ( (NumSnowLayerNeg(I,J) == 0) .and. (SnowWaterEquiv(I,J) > 0.0) .and. (WaterPhaseChg(I,1,J) > 0.0) ) then
       SnowWaterPrev  = SnowWaterEquiv(I,J)
       SnowWaterEquiv(I,J) = max(0.0, SnowWaterPrev-WaterPhaseChg(I,1,J))
       SnowWaterRatio = SnowWaterEquiv(I,J) / SnowWaterPrev
       SnowDepth(I,J)      = max(0.0, SnowWaterRatio*SnowDepth(I,J) )
       SnowDepth(I,J)      = min(max(SnowDepth(I,J),SnowWaterEquiv(I,J)/500.0), SnowWaterEquiv(I,J)/50.0)      ! limit adjustment to a reasonable density
       EnergyResLeft  = EnergyRes(I,1,J) - ConstLatHeatFusion * (SnowWaterPrev - SnowWaterEquiv(I,J)) / MainTimeStep
       if ( EnergyResLeft > 0.0 ) then
          WaterPhaseChg(I,1,J) = EnergyResLeft * MainTimeStep / ConstLatHeatFusion
          EnergyRes(I,1,J)     = EnergyResLeft
       else
          WaterPhaseChg(I,1,J) = 0.0
          EnergyRes(I,1,J)     = 0.0
       endif
       MeltGroundSnow(I,J)     = max(0.0, (SnowWaterPrev-SnowWaterEquiv(I,J))) / MainTimeStep
       HeatLhTotPhsChg    = ConstLatHeatFusion * MeltGroundSnow(I,J)
       PondSfcThinSnwMelt(I,J) = SnowWaterPrev - SnowWaterEquiv(I,J)
    endif

    ! The rate of melting and freezing for multi-layer snow and soil
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       if ( (IndexPhaseChange(I,LoopInd,J) > 0) .and. (abs(EnergyRes(I,LoopInd,J)) > 0.0) ) then
          EnergyResLeft = 0.0
          if ( WaterPhaseChg(I,LoopInd,J) > 0.0 ) then
             MassWatIceTmp(I,LoopInd,J) = max(0.0, MassWatIceInit(I,LoopInd,J)-WaterPhaseChg(I,LoopInd,J))
             EnergyResLeft          = EnergyRes(I,LoopInd,J) - ConstLatHeatFusion * &
                                      (MassWatIceInit(I,LoopInd,J) - MassWatIceTmp(I,LoopInd,J)) / MainTimeStep
          elseif ( WaterPhaseChg(I,LoopInd,J) < 0.0 ) then
             if ( LoopInd <= 0 ) then  ! snow layer
                MassWatIceTmp(I,LoopInd,J) = min(MassWatTotInit(I,LoopInd,J), MassWatIceInit(I,LoopInd,J)-WaterPhaseChg(I,LoopInd,J))
             else                      ! soil layer
                if ( MassWatTotInit(I,LoopInd,J) < SoilSupercoolWater(I,LoopInd,J) ) then
                   MassWatIceTmp(I,LoopInd,J) = 0.0
                else
                   MassWatIceTmp(I,LoopInd,J) = min(MassWatTotInit(I,LoopInd,J)-SoilSupercoolWater(I,LoopInd,J), &
                                                MassWatIceInit(I,LoopInd,J)-WaterPhaseChg(I,LoopInd,J))
                   MassWatIceTmp(I,LoopInd,J) = max(MassWatIceTmp(I,LoopInd,J), 0.0)
                endif
             endif
             EnergyResLeft = EnergyRes(I,LoopInd,J) - ConstLatHeatFusion * (MassWatIceInit(I,LoopInd,J) - &
                                                                        MassWatIceTmp(I,LoopInd,J)) / MainTimeStep
          endif
          MassWatLiqTmp(I,LoopInd,J) = max(0.0, MassWatTotInit(I,LoopInd,J)-MassWatIceTmp(I,LoopInd,J)) ! update liquid water mass

          ! update soil/snow temperature and energy surplus/loss
          if ( abs(EnergyResLeft) > 0.0 ) then
             TemperatureSoilSnow(I,LoopInd,J) = TemperatureSoilSnow(I,LoopInd,J) + PhaseChgFacSoilSnow(I,LoopInd,J) * EnergyResLeft
             if ( LoopInd <= 0 ) then  ! snow
                if ( (MassWatLiqTmp(I,LoopInd,J)*MassWatIceTmp(I,LoopInd,J)) > 0.0 ) &
                   TemperatureSoilSnow(I,LoopInd,J) = ConstFreezePoint
                if ( MassWatIceTmp(I,LoopInd,J) == 0.0 ) then         ! BARLAGE
                   TemperatureSoilSnow(I,LoopInd,J) = ConstFreezePoint
                   EnergyRes(I,LoopInd+1,J)         = EnergyRes(I,LoopInd+1,J) + EnergyResLeft
                   WaterPhaseChg(I,LoopInd+1,J)     = EnergyRes(I,LoopInd+1,J) * MainTimeStep / ConstLatHeatFusion
                endif
             endif
          endif
          HeatLhTotPhsChg = HeatLhTotPhsChg + ConstLatHeatFusion * &
                            (MassWatIceInit(I,LoopInd,J) - MassWatIceTmp(I,LoopInd,J)) / MainTimeStep
          ! snow melting rate
          if ( LoopInd < 1 ) then
             MeltGroundSnow(I,J) = MeltGroundSnow(I,J) + max(0.0, (MassWatIceInit(I,LoopInd,J)-MassWatIceTmp(I,LoopInd,J))) / MainTimeStep
             if (OptSnowAlbedo == 3) then
                SnowFreezeRate(I,LoopInd,J) = max(0.0, (MassWatIceTmp(I,LoopInd,J)-MassWatIceInit(I,LoopInd,J))) / MainTimeStep
             endif
          endif
       endif
    enddo

    !--- update snow and soil ice and liquid content
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, 0     ! snow
       SnowLiqWater(I,LoopInd,J) = MassWatLiqTmp(I,LoopInd,J)
       SnowIce(I,LoopInd,J)      = MassWatIceTmp(I,LoopInd,J)
    enddo
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer       ! soil
       SoilLiqWater(I,LoopInd,J) = MassWatLiqTmp(I,LoopInd,J) / (1000.0 * ThicknessSnowSoilLayer(I,LoopInd,J))
       SoilMoisture(I,LoopInd,J) = (MassWatLiqTmp(I,LoopInd,J)+MassWatIceTmp(I,LoopInd,J)) / (1000.0*ThicknessSnowSoilLayer(I,LoopInd,J))
    enddo

   end do
   end do
    !$acc end data
    deallocate(EnergyRes, WaterPhaseChg, MassWatTotInit, MassWatIceInit, MassWatLiqInit)
    deallocate(MassWatIceTmp, MassWatLiqTmp)

    end associate

  end subroutine SoilSnowWaterPhaseChange

end module SoilSnowWaterPhaseChangeMod
