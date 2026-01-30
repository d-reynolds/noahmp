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
    real(kind=kind_noahmp) :: EnergyRes(-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)          ! energy residual [w/m2]
    real(kind=kind_noahmp) :: WaterPhaseChg(-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)      ! melting or freezing water [kg/m2]
    real(kind=kind_noahmp) :: MassWatTotInit(-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)     ! initial total water (ice + liq) mass
    real(kind=kind_noahmp) :: MassWatIceInit(-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)     ! initial ice content
    real(kind=kind_noahmp) :: MassWatLiqInit(-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)     ! initial liquid content
    real(kind=kind_noahmp) :: MassWatIceTmp(-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)      ! soil/snow ice mass [mm]
    real(kind=kind_noahmp) :: MassWatLiqTmp(-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)      ! soil/snow liquid water mass [mm]
    integer                               :: I, J                           ! grid indices
    !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd, EnergyResLeft, SnowWaterPrev, SnowWaterRatio, HeatLhTotPhsChg) &
    !$acc                                                       private(EnergyRes, WaterPhaseChg, MassWatTotInit, MassWatIceInit, MassWatLiqInit, MassWatIceTmp, MassWatLiqTmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
! --------------------------------------------------------------------
    associate(                                                                       &
              OptSnowAlbedo          => noahmp%config%nmlist%OptSnowAlbedo          ,& ! in,    options for ground snow surface albedo
              OptSoilSupercoolWater  => noahmp%config%nmlist%OptSoilSupercoolWater  ,& ! in,    options for soil supercooled liquid water
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax        ,& ! in,    maximum number of snow layers
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg(I,J)        ,& ! in,    actual number of snow layers (negative)
              MainTimeStep           => noahmp%config%domain%MainTimeStep           ,& ! in,    main noahmp timestep [s]
              SurfaceType            => noahmp%config%domain%SurfaceType(I,J)            ,& ! in,    surface type 1-soil; 2-lake
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
              PhaseChgFacSoilSnow    => noahmp%energy%state%PhaseChgFacSoilSnow     ,& ! in,    energy factor for soil & snow phase change
              TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow     ,& ! inout, snow and soil layer temperature [K]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! inout, soil water content [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! inout, total soil moisture [m3/m3]
              SnowIce                => noahmp%water%state%SnowIce                  ,& ! inout, snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater             ,& ! inout, snow layer liquid water [mm]
              SnowDepth              => noahmp%water%state%SnowDepth(I,J)                ,& ! inout, snow depth [m]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv(I,J)           ,& ! inout, snow water equivalent [mm]
              IndexPhaseChange       => noahmp%water%state%IndexPhaseChange         ,& ! out,   phase change index [0-none;1-melt;2-refreeze]
              SoilSupercoolWater     => noahmp%water%state%SoilSupercoolWater       ,& ! out,   supercooled water in soil [kg/m2]
              PondSfcThinSnwMelt     => noahmp%water%state%PondSfcThinSnwMelt(I,J)       ,& ! out,   surface ponding [mm] from melt when thin snow w/o layer
              MeltGroundSnow         => noahmp%water%flux%MeltGroundSnow(I,J)            ,& ! out,   ground snowmelt rate [mm/s]
              SnowFreezeRate         => noahmp%water%flux%SnowFreezeRate             & ! out,   rate of snow freezing [mm/s]
             )
! ----------------------------------------------------------------------

    !--- Initialization
    MeltGroundSnow     = 0.0
    PondSfcThinSnwMelt = 0.0
    HeatLhTotPhsChg    = 0.0

    ! supercooled water content
    !$acc loop seq
    do LoopInd = -NumSnowLayerMax+1, NumSoilLayer 
         EnergyRes(LoopInd)          = 0.0
         WaterPhaseChg(LoopInd)      = 0.0
         MassWatTotInit(LoopInd)     = 0.0
         MassWatIceInit(LoopInd)     = 0.0
         MassWatLiqInit(LoopInd)     = 0.0
         MassWatIceTmp(LoopInd)      = 0.0
         MassWatLiqTmp(LoopInd)      = 0.0
         SnowFreezeRate(I,LoopInd,J) = 0.0
         SoilSupercoolWater(I,LoopInd,J) = 0.0
    enddo

    ! snow layer water mass
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg+1, 0
       MassWatIceTmp(LoopInd) = SnowIce(I,LoopInd,J)
       MassWatLiqTmp(LoopInd) = SnowLiqWater(I,LoopInd,J)
    enddo

    ! soil layer water mass
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       MassWatLiqTmp(LoopInd) = SoilLiqWater(I,LoopInd,J) * ThicknessSnowSoilLayer(I,LoopInd,J) * 1000.0
       MassWatIceTmp(LoopInd) = (SoilMoisture(I,LoopInd,J) - SoilLiqWater(I,LoopInd,J)) * ThicknessSnowSoilLayer(I,LoopInd,J) * 1000.0
    enddo

    ! other required variables
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       IndexPhaseChange(I,LoopInd,J) = 0
       EnergyRes(LoopInd)        = 0.0
       WaterPhaseChg(LoopInd)    = 0.0
       MassWatIceInit(LoopInd)   = MassWatIceTmp(LoopInd)
       MassWatLiqInit(LoopInd)   = MassWatLiqTmp(LoopInd)
       MassWatTotInit(LoopInd)   = MassWatIceTmp(LoopInd) + MassWatLiqTmp(LoopInd)
    enddo

    !--- compute soil supercool water content
    if ( SurfaceType == 1 ) then ! land points
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
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       if ( (MassWatIceTmp(LoopInd) > 0.0) .and. (TemperatureSoilSnow(I,LoopInd,J) >= ConstFreezePoint) ) then
          IndexPhaseChange(I,LoopInd,J) = 1  ! melting
       endif
       if ( (MassWatLiqTmp(LoopInd) > SoilSupercoolWater(I,LoopInd,J)) .and. &
            (TemperatureSoilSnow(I,LoopInd,J) < ConstFreezePoint) ) then
          IndexPhaseChange(I,LoopInd,J) = 2  ! freezing
       endif
       ! If snow exists, but its thickness is not enough to create a layer
       if ( (NumSnowLayerNeg == 0) .and. (SnowWaterEquiv > 0.0) .and. (LoopInd == 1) ) then
          if ( TemperatureSoilSnow(I,LoopInd,J) >= ConstFreezePoint ) then
             IndexPhaseChange(I,LoopInd,J) = 1
          endif
       endif
    enddo

    !--- Calculate the energy surplus and loss for melting and freezing
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       if ( IndexPhaseChange(I,LoopInd,J) > 0 ) then
          EnergyRes(LoopInd)           = (TemperatureSoilSnow(I,LoopInd,J)-ConstFreezePoint) / PhaseChgFacSoilSnow(I,LoopInd,J)
          TemperatureSoilSnow(I,LoopInd,J) = ConstFreezePoint
       endif
       if ( (IndexPhaseChange(I,LoopInd,J) == 1) .and. (EnergyRes(LoopInd) < 0.0) ) then
          EnergyRes(LoopInd)        = 0.0
          IndexPhaseChange(I,LoopInd,J) = 0
       endif
       if ( (IndexPhaseChange(I,LoopInd,J) == 2) .and. (EnergyRes(LoopInd) > 0.0) ) then
          EnergyRes(LoopInd)        = 0.0
          IndexPhaseChange(I,LoopInd,J) = 0
       endif
       WaterPhaseChg(LoopInd) = EnergyRes(LoopInd) * MainTimeStep / ConstLatHeatFusion
    enddo

    !--- The rate of melting for snow without a layer, needs more work.
    if ( (NumSnowLayerNeg == 0) .and. (SnowWaterEquiv > 0.0) .and. (WaterPhaseChg(1) > 0.0) ) then
       SnowWaterPrev  = SnowWaterEquiv
       SnowWaterEquiv = max(0.0, SnowWaterPrev-WaterPhaseChg(1))
       SnowWaterRatio = SnowWaterEquiv / SnowWaterPrev
       SnowDepth      = max(0.0, SnowWaterRatio*SnowDepth )
       SnowDepth      = min(max(SnowDepth,SnowWaterEquiv/500.0), SnowWaterEquiv/50.0)      ! limit adjustment to a reasonable density
       EnergyResLeft  = EnergyRes(1) - ConstLatHeatFusion * (SnowWaterPrev - SnowWaterEquiv) / MainTimeStep
       if ( EnergyResLeft > 0.0 ) then
          WaterPhaseChg(1) = EnergyResLeft * MainTimeStep / ConstLatHeatFusion
          EnergyRes(1)     = EnergyResLeft
       else
          WaterPhaseChg(1) = 0.0
          EnergyRes(1)     = 0.0
       endif
       MeltGroundSnow     = max(0.0, (SnowWaterPrev-SnowWaterEquiv)) / MainTimeStep
       HeatLhTotPhsChg    = ConstLatHeatFusion * MeltGroundSnow
       PondSfcThinSnwMelt = SnowWaterPrev - SnowWaterEquiv
    endif

    ! The rate of melting and freezing for multi-layer snow and soil
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       if ( (IndexPhaseChange(I,LoopInd,J) > 0) .and. (abs(EnergyRes(LoopInd)) > 0.0) ) then
          EnergyResLeft = 0.0
          if ( WaterPhaseChg(LoopInd) > 0.0 ) then
             MassWatIceTmp(LoopInd) = max(0.0, MassWatIceInit(LoopInd)-WaterPhaseChg(LoopInd))
             EnergyResLeft          = EnergyRes(LoopInd) - ConstLatHeatFusion * &
                                      (MassWatIceInit(LoopInd) - MassWatIceTmp(LoopInd)) / MainTimeStep
          elseif ( WaterPhaseChg(LoopInd) < 0.0 ) then
             if ( LoopInd <= 0 ) then  ! snow layer
                MassWatIceTmp(LoopInd) = min(MassWatTotInit(LoopInd), MassWatIceInit(LoopInd)-WaterPhaseChg(LoopInd))
             else                      ! soil layer
                if ( MassWatTotInit(LoopInd) < SoilSupercoolWater(I,LoopInd,J) ) then
                   MassWatIceTmp(LoopInd) = 0.0
                else
                   MassWatIceTmp(LoopInd) = min(MassWatTotInit(LoopInd)-SoilSupercoolWater(I,LoopInd,J), &
                                                MassWatIceInit(LoopInd)-WaterPhaseChg(LoopInd))
                   MassWatIceTmp(LoopInd) = max(MassWatIceTmp(LoopInd), 0.0)
                endif
             endif
             EnergyResLeft = EnergyRes(LoopInd) - ConstLatHeatFusion * (MassWatIceInit(LoopInd) - &
                                                                        MassWatIceTmp(LoopInd)) / MainTimeStep
          endif
          MassWatLiqTmp(LoopInd) = max(0.0, MassWatTotInit(LoopInd)-MassWatIceTmp(LoopInd)) ! update liquid water mass

          ! update soil/snow temperature and energy surplus/loss
          if ( abs(EnergyResLeft) > 0.0 ) then
             TemperatureSoilSnow(I,LoopInd,J) = TemperatureSoilSnow(I,LoopInd,J) + PhaseChgFacSoilSnow(I,LoopInd,J) * EnergyResLeft
             if ( LoopInd <= 0 ) then  ! snow
                if ( (MassWatLiqTmp(LoopInd)*MassWatIceTmp(LoopInd)) > 0.0 ) &
                   TemperatureSoilSnow(I,LoopInd,J) = ConstFreezePoint
                if ( MassWatIceTmp(LoopInd) == 0.0 ) then         ! BARLAGE
                   TemperatureSoilSnow(I,LoopInd,J) = ConstFreezePoint
                   EnergyRes(LoopInd+1)         = EnergyRes(LoopInd+1) + EnergyResLeft
                   WaterPhaseChg(LoopInd+1)     = EnergyRes(LoopInd+1) * MainTimeStep / ConstLatHeatFusion
                endif
             endif
          endif
          HeatLhTotPhsChg = HeatLhTotPhsChg + ConstLatHeatFusion * &
                            (MassWatIceInit(LoopInd) - MassWatIceTmp(LoopInd)) / MainTimeStep
          ! snow melting rate
          if ( LoopInd < 1 ) then
             MeltGroundSnow = MeltGroundSnow + max(0.0, (MassWatIceInit(LoopInd)-MassWatIceTmp(LoopInd))) / MainTimeStep
             if (OptSnowAlbedo == 3) then
                SnowFreezeRate(I,LoopInd,J) = max(0.0, (MassWatIceTmp(LoopInd)-MassWatIceInit(LoopInd))) / MainTimeStep
             endif
          endif
       endif
    enddo

    !--- update snow and soil ice and liquid content
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg+1, 0     ! snow
       SnowLiqWater(I,LoopInd,J) = MassWatLiqTmp(LoopInd)
       SnowIce(I,LoopInd,J)      = MassWatIceTmp(LoopInd)
    enddo
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer       ! soil
       SoilLiqWater(I,LoopInd,J) = MassWatLiqTmp(LoopInd) / (1000.0 * ThicknessSnowSoilLayer(I,LoopInd,J))
       SoilMoisture(I,LoopInd,J) = (MassWatLiqTmp(LoopInd)+MassWatIceTmp(LoopInd)) / (1000.0*ThicknessSnowSoilLayer(I,LoopInd,J))
    enddo

    end associate
   end do
   end do


  end subroutine SoilSnowWaterPhaseChange

end module SoilSnowWaterPhaseChangeMod
