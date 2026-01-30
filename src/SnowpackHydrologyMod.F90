module SnowpackHydrologyMod

!!! Snowpack hydrology processes (sublimation/frost, evaporation/dew, meltwater) (2D GPU-optimized)
!!! Update snowpack ice and liquid water content

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowLayerCombineMod, only : SnowLayerCombine

  implicit none

contains

  subroutine SnowpackHydrology(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SNOWH2O
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                          ! grid indices
    integer                          :: LoopInd                       ! do loop/array indices
    real(kind=kind_noahmp)           :: InflowSnowLayer               ! water flow into each snow layer [mm]
    real(kind=kind_noahmp)           :: SnowIceTmp                    ! ice mass after minus sublimation
    real(kind=kind_noahmp)           :: SnowWaterRatio                ! ratio of SWE after frost & sublimation to original SWE
    real(kind=kind_noahmp)           :: SnowWaterTmp                  ! temporary SWE

! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp) &
    !$acc private(LoopInd, InflowSnowLayer, SnowIceTmp, SnowWaterRatio, SnowWaterTmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                       &
                  NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax        ,& ! in,    maximum number of snow layers
                  MainTimeStep           => noahmp%config%domain%MainTimeStep           ,& ! in,    noahmp main time step [s]
                  FrostSnowSfcIce        => noahmp%water%flux%FrostSnowSfcIce(I,J)      ,& ! in,    snow surface frost rate [mm/s]
                  SublimSnowSfcIce       => noahmp%water%flux%SublimSnowSfcIce(I,J)     ,& ! in,    snow surface sublimation rate [mm/s]
                  RainfallGround         => noahmp%water%flux%RainfallGround(I,J)       ,& ! in,    ground surface rain rate [mm/s]
                  SnowLiqFracMax         => noahmp%water%param%SnowLiqFracMax(I,J)           ,& ! in,    maximum liquid water fraction in snow
                  SnowLiqHoldCap         => noahmp%water%param%SnowLiqHoldCap(I,J)           ,& ! in,    liquid water holding capacity for snowpack [m3/m3]
                  SnowLiqReleaseFac      => noahmp%water%param%SnowLiqReleaseFac(I,J)        ,& ! in,    snowpack water release timescale factor [1/s]
                  NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg(I,J)   ,& ! inout, actual number of snow layers (negative)
                  ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m] - indexed as (I,:,J)
                  SnowDepth              => noahmp%water%state%SnowDepth(I,J)           ,& ! inout, snow depth [m]
                  SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv(I,J)      ,& ! inout, snow water equivalent [mm]
                  SnowIce                => noahmp%water%state%SnowIce                  ,& ! inout, snow layer ice [mm] - indexed as (I,:,J)
                  SnowLiqWater           => noahmp%water%state%SnowLiqWater             ,& ! inout, snow layer liquid water [mm] - indexed as (I,:,J)
                  SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! inout, soil liquid moisture [m3/m3] - indexed as (I,:,J)
                  SoilIce                => noahmp%water%state%SoilIce                  ,& ! inout, soil ice moisture [m3/m3] - indexed as (I,:,J)
                  SnowIceVol             => noahmp%water%state%SnowIceVol               ,& ! inout, partial volume of snow ice [m3/m3] - indexed as (I,:,J)
                  SnowLiqWaterVol        => noahmp%water%state%SnowLiqWaterVol          ,& ! inout, partial volume of snow liquid water [m3/m3] - indexed as (I,:,J)
                  SnowEffPorosity        => noahmp%water%state%SnowEffPorosity          ,& ! out,   snow effective porosity [m3/m3] - indexed as (I,:,J)
                  SnowBotOutflow         => noahmp%water%flux%SnowBotOutflow(I,J)       ,& ! out,   total water (snowmelt + rain through pack) out of snowpack bottom [mm/s]
                  OutflowSnowLayer       => noahmp%water%flux%OutflowSnowLayer           & ! out,   water flow out of each snow layer [mm/s] - indexed as (I,:,J)
                 )
! ----------------------------------------------------------------------

        ! initialization
        !$acc loop seq
        do LoopInd = -NumSnowLayerMax+1, 0
           SnowEffPorosity(I,LoopInd,J) = 0.0
           OutflowSnowLayer(I,LoopInd,J) = 0.0
        enddo
        SnowBotOutflow  = 0.0
        InflowSnowLayer = 0.0

        ! for the case when SnowWaterEquiv becomes '0' after 'COMBINE'
        if ( SnowWaterEquiv == 0.0 ) then
           SoilIce(I,1,J) = SoilIce(I,1,J) + (FrostSnowSfcIce-SublimSnowSfcIce) * MainTimeStep / &
                                     (ThicknessSnowSoilLayer(I,1,J)*1000.0)  ! Barlage: SoilLiqWater->SoilIce v3.6
           if ( SoilIce(I,1,J) < 0.0 ) then
              SoilLiqWater(I,1,J) = SoilLiqWater(I,1,J) + SoilIce(I,1,J)
              SoilIce(I,1,J)      = 0.0
           endif
        endif

        ! for shallow snow without a layer
        ! snow surface sublimation may be larger than existing snow mass. To conserve water,
        ! excessive sublimation is used to reduce soil water. Smaller time steps would tend to aviod this problem.
        if ( (NumSnowLayerNeg == 0) .and. (SnowWaterEquiv > 0.0) ) then
           SnowWaterTmp   = SnowWaterEquiv
           SnowWaterEquiv = SnowWaterEquiv - SublimSnowSfcIce*MainTimeStep + FrostSnowSfcIce*MainTimeStep
           SnowWaterRatio = SnowWaterEquiv / SnowWaterTmp
           SnowDepth      = max(0.0, SnowWaterRatio*SnowDepth )
           SnowDepth      = min(max(SnowDepth,SnowWaterEquiv/500.0), SnowWaterEquiv/50.0)    ! limit adjustment to a reasonable density
           if ( SnowWaterEquiv < 0.0 ) then
              SoilIce(I,1,J)     = SoilIce(I,1,J) + SnowWaterEquiv / (ThicknessSnowSoilLayer(I,1,J)*1000.0)
              SnowWaterEquiv = 0.0
              SnowDepth      = 0.0
           endif
           if ( SoilIce(I,1,J) < 0.0 ) then
              SoilLiqWater(I,1,J) = SoilLiqWater(I,1,J) + SoilIce(I,1,J)
              SoilIce(I,1,J)      = 0.0
           endif
        endif

        if ( (SnowDepth <= 1.0e-8) .or. (SnowWaterEquiv <= 1.0e-6) ) then
           SnowDepth      = 0.0
           SnowWaterEquiv = 0.0
        endif

        ! for multi-layer (>=1) snow
        if ( NumSnowLayerNeg < 0 ) then
          SnowWaterTmp = SnowIce(I,NumSnowLayerNeg+1,J) + SnowLiqWater(I,NumSnowLayerNeg+1,J) ! top layer total snow water before sublimation
          SnowIceTmp = SnowIce(I,NumSnowLayerNeg+1,J) - SublimSnowSfcIce*MainTimeStep + FrostSnowSfcIce*MainTimeStep
          SnowIce(I,NumSnowLayerNeg+1,J) = SnowIceTmp
      if ( (SnowIceTmp < 1.0e-6) .and. (NumSnowLayerNeg < 0) ) call SnowLayerCombine(noahmp, I, J)
          if ( (SnowIceTmp >= 1.0e-6) .and. (NumSnowLayerNeg < 0) ) then ! re-adjust snow layer thickness
             ThicknessSnowSoilLayer(I,NumSnowLayerNeg+1,J) = ThicknessSnowSoilLayer(I,NumSnowLayerNeg+1,J) * &
                           (SnowIce(I,NumSnowLayerNeg+1,J) + SnowLiqWater(I,NumSnowLayerNeg+1,J)) / SnowWaterTmp ! assuming same snow density
          endif
          if ( NumSnowLayerNeg < 0 ) then
             SnowLiqWater(I,NumSnowLayerNeg+1,J) = SnowLiqWater(I,NumSnowLayerNeg+1,J) + RainfallGround * MainTimeStep
             SnowLiqWater(I,NumSnowLayerNeg+1,J) = max(0.0, SnowLiqWater(I,NumSnowLayerNeg+1,J))
          endif
        endif

        ! Porosity and partial volume
        !$acc loop seq
        do LoopInd = NumSnowLayerNeg+1, 0
           SnowIceVol(I,LoopInd,J)      = min(1.0, SnowIce(I,LoopInd,J)/(ThicknessSnowSoilLayer(I,LoopInd,J)*ConstDensityIce))
           SnowEffPorosity(I,LoopInd,J) = 1.0 - SnowIceVol(I,LoopInd,J)
        enddo

        ! compute inter-layer snow water flow
        !$acc loop seq
        do LoopInd = NumSnowLayerNeg+1, 0
           SnowLiqWater(I,LoopInd,J)     = SnowLiqWater(I,LoopInd,J) + InflowSnowLayer
           SnowLiqWaterVol(I,LoopInd,J)  = SnowLiqWater(I,LoopInd,J) / (ThicknessSnowSoilLayer(I,LoopInd,J)*ConstDensityWater)
           OutflowSnowLayer(I,LoopInd,J) = max(0.0, (SnowLiqWaterVol(I,LoopInd,J)-SnowLiqHoldCap*SnowEffPorosity(I,LoopInd,J)) * &
                                           ThicknessSnowSoilLayer(I,LoopInd,J))
           if ( LoopInd == 0 ) then
              OutflowSnowLayer(I,LoopInd,J) = max((SnowLiqWaterVol(I,LoopInd,J)-SnowEffPorosity(I,LoopInd,J)) * ThicknessSnowSoilLayer(I,LoopInd,J), &
                                              SnowLiqReleaseFac * MainTimeStep * OutflowSnowLayer(I,LoopInd,J))
           endif
           OutflowSnowLayer(I,LoopInd,J)    = OutflowSnowLayer(I,LoopInd,J) * ConstDensityWater
           SnowLiqWater(I,LoopInd,J) = SnowLiqWater(I,LoopInd,J) - OutflowSnowLayer(I,LoopInd,J)
           if ( (SnowLiqWater(I,LoopInd,J)/(SnowIce(I,LoopInd,J)+SnowLiqWater(I,LoopInd,J))) > SnowLiqFracMax ) then
              OutflowSnowLayer(I,LoopInd,J) = OutflowSnowLayer(I,LoopInd,J) + (SnowLiqWater(I,LoopInd,J) - &
                                          SnowLiqFracMax / (1.0-SnowLiqFracMax) * SnowIce(I,LoopInd,J))
              SnowLiqWater(I,LoopInd,J) = SnowLiqFracMax / (1.0 - SnowLiqFracMax) * SnowIce(I,LoopInd,J)
           endif
           InflowSnowLayer = OutflowSnowLayer(I,LoopInd,J)
           SnowLiqWaterVol(I,LoopInd,J)  = SnowLiqWater(I,LoopInd,J) / (ThicknessSnowSoilLayer(I,LoopInd,J)*ConstDensityWater) ! update SnowLiqWaterVol
        enddo

        ! update snow depth
        !$acc loop seq
        do LoopInd = NumSnowLayerNeg+1, 0
           ThicknessSnowSoilLayer(I,LoopInd,J) = max(ThicknessSnowSoilLayer(I,LoopInd,J), &
                                                 SnowLiqWater(I,LoopInd,J)/ConstDensityWater+SnowIce(I,LoopInd,J)/ConstDensityIce)
        enddo

        ! Liquid water from snow bottom to soil [mm/s]
        SnowBotOutflow = OutflowSnowLayer(I,0,J) / MainTimeStep
        !$acc loop seq
        do LoopInd = -NumSnowLayerMax+1, 0
           OutflowSnowLayer(I,LoopInd,J) = OutflowSnowLayer(I,LoopInd,J) / MainTimeStep
        enddo

        end associate

      enddo
    enddo

  end subroutine SnowpackHydrology

end module SnowpackHydrologyMod
