module SnowfallBelowCanopyMod

!!! Snowfall process after canopy interception (2D GPU-optimized)
!!! Update snow water equivalent and snow depth

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowfallAfterCanopyIntercept(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SNOWFALL
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J             ! grid indices
    integer                          :: IndNewSnowLayer  ! 0-no new layers, 1-creating new layers

! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp) private(IndNewSnowLayer)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                        &
                  MainTimeStep            => noahmp%config%domain%MainTimeStep           ,& ! in,    noahmp main time step [s]
                  SnowfallGround          => noahmp%water%flux%SnowfallGround(I,J)       ,& ! in,    snowfall rate at ground [mm/s]
                  SnowDepthIncr           => noahmp%water%flux%SnowDepthIncr(I,J)        ,& ! in,    snow depth increasing rate [m/s] due to snowfall
                  TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight(I,J) ,& ! in,    air temperature [K] at reference height
                  NumSnowLayerNeg         => noahmp%config%domain%NumSnowLayerNeg(I,J)   ,& ! inout, actual number of snow layers (negative)
                  ThicknessSnowSoilLayer  => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m]
                  SnowDepth               => noahmp%water%state%SnowDepth(I,J)           ,& ! inout, snow depth [m]
                  SnowWaterEquiv          => noahmp%water%state%SnowWaterEquiv(I,J)      ,& ! inout, snow water equivalent [mm]
                  SnowIce                 => noahmp%water%state%SnowIce                  ,& ! inout, snow layer ice [mm]
                  SnowLiqWater            => noahmp%water%state%SnowLiqWater             ,& ! inout, snow layer liquid water [mm]
                  TemperatureSoilSnow     => noahmp%energy%state%TemperatureSoilSnow      & ! inout, snow and soil layer temperature [K]

                 )
! ----------------------------------------------------------------------

    IndNewSnowLayer = 0

    ! shallow snow / no layer
    if ( (NumSnowLayerNeg == 0) .and. (SnowfallGround > 0.0) ) then
       SnowDepth      = SnowDepth + SnowDepthIncr * MainTimeStep
       SnowWaterEquiv = SnowWaterEquiv + SnowfallGround * MainTimeStep
    endif

    ! creating a new layer
    !if ( (NumSnowLayerNeg == 0)  .and. (SnowfallGround > 0.0) .and. (SnowDepth >= 0.05) ) then
    !if ( (NumSnowLayerNeg == 0)  .and. (SnowfallGround > 0.0) .and. (SnowDepth >= 0.025) ) then !MB: change limit
    ! C.He: remove SnowfallGround > 0.0 to allow adjusting snow layer number based on SnowDepth when no snowfall
    if ( (NumSnowLayerNeg == 0) .and. (SnowDepth >= 0.025) ) then
       NumSnowLayerNeg                                 = -1
       IndNewSnowLayer                                 =  1
       ThicknessSnowSoilLayer(I,0,J) = SnowDepth
       SnowDepth                                       = 0.0
       TemperatureSoilSnow(I,0,J)  = min(273.16, TemperatureAirRefHeight)   ! temporary setup
       SnowIce(I,0,J)               = SnowWaterEquiv
       SnowLiqWater(I,0,J)          = 0.0
    endif

    ! snow with layers
    if ( (NumSnowLayerNeg < 0) .and. (IndNewSnowLayer == 0) .and. (SnowfallGround > 0.0) ) then
       SnowIce(I,NumSnowLayerNeg+1,J) = SnowIce(I,NumSnowLayerNeg+1,J) + &
                                                             SnowfallGround * MainTimeStep
       ThicknessSnowSoilLayer(I,NumSnowLayerNeg+1,J) = &
           ThicknessSnowSoilLayer(I,NumSnowLayerNeg+1,J) + SnowDepthIncr * MainTimeStep
    endif

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine SnowfallAfterCanopyIntercept

end module SnowfallBelowCanopyMod
