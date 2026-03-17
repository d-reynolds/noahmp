module GlacierThermalDiffusionMod

!!! Solve glacier ice and snow layer thermal diffusion
!!! Calculate the right hand side of the time tendency term of the glacier
!!! and snow thermal diffusion equation. Currently snow and glacier ice layers
!!! are coupled in solving the equations. Also compute/prepare the matrix
!!! coefficients for the tri-diagonal matrix of the implicit time scheme.
!!! Then solve the tridiagonal matrix system and update temperatures.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use MatrixSolverTriDiagonalMod, only : MatrixSolverTriDiagonal

  implicit none

contains

  subroutine GlacierThermalDiffusion(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: HRT_GLACIER and HSTEP_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)     , intent(inout) :: noahmp
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatRight  ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft1  ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft2  ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft3  ! left-hand side term of the matrix

! local variable
    integer                               :: I, J                     ! grid indices
    integer                               :: LoopInd                  ! loop index
    real(kind=kind_noahmp)                :: DepthSnowSoilTmp         ! temporary snow/soil layer depth [m]
    real(kind=kind_noahmp)                :: DepthSnowSoilInv (-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)  ! inverse of snow/soil layer depth [1/m]
    real(kind=kind_noahmp)                :: HeatCapacPerArea (-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)  ! Heat capacity of soil/snow per area [J/m2/K]
    real(kind=kind_noahmp)                :: TempGradDepth    (-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)  ! temperature gradient (derivative) with soil/snow depth [K/m]
    real(kind=kind_noahmp)                :: EnergyExcess     (-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)  ! energy flux excess in soil/snow [W/m2]

! --------------------------------------------------------------------
        associate(                                                                                &
                  NumSoilLayer             => noahmp%config%domain%NumSoilLayer                  ,& ! in,  number of soil layers
                  NumSnowLayerMax          => noahmp%config%domain%NumSnowLayerMax               ,& ! in,  maximum number of snow layers
                  NumSnowLayerNeg          => noahmp%config%domain%NumSnowLayerNeg          ,& ! in,  actual number of snow layers (negative)
                  DepthSnowSoilLayer       => noahmp%config%domain%DepthSnowSoilLayer            ,& ! in,  depth of snow/soil layer-bottom [m]
                  OptSoilTemperatureBottom => noahmp%config%nmlist%OptSoilTemperatureBottom      ,& ! in,  options for lower boundary condition of soil temperature
                  OptSnowSoilTempTime      => noahmp%config%nmlist%OptSnowSoilTempTime           ,& ! in,  options for snow/soil temperature time scheme
                  TemperatureSoilBottom    => noahmp%forcing%TemperatureSoilBottom          ,& ! in,  bottom boundary soil temperature [K]
                  DepthSoilTempBotToSno    => noahmp%energy%state%DepthSoilTempBotToSno     ,& ! in,  depth of lower boundary condition [m] from snow surface
                  TemperatureSoilSnow      => noahmp%energy%state%TemperatureSoilSnow            ,& ! inout, snow and soil layer temperature [K]
                  ThermConductSoilSnow     => noahmp%energy%state%ThermConductSoilSnow           ,& ! in,  thermal conductivity [W/m/K] for all soil & snow
                  HeatCapacSoilSnow        => noahmp%energy%state%HeatCapacSoilSnow              ,& ! in,  heat capacity [J/m3/K] for all soil & snow
                  HeatGroundTot            => noahmp%energy%flux%HeatGroundTot              ,& ! in,  total ground heat flux [W/m2] (+ to soil/snow)
                  RadSwPenetrateGrd        => noahmp%energy%flux%RadSwPenetrateGrd               ,& ! in,  light penetrating through soil/snow water [W/m2]
                  HeatFromSoilBot          => noahmp%energy%flux%HeatFromSoilBot             & ! out, energy influx from soil bottom [W/m2]
                 )

    ! !$acc parallel loop collapse(2) gang vector default(present) &
    ! !$acc private(LoopInd,DepthSnowSoilTmp,DepthSnowSoilInv,HeatCapacPerArea) &
    ! !$acc private(TempGradDepth,EnergyExcess) &
    ! !$acc private(MatLeft3Tmp,MatRightTmp,MatSolution)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if (noahmp%config%domain%IndicatorIceSfc(I,J) /= -1) cycle  ! only process glacier points

        ! initialization
        ! !$acc loop seq
        do LoopInd = -NumSnowLayerMax+1, NumSoilLayer
           MatRight(I,LoopInd,J)         = 0.0
           MatLeft1(I,LoopInd,J)         = 0.0
           MatLeft2(I,LoopInd,J)         = 0.0
           MatLeft3(I,LoopInd,J)         = 0.0
           DepthSnowSoilInv(LoopInd) = 0.0
           HeatCapacPerArea(LoopInd) = 0.0
           TempGradDepth(LoopInd)    = 0.0
           EnergyExcess(LoopInd)     = 0.0
        enddo

        ! compute gradient and flux of glacier/snow thermal diffusion
        ! !$acc loop seq
        do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
           if ( LoopInd == (NumSnowLayerNeg(I,J)+1) ) then
              HeatCapacPerArea(LoopInd) = - DepthSnowSoilLayer(I,LoopInd,J) * HeatCapacSoilSnow(I,LoopInd,J)
              DepthSnowSoilTmp          = - DepthSnowSoilLayer(I,LoopInd+1,J)
              DepthSnowSoilInv(LoopInd) = 2.0 / DepthSnowSoilTmp
              TempGradDepth(LoopInd)    = 2.0 * (TemperatureSoilSnow(I,LoopInd,J) - TemperatureSoilSnow(I,LoopInd+1,J)) / DepthSnowSoilTmp
              EnergyExcess(LoopInd)     = ThermConductSoilSnow(I,LoopInd,J) * TempGradDepth(LoopInd) - &
                                          HeatGroundTot(I,J) - RadSwPenetrateGrd(I,LoopInd,J)
           elseif ( LoopInd < NumSoilLayer ) then
              HeatCapacPerArea(LoopInd) = (DepthSnowSoilLayer(I,LoopInd-1,J) - DepthSnowSoilLayer(I,LoopInd,J)) * HeatCapacSoilSnow(I,LoopInd,J)
              DepthSnowSoilTmp          = DepthSnowSoilLayer(I,LoopInd-1,J) - DepthSnowSoilLayer(I,LoopInd+1,J)
              DepthSnowSoilInv(LoopInd) = 2.0 / DepthSnowSoilTmp
              TempGradDepth(LoopInd)    = 2.0 * (TemperatureSoilSnow(I,LoopInd,J) - TemperatureSoilSnow(I,LoopInd+1,J)) / DepthSnowSoilTmp
              EnergyExcess(LoopInd)     = (ThermConductSoilSnow(I,LoopInd,J)*TempGradDepth(LoopInd) - &
                                          ThermConductSoilSnow(I,LoopInd-1,J)*TempGradDepth(LoopInd-1) ) - RadSwPenetrateGrd(I,LoopInd,J)
           elseif ( LoopInd == NumSoilLayer ) then
              HeatCapacPerArea(LoopInd) = (DepthSnowSoilLayer(I,LoopInd-1,J) - DepthSnowSoilLayer(I,LoopInd,J)) * HeatCapacSoilSnow(I,LoopInd,J)
              DepthSnowSoilTmp          =  DepthSnowSoilLayer(I,LoopInd-1,J) - DepthSnowSoilLayer(I,LoopInd,J)
              if ( OptSoilTemperatureBottom == 1 ) then
                 HeatFromSoilBot(I,J)        = 0.0
              endif
              if ( OptSoilTemperatureBottom == 2 ) then
                 TempGradDepth(LoopInd) = (TemperatureSoilSnow(I,LoopInd,J) - TemperatureSoilBottom(I,J)) / &
                                          (0.5 * (DepthSnowSoilLayer(I,LoopInd-1,J)+DepthSnowSoilLayer(I,LoopInd,J)) - DepthSoilTempBotToSno(I,J))
                 HeatFromSoilBot(I,J)        = -ThermConductSoilSnow(I,LoopInd,J) * TempGradDepth(LoopInd)
              endif
              EnergyExcess(LoopInd)     = (-HeatFromSoilBot(I,J) - ThermConductSoilSnow(I,LoopInd-1,J)*TempGradDepth(LoopInd-1)) - &
                                          RadSwPenetrateGrd(I,LoopInd,J)
           endif
        enddo

        ! prepare the matrix coefficients for the tri-diagonal matrix
        ! !$acc loop seq
        do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
           if ( LoopInd == (NumSnowLayerNeg(I,J)+1) ) then
              MatLeft1(I,LoopInd,J)    = 0.0
              MatLeft3(I,LoopInd,J)    = - ThermConductSoilSnow(I,LoopInd,J) * DepthSnowSoilInv(LoopInd) / HeatCapacPerArea(LoopInd)
              if ( (OptSnowSoilTempTime == 1) .or. (OptSnowSoilTempTime == 3) ) then
                 MatLeft2(I,LoopInd,J) = - MatLeft3(I,LoopInd,J)
              endif
              if ( OptSnowSoilTempTime == 2 ) then
                 MatLeft2(I,LoopInd,J) = - MatLeft3(I,LoopInd,J) + ThermConductSoilSnow(I,LoopInd,J) / &
                                    (0.5*DepthSnowSoilLayer(I,LoopInd,J)*DepthSnowSoilLayer(I,LoopInd,J)*HeatCapacSoilSnow(I,LoopInd,J))
              endif
           elseif ( LoopInd < NumSoilLayer ) then
              MatLeft1(I,LoopInd,J)    = - ThermConductSoilSnow(I,LoopInd-1,J) * DepthSnowSoilInv(LoopInd-1) / HeatCapacPerArea(LoopInd)
              MatLeft3(I,LoopInd,J)    = - ThermConductSoilSnow(I,LoopInd,J  ) * DepthSnowSoilInv(LoopInd  ) / HeatCapacPerArea(LoopInd)
              MatLeft2(I,LoopInd,J)    = - (MatLeft1(I,LoopInd,J) + MatLeft3 (I,LoopInd,J))
           elseif ( LoopInd == NumSoilLayer ) then
              MatLeft1(I,LoopInd,J)    = - ThermConductSoilSnow(I,LoopInd-1,J) * DepthSnowSoilInv(LoopInd-1) / HeatCapacPerArea(LoopInd)
              MatLeft3(I,LoopInd,J)    = 0.0
              MatLeft2(I,LoopInd,J)    = - (MatLeft1(I,LoopInd,J) + MatLeft3(I,LoopInd,J))
           endif
              MatRight(I,LoopInd,J)    = EnergyExcess(LoopInd) / (-HeatCapacPerArea(LoopInd))
        enddo



      end do
    end do
    ! !$acc end parallel loop


        end associate

  end subroutine GlacierThermalDiffusion

end module GlacierThermalDiffusionMod
