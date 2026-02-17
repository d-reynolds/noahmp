module GlacierTemperatureSolverMod

!!! Compute Glacier and snow layer temperature using tri-diagonal matrix solution
!!! Dependent on the output from GlacierThermalDiffusion module

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use MatrixSolverTriDiagonalMod, only : MatrixSolverTriDiagonal

  implicit none

contains

  subroutine GlacierTemperatureSolver(noahmp, TimeStep, MatLeft1, MatLeft2, MatLeft3, MatRight)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: HSTEP_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)     , intent(inout) :: noahmp
    real(kind=kind_noahmp), intent(in)    :: TimeStep                             ! timestep (may not be the same as model timestep)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatRight  ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft1  ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft2  ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft3  ! left-hand side term of the matrix

! local variable
    integer                               :: I, J              ! grid indices
    integer                               :: LoopInd           ! layer loop index
    real(kind=kind_noahmp) :: MatRightTmp(-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)              ! temporary MatRight matrix coefficient
    real(kind=kind_noahmp) :: MatLeft3Tmp(-noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer)              ! temporary MatLeft3 matrix coefficient


    associate(                                                                &
              NumSoilLayer        => noahmp%config%domain%NumSoilLayer       ,& ! in,    number of glacier/soil layers
              NumSnowLayerMax     => noahmp%config%domain%NumSnowLayerMax    ,& ! in,    maximum number of snow layers
              NumSnowLayerNeg     => noahmp%config%domain%NumSnowLayerNeg    ,& ! in,    actual number of snow layers (negative)
              TemperatureSoilSnow => noahmp%energy%state%TemperatureSoilSnow  & ! inout, snow and glacier layer temperature [K]
             )

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd, MatRightTmp, MatLeft3Tmp) &
    !$acc firstprivate(TimeStep)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! initialization
      !$acc loop seq
      do LoopInd = -NumSnowLayerMax+1, NumSoilLayer
        MatRightTmp(LoopInd) = 0.0
        MatLeft3Tmp(LoopInd) = 0.0
      enddo

    ! update tri-diagonal matrix elements
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       MatRight(I,LoopInd,J) =       MatRight(I,LoopInd,J) * TimeStep
       MatLeft1(I,LoopInd,J) =       MatLeft1(I,LoopInd,J) * TimeStep
       MatLeft2(I,LoopInd,J) = 1.0 + MatLeft2(I,LoopInd,J) * TimeStep
       MatLeft3(I,LoopInd,J) =       MatLeft3(I,LoopInd,J) * TimeStep
       MatRightTmp(LoopInd) = MatRight(I,LoopInd,J)
       MatLeft3Tmp(LoopInd) = MatLeft3(I,LoopInd,J)

    enddo

    ! solve the tri-diagonal matrix equation
    call MatrixSolverTriDiagonal(MatLeft3(I,:,J),MatLeft1(I,:,J),MatLeft2(I,:,J),MatLeft3Tmp,MatRightTmp,MatRight(I,:,J),&
                                 NumSnowLayerNeg(I,J)+1,NumSoilLayer,NumSnowLayerMax)

    ! update snow & glacier temperature
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       TemperatureSoilSnow(I,LoopInd,J) = TemperatureSoilSnow(I,LoopInd,J) + MatLeft3(I,LoopInd,J)
    enddo


      end do
    end do
    !$acc end parallel loop



    end associate

  end subroutine GlacierTemperatureSolver

end module GlacierTemperatureSolverMod
