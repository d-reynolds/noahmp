module GlacierTemperatureMainMod

!!! Main module to compute snow (if exists) and glacier ice temperature.
!!! Note that snow temperatures during melting season may exceed melting
!!! point but later in GlacierPhaseChange subroutine the snow
!!! temperatures are reset to melting point for melting snow.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use GlacierThermalDiffusionMod,  only : GlacierThermalDiffusion
  use GlacierTemperatureSolverMod, only : GlacierTemperatureSolver
  implicit none

contains

  subroutine GlacierTemperatureMain(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: TSNOSOI_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)     , intent(inout) :: noahmp

! local variable
    integer                               :: I, J         ! grid indices
    integer                               :: IndLoop      ! snow and soil layer loop
    real(kind=kind_noahmp)                :: RadSwAbsSum  ! sum of RadSwAbsSnowSoilLayer
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatRight  ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatLeft1  ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatLeft2  ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatLeft3  ! left-hand side term of the matrix

! --------------------------------------------------------------------

    ! initialization and compute solar penetration
    !$acc parallel loop collapse(2) gang vector present(noahmp) private(IndLoop, RadSwAbsSum)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                         &
                  NumSoilLayer          => noahmp%config%domain%NumSoilLayer              ,& ! in,  number of glacier/soil layers
                  NumSnowLayerMax       => noahmp%config%domain%NumSnowLayerMax           ,& ! in,  maximum number of snow layers
                  NumSnowLayerNeg       => noahmp%config%domain%NumSnowLayerNeg(I,J)      ,& ! in,  actual number of snow layers (negative)
                  DepthSoilTempBottom   => noahmp%config%domain%DepthSoilTempBottom(I,J)  ,& ! in,  depth [m] from glacier surface for lower soil temperature boundary
                  OptSnowAlbedo         => noahmp%config%nmlist%OptSnowAlbedo             ,& ! in,  options for ground snow surface albedo
                  SnowDepth             => noahmp%water%state%SnowDepth(I,J)              ,& ! in,  snow depth [m]
                  RadSwAbsSnowSoilLayer => noahmp%energy%flux%RadSwAbsSnowSoilLayer       ,& ! in,  total absorbed solar radiation by snow for each layer [W/m2]
                  RadSwAbsGrd           => noahmp%energy%flux%RadSwAbsGrd(I,J)            ,& ! in,  solar radiation absorbed by ground [W/m2]
                  DepthSoilTempBotToSno => noahmp%energy%state%DepthSoilTempBotToSno(I,J) ,& ! out, depth of lower boundary condition [m] from snow surface
                  RadSwPenetrateGrd     => noahmp%energy%flux%RadSwPenetrateGrd            & ! out, light penetrating through snow/ice [W/m2]
                 )
! ----------------------------------------------------------------------

        ! initialize RadSwPenetrateGrd
        !$acc loop seq
        do IndLoop = -NumSnowLayerMax+1, NumSoilLayer
           RadSwPenetrateGrd(I,IndLoop,J) = 0.0
        enddo

        ! compute solar penetration through snowpack and glacier ice
        if (OptSnowAlbedo == 3 .and. NumSnowLayerNeg < 0) then
           ! Check if sum > 0 (replacing sum() intrinsic)
           RadSwAbsSum = 0.0
           !$acc loop seq
           do IndLoop = lbound(RadSwAbsSnowSoilLayer,2), ubound(RadSwAbsSnowSoilLayer,2)
              RadSwAbsSum = RadSwAbsSum + RadSwAbsSnowSoilLayer(I,IndLoop,J)
           enddo
           if (RadSwAbsSum > 0.0) then
              !$acc loop seq
              do IndLoop = NumSnowLayerNeg+1, 1, 1
                 if (IndLoop == NumSnowLayerNeg+1) then
                    RadSwPenetrateGrd(I,IndLoop,J) = RadSwAbsSnowSoilLayer(I,IndLoop,J) - RadSwAbsGrd
                 else
                    RadSwPenetrateGrd(I,IndLoop,J) = RadSwAbsSnowSoilLayer(I,IndLoop,J)
                 endif
              enddo
           endif
        endif

        ! adjust DepthSoilTempBottom from glacier ice surface to DepthSoilTempBotToSno from snow surface
        DepthSoilTempBotToSno = DepthSoilTempBottom - SnowDepth

        end associate

      end do
    end do
    !$acc end parallel loop

    ! initialization
    if (.not. allocated(MatRight)) allocate(MatRight(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                                                     -noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer,                   &
                                                     noahmp%config%domain%JTS:noahmp%config%domain%JTE))

    if (.not. allocated(MatLeft1)) allocate(MatLeft1(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                                                     -noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer,                   &
                                                     noahmp%config%domain%JTS:noahmp%config%domain%JTE))

    if (.not. allocated(MatLeft2)) allocate(MatLeft2(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                                                     -noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer,                   &
                                                     noahmp%config%domain%JTS:noahmp%config%domain%JTE))

    if (.not. allocated(MatLeft3)) allocate(MatLeft3(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                                                     -noahmp%config%domain%NumSnowLayerMax+1:noahmp%config%domain%NumSoilLayer,                   &
                                                     noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    MatRight(:,:,:) = 0.0
    MatLeft1(:,:,:) = 0.0
    MatLeft2(:,:,:) = 0.0
    MatLeft3(:,:,:) = 0.0

    !$acc data copyin(MatRight, MatLeft1, MatLeft2, MatLeft3)

    ! compute soil temperatures
    call GlacierThermalDiffusion(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)
    call GlacierTemperatureSolver(noahmp, noahmp%config%domain%MainTimeStep, MatLeft1, MatLeft2, MatLeft3, MatRight)
    !$acc end data

    ! deallocate local arrays to avoid memory leaks
    deallocate(MatRight)
    deallocate(MatLeft1)
    deallocate(MatLeft2)
    deallocate(MatLeft3)

  end subroutine GlacierTemperatureMain

end module GlacierTemperatureMainMod
