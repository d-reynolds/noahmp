module SoilSnowTemperatureMainMod

!!! Main module to compute snow (if exists) and soil layer temperature. 
!!! Note that snow temperatures during melting season may exceed melting 
!!! point but later in SoilSnowPhaseChange subroutine the snow
!!! temperatures are reset to melting point for melting snow.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SoilSnowTemperatureSolverMod, only : SoilSnowTemperatureSolver
  use SoilSnowThermalDiffusionMod,  only : SoilSnowThermalDiffusion 

  implicit none

contains

  subroutine SoilSnowTemperatureMain(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: TSNOSOI
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)     , intent(inout) :: noahmp

! local variable
    integer                                           :: IndLoop      ! snow and soil layer loop
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatRight     ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatLeft1     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatLeft2     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatLeft3     ! left-hand side term
    integer                                          :: I, J         ! grid indices
    real(kind=kind_noahmp)                           :: RadSwAbsSum  ! temporary variable for sum of absorbed solar radiation by snow and soil layers
    associate(                                                                    &
              NumSoilLayer          => noahmp%config%domain%NumSoilLayer         ,& ! in,  number of soil layers
              NumSnowLayerMax       => noahmp%config%domain%NumSnowLayerMax      ,& ! in,  maximum number of snow layers
              NumSnowLayerNeg       => noahmp%config%domain%NumSnowLayerNeg      ,& ! in,  actual number of snow layers (negative)
              DepthSoilTempBottom   => noahmp%config%domain%DepthSoilTempBottom  ,& ! in,  depth [m] from soil surface for soil temp. lower boundary
              OptSnowAlbedo         => noahmp%config%nmlist%OptSnowAlbedo        ,& ! in,  options for ground snow surface albedo
              SnowDepth             => noahmp%water%state%SnowDepth              ,& ! in,  snow depth [m]
              RadSwAbsSnowSoilLayer => noahmp%energy%flux%RadSwAbsSnowSoilLayer  ,& ! in,  total absorbed solar radiation by snow for each layer [W/m2]
              RadSwAbsGrd           => noahmp%energy%flux%RadSwAbsGrd            ,& ! in,  solar radiation absorbed by ground [W/m2]
              DepthSoilTempBotToSno => noahmp%energy%state%DepthSoilTempBotToSno ,& ! out, depth [m] of soil temp. lower boundary from snow surface
              HeatFromSoilBot       => noahmp%energy%flux%HeatFromSoilBot        ,& ! out, energy influx from soil bottom during soil timestep [J/m2]
              RadSwPenetrateGrd     => noahmp%energy%flux%RadSwPenetrateGrd       & ! out, light penetrating through soil/snow water [W/m2]
             )

    !$acc parallel loop collapse(2) gang vector default(present) private(IndLoop, RadSwAbsSum)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! compute solar penetration through snowpack and soil
    !$acc loop seq
    do IndLoop = -NumSnowLayerMax+1, NumSoilLayer
      RadSwPenetrateGrd(I,IndLoop,J) = 0.0
    enddo

    if (OptSnowAlbedo == 3 .and. NumSnowLayerNeg(I,J) < 0) then
      ! Check if sum > 0 (replacing sum() intrinsic)
      RadSwAbsSum = 0.0
      !$acc loop seq
      do IndLoop = lbound(RadSwAbsSnowSoilLayer,2), ubound(RadSwAbsSnowSoilLayer,2)
         RadSwAbsSum = RadSwAbsSum + RadSwAbsSnowSoilLayer(I,IndLoop,J)
      enddo
          if (RadSwAbsSum > 0.0) then
          !$acc loop seq
          do IndLoop = NumSnowLayerNeg(I,J)+1, 1, 1
             if (IndLoop == NumSnowLayerNeg(I,J)+1) then
                RadSwPenetrateGrd(I,IndLoop,J) = RadSwAbsSnowSoilLayer(I,IndLoop,J) - RadSwAbsGrd(I,J) 
             else
                RadSwPenetrateGrd(I,IndLoop,J) = RadSwAbsSnowSoilLayer(I,IndLoop,J)
             endif
          enddo
       endif
    endif

    ! adjust DepthSoilTempBottom from soil surface to DepthSoilTempBotToSno from snow surface
    DepthSoilTempBotToSno(I,J) = DepthSoilTempBottom(I,J) - SnowDepth(I,J)


      enddo
   enddo

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
    call SoilSnowThermalDiffusion(noahmp, noahmp%config%domain%SoilTimeStep, MatLeft1, MatLeft2, MatLeft3, MatRight)
    call SoilSnowTemperatureSolver(noahmp, noahmp%config%domain%SoilTimeStep, MatLeft1, MatLeft2, MatLeft3, MatRight)

    !$acc end data

    ! deallocate local arrays to avoid memory leaks
    deallocate(MatRight)
    deallocate(MatLeft1)
    deallocate(MatLeft2)
    deallocate(MatLeft3)



    end associate

  end subroutine SoilSnowTemperatureMain

end module SoilSnowTemperatureMainMod
