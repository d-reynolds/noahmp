module RunoffSubSurfaceDrainageMod

!!! Calculate subsurface runoff using derived soil water drainage rate

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine RunoffSubSurfaceDrainage(noahmp)

! ------------------------ Code history --------------------------------------------------
! Originally embeded in WATER subroutine instead of as a separate subroutine
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J      ! grid indices

! --------------------------------------------------------------------
   !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                           &
              DrainSoilBot     => noahmp%water%flux%DrainSoilBot(I,J)    ,& ! in,    soil bottom drainage [mm/s]
              RunoffSubsurface => noahmp%water%flux%RunoffSubsurface(I,J) & ! inout, subsurface runoff [mm/s]
             )
! ----------------------------------------------------------------------

    ! compuate subsurface runoff mm/s
    RunoffSubsurface = RunoffSubsurface + DrainSoilBot

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine RunoffSubSurfaceDrainage

end module RunoffSubSurfaceDrainageMod
