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
    associate(                                                           &
              DrainSoilBot     => noahmp%water%flux%DrainSoilBot    ,& ! in,    soil bottom drainage [mm/s]
              RunoffSubsurface => noahmp%water%flux%RunoffSubsurface & ! inout, subsurface runoff [mm/s]
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points


    ! compuate subsurface runoff mm/s
    RunoffSubsurface(I,J) = RunoffSubsurface(I,J) + DrainSoilBot(I,J)


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine RunoffSubSurfaceDrainage

end module RunoffSubSurfaceDrainageMod
