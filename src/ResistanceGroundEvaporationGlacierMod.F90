module ResistanceGroundEvaporationGlacierMod

!!! Compute glacier surface resistance to ground evaporation/sublimation
!!! It represents the resistance imposed by the molecular diffusion in
!!! surface (as opposed to aerodynamic resistance computed elsewhere in the model)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine ResistanceGroundEvaporationGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in ENERGY_GLACIER subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J      ! grid indices

! --------------------------------------------------------------------
   !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                               &
              ResistanceGrdEvap => noahmp%energy%state%ResistanceGrdEvap(I,J),& ! out, ground surface resistance [s/m] to evaporation
              RelHumidityGrd    => noahmp%energy%state%RelHumidityGrd(I,J)    & ! out, raltive humidity in surface glacier/snow air space
             )
! ----------------------------------------------------------------------

    ResistanceGrdEvap = 1.0
    RelHumidityGrd    = 1.0

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine ResistanceGroundEvaporationGlacier

end module ResistanceGroundEvaporationGlacierMod
