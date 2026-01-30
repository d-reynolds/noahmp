module PsychrometricVariableGlacierMod

!!! Compute psychrometric variables for glacier ground

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine PsychrometricVariableGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in ENERGY_GLACIER subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J      ! grid indices

! --------------------------------------------------------------------
   !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                                &
              PressureAirRefHeight => noahmp%forcing%PressureAirRefHeight(I,J),& ! in,  air pressure [Pa] at reference height
              LatHeatVapGrd        => noahmp%energy%state%LatHeatVapGrd(I,J)  ,& ! out, latent heat of vaporization/subli [J/kg], ground
              PsychConstGrd        => noahmp%energy%state%PsychConstGrd(I,J)   & ! out, psychrometric constant [Pa/K], ground
             )
! ----------------------------------------------------------------------

    LatHeatVapGrd = ConstLatHeatSublim
    PsychConstGrd = ConstHeatCapacAir * PressureAirRefHeight / (0.622 * LatHeatVapGrd)

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine PsychrometricVariableGlacier

end module PsychrometricVariableGlacierMod
