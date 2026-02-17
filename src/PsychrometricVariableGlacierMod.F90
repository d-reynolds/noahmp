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
    associate(                                                                &
              PressureAirRefHeight => noahmp%forcing%PressureAirRefHeight,& ! in,  air pressure [Pa] at reference height
              LatHeatVapGrd        => noahmp%energy%state%LatHeatVapGrd  ,& ! out, latent heat of vaporization/subli [J/kg], ground
              PsychConstGrd        => noahmp%energy%state%PsychConstGrd   & ! out, psychrometric constant [Pa/K], ground
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    LatHeatVapGrd(I,J) = ConstLatHeatSublim
    PsychConstGrd(I,J) = ConstHeatCapacAir * PressureAirRefHeight(I,J) / (0.622 * LatHeatVapGrd(I,J))


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine PsychrometricVariableGlacier

end module PsychrometricVariableGlacierMod
