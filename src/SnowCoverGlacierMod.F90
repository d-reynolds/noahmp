module SnowCoverGlacierMod

!!! Compute glacier ground snow cover fraction

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowCoverGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in RADIATION_GLACIER subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J      ! grid indices

! --------------------------------------------------------------------
    associate(                                                        &
              SnowWaterEquiv => noahmp%water%state%SnowWaterEquiv,& ! in,  snow water equivalent [mm]
              SnowCoverFrac  => noahmp%water%state%SnowCoverFrac  & ! out, snow cover fraction
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if (noahmp%config%domain%IndicatorIceSfc(I,J) /= -1) cycle

    SnowCoverFrac(I,J) = 0.0
    if ( SnowWaterEquiv(I,J) > 0.0 ) SnowCoverFrac(I,J) = 1.0


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine SnowCoverGlacier

end module SnowCoverGlacierMod
