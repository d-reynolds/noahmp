module SurfaceEmissivityGlacierMod

!!! Compute glacier surface longwave emissivity

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SurfaceEmissivityGlacier(noahmp)

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

    associate(                                                             &
              EmissivitySnow   => noahmp%energy%param%EmissivitySnow(I,J)  ,& ! in,  snow emissivity
              EmissivityIceSfc => noahmp%energy%param%EmissivityIceSfc(I,J),& ! in,  emissivity ice surface
              SnowCoverFrac    => noahmp%water%state%SnowCoverFrac(I,J)    ,& ! in,  snow cover fraction
              EmissivityGrd    => noahmp%energy%state%EmissivityGrd(I,J)   ,& ! out, ground emissivity
              EmissivitySfc    => noahmp%energy%state%EmissivitySfc(I,J)    & ! out, surface emissivity
             )
! ----------------------------------------------------------------------

    ! ground emissivity
    EmissivityGrd = EmissivityIceSfc * (1.0 - SnowCoverFrac) + EmissivitySnow * SnowCoverFrac

    ! surface emissivity
    EmissivitySfc = EmissivityGrd

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine SurfaceEmissivityGlacier

end module SurfaceEmissivityGlacierMod
