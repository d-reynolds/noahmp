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
    associate(                                                             &
              EmissivitySnow   => noahmp%energy%param%EmissivitySnow  ,& ! in,  snow emissivity
              EmissivityIceSfc => noahmp%energy%param%EmissivityIceSfc,& ! in,  emissivity ice surface
              SnowCoverFrac    => noahmp%water%state%SnowCoverFrac    ,& ! in,  snow cover fraction
              EmissivityGrd    => noahmp%energy%state%EmissivityGrd   ,& ! out, ground emissivity
              EmissivitySfc    => noahmp%energy%state%EmissivitySfc    & ! out, surface emissivity
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! ground emissivity
    EmissivityGrd(I,J) = EmissivityIceSfc(I,J) * (1.0 - SnowCoverFrac(I,J)) + EmissivitySnow(I,J) * SnowCoverFrac(I,J)

    ! surface emissivity
    EmissivitySfc(I,J) = EmissivityGrd(I,J)


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine SurfaceEmissivityGlacier

end module SurfaceEmissivityGlacierMod
