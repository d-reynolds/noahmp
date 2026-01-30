module GroundRoughnessPropertyGlacierMod

!!! Compute glacier ground roughness length, displacement height, and surface reference height

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine GroundRoughnessPropertyGlacier(noahmp)

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

    associate(                                                                &
              RefHeightAboveSfc => noahmp%config%domain%RefHeightAboveSfc(I,J),& ! in,  reference height [m] above surface zero plane
              SnowDepth         => noahmp%water%state%SnowDepth(I,J)          ,& ! in,  snow depth [m]
              RoughLenMomSnow   => noahmp%energy%param%RoughLenMomSnow(I,J)   ,& ! in,  snow surface roughness length [m]
              RoughLenMomSfc    => noahmp%energy%state%RoughLenMomSfc(I,J)    ,& ! out, roughness length [m], momentum, surface
              RoughLenMomGrd    => noahmp%energy%state%RoughLenMomGrd(I,J)    ,& ! out, roughness length [m], momentum, ground
              ZeroPlaneDispSfc  => noahmp%energy%state%ZeroPlaneDispSfc(I,J)  ,& ! out, surface zero plane displacement [m]
              ZeroPlaneDispGrd  => noahmp%energy%state%ZeroPlaneDispGrd(I,J)  ,& ! out, ground zero plane displacement [m]
              RefHeightAboveGrd => noahmp%energy%state%RefHeightAboveGrd(I,J)  & ! out, reference height [m] above ground
             )
! ----------------------------------------------------------------------

    ! ground roughness length
    RoughLenMomGrd    = RoughLenMomSnow
    RoughLenMomSfc    = RoughLenMomGrd

    ! surface roughness length and displacement height
    ZeroPlaneDispGrd  = SnowDepth
    ZeroPlaneDispSfc  = ZeroPlaneDispGrd

    ! reference height above ground
    RefHeightAboveGrd = ZeroPlaneDispSfc + RefHeightAboveSfc

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine GroundRoughnessPropertyGlacier

end module GroundRoughnessPropertyGlacierMod
