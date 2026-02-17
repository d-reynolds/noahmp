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
    associate(                                                                &
              RefHeightAboveSfc => noahmp%config%domain%RefHeightAboveSfc,& ! in,  reference height [m] above surface zero plane
              SnowDepth         => noahmp%water%state%SnowDepth          ,& ! in,  snow depth [m]
              RoughLenMomSnow   => noahmp%energy%param%RoughLenMomSnow   ,& ! in,  snow surface roughness length [m]
              RoughLenMomSfc    => noahmp%energy%state%RoughLenMomSfc    ,& ! out, roughness length [m], momentum, surface
              RoughLenMomGrd    => noahmp%energy%state%RoughLenMomGrd    ,& ! out, roughness length [m], momentum, ground
              ZeroPlaneDispSfc  => noahmp%energy%state%ZeroPlaneDispSfc  ,& ! out, surface zero plane displacement [m]
              ZeroPlaneDispGrd  => noahmp%energy%state%ZeroPlaneDispGrd  ,& ! out, ground zero plane displacement [m]
              RefHeightAboveGrd => noahmp%energy%state%RefHeightAboveGrd  & ! out, reference height [m] above ground
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if (noahmp%config%domain%IndicatorIceSfc(I,J) /= -1) cycle


    ! ground roughness length
    RoughLenMomGrd(I,J)    = RoughLenMomSnow(I,J)
    RoughLenMomSfc(I,J)    = RoughLenMomGrd(I,J)

    ! surface roughness length and displacement height
    ZeroPlaneDispGrd(I,J)  = SnowDepth(I,J)
    ZeroPlaneDispSfc(I,J)  = ZeroPlaneDispGrd(I,J)

    ! reference height above ground
    RefHeightAboveGrd(I,J) = ZeroPlaneDispSfc(I,J) + RefHeightAboveSfc(I,J)


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine GroundRoughnessPropertyGlacier

end module GroundRoughnessPropertyGlacierMod
