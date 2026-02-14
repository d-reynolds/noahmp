module GroundRoughnessPropertyMod

!!! Compute ground roughness length, displacement height, and surface reference height (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine GroundRoughnessProperty(noahmp, FlagVegSfc)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in ENERGY subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)              , intent(inout) :: noahmp
    logical, dimension(noahmp%config%domain%ITS:noahmp%config%domain%ITE,noahmp%config%domain%JTS:noahmp%config%domain%JTE), intent(in   ) :: FlagVegSfc          ! flag: true if vegetated surface

! local variables
    integer :: I, J  ! grid indices

! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp, FlagVegSfc)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                 &
                  SurfaceType       => noahmp%config%domain%SurfaceType(I,J)           ,& ! in,  surface type 1-soil; 2-lake
                  RefHeightAboveSfc => noahmp%config%domain%RefHeightAboveSfc(I,J)     ,& ! in,  reference height [m] above surface zero plane
                  FlagUrban         => noahmp%config%domain%FlagUrban(I,J)        ,& ! in,  logical flag for urban grid
                  SnowCoverFrac     => noahmp%water%state%SnowCoverFrac(I,J)      ,& ! in,  snow cover fraction
                  SnowDepth         => noahmp%water%state%SnowDepth(I,J)          ,& ! in,  snow depth [m]
                  HeightCanopyTop   => noahmp%energy%param%HeightCanopyTop(I,J)        ,& ! in,  top of canopy [m]
                  RoughLenMomVeg    => noahmp%energy%param%RoughLenMomVeg(I,J)         ,& ! in,  momentum roughness length vegetated [m]
                  RoughLenMomSnow   => noahmp%energy%param%RoughLenMomSnow(I,J)        ,& ! in,  snow surface roughness length [m]
                  RoughLenMomSoil   => noahmp%energy%param%RoughLenMomSoil(I,J)        ,& ! in,  bare-soil roughness length [m]
                  RoughLenMomLake   => noahmp%energy%param%RoughLenMomLake(I,J)        ,& ! in,  lake surface roughness length [m]
                  TemperatureGrd    => noahmp%energy%state%TemperatureGrd(I,J)    ,& ! in,  ground temperature [K]
                  RoughLenMomSfc    => noahmp%energy%state%RoughLenMomSfc(I,J)    ,& ! out, roughness length [m], momentum, surface
                  RoughLenMomGrd    => noahmp%energy%state%RoughLenMomGrd(I,J)    ,& ! out, roughness length [m], momentum, ground
                  ZeroPlaneDispSfc  => noahmp%energy%state%ZeroPlaneDispSfc(I,J)  ,& ! out, surface zero plane displacement [m]
                  ZeroPlaneDispGrd  => noahmp%energy%state%ZeroPlaneDispGrd(I,J)  ,& ! out, ground zero plane displacement [m]
                  RefHeightAboveGrd => noahmp%energy%state%RefHeightAboveGrd(I,J)  & ! out, reference height [m] above ground
                 )
! ----------------------------------------------------------------------

        ! ground roughness length
        if ( SurfaceType == 2 ) then ! Lake
           if ( TemperatureGrd <= ConstFreezePoint ) then
              RoughLenMomGrd = RoughLenMomLake * (1.0-SnowCoverFrac) + SnowCoverFrac * RoughLenMomSnow
           else
              RoughLenMomGrd = RoughLenMomLake
           endif
        else                         ! soil
           RoughLenMomGrd    = RoughLenMomSoil * (1.0-SnowCoverFrac) + SnowCoverFrac * RoughLenMomSnow
        endif

        ! surface roughness length and displacement height
        ZeroPlaneDispGrd     = SnowDepth
        if ( FlagVegSfc(I,J) .eqv. .true. ) then
           RoughLenMomSfc    = RoughLenMomVeg
           ZeroPlaneDispSfc  = 0.65 * HeightCanopyTop
           if ( SnowDepth > ZeroPlaneDispSfc ) ZeroPlaneDispSfc = SnowDepth
        else
           RoughLenMomSfc    = RoughLenMomGrd
           ZeroPlaneDispSfc  = ZeroPlaneDispGrd
        endif

        ! special case for urban
        if ( FlagUrban .eqv. .true. ) then
           RoughLenMomGrd    = RoughLenMomVeg
           ZeroPlaneDispGrd  = 0.65 * HeightCanopyTop
           RoughLenMomSfc    = RoughLenMomGrd
           ZeroPlaneDispSfc  = ZeroPlaneDispGrd
        endif

        ! reference height above ground
        RefHeightAboveGrd    = ZeroPlaneDispSfc + RefHeightAboveSfc
        if ( ZeroPlaneDispGrd >= RefHeightAboveGrd ) RefHeightAboveGrd = ZeroPlaneDispGrd + RefHeightAboveSfc

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine GroundRoughnessProperty

end module GroundRoughnessPropertyMod
