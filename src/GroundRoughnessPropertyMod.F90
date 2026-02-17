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
        associate(                                                                 &
                  SurfaceType       => noahmp%config%domain%SurfaceType           ,& ! in,  surface type 1-soil; 2-lake
                  RefHeightAboveSfc => noahmp%config%domain%RefHeightAboveSfc     ,& ! in,  reference height [m] above surface zero plane
                  FlagUrban         => noahmp%config%domain%FlagUrban        ,& ! in,  logical flag for urban grid
                  SnowCoverFrac     => noahmp%water%state%SnowCoverFrac      ,& ! in,  snow cover fraction
                  SnowDepth         => noahmp%water%state%SnowDepth          ,& ! in,  snow depth [m]
                  HeightCanopyTop   => noahmp%energy%param%HeightCanopyTop        ,& ! in,  top of canopy [m]
                  RoughLenMomVeg    => noahmp%energy%param%RoughLenMomVeg         ,& ! in,  momentum roughness length vegetated [m]
                  RoughLenMomSnow   => noahmp%energy%param%RoughLenMomSnow        ,& ! in,  snow surface roughness length [m]
                  RoughLenMomSoil   => noahmp%energy%param%RoughLenMomSoil        ,& ! in,  bare-soil roughness length [m]
                  RoughLenMomLake   => noahmp%energy%param%RoughLenMomLake        ,& ! in,  lake surface roughness length [m]
                  TemperatureGrd    => noahmp%energy%state%TemperatureGrd    ,& ! in,  ground temperature [K]
                  RoughLenMomSfc    => noahmp%energy%state%RoughLenMomSfc    ,& ! out, roughness length [m], momentum, surface
                  RoughLenMomGrd    => noahmp%energy%state%RoughLenMomGrd    ,& ! out, roughness length [m], momentum, ground
                  ZeroPlaneDispSfc  => noahmp%energy%state%ZeroPlaneDispSfc  ,& ! out, surface zero plane displacement [m]
                  ZeroPlaneDispGrd  => noahmp%energy%state%ZeroPlaneDispGrd  ,& ! out, ground zero plane displacement [m]
                  RefHeightAboveGrd => noahmp%energy%state%RefHeightAboveGrd  & ! out, reference height [m] above ground
                 )

    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


        ! ground roughness length
        if ( SurfaceType(I,J) == 2 ) then ! Lake
           if ( TemperatureGrd(I,J) <= ConstFreezePoint ) then
              RoughLenMomGrd(I,J) = RoughLenMomLake(I,J) * (1.0-SnowCoverFrac(I,J)) + SnowCoverFrac(I,J) * RoughLenMomSnow(I,J)
           else
              RoughLenMomGrd(I,J) = RoughLenMomLake(I,J)
           endif
        else                         ! soil
           RoughLenMomGrd(I,J)    = RoughLenMomSoil(I,J) * (1.0-SnowCoverFrac(I,J)) + SnowCoverFrac(I,J) * RoughLenMomSnow(I,J)
        endif

        ! surface roughness length and displacement height
        ZeroPlaneDispGrd(I,J)     = SnowDepth(I,J)
        if ( FlagVegSfc(I,J) .eqv. .true. ) then
           RoughLenMomSfc(I,J)    = RoughLenMomVeg(I,J)
           ZeroPlaneDispSfc(I,J)  = 0.65 * HeightCanopyTop(I,J)
           if ( SnowDepth(I,J) > ZeroPlaneDispSfc(I,J) ) ZeroPlaneDispSfc(I,J) = SnowDepth(I,J)
        else
           RoughLenMomSfc(I,J)    = RoughLenMomGrd(I,J)
           ZeroPlaneDispSfc(I,J)  = ZeroPlaneDispGrd(I,J)
        endif

        ! special case for urban
        if ( FlagUrban(I,J) .eqv. .true. ) then
           RoughLenMomGrd(I,J)    = RoughLenMomVeg(I,J)
           ZeroPlaneDispGrd(I,J)  = 0.65 * HeightCanopyTop(I,J)
           RoughLenMomSfc(I,J)    = RoughLenMomGrd(I,J)
           ZeroPlaneDispSfc(I,J)  = ZeroPlaneDispGrd(I,J)
        endif

        ! reference height above ground
        RefHeightAboveGrd(I,J)    = ZeroPlaneDispSfc(I,J) + RefHeightAboveSfc(I,J)
        if ( ZeroPlaneDispGrd(I,J) >= RefHeightAboveGrd(I,J) ) RefHeightAboveGrd(I,J) = ZeroPlaneDispGrd(I,J) + RefHeightAboveSfc(I,J)


      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine GroundRoughnessProperty

end module GroundRoughnessPropertyMod
