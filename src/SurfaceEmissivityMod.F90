module SurfaceEmissivityMod

!!! Compute ground, vegetation, and total surface longwave emissivity

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SurfaceEmissivity(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in ENERGY subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: I, J         ! grid indices

! --------------------------------------------------------------------
   !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                              &
              IndicatorIceSfc    => noahmp%config%domain%IndicatorIceSfc(I,J)   ,& ! in,  indicator for ice point: 1->seaice; -1->land ice; 0->soil
              SurfaceType        => noahmp%config%domain%SurfaceType(I,J)       ,& ! in,  surface type 1-soil; 2-lake
              EmissivitySnow     => noahmp%energy%param%EmissivitySnow(I,J)     ,& ! in,  snow emissivity
              EmissivitySoilLake => noahmp%energy%param%EmissivitySoilLake ,& ! in,  emissivity soil surface
              EmissivityIceSfc   => noahmp%energy%param%EmissivityIceSfc(I,J)   ,& ! in,  emissivity ice surface
              SnowCoverFrac      => noahmp%water%state%SnowCoverFrac(I,J)       ,& ! in,  snow cover fraction
              LeafAreaIndEff     => noahmp%energy%state%LeafAreaIndEff(I,J)     ,& ! in,  leaf area index, after burying by snow
              StemAreaIndEff     => noahmp%energy%state%StemAreaIndEff(I,J)     ,& ! in,  stem area index, after burying by snow
              VegFrac            => noahmp%energy%state%VegFrac(I,J)            ,& ! in,  greeness vegetation fraction
              EmissivityVeg      => noahmp%energy%state%EmissivityVeg(I,J)      ,& ! out, vegetation emissivity
              EmissivityGrd      => noahmp%energy%state%EmissivityGrd(I,J)      ,& ! out, ground emissivity
              EmissivitySfc      => noahmp%energy%state%EmissivitySfc(I,J)       & ! out, surface emissivity
             )
! ----------------------------------------------------------------------

    if ( IndicatorIceSfc == 0 .or. IndicatorIceSfc == 1 ) then
      ! vegetation emissivity
      EmissivityVeg = 1.0 - exp(-(LeafAreaIndEff + StemAreaIndEff) / 1.0)

      ! ground emissivity
      if ( IndicatorIceSfc == 1 ) then
        EmissivityGrd = EmissivityIceSfc * (1.0-SnowCoverFrac) + EmissivitySnow * SnowCoverFrac
      else
        EmissivityGrd = EmissivitySoilLake(I,SurfaceType,J) * (1.0-SnowCoverFrac) + EmissivitySnow * SnowCoverFrac
      endif

      ! net surface emissivity
      EmissivitySfc = VegFrac * (EmissivityGrd*(1-EmissivityVeg) + EmissivityVeg + &
                      EmissivityVeg*(1-EmissivityVeg)*(1-EmissivityGrd)) + (1-VegFrac) * EmissivityGrd
    else if ( IndicatorIceSfc == -1 ) then ! glacier ice point
      ! ground emissivity
      EmissivityGrd = EmissivityIceSfc * (1.0 - SnowCoverFrac) + EmissivitySnow * SnowCoverFrac

      ! surface emissivity
      EmissivitySfc = EmissivityGrd
    endif
    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine SurfaceEmissivity

end module SurfaceEmissivityMod
