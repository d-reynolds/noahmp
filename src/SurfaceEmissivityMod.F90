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
    associate(                                                              &
              IndicatorIceSfc    => noahmp%config%domain%IndicatorIceSfc   ,& ! in,  indicator for ice point: 1->seaice; -1->land ice; 0->soil
              SurfaceType        => noahmp%config%domain%SurfaceType       ,& ! in,  surface type 1-soil; 2-lake
              EmissivitySnow     => noahmp%energy%param%EmissivitySnow     ,& ! in,  snow emissivity
              EmissivitySoilLake => noahmp%energy%param%EmissivitySoilLake ,& ! in,  emissivity soil surface
              EmissivityIceSfc   => noahmp%energy%param%EmissivityIceSfc   ,& ! in,  emissivity ice surface
              SnowCoverFrac      => noahmp%water%state%SnowCoverFrac       ,& ! in,  snow cover fraction
              LeafAreaIndEff     => noahmp%energy%state%LeafAreaIndEff     ,& ! in,  leaf area index, after burying by snow
              StemAreaIndEff     => noahmp%energy%state%StemAreaIndEff     ,& ! in,  stem area index, after burying by snow
              VegFrac            => noahmp%energy%state%VegFrac            ,& ! in,  greeness vegetation fraction
              EmissivityVeg      => noahmp%energy%state%EmissivityVeg      ,& ! out, vegetation emissivity
              EmissivityGrd      => noahmp%energy%state%EmissivityGrd      ,& ! out, ground emissivity
              EmissivitySfc      => noahmp%energy%state%EmissivitySfc       & ! out, surface emissivity
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    if ( IndicatorIceSfc(I,J) == 0 .or. IndicatorIceSfc(I,J) == 1 ) then
      ! vegetation emissivity
      EmissivityVeg(I,J) = 1.0 - exp(-(LeafAreaIndEff(I,J) + StemAreaIndEff(I,J)) / 1.0)

      ! ground emissivity
      if ( IndicatorIceSfc(I,J) == 1 ) then
        EmissivityGrd(I,J) = EmissivityIceSfc(I,J) * (1.0-SnowCoverFrac(I,J)) + EmissivitySnow(I,J) * SnowCoverFrac(I,J)
      else
        EmissivityGrd(I,J) = EmissivitySoilLake(I,SurfaceType(I,J),J) * (1.0-SnowCoverFrac(I,J)) + EmissivitySnow(I,J) * SnowCoverFrac(I,J)
      endif

      ! net surface emissivity
      EmissivitySfc(I,J) = VegFrac(I,J) * (EmissivityGrd(I,J)*(1-EmissivityVeg(I,J)) + EmissivityVeg(I,J) + &
                      EmissivityVeg(I,J)*(1-EmissivityVeg(I,J))*(1-EmissivityGrd(I,J))) + (1-VegFrac(I,J)) * EmissivityGrd(I,J)
    else if ( IndicatorIceSfc(I,J) == -1 ) then ! glacier ice point
      ! ground emissivity
      EmissivityGrd(I,J) = EmissivityIceSfc(I,J) * (1.0 - SnowCoverFrac(I,J)) + EmissivitySnow(I,J) * SnowCoverFrac(I,J)

      ! surface emissivity
      EmissivitySfc(I,J) = EmissivityGrd(I,J)
    endif

      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine SurfaceEmissivity

end module SurfaceEmissivityMod
