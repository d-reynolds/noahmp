module GroundAlbedoGlacierMod

!!! Compute glacier ground albedo based on snow and ice albedo

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine GroundAlbedoGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: RADIATION_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J          ! grid indices
    integer                          :: IndSwBnd      ! solar radiation band index

! --------------------------------------------------------------------
    associate(                                                           &
              NumSwRadBand  => noahmp%config%domain%NumSwRadBand         ,& ! in,  number of solar radiation wave bands
              SnowCoverFrac => noahmp%water%state%SnowCoverFrac     ,& ! in,  snow cover fraction
              AlbedoLandIce => noahmp%energy%param%AlbedoLandIce         ,& ! in,  albedo land ice: 1=vis, 2=nir
              AlbedoSnowDir => noahmp%energy%state%AlbedoSnowDir         ,& ! in,  snow albedo for direct(1=vis, 2=nir)
              AlbedoSnowDif => noahmp%energy%state%AlbedoSnowDif         ,& ! in,  snow albedo for diffuse(1=vis, 2=nir)
              AlbedoGrdDir  => noahmp%energy%state%AlbedoGrdDir          ,& ! out, ground albedo (direct beam: vis, nir)
              AlbedoGrdDif  => noahmp%energy%state%AlbedoGrdDif          ,& ! out, ground albedo (diffuse: vis, nir)
              AlbedoSfcDir  => noahmp%energy%state%AlbedoSfcDir          ,& ! out, surface albedo (direct beam: vis, nir)
              AlbedoSfcDif  => noahmp%energy%state%AlbedoSfcDif           & ! out, surface albedo (diffuse: vis, nir)
             )

   !$acc parallel loop collapse(2) gang vector default(present) private(IndSwBnd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if (noahmp%config%domain%IndicatorIceSfc(I,J) /= -1) cycle

      ! solar radiation process is only done if there is light
      if ( noahmp%config%domain%CosSolarZenithAngle(I,J) <= 0 ) cycle


    !$acc loop seq
    do IndSwBnd = 1, NumSwRadBand

       AlbedoGrdDir(I,IndSwBnd,J) = AlbedoLandIce(I,IndSwBnd,J)*(1.0-SnowCoverFrac(I,J)) + AlbedoSnowDir(I,IndSwBnd,J)*SnowCoverFrac(I,J)
       AlbedoGrdDif(I,IndSwBnd,J) = AlbedoLandIce(I,IndSwBnd,J)*(1.0-SnowCoverFrac(I,J)) + AlbedoSnowDif(I,IndSwBnd,J)*SnowCoverFrac(I,J)
       AlbedoSfcDir(I,IndSwBnd,J) = AlbedoGrdDir(I,IndSwBnd,J)
       AlbedoSfcDif(I,IndSwBnd,J) = AlbedoGrdDif(I,IndSwBnd,J)

    enddo


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine GroundAlbedoGlacier

end module GroundAlbedoGlacierMod
