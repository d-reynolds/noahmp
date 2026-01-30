module SurfaceAlbedoGlacierMod

!!! Compute glacier surface albedo

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowAgingBatsMod,       only : SnowAgingBats
  use SnowAlbedoBatsMod,      only : SnowAlbedoBats
  use SnowAlbedoClassMod,     only : SnowAlbedoClass
  use GroundAlbedoGlacierMod, only : GroundAlbedoGlacier
  use SnowAlbedoSnicarMod,    only : SnowAlbedoSnicar
  use SnowFreshRadiusMod,     only : SnowFreshRadius
  use SnowAgingSnicarMod,     only : SnowAgingSnicar

  implicit none

contains

  subroutine SurfaceAlbedoGlacier(noahmp)

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
    integer                          :: IndBand       ! solar band index
    integer                          :: IndSnow       ! snow layer index
! --------------------------------------------------------------------
   !$acc parallel loop collapse(2) gang vector present(noahmp) private(IndBand)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                                    &
              NumSwRadBand        => noahmp%config%domain%NumSwRadBand            ,& ! in,  number of solar radiation wave bands
              CosSolarZenithAngle => noahmp%config%domain%CosSolarZenithAngle(I,J),& ! in,  cosine solar zenith angle
              OptSnowAlbedo       => noahmp%config%nmlist%OptSnowAlbedo           ,& ! in,  options for ground snow surface albedo
              AlbedoGrdDir        => noahmp%energy%state%AlbedoGrdDir             ,& ! out, ground albedo (direct beam: vis, nir)
              AlbedoGrdDif        => noahmp%energy%state%AlbedoGrdDif             ,& ! out, ground albedo (diffuse: vis, nir)
              AlbedoSnowDir       => noahmp%energy%state%AlbedoSnowDir            ,& ! out, snow albedo for direct(1=vis, 2=nir)
              AlbedoSnowDif       => noahmp%energy%state%AlbedoSnowDif            ,& ! out, snow albedo for diffuse(1=vis, 2=nir)
              AlbedoSfcDir        => noahmp%energy%state%AlbedoSfcDir             ,& ! out, surface albedo (direct)
              AlbedoSfcDif        => noahmp%energy%state%AlbedoSfcDif             ,& ! out, surface albedo (diffuse)
              FracRadSwAbsSnowDir => noahmp%energy%flux%FracRadSwAbsSnowDir       ,& ! out, direct solar flux factor absorbed by snow [frc]
              FracRadSwAbsSnowDif => noahmp%energy%flux%FracRadSwAbsSnowDif        & ! out, diffuse solar flux factor absorbed by snow [frc]
             )
! ----------------------------------------------------------------------

    ! initialization
    !$acc loop seq
    do IndBand = 1, NumSwRadBand
       AlbedoSfcDir (I,IndBand,J) = 0.0
       AlbedoSfcDif (I,IndBand,J) = 0.0
       AlbedoGrdDir (I,IndBand,J) = 0.0
       AlbedoGrdDif (I,IndBand,J) = 0.0
       AlbedoSnowDir(I,IndBand,J) = 0.0
       AlbedoSnowDif(I,IndBand,J) = 0.0
       do IndSnow = -noahmp%config%domain%NumSnowLayerMax+1,1
         FracRadSwAbsSnowDir(I,IndSnow,IndBand,J) = 0.0
         FracRadSwAbsSnowDif(I,IndSnow,IndBand,J) = 0.0
       enddo
    enddo
    end associate

      end do
    end do
   !$acc end parallel loop

    ! snow aging (allow nighttime BATS snow albedo aging)
    if ( noahmp%config%nmlist%OptSnowAlbedo == 1 ) call SnowAgingBats(noahmp)

    ! snow grain size and aging for SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       call SnowFreshRadius(noahmp)
       call SnowAgingSnicar(noahmp)
    endif

      ! snow albedo
      if ( noahmp%config%nmlist%OptSnowAlbedo == 1 ) call SnowAlbedoBats(noahmp)
      if ( noahmp%config%nmlist%OptSnowAlbedo == 2 ) call SnowAlbedoClass(noahmp)
      if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) call SnowAlbedoSnicar(noahmp)

      ! ground albedo
      call GroundAlbedoGlacier(noahmp)

   !$acc parallel loop collapse(2) gang vector present(noahmp) private(IndBand)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                                    &
              NumSwRadBand        => noahmp%config%domain%NumSwRadBand            ,& ! in,  number of solar radiation wave bands
              CosSolarZenithAngle => noahmp%config%domain%CosSolarZenithAngle(I,J),& ! in,  cosine solar zenith angle
              AlbedoGrdDir        => noahmp%energy%state%AlbedoGrdDir             ,& ! out, ground albedo (direct beam: vis, nir)
              AlbedoGrdDif        => noahmp%energy%state%AlbedoGrdDif             ,& ! out, ground albedo (diffuse: vis, nir)
              AlbedoSfcDir        => noahmp%energy%state%AlbedoSfcDir             ,& ! out, surface albedo (direct)
              AlbedoSfcDif        => noahmp%energy%state%AlbedoSfcDif              & ! out, surface albedo (diffuse)
             )
! ----------------------------------------------------------------------
    ! solar radiation process is only done if there is light
    if ( CosSolarZenithAngle > 0 ) then
       ! surface albedo
       !$acc loop seq
       do IndBand = 1, NumSwRadBand
          AlbedoSfcDir(I,IndBand,J) = AlbedoGrdDir(I,IndBand,J)
          AlbedoSfcDif(I,IndBand,J) = AlbedoGrdDif(I,IndBand,J)
       enddo
    endif  ! CosSolarZenithAngle > 0

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine SurfaceAlbedoGlacier

end module SurfaceAlbedoGlacierMod
