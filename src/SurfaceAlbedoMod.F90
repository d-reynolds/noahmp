module SurfaceAlbedoMod

!!! Compute total surface albedo and vegetation radiative fluxes 
!!! per unit incoming direct and diffuse radiation and sunlit fraction of canopy

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowAgingBatsMod,            only : SnowAgingBats
  use SnowAlbedoBatsMod,           only : SnowAlbedoBats
  use SnowAlbedoClassMod,          only : SnowAlbedoClass
  use SnowAlbedoSnicarMod,         only : SnowAlbedoSnicar
  use GroundAlbedoMod,             only : GroundAlbedo
  use CanopyRadiationTwoStreamMod, only : CanopyRadiationTwoStream
  use SnowFreshRadiusMod,          only : SnowFreshRadius
  use SnowAgingSnicarMod,          only : SnowAgingSnicar

  implicit none

contains

  subroutine SurfaceAlbedo(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ALBEDO
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J            ! grid indices
    integer                          :: IndSnow         ! snow layer index
    integer                          :: IndBand         ! waveband indices
    integer                          :: IndDif          ! direct beam: IndDif=0; diffuse: IndDif=1
    real(kind=kind_noahmp)           :: LeafWgt         ! fraction of LeafAreaIndex+StemAreaIndex that is LeafAreaIndex
    real(kind=kind_noahmp)           :: StemWgt         ! fraction of LeafAreaIndex+StemAreaIndex that is StemAreaIndex
    real(kind=kind_noahmp)           :: MinThr          ! prevents overflow for division by zero
    real(kind=kind_noahmp)           :: LightExtDir     ! optical depth direct beam per unit leaf + stem area

    MinThr           = 1.0e-06

! --------------------------------------------------------------------
   !$acc parallel loop collapse(2) gang vector present(noahmp) &
   !$acc private(IndBand, IndDif, LeafWgt, StemWgt, MinThr, LightExtDir)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                                 &
              NumSwRadBand        => noahmp%config%domain%NumSwRadBand             ,& ! in,  number of solar radiation wave bands
              CosSolarZenithAngle => noahmp%config%domain%CosSolarZenithAngle(I,J) ,& ! in,  cosine solar zenith angle
              OptSnowAlbedo       => noahmp%config%nmlist%OptSnowAlbedo       ,& ! in,  options for ground snow surface albedo
              ReflectanceLeaf     => noahmp%energy%param%ReflectanceLeaf      ,& ! in,  leaf reflectance: 1=vis, 2=nir
              ReflectanceStem     => noahmp%energy%param%ReflectanceStem      ,& ! in,  stem reflectance: 1=vis, 2=nir
              TransmittanceLeaf   => noahmp%energy%param%TransmittanceLeaf    ,& ! in,  leaf transmittance: 1=vis, 2=nir
              TransmittanceStem   => noahmp%energy%param%TransmittanceStem    ,& ! in,  stem transmittance: 1=vis, 2=nir
              LeafAreaIndEff      => noahmp%energy%state%LeafAreaIndEff(I,J)       ,& ! in,  leaf area index, after burying by snow
              StemAreaIndEff      => noahmp%energy%state%StemAreaIndEff(I,J)       ,& ! in,  stem area index, after burying by snow
              AlbedoGrdDir        => noahmp%energy%state%AlbedoGrdDir         ,& ! out, ground albedo (direct beam: vis, nir)
              AlbedoGrdDif        => noahmp%energy%state%AlbedoGrdDif         ,& ! out, ground albedo (diffuse: vis, nir)
              AlbedoSnowDir       => noahmp%energy%state%AlbedoSnowDir        ,& ! out, snow albedo for direct(1=vis, 2=nir)
              AlbedoSnowDif       => noahmp%energy%state%AlbedoSnowDif        ,& ! out, snow albedo for diffuse(1=vis, 2=nir)
              AlbedoSfcDir        => noahmp%energy%state%AlbedoSfcDir         ,& ! out, surface albedo (direct)
              AlbedoSfcDif        => noahmp%energy%state%AlbedoSfcDif         ,& ! out, surface albedo (diffuse)
              CanopySunlitFrac    => noahmp%energy%state%CanopySunlitFrac(I,J)     ,& ! out, sunlit fraction of canopy
              CanopyShadeFrac     => noahmp%energy%state%CanopyShadeFrac(I,J)      ,& ! out, shaded fraction of canopy
              LeafAreaIndSunlit   => noahmp%energy%state%LeafAreaIndSunlit(I,J)    ,& ! out, sunlit leaf area
              LeafAreaIndShade    => noahmp%energy%state%LeafAreaIndShade(I,J)     ,& ! out, shaded leaf area
              GapBtwCanopy        => noahmp%energy%state%GapBtwCanopy(I,J)         ,& ! out, between canopy gap fraction for beam
              GapInCanopy         => noahmp%energy%state%GapInCanopy(I,J)          ,& ! out, within canopy gap fraction for beam
              ReflectanceVeg      => noahmp%energy%state%ReflectanceVeg       ,& ! out, leaf/stem reflectance weighted by fraction LAI and SAI
              TransmittanceVeg    => noahmp%energy%state%TransmittanceVeg     ,& ! out, leaf/stem transmittance weighted by fraction LAI and SAI
              VegAreaIndEff       => noahmp%energy%state%VegAreaIndEff(I,J)        ,& ! out, one-sided leaf+stem area index [m2/m2]
              VegAreaProjDir      => noahmp%energy%state%VegAreaProjDir(I,J)       ,& ! out, projected leaf+stem area in solar direction
              FracRadSwAbsSnowDir => noahmp%energy%flux%FracRadSwAbsSnowDir   ,& ! out, direct solar flux factor absorbed by snow [frc]
              FracRadSwAbsSnowDif => noahmp%energy%flux%FracRadSwAbsSnowDif   ,& ! out, diffuse solar flux factor absorbed by snow [frc]
              RadSwAbsVegDir      => noahmp%energy%flux%RadSwAbsVegDir        ,& ! out, flux abs by veg (per unit direct flux)
              RadSwAbsVegDif      => noahmp%energy%flux%RadSwAbsVegDif        ,& ! out, flux abs by veg (per unit diffuse flux)
              RadSwDirTranGrdDir  => noahmp%energy%flux%RadSwDirTranGrdDir    ,& ! out, down direct flux below veg (per unit dir flux)
              RadSwDifTranGrdDir  => noahmp%energy%flux%RadSwDifTranGrdDir    ,& ! out, down diffuse flux below veg (per unit dir flux)
              RadSwDifTranGrdDif  => noahmp%energy%flux%RadSwDifTranGrdDif    ,& ! out, down diffuse flux below veg (per unit dif flux)
              RadSwDirTranGrdDif  => noahmp%energy%flux%RadSwDirTranGrdDif    ,& ! out, down direct flux below veg per unit dif flux (= 0)
              RadSwReflVegDir     => noahmp%energy%flux%RadSwReflVegDir       ,& ! out, flux reflected by veg layer (per unit direct flux)
              RadSwReflVegDif     => noahmp%energy%flux%RadSwReflVegDif       ,& ! out, flux reflected by veg layer (per unit diffuse flux)
              RadSwReflGrdDir     => noahmp%energy%flux%RadSwReflGrdDir       ,& ! out, flux reflected by ground (per unit direct flux)
              RadSwReflGrdDif     => noahmp%energy%flux%RadSwReflGrdDif        & ! out, flux reflected by ground (per unit diffuse flux)
             )
! ----------------------------------------------------------------------

    ! initialization
    GapBtwCanopy     = 0.0
    GapInCanopy      = 0.0
    VegAreaProjDir   = 0.0
    ReflectanceVeg   = 0.0
    TransmittanceVeg = 0.0
    CanopySunlitFrac = 0.0
    VegAreaIndEff = LeafAreaIndEff + StemAreaIndEff

    !$acc loop seq
    do IndBand = 1, NumSwRadBand
       AlbedoSfcDir      (I,IndBand,J) = 0.0
       AlbedoSfcDif      (I,IndBand,J) = 0.0
       AlbedoGrdDir      (I,IndBand,J) = 0.0
       AlbedoGrdDif      (I,IndBand,J) = 0.0
       AlbedoSnowDir     (I,IndBand,J) = 0.0
       AlbedoSnowDif     (I,IndBand,J) = 0.0
       RadSwAbsVegDir    (I,IndBand,J) = 0.0
       RadSwAbsVegDif    (I,IndBand,J) = 0.0
       RadSwDirTranGrdDir(I,IndBand,J) = 0.0
       RadSwDirTranGrdDif(I,IndBand,J) = 0.0
       RadSwDifTranGrdDir(I,IndBand,J) = 0.0
       RadSwDifTranGrdDif(I,IndBand,J) = 0.0
       RadSwReflVegDir   (I,IndBand,J) = 0.0
       RadSwReflVegDif   (I,IndBand,J) = 0.0
       RadSwReflGrdDir   (I,IndBand,J) = 0.0
       RadSwReflGrdDif   (I,IndBand,J) = 0.0
       do IndSnow = -2, 1
         FracRadSwAbsSnowDir(I,IndSnow,IndBand,J) = 0.0
         FracRadSwAbsSnowDif(I,IndSnow,IndBand,J) = 0.0
       enddo
    enddo

    ! solar radiation process is only done if there is light
    if ( CosSolarZenithAngle > 0 ) then

       ! weight reflectance/transmittance by LeafAreaIndex and StemAreaIndex
       LeafWgt = LeafAreaIndEff / max(VegAreaIndEff, MinThr)
       StemWgt = StemAreaIndEff / max(VegAreaIndEff, MinThr)
       !$acc loop seq
       do IndBand = 1, NumSwRadBand
          ReflectanceVeg(I,IndBand,J)   = max(ReflectanceLeaf(I,IndBand,J)*LeafWgt+ReflectanceStem(I,IndBand,J)*StemWgt, MinThr)
          TransmittanceVeg(I,IndBand,J) = max(TransmittanceLeaf(I,IndBand,J)*LeafWgt+TransmittanceStem(I,IndBand,J)*StemWgt, MinThr)
       enddo

    endif
      end associate
      enddo
   enddo
   !$acc end parallel loop

    ! snow aging (allow nighttime BATS snow albedo aging)
    if ( noahmp%config%nmlist%OptSnowAlbedo == 1 ) call SnowAgingBats(noahmp)

    ! snow grain size and aging for SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       call SnowFreshRadius(noahmp)
       call SnowAgingSnicar(noahmp)
    endif

      ! snow albedos
      if ( noahmp%config%nmlist%OptSnowAlbedo == 1 ) call SnowAlbedoBats(noahmp)
      if ( noahmp%config%nmlist%OptSnowAlbedo == 2 ) call SnowAlbedoClass(noahmp)
      if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) call SnowAlbedoSnicar(noahmp)

      ! ground surface albedo
      call GroundAlbedo(noahmp)

      ! loop over shortwave bands to calculate surface albedos and solar
      ! fluxes for unit incoming direct (IndDif=0) and diffuse flux (IndDif=1)
      do IndBand = 1, noahmp%config%domain%NumSwRadBand
         IndDif = 0      ! direct
         call CanopyRadiationTwoStream(noahmp, IndBand, IndDif)
         IndDif = 1      ! diffuse
         call CanopyRadiationTwoStream(noahmp, IndBand, IndDif)
      enddo


   !$acc parallel loop collapse(2) gang vector present(noahmp) &
   !$acc private(IndBand, IndDif, LeafWgt, StemWgt, MinThr, LightExtDir)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                                 &
              CosSolarZenithAngle => noahmp%config%domain%CosSolarZenithAngle(I,J) ,& ! in,  cosine solar zenith angle
              LeafAreaIndEff      => noahmp%energy%state%LeafAreaIndEff(I,J)       ,& ! in,  leaf area index, after burying by snow
              CanopySunlitFrac    => noahmp%energy%state%CanopySunlitFrac(I,J)     ,& ! out, sunlit fraction of canopy
              CanopyShadeFrac     => noahmp%energy%state%CanopyShadeFrac(I,J)      ,& ! out, shaded fraction of canopy
              LeafAreaIndSunlit   => noahmp%energy%state%LeafAreaIndSunlit(I,J)    ,& ! out, sunlit leaf area
              LeafAreaIndShade    => noahmp%energy%state%LeafAreaIndShade(I,J)     ,& ! out, shaded leaf area
              ReflectanceVeg      => noahmp%energy%state%ReflectanceVeg       ,& ! out, leaf/stem reflectance weighted by fraction LAI and SAI
              TransmittanceVeg    => noahmp%energy%state%TransmittanceVeg     ,& ! out, leaf/stem transmittance weighted by fraction LAI and SAI
              VegAreaIndEff       => noahmp%energy%state%VegAreaIndEff(I,J)        ,& ! out, one-sided leaf+stem area index [m2/m2]
              VegAreaProjDir      => noahmp%energy%state%VegAreaProjDir(I,J)        & ! out, projected leaf+stem area in solar direction
             )

    ! solar radiation process is only done if there is light
    if ( CosSolarZenithAngle > 0 ) then
       ! sunlit fraction of canopy. set CanopySunlitFrac = 0 if CanopySunlitFrac < 0.01.
       LightExtDir      = VegAreaProjDir / CosSolarZenithAngle * sqrt(1.0-ReflectanceVeg(I,1,J)-TransmittanceVeg(I,1,J))
       CanopySunlitFrac = (1.0 - exp(-LightExtDir*VegAreaIndEff)) / max(LightExtDir*VegAreaIndEff, MinThr)
       LightExtDir      = CanopySunlitFrac
       if ( LightExtDir < 0.01 ) then
          LeafWgt = 0.0
       else
          LeafWgt = LightExtDir
       endif
       CanopySunlitFrac = LeafWgt

    endif  ! CosSolarZenithAngle > 0

    ! shaded canopy fraction
    CanopyShadeFrac   = 1.0 - CanopySunlitFrac
    LeafAreaIndSunlit = LeafAreaIndEff * CanopySunlitFrac
    LeafAreaIndShade  = LeafAreaIndEff * CanopyShadeFrac

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine SurfaceAlbedo

end module SurfaceAlbedoMod
