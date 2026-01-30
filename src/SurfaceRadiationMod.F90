module SurfaceRadiationMod

!!! Compute surface (ground and vegetation) radiative fluxes (absorption and reflection)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SurfaceRadiation(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SURRAD
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                    ! grid indices
    integer                          :: IndBand                 ! waveband indices (1=vis, 2=nir)
    integer                                             :: IndLoop                 ! snow and soil layer loop
    real(kind=kind_noahmp)                              :: MinThr                  ! prevents overflow for division by zero
    real(kind=kind_noahmp)                              :: RadSwAbsGrdTmp          ! ground absorbed solar radiation [W/m2]
    real(kind=kind_noahmp)                              :: RadSwReflSfcNir         ! surface reflected solar radiation NIR [W/m2]
    real(kind=kind_noahmp)                              :: RadSwReflSfcVis         ! surface reflected solar radiation VIS [W/m2]
    real(kind=kind_noahmp)                              :: LeafAreaIndFrac         ! leaf area fraction of canopy
    real(kind=kind_noahmp)                              :: RadSwTranGrdDir         ! transmitted solar radiation at ground: direct [W/m2]
    real(kind=kind_noahmp)                              :: RadSwTranGrdDif         ! transmitted solar radiation at ground: diffuse [W/m2]
    real(kind=kind_noahmp)                              :: RadSwAbsCanDir(1:noahmp%config%domain%NumSwRadBand)          ! direct beam absorbed by canopy [W/m2]
    real(kind=kind_noahmp)                              :: RadSwAbsCanDif(1:noahmp%config%domain%NumSwRadBand)          ! diffuse radiation absorbed by canopy [W/m2]
    real(kind=kind_noahmp)                              :: FracRadSwAbsSnowDirMean(-noahmp%config%domain%NumSnowLayerMax+1:1,1:noahmp%config%domain%NumSwRadBand)  ! direct solar flux factor absorbed by snow [frc] scaling
    real(kind=kind_noahmp)                              :: FracRadSwAbsSnowDifMean(-noahmp%config%domain%NumSnowLayerMax+1:1,1:noahmp%config%domain%NumSwRadBand)  ! diffuse solar flux factor absorbed by snow [frc] scaling


    !$acc parallel loop collapse(2) gang vector present(noahmp) &
    !$acc private(IndBand,IndLoop,MinThr,RadSwAbsGrdTmp,RadSwReflSfcNir,RadSwReflSfcVis) &
    !$acc private(LeafAreaIndFrac,RadSwTranGrdDir,RadSwTranGrdDif,RadSwAbsCanDir) &
    !$acc private(RadSwAbsCanDif) &
    !$acc private(FracRadSwAbsSnowDirMean) &
    !$acc private(FracRadSwAbsSnowDifMean)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
! --------------------------------------------------------------------
    associate(                                                                   &
              NumSwRadBand          => noahmp%config%domain%NumSwRadBand        ,& ! in,  number of solar radiation wave bands
              OptSnowAlbedo         => noahmp%config%nmlist%OptSnowAlbedo       ,& ! in,  options for ground snow surface albedo
              NumSnowLayerMax       => noahmp%config%domain%NumSnowLayerMax     ,& ! in,  maximum number of snow layers 
              NumSnowLayerNeg       => noahmp%config%domain%NumSnowLayerNeg(I,J)     ,& ! in,  actual number of snow layers (negative)
              SnowCoverFrac         => noahmp%water%state%SnowCoverFrac(I,J)         ,& ! in,  snow cover fraction
              LeafAreaIndEff        => noahmp%energy%state%LeafAreaIndEff(I,J)       ,& ! in,  leaf area index, after burying by snow
              VegAreaIndEff         => noahmp%energy%state%VegAreaIndEff(I,J)        ,& ! in,  one-sided leaf+stem area index [m2/m2]
              CanopySunlitFrac      => noahmp%energy%state%CanopySunlitFrac(I,J)     ,& ! in,  sunlit fraction of canopy
              CanopyShadeFrac       => noahmp%energy%state%CanopyShadeFrac(I,J)      ,& ! in,  shaded fraction of canopy
              LeafAreaIndSunlit     => noahmp%energy%state%LeafAreaIndSunlit(I,J)    ,& ! in,  sunlit leaf area
              LeafAreaIndShade      => noahmp%energy%state%LeafAreaIndShade(I,J)     ,& ! in,  shaded leaf area
              AlbedoGrdDir          => noahmp%energy%state%AlbedoGrdDir         ,& ! in,  ground albedo (direct beam: vis, nir)
              AlbedoGrdDif          => noahmp%energy%state%AlbedoGrdDif         ,& ! in,  ground albedo (diffuse: vis, nir)
              AlbedoSfcDir          => noahmp%energy%state%AlbedoSfcDir         ,& ! in,  surface albedo (direct)
              AlbedoSfcDif          => noahmp%energy%state%AlbedoSfcDif         ,& ! in,  surface albedo (diffuse)
              AlbedoSnowDir         => noahmp%energy%state%AlbedoSnowDir        ,& ! in,  snow albedo for direct(1=vis, 2=nir)
              AlbedoSnowDif         => noahmp%energy%state%AlbedoSnowDif        ,& ! in,  snow albedo for diffuse(1=vis, 2=nir)
              AlbedoSoilDir         => noahmp%energy%state%AlbedoSoilDir        ,& ! in,  soil albedo (direct)
              AlbedoSoilDif         => noahmp%energy%state%AlbedoSoilDif        ,& ! in,  soil albedo (diffuse)
              RadSwDownDir          => noahmp%energy%flux%RadSwDownDir          ,& ! in,  incoming direct solar radiation [W/m2]
              RadSwDownDif          => noahmp%energy%flux%RadSwDownDif          ,& ! in,  incoming diffuse solar radiation [W/m2]
              RadSwAbsVegDir        => noahmp%energy%flux%RadSwAbsVegDir        ,& ! in,  flux abs by veg (per unit direct flux)
              RadSwAbsVegDif        => noahmp%energy%flux%RadSwAbsVegDif        ,& ! in,  flux abs by veg (per unit diffuse flux)
              RadSwDirTranGrdDir    => noahmp%energy%flux%RadSwDirTranGrdDir    ,& ! in,  down direct flux below veg (per unit dir flux)
              RadSwDifTranGrdDir    => noahmp%energy%flux%RadSwDifTranGrdDir    ,& ! in,  down diffuse flux below veg (per unit dir flux)
              RadSwDifTranGrdDif    => noahmp%energy%flux%RadSwDifTranGrdDif    ,& ! in,  down diffuse flux below veg (per unit dif flux)
              RadSwReflVegDir       => noahmp%energy%flux%RadSwReflVegDir       ,& ! in,  flux reflected by veg layer (per unit direct flux)
              RadSwReflVegDif       => noahmp%energy%flux%RadSwReflVegDif       ,& ! in,  flux reflected by veg layer (per unit diffuse flux)
              RadSwReflGrdDir       => noahmp%energy%flux%RadSwReflGrdDir       ,& ! in,  flux reflected by ground (per unit direct flux)
              RadSwReflGrdDif       => noahmp%energy%flux%RadSwReflGrdDif       ,& ! in,  flux reflected by ground (per unit diffuse flux)
              FracRadSwAbsSnowDir   => noahmp%energy%flux%FracRadSwAbsSnowDir   ,& ! in,  direct solar flux factor absorbed by snow [frc]
              FracRadSwAbsSnowDif   => noahmp%energy%flux%FracRadSwAbsSnowDif   ,& ! in,  diffuse solar flux factor absorbed by snow [frc]
              RadPhotoActAbsSunlit  => noahmp%energy%flux%RadPhotoActAbsSunlit(I,J)  ,& ! out, average absorbed par for sunlit leaves [W/m2]
              RadPhotoActAbsShade   => noahmp%energy%flux%RadPhotoActAbsShade(I,J)   ,& ! out, average absorbed par for shaded leaves [W/m2]
              RadSwAbsVeg           => noahmp%energy%flux%RadSwAbsVeg(I,J)           ,& ! out, solar radiation absorbed by vegetation [W/m2]
              RadSwAbsGrd           => noahmp%energy%flux%RadSwAbsGrd(I,J)           ,& ! out, solar radiation absorbed by ground [W/m2]
              RadSwAbsSfc           => noahmp%energy%flux%RadSwAbsSfc(I,J)           ,& ! out, total absorbed solar radiation [W/m2]
              RadSwAbsSnowSoilLayer => noahmp%energy%flux%RadSwAbsSnowSoilLayer ,& ! out, total absorbed solar radiation by snow for each layer [W/m2]
              RadSwReflSfc          => noahmp%energy%flux%RadSwReflSfc(I,J)          ,& ! out, total reflected solar radiation [W/m2]
              RadSwReflVeg          => noahmp%energy%flux%RadSwReflVeg(I,J)          ,& ! out, reflected solar radiation by vegetation [W/m2]
              RadSwReflGrd          => noahmp%energy%flux%RadSwReflGrd(I,J)           & ! out, reflected solar radiation by ground [W/m2]
             )
! ----------------------------------------------------------------------

    ! initialization
    MinThr                 = 1.0e-6
    RadSwAbsGrd            = 0.0
    RadSwAbsVeg            = 0.0
    RadSwAbsSfc            = 0.0
    RadSwReflSfc           = 0.0
    RadSwReflVeg           = 0.0
    RadSwReflGrd           = 0.0
    RadPhotoActAbsSunlit   = 0.0
    RadPhotoActAbsShade    = 0.0

    !$acc loop seq
    do IndBand = 1, NumSwRadBand
      RadSwAbsCanDir(IndBand)       = 0.0
      RadSwAbsCanDif(IndBand)       = 0.0
    enddo
    if ( OptSnowAlbedo == 3 ) then
       !$acc loop seq
       do IndLoop = -NumSnowLayerMax+1, 1, 1
          RadSwAbsSnowSoilLayer(I,IndLoop,J) = 0.0
          do IndBand = 1, NumSwRadBand
             FracRadSwAbsSnowDirMean(IndLoop,IndBand) = 0.0
             FracRadSwAbsSnowDifMean(IndLoop,IndBand) = 0.0
          enddo
       enddo
    endif

    ! aggregate radiative flux
    do IndBand = 1, NumSwRadBand
       ! absorbed by canopy
       RadSwAbsCanDir(IndBand) = RadSwDownDir(I,IndBand,J) * RadSwAbsVegDir(I,IndBand,J)
       RadSwAbsCanDif(IndBand) = RadSwDownDif(I,IndBand,J) * RadSwAbsVegDif(I,IndBand,J)
       RadSwAbsVeg             = RadSwAbsVeg + RadSwAbsCanDir(IndBand) + RadSwAbsCanDif(IndBand)
       RadSwAbsSfc             = RadSwAbsSfc + RadSwAbsCanDir(IndBand) + RadSwAbsCanDif(IndBand)
       ! transmitted solar fluxes incident on ground
       RadSwTranGrdDir         = RadSwDownDir(I,IndBand,J) * RadSwDirTranGrdDir(I,IndBand,J)
       RadSwTranGrdDif         = RadSwDownDir(I,IndBand,J) * RadSwDifTranGrdDir(I,IndBand,J) + &
                                 RadSwDownDif(I,IndBand,J) * RadSwDifTranGrdDif(I,IndBand,J)

       ! solar radiation absorbed by ground surface
       RadSwAbsGrdTmp          = RadSwTranGrdDir * (1.0 - AlbedoGrdDir(I,IndBand,J)) + &
                                 RadSwTranGrdDif * (1.0 - AlbedoGrdDif(I,IndBand,J))
       RadSwAbsGrd             = RadSwAbsGrd + RadSwAbsGrdTmp
       RadSwAbsSfc             = RadSwAbsSfc + RadSwAbsGrdTmp

       ! SNICAR snow layer absorption
       if ( OptSnowAlbedo == 3 ) then
          !$acc loop seq
          do IndLoop = -NumSnowLayerMax+1, 1, 1
             FracRadSwAbsSnowDirMean(IndLoop,IndBand) = FracRadSwAbsSnowDir(I,IndLoop,IndBand,J) * SnowCoverFrac +    &
                                                        ((1.0 - SnowCoverFrac) * (1.0 - AlbedoSoilDir(I,IndBand,J)) * &
                                                        (FracRadSwAbsSnowDir(I,IndLoop,IndBand,J)/(1.0 - AlbedoSnowDir(I,IndBand,J))))
             FracRadSwAbsSnowDifMean(IndLoop,IndBand) = FracRadSwAbsSnowDif(I,IndLoop,IndBand,J) * SnowCoverFrac +    &
                                                        ((1.0 - SnowCoverFrac) * (1.0 - AlbedoSoilDif(I,IndBand,J))*  &
                                                        (FracRadSwAbsSnowDif(I,IndLoop,IndBand,J)/(1.0 - AlbedoSnowDif(I,IndBand,J))))
             RadSwAbsSnowSoilLayer(I,IndLoop,J) = RadSwAbsSnowSoilLayer(I,IndLoop,J) + &
                                              RadSwTranGrdDir * FracRadSwAbsSnowDirMean(IndLoop,IndBand) + &
                                              RadSwTranGrdDif * FracRadSwAbsSnowDifMean(IndLoop,IndBand) 
          enddo
       endif
    enddo

    if (OptSnowAlbedo == 3 .and. NumSnowLayerNeg == 0) then
       !$acc loop seq
       do IndLoop = -NumSnowLayerMax+1, 1, 1
          RadSwAbsSnowSoilLayer(I,IndLoop,J) = 0.0
       enddo
       RadSwAbsSnowSoilLayer(I,1,J) = RadSwAbsGrd
    endif

    ! partition visible canopy absorption to sunlit and shaded fractions
    ! to get average absorbed par for sunlit and shaded leaves
    LeafAreaIndFrac = LeafAreaIndEff / max(VegAreaIndEff, MinThr)
    if ( CanopySunlitFrac > 0.0 ) then
       RadPhotoActAbsSunlit = (RadSwAbsCanDir(1) + CanopySunlitFrac * RadSwAbsCanDif(1)) * &
                              LeafAreaIndFrac / max(LeafAreaIndSunlit, MinThr)
       RadPhotoActAbsShade = (CanopyShadeFrac * RadSwAbsCanDif(1)) * &
                              LeafAreaIndFrac / max(LeafAreaIndShade, MinThr)
    else
       RadPhotoActAbsSunlit = 0.0
       RadPhotoActAbsShade  = (RadSwAbsCanDir(1) + RadSwAbsCanDif(1)) * &
                              LeafAreaIndFrac / max(LeafAreaIndShade, MinThr)
    endif

    ! reflected solar radiation
    RadSwReflSfcVis = AlbedoSfcDir(I,1,J) * RadSwDownDir(I,1,J) + AlbedoSfcDif(I,1,J) * RadSwDownDif(I,1,J)
    RadSwReflSfcNir = AlbedoSfcDir(I,2,J) * RadSwDownDir(I,2,J) + AlbedoSfcDif(I,2,J) * RadSwDownDif(I,2,J)
    RadSwReflSfc    = RadSwReflSfcVis + RadSwReflSfcNir

    ! reflected solar radiation of veg. and ground (combined ground)
    RadSwReflVeg = RadSwReflVegDir(I,1,J)*RadSwDownDir(I,1,J) + RadSwReflVegDif(I,1,J)*RadSwDownDif(I,1,J) + &                   
                        RadSwReflVegDir(I,2,J)*RadSwDownDir(I,2,J) + RadSwReflVegDif(I,2,J)*RadSwDownDif(I,2,J)
    RadSwReflGrd = RadSwReflGrdDir(I,1,J)*RadSwDownDir(I,1,J) + RadSwReflGrdDif(I,1,J)*RadSwDownDif(I,1,J) + &
                        RadSwReflGrdDir(I,2,J)*RadSwDownDir(I,2,J) + RadSwReflGrdDif(I,2,J)*RadSwDownDif(I,2,J)

    end associate
      end do
    end do
    !$acc end parallel loop


  end subroutine SurfaceRadiation

end module SurfaceRadiationMod
