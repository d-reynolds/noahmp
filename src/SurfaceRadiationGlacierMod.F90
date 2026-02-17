module SurfaceRadiationGlacierMod

!!! Compute glacier surface radiative fluxes (absorption and reflection)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SurfaceRadiationGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: RADIATION_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                    ! grid indices
    integer                          :: IndBand                 ! waveband indices (1=vis, 2=nir)
    integer                          :: IndLoop                 ! snow and soil layer loop
    real(kind=kind_noahmp)           :: RadSwAbsGrdTmp          ! ground absorbed solar radiation [W/m2]
    real(kind=kind_noahmp)           :: RadSwReflGrdTmp         ! ground reflected solar radiation [W/m2]
    real(kind=kind_noahmp)           :: FracRadSwAbsSnowDirMean ! direct solar flux factor absorbed by snow [frc] scaling
    real(kind=kind_noahmp)           :: FracRadSwAbsSnowDifMean ! diffuse solar flux factor absorbed by snow [frc] scaling

! -----------------------------------------------------------------
    associate(                                                                     &
              NumSwRadBand          => noahmp%config%domain%NumSwRadBand           ,& ! in,  number of solar radiation wave bands
              NumSnowLayerMax       => noahmp%config%domain%NumSnowLayerMax        ,& ! in,  maximum number of snow layers
              NumSnowLayerNeg       => noahmp%config%domain%NumSnowLayerNeg   ,& ! in,  actual number of snow layers (negative)
              OptSnowAlbedo         => noahmp%config%nmlist%OptSnowAlbedo          ,& ! in,  options for ground snow surface albedo
              RadSwDownDir          => noahmp%energy%flux%RadSwDownDir             ,& ! in,  incoming direct solar radiation [W/m2]
              RadSwDownDif          => noahmp%energy%flux%RadSwDownDif             ,& ! in,  incoming diffuse solar radiation [W/m2]
              FracRadSwAbsSnowDir   => noahmp%energy%flux%FracRadSwAbsSnowDir      ,& ! in,  direct solar flux factor absorbed by snow [frc]
              FracRadSwAbsSnowDif   => noahmp%energy%flux%FracRadSwAbsSnowDif      ,& ! in,  diffuse solar flux factor absorbed by snow [frc]
              AlbedoGrdDir          => noahmp%energy%state%AlbedoGrdDir            ,& ! in,  ground albedo (direct beam: vis, nir)
              AlbedoGrdDif          => noahmp%energy%state%AlbedoGrdDif            ,& ! in,  ground albedo (diffuse: vis, nir)
              AlbedoSnowDir         => noahmp%energy%state%AlbedoSnowDir           ,& ! in,  snow albedo for direct(1=vis, 2=nir)
              AlbedoSnowDif         => noahmp%energy%state%AlbedoSnowDif           ,& ! in,  snow albedo for diffuse(1=vis, 2=nir)
              SnowCoverFrac         => noahmp%water%state%SnowCoverFrac       ,& ! in,  snow cover fraction
              AlbedoLandIce         => noahmp%energy%param%AlbedoLandIce           ,& ! in,  albedo land ice: 1=vis, 2=nir
              RadSwAbsGrd           => noahmp%energy%flux%RadSwAbsGrd         ,& ! out, solar radiation absorbed by ground [W/m2]
              RadSwAbsSfc           => noahmp%energy%flux%RadSwAbsSfc         ,& ! out, total absorbed solar radiation [W/m2]
              RadSwReflSfc          => noahmp%energy%flux%RadSwReflSfc        ,& ! out, total reflected solar radiation [W/m2]
              RadSwAbsSnowSoilLayer => noahmp%energy%flux%RadSwAbsSnowSoilLayer     & ! out, total absorbed solar radiation by snow for each layer [W/m2]
             )

   !$acc parallel loop collapse(2) gang vector default(present) private(IndBand,IndLoop,RadSwAbsGrdTmp,RadSwReflGrdTmp,FracRadSwAbsSnowDirMean,FracRadSwAbsSnowDifMean)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! initialization
    RadSwAbsGrd(I,J)  = 0.0
    RadSwAbsSfc(I,J)  = 0.0
    RadSwReflSfc(I,J) = 0.0
    if ( OptSnowAlbedo == 3 ) then
       !$acc loop seq
       do IndLoop = -NumSnowLayerMax+1, 1
          RadSwAbsSnowSoilLayer(I,IndLoop,J) = 0.0
       enddo
    endif

    ! aggretate radiative flux
    !$acc loop seq
    do IndBand = 1, NumSwRadBand
       ! solar radiation absorbed by glacier surface
       RadSwAbsGrdTmp  = RadSwDownDir(I,IndBand,J) * (1.0 - AlbedoGrdDir(I,IndBand,J)) + &
                         RadSwDownDif(I,IndBand,J) * (1.0 - AlbedoGrdDif(I,IndBand,J))
       RadSwAbsGrd(I,J)     = RadSwAbsGrd(I,J) + RadSwAbsGrdTmp
       RadSwAbsSfc(I,J)     = RadSwAbsSfc(I,J) + RadSwAbsGrdTmp
      
       ! solar radiation reflected by glacier surface
       RadSwReflGrdTmp = RadSwDownDir(I,IndBand,J) * AlbedoGrdDir(I,IndBand,J) + &
                         RadSwDownDif(I,IndBand,J) * AlbedoGrdDif(I,IndBand,J)
       RadSwReflSfc(I,J)    = RadSwReflSfc(I,J) + RadSwReflGrdTmp

       if ( OptSnowAlbedo == 3 ) then
          !$acc loop seq
          do IndLoop = -NumSnowLayerMax+1, 1
             FracRadSwAbsSnowDirMean = FracRadSwAbsSnowDir(I,IndLoop,IndBand,J) * SnowCoverFrac(I,J) +    &
                                       ((1.0 - SnowCoverFrac(I,J)) * (1.0 - AlbedoLandIce(I,IndBand,J)) * &
                                       (FracRadSwAbsSnowDir(I,IndLoop,IndBand,J)/(1.0-AlbedoSnowDir(I,IndBand,J))))
             FracRadSwAbsSnowDifMean = FracRadSwAbsSnowDif(I,IndLoop,IndBand,J) * SnowCoverFrac(I,J) +    &
                                       ((1.0 - SnowCoverFrac(I,J)) * (1.0 - AlbedoLandIce(I,IndBand,J)) * &
                                       (FracRadSwAbsSnowDif(I,IndLoop,IndBand,J)/(1.0-AlbedoSnowDif(I,IndBand,J))))
             RadSwAbsSnowSoilLayer(I,IndLoop,J) = RadSwAbsSnowSoilLayer(I,IndLoop,J) + &
                                                  RadSwDownDir(I,IndBand,J) * FracRadSwAbsSnowDirMean +  &
                                                  RadSwDownDif(I,IndBand,J) * FracRadSwAbsSnowDifMean
          enddo
       endif
    enddo

    if (OptSnowAlbedo == 3 .and. NumSnowLayerNeg(I,J) == 0) then
       !$acc loop seq
       do IndLoop = -NumSnowLayerMax+1, 1
          RadSwAbsSnowSoilLayer(I,IndLoop,J) = 0.0
       enddo
       RadSwAbsSnowSoilLayer(I,1,J) = RadSwAbsGrd(I,J)
    endif


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine SurfaceRadiationGlacier

end module SurfaceRadiationGlacierMod
