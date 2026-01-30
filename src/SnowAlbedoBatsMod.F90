module SnowAlbedoBatsMod

!!! Compute snow albedo based on BATS scheme (Yang et al. (1997) J.of Climate) (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowAlbedoBats(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SNOWALB_BATS
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                ! grid indices
    real(kind=kind_noahmp)           :: ZenithAngFac        ! solar zenith angle correction factor
    real(kind=kind_noahmp)           :: ZenithAngFacTmp     ! temperary zenith angle correction factor
    real(kind=kind_noahmp)           :: SolarAngleFac2      ! 2.0 * SolarAngleFac
    real(kind=kind_noahmp)           :: SolarAngleFac1      ! 1 / SolarAngleFac
    real(kind=kind_noahmp)           :: SolarAngleFac       ! adjustable solar zenith angle factor
    integer                          :: LoopInd             ! loop index
! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        ! solar radiation process is only done if there is light
        if ( noahmp%config%domain%CosSolarZenithAngle(I,J) < 0 ) cycle

        associate(                                                                     &
                  NumSwRadBand        => noahmp%config%domain%NumSwRadBand            ,& ! in,  number of solar radiation wave bands
                  CosSolarZenithAngle => noahmp%config%domain%CosSolarZenithAngle(I,J),& ! in,  cosine solar zenith angle
                  SolarZenithAdjBats  => noahmp%energy%param%SolarZenithAdjBats(I,J)       ,& ! in,  zenith angle snow albedo adjustment
                  FreshSnoAlbVisBats  => noahmp%energy%param%FreshSnoAlbVisBats(I,J)       ,& ! in,  new snow visible albedo
                  FreshSnoAlbNirBats  => noahmp%energy%param%FreshSnoAlbNirBats(I,J)       ,& ! in,  new snow NIR albedo
                  SnoAgeFacDifVisBats => noahmp%energy%param%SnoAgeFacDifVisBats(I,J)      ,& ! in,  age factor for diffuse visible snow albedo
                  SnoAgeFacDifNirBats => noahmp%energy%param%SnoAgeFacDifNirBats(I,J)      ,& ! in,  age factor for diffuse NIR snow albedo
                  SzaFacDirVisBats    => noahmp%energy%param%SzaFacDirVisBats(I,J)         ,& ! in,  cosz factor for direct visible snow albedo
                  SzaFacDirNirBats    => noahmp%energy%param%SzaFacDirNirBats(I,J)         ,& ! in,  cosz factor for direct NIR snow albedo
                  SnowAgeFac          => noahmp%energy%state%SnowAgeFac(I,J)          ,& ! in,  snow age factor
                  AlbedoSnowDir       => noahmp%energy%state%AlbedoSnowDir            ,& ! out, snow albedo for direct(1=vis, 2=nir) (3D)
                  AlbedoSnowDif       => noahmp%energy%state%AlbedoSnowDif             & ! out, snow albedo for diffuse(1=vis, 2=nir) (3D)
                 )
! ----------------------------------------------------------------------

        ! initialization
        !$acc loop seq
        do LoopInd = 1, NumSwRadBand
          AlbedoSnowDir(I,LoopInd,J) = 0.0
          AlbedoSnowDif(I,LoopInd,J) = 0.0
        enddo

        ! when CosSolarZenithAngle > 0
        SolarAngleFac    = SolarZenithAdjBats
        SolarAngleFac1   = 1.0 / SolarAngleFac
        SolarAngleFac2   = 2.0 * SolarAngleFac
        ZenithAngFacTmp  = (1.0 + SolarAngleFac1) / (1.0 + SolarAngleFac2*CosSolarZenithAngle) - SolarAngleFac1
        ZenithAngFac     = amax1(ZenithAngFacTmp, 0.0)
        AlbedoSnowDif(I,1,J) = FreshSnoAlbVisBats * (1.0 - SnoAgeFacDifVisBats * SnowAgeFac)
        AlbedoSnowDif(I,2,J) = FreshSnoAlbNirBats * (1.0 - SnoAgeFacDifNirBats * SnowAgeFac)
        AlbedoSnowDir(I,1,J) = AlbedoSnowDif(I,1,J) + SzaFacDirVisBats * ZenithAngFac * (1.0 - AlbedoSnowDif(I,1,J))
        AlbedoSnowDir(I,2,J) = AlbedoSnowDif(I,2,J) + SzaFacDirNirBats * ZenithAngFac * (1.0 - AlbedoSnowDif(I,2,J))

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine SnowAlbedoBats

end module SnowAlbedoBatsMod
