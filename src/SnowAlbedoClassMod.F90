module SnowAlbedoClassMod

!!! Compute snow albedo based on the CLASS scheme (Verseghy, 1991) (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowAlbedoClass(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SNOWALB_CLASS
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                  ! grid indices
    real(kind=kind_noahmp)           :: SnowAlbedoTmp         ! temporary snow albedo
    integer                          :: LoopInd    ! loop index for array section expansion

! --------------------------------------------------------------------
        associate(                                                                     &
                  NumSwRadBand         => noahmp%config%domain%NumSwRadBand           ,& ! in,  number of solar radiation wave bands
                  MainTimeStep         => noahmp%config%domain%MainTimeStep           ,& ! in,  noahmp main time step [s]
                  SnowfallGround       => noahmp%water%flux%SnowfallGround       ,& ! in,  snowfall at ground [mm/s]
                  SnowMassFullCoverOld => noahmp%water%param%SnowMassFullCoverOld     ,& ! in,  new snow mass to fully cover old snow [mm]
                  SnowAlbRefClass      => noahmp%energy%param%SnowAlbRefClass         ,& ! in,  reference snow albedo in CLASS scheme
                  SnowAgeFacClass      => noahmp%energy%param%SnowAgeFacClass         ,& ! in,  snow aging e-folding time [s]
                  SnowAlbFreshClass    => noahmp%energy%param%SnowAlbFreshClass       ,& ! in,  fresh snow albedo
                  AlbedoSnowPrev       => noahmp%energy%state%AlbedoSnowPrev     ,& ! in,  snow albedo at last time step
                  AlbedoSnowDir        => noahmp%energy%state%AlbedoSnowDir           ,& ! out, snow albedo for direct (1=vis, 2=nir) (3D)
                  AlbedoSnowDif        => noahmp%energy%state%AlbedoSnowDif            & ! out, snow albedo for diffuse (1=vis, 2=nir) (3D)
                 )

    !$acc parallel loop collapse(2) gang vector default(present) private(SnowAlbedoTmp) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        ! solar radiation process is only done if there is light
        if ( noahmp%config%domain%CosSolarZenithAngle(I,J) <= 0 ) cycle


        ! initialization
        !$acc loop seq
        do LoopInd = 1, NumSwRadBand
           AlbedoSnowDir(I,LoopInd,J) = 0.0
           AlbedoSnowDif(I,LoopInd,J) = 0.0
        enddo

        ! when CosSolarZenithAngle > 0
        SnowAlbedoTmp = SnowAlbRefClass(I,J) + (AlbedoSnowPrev(I,J)-SnowAlbRefClass(I,J)) * exp(-0.01*MainTimeStep/SnowAgeFacClass(I,J))

        ! 1 mm fresh snow(SWE) -- 10mm snow depth, assumed the fresh snow density 100kg/m3
        ! here assume 1cm snow depth will fully cover the old snow
        if ( SnowfallGround(I,J) > 0.0 ) then
           SnowAlbedoTmp = SnowAlbedoTmp + min(SnowfallGround(I,J), SnowMassFullCoverOld(I,J)/MainTimeStep) * &
                                           (SnowAlbFreshClass(I,J)-SnowAlbedoTmp) / (SnowMassFullCoverOld(I,J)/MainTimeStep)
        endif

        AlbedoSnowDif(I,1,J) = SnowAlbedoTmp
        AlbedoSnowDif(I,2,J) = SnowAlbedoTmp
        AlbedoSnowDir(I,1,J) = SnowAlbedoTmp
        AlbedoSnowDir(I,2,J) = SnowAlbedoTmp

        AlbedoSnowPrev(I,J)   = SnowAlbedoTmp


      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine SnowAlbedoClass

end module SnowAlbedoClassMod
