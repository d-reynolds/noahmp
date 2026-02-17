module SnowAgingBatsMod

!!! Estimate snow age based on BATS snow albedo scheme for use in BATS snow albedo calculation (2D GPU-optimized)
!!! Reference: Yang et al. (1997) J.of Climate

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowAgingBats(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SNOW_AGE
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                ! grid indices
    real(kind=kind_noahmp)           :: SnowAgeFacTot       ! total aging effects
    real(kind=kind_noahmp)           :: SnowAgeVapEff       ! effects of grain growth due to vapor diffusion
    real(kind=kind_noahmp)           :: SnowAgeFrzEff       ! effects of grain growth at freezing of melt water
    real(kind=kind_noahmp)           :: SnowAgeSootEff      ! effects of soot
    real(kind=kind_noahmp)           :: SnowAgeChg          ! nondimensional snow age change
    real(kind=kind_noahmp)           :: SnowAgeTmp          ! temporary nondimensional snow age
    real(kind=kind_noahmp)           :: SnowFreshFac        ! fresh snowfall factor
    real(kind=kind_noahmp)           :: SnowAgeTimeFac      ! snow aging time factor
    real(kind=kind_noahmp)           :: SnowGrowVapExp      ! snow vapor diffusion growth exponential factor

! --------------------------------------------------------------------
        associate(                                                                     &
                  MainTimeStep         => noahmp%config%domain%MainTimeStep           ,& ! in,    main noahmp timestep [s]
                  SnowMassFullCoverOld => noahmp%water%param%SnowMassFullCoverOld ,& ! in,    new snow mass to fully cover old snow [mm]
                  SnowAgeFacBats       => noahmp%energy%param%SnowAgeFacBats     ,& ! in,    snow aging parameter
                  SnowGrowVapFacBats   => noahmp%energy%param%SnowGrowVapFacBats ,& ! in,    vapor diffusion snow growth factor
                  SnowGrowFrzFacBats   => noahmp%energy%param%SnowGrowFrzFacBats ,& ! in,    extra snow growth factor near freezing
                  SnowSootFacBats      => noahmp%energy%param%SnowSootFacBats    ,& ! in,    dirt and soot effect factor
                  TemperatureGrd       => noahmp%energy%state%TemperatureGrd     ,& ! in,    ground temperature [K]
                  SnowWaterEquiv       => noahmp%water%state%SnowWaterEquiv      ,& ! in,    snow water equivalent [mm]
                  SnowWaterEquivPrev   => noahmp%water%state%SnowWaterEquivPrev  ,& ! in,    snow water equivalent at previous time step [mm]
                  SnowAgeNondim        => noahmp%energy%state%SnowAgeNondim      ,& ! inout, non-dimensional snow age
                  SnowAgeFac           => noahmp%energy%state%SnowAgeFac          & ! out,   snow age factor
                 )

    !$acc parallel loop collapse(2) gang vector default(present) private(SnowAgeChg, SnowAgeFacTot, SnowAgeFrzEff, &
    !$acc SnowAgeSootEff, SnowAgeTimeFac, SnowAgeTmp, SnowAgeVapEff, SnowFreshFac, SnowGrowVapExp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        ! solar radiation process is only done if there is light
        if ( noahmp%config%domain%CosSolarZenithAngle(I,J) <= 0 ) cycle


        if ( SnowWaterEquiv(I,J) <= 0.0 ) then
           SnowAgeNondim(I,J)  = 0.0
        else
           SnowAgeTimeFac = MainTimeStep / SnowAgeFacBats(I,J)
           SnowGrowVapExp = SnowGrowVapFacBats(I,J) * (1.0/ConstFreezePoint - 1.0/TemperatureGrd(I,J))
           SnowAgeVapEff  = exp(SnowGrowVapExp)
           SnowAgeFrzEff  = exp(amin1(0.0, SnowGrowFrzFacBats(I,J)*SnowGrowVapExp))
           SnowAgeSootEff = SnowSootFacBats(I,J)
           SnowAgeFacTot  = SnowAgeVapEff + SnowAgeFrzEff + SnowAgeSootEff
           SnowAgeChg     = SnowAgeTimeFac * SnowAgeFacTot
           SnowFreshFac   = amax1(0.0, SnowWaterEquiv(I,J)-SnowWaterEquivPrev(I,J)) / SnowMassFullCoverOld(I,J)
           SnowAgeTmp     = (SnowAgeNondim(I,J) + SnowAgeChg) * (1.0 - SnowFreshFac)
           SnowAgeNondim(I,J)  = amax1(0.0, SnowAgeTmp)
        endif

        SnowAgeFac(I,J) = SnowAgeNondim(I,J) / (SnowAgeNondim(I,J) + 1.0)


      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine SnowAgingBats

end module SnowAgingBatsMod
