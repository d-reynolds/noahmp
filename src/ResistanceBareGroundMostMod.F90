module ResistanceBareGroundMostMod

!!! Compute bare ground resistance and drag coefficient for momentum and heat
!!! based on Monin-Obukhov (M-O) Similarity Theory (MOST) (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine ResistanceBareGroundMOST(noahmp, IndIter, HeatSensibleTmp, MoStabParaSgn)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SFCDIF1 for bare ground portion
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

! in & out variables
    integer               , intent(in   ) :: IndIter                     ! iteration index
    integer, allocatable, intent(inout   ) :: MoStabParaSgn(:,:)   ! number of times moz changes sign
    real(kind=kind_noahmp), allocatable, intent(in   ) :: HeatSensibleTmp(:,:)   ! temporary effective vegetation area index with constraint (<=6.0)
    type(noahmp_type)     , intent(inout) :: noahmp

! local variable
    real(kind=kind_noahmp)                :: MPE                         ! prevents overflow for division by zero
    real(kind=kind_noahmp)                :: TMPCM                       ! temporary calculation for CM
    real(kind=kind_noahmp)                :: TMPCH                       ! temporary calculation for CH
    real(kind=kind_noahmp)                :: FMNEW                       ! stability correction factor, momentum, for current moz
    real(kind=kind_noahmp)                :: FHNEW                       ! stability correction factor, sen heat, for current moz
    real(kind=kind_noahmp)                :: MOZOLD                      ! Monin-Obukhov stability parameter from prior iteration
    real(kind=kind_noahmp)                :: TMP1,TMP2,TMP3,TMP4,TMP5    ! temporary calculation
    real(kind=kind_noahmp)                :: TVIR                        ! temporary virtual temperature (k)
    real(kind=kind_noahmp)                :: TMPCM2                      ! temporary calculation for CM2
    real(kind=kind_noahmp)                :: TMPCH2                      ! temporary calculation for CH2
    real(kind=kind_noahmp)                :: FM2NEW                      ! stability correction factor, momentum, for current moz
    real(kind=kind_noahmp)                :: FH2NEW                      ! stability correction factor, sen heat, for current moz
    real(kind=kind_noahmp)                :: TMP12,TMP22,TMP32           ! temporary calculation
    real(kind=kind_noahmp)                :: CMFM, CHFH, CM2FM2, CH2FH2  ! temporary calculation
    integer                               :: I, J                         ! grid indices

        associate(                                                                     &
                  TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight ,& ! in,    air temperature [K] at reference height
                  SpecHumidityRefHeight   => noahmp%forcing%SpecHumidityRefHeight ,& ! in,    specific humidity [kg/kg] at reference height
                  RefHeightAboveGrd       => noahmp%energy%state%RefHeightAboveGrd ,& ! in,    reference height [m] above ground
                  DensityAirRefHeight     => noahmp%energy%state%DensityAirRefHeight ,& ! in,    density air [kg/m3]
                  WindSpdRefHeight        => noahmp%energy%state%WindSpdRefHeight ,& ! in,    wind speed [m/s] at reference height
                  ZeroPlaneDispGrd        => noahmp%energy%state%ZeroPlaneDispGrd ,& ! in,    ground zero plane displacement [m]
                  RoughLenShBareGrd       => noahmp%energy%state%RoughLenShBareGrd ,& ! in,    roughness length [m], sensible heat, bare ground
                  RoughLenMomGrd          => noahmp%energy%state%RoughLenMomGrd   ,& ! in,    roughness length [m], momentum, ground
                  MoStabCorrMomBare       => noahmp%energy%state%MoStabCorrMomBare ,& ! inout, M-O momentum stability correction, above ZeroPlaneDisp, bare ground
                  MoStabCorrShBare        => noahmp%energy%state%MoStabCorrShBare  ,& ! inout, M-O sen heat stability correction, above ZeroPlaneDisp, bare ground
                  MoStabCorrMomBare2m     => noahmp%energy%state%MoStabCorrMomBare2m ,& ! inout, M-O momentum stability correction, 2m, bare ground
                  MoStabCorrShBare2m      => noahmp%energy%state%MoStabCorrShBare2m ,& ! inout, M-O sen heat stability correction, 2m, bare ground
                  FrictionVelBare         => noahmp%energy%state%FrictionVelBare   ,& ! inout, friction velocity [m/s], bare ground
                  MoStabParaBare          => noahmp%energy%state%MoStabParaBare    ,& ! inout, Monin-Obukhov stability (z/L), above ZeroPlaneDisp, bare ground
                  MoStabParaBare2m        => noahmp%energy%state%MoStabParaBare2m  ,& ! out,   Monin-Obukhov stability (z/L), 2m, bare ground
                  MoLengthBare            => noahmp%energy%state%MoLengthBare      ,& ! out,   Monin-Obukhov length [m], above ZeroPlaneDisp, bare ground
                  ExchCoeffMomBare        => noahmp%energy%state%ExchCoeffMomBare  ,& ! out,   exchange coeff [m/s] for momentum, above ZeroPlaneDisp, bare ground
                  ExchCoeffShBare         => noahmp%energy%state%ExchCoeffShBare   ,& ! out,   exchange coeff [m/s]  for heat, above ZeroPlaneDisp, bare ground
                  ExchCoeffSh2mBareMo     => noahmp%energy%state%ExchCoeffSh2mBareMo ,& ! out,   exchange coeff [m/s] for heat, 2m, bare ground
                  ResistanceMomBareGrd    => noahmp%energy%state%ResistanceMomBareGrd ,& ! out,   aerodynamic resistance for momentum [s/m], bare ground
                  ResistanceShBareGrd     => noahmp%energy%state%ResistanceShBareGrd ,& ! out,   aerodynamic resistance for sensible heat [s/m], bare ground
                  ResistanceLhBareGrd     => noahmp%energy%state%ResistanceLhBareGrd  & ! out,   aerodynamic resistance for water vapor [s/m], bare ground
                 )

   !$acc parallel loop collapse(2) gang vector default(present) private(MPE, TMPCM, TMPCH, FMNEW, FHNEW, MOZOLD, &
   !$acc TMP1,TMP2,TMP3,TMP4,TMP5, TVIR) private(TMPCM2, TMPCH2, FM2NEW, FH2NEW, TMP12,TMP22,TMP32, CMFM, CHFH, &
   !$acc CM2FM2, CH2FH2 ) firstprivate(IndIter)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! initialization
    MPE    = 1.0e-6
    MOZOLD = MoStabParaBare(I,J)  ! M-O stability parameter for next iteration
    if ( RefHeightAboveGrd(I,J) <= ZeroPlaneDispGrd(I,J) ) then
#ifndef _OPENACC        
         write(*,*) "WARNING: critical problem: RefHeightAboveGrd(I,J) <= ZeroPlaneDispGrd(I,J); model stops"
#endif
       stop "Error in ResistanceBareGroundMostMod.F90"
    endif

    ! temporary drag coefficients
    TMPCM  = log((RefHeightAboveGrd(I,J) - ZeroPlaneDispGrd(I,J)) / RoughLenMomGrd(I,J))
    TMPCH  = log((RefHeightAboveGrd(I,J) - ZeroPlaneDispGrd(I,J)) / RoughLenShBareGrd(I,J))
    TMPCM2 = log((2.0 + RoughLenMomGrd(I,J)) / RoughLenMomGrd(I,J))
    TMPCH2 = log((2.0 + RoughLenShBareGrd(I,J)) / RoughLenShBareGrd(I,J))

    ! compute M-O stability parameter
    if ( IndIter == 1 ) then
       FrictionVelBare(I,J)  = 0.0
       MoStabParaBare(I,J)   = 0.0
       MoLengthBare(I,J)     = 0.0
       MoStabParaBare2m(I,J) = 0.0
    else
       TVIR = (1.0 + 0.61*SpecHumidityRefHeight(I,J)) * TemperatureAirRefHeight(I,J)
       TMP1 = ConstVonKarman * (ConstGravityAcc/TVIR) * HeatSensibleTmp(I,J) / (DensityAirRefHeight(I,J)*ConstHeatCapacAir)
       if ( abs(TMP1) <= MPE ) TMP1 = MPE
       MoLengthBare(I,J)     = -1.0 * FrictionVelBare(I,J)**3 / TMP1
       MoStabParaBare(I,J)   = min((RefHeightAboveGrd(I,J) - ZeroPlaneDispGrd(I,J)) / MoLengthBare(I,J), 1.0)
       MoStabParaBare2m(I,J) = min((2.0 + RoughLenShBareGrd(I,J)) / MoLengthBare(I,J), 1.0)
    endif

    ! accumulate number of times moz changes sign.
    if ( MOZOLD*MoStabParaBare(I,J) < 0.0 ) MoStabParaSgn(I,J) = MoStabParaSgn(I,J) + 1
    if ( MoStabParaSgn(I,J) >= 2 ) then
       MoStabParaBare(I,J)      = 0.0
       MoStabCorrMomBare(I,J)   = 0.0
       MoStabCorrShBare(I,J)    = 0.0
       MoStabParaBare2m(I,J)    = 0.0
       MoStabCorrMomBare2m(I,J) = 0.0
       MoStabCorrShBare2m(I,J)  = 0.0
    endif

    ! evaluate stability-dependent variables using moz from prior iteration
    if ( MoStabParaBare(I,J) < 0.0 ) then
       TMP1   = (1.0 - 16.0 * MoStabParaBare(I,J))**0.25
       TMP2   = log((1.0 + TMP1*TMP1) / 2.0)
       TMP3   = log((1.0 + TMP1) / 2.0)
       FMNEW  = 2.0 * TMP3 + TMP2 - 2.0 * atan(TMP1) + 1.5707963
       FHNEW  = 2 * TMP2
       ! 2-meter quantities
       TMP12  = (1.0 - 16.0 * MoStabParaBare2m(I,J))**0.25
       TMP22  = log((1.0 + TMP12*TMP12) / 2.0)
       TMP32  = log((1.0 + TMP12) / 2.0)
       FM2NEW = 2.0 * TMP32 + TMP22 - 2.0 * atan(TMP12) + 1.5707963
       FH2NEW = 2 * TMP22
    else
       FMNEW  = -5.0 * MoStabParaBare(I,J)
       FHNEW  = FMNEW
       FM2NEW = -5.0 * MoStabParaBare2m(I,J)
       FH2NEW = FM2NEW
    endif

    ! except for first iteration, weight stability factors for previous
    ! iteration to help avoid flip-flops from one iteration to the next
    if ( IndIter == 1 ) then
       MoStabCorrMomBare(I,J)   = FMNEW
       MoStabCorrShBare(I,J)    = FHNEW
       MoStabCorrMomBare2m(I,J) = FM2NEW
       MoStabCorrShBare2m(I,J)  = FH2NEW
    else
       MoStabCorrMomBare(I,J)   = 0.5 * (MoStabCorrMomBare(I,J)   + FMNEW)
       MoStabCorrShBare(I,J)    = 0.5 * (MoStabCorrShBare(I,J)    + FHNEW)
       MoStabCorrMomBare2m(I,J) = 0.5 * (MoStabCorrMomBare2m(I,J) + FM2NEW)
       MoStabCorrShBare2m(I,J)  = 0.5 * (MoStabCorrShBare2m(I,J)  + FH2NEW)
    endif

    ! exchange coefficients
    MoStabCorrShBare(I,J)    = min(MoStabCorrShBare(I,J)   , 0.9*TMPCH )
    MoStabCorrMomBare(I,J)   = min(MoStabCorrMomBare(I,J)  , 0.9*TMPCM )
    MoStabCorrShBare2m(I,J)  = min(MoStabCorrShBare2m(I,J) , 0.9*TMPCH2)
    MoStabCorrMomBare2m(I,J) = min(MoStabCorrMomBare2m(I,J), 0.9*TMPCM2)
    CMFM   = TMPCM  - MoStabCorrMomBare(I,J)
    CHFH   = TMPCH  - MoStabCorrShBare(I,J)
    CM2FM2 = TMPCM2 - MoStabCorrMomBare2m(I,J)
    CH2FH2 = TMPCH2 - MoStabCorrShBare2m(I,J)
    if ( abs(CMFM) <= MPE )   CMFM   = MPE
    if ( abs(CHFH) <= MPE )   CHFH   = MPE
    if ( abs(CM2FM2) <= MPE ) CM2FM2 = MPE
    if ( abs(CH2FH2) <= MPE ) CH2FH2 = MPE
    ExchCoeffMomBare(I,J)    = ConstVonKarman * ConstVonKarman / (CMFM * CMFM)
    ExchCoeffShBare(I,J)     = ConstVonKarman * ConstVonKarman / (CMFM * CHFH)
    !ExchCoeffSh2mBareMo = ConstVonKarman * ConstVonKarman / (CM2FM2 * CH2FH2)

    ! friction velocity
    FrictionVelBare(I,J)     = WindSpdRefHeight(I,J) * sqrt(ExchCoeffMomBare(I,J))
    ExchCoeffSh2mBareMo(I,J) = ConstVonKarman * FrictionVelBare(I,J) / CH2FH2

    ! aerodynamic resistance
    ResistanceMomBareGrd(I,J) = max(1.0, 1.0/(ExchCoeffMomBare(I,J)*WindSpdRefHeight(I,J)))
    ResistanceShBareGrd(I,J)  = max(1.0, 1.0/(ExchCoeffShBare(I,J)*WindSpdRefHeight(I,J)))
    ResistanceLhBareGrd(I,J)  = ResistanceShBareGrd(I,J)

      end do
    end do


        end associate

  end subroutine ResistanceBareGroundMOST

end module ResistanceBareGroundMostMod
