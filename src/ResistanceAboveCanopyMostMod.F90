module ResistanceAboveCanopyMostMod

!!! Compute surface resistance and drag coefficient for momentum and heat
!!! based on Monin-Obukhov (M-O) Similarity Theory (MOST) (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine ResistanceAboveCanopyMOST(noahmp, IterationInd, HeatSensibleTmp, MoStabParaSgn)
   !$acc routine gang
! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SFCDIF1 for vegetated portion
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    integer               , intent(in   ) :: IterationInd                ! iteration index
    integer, allocatable, intent(inout   ) :: MoStabParaSgn(:,:)   ! number of times moz changes sign
    real(kind=kind_noahmp), allocatable, intent(in   ) :: HeatSensibleTmp(:,:)   ! temporary effective vegetation area index with constraint (<=6.0)
    type(noahmp_type)     , intent(inout) :: noahmp

! local variable
    integer                               :: I, J            ! grid indices
    real(kind=kind_noahmp)                :: MPE             ! prevents overflow for division by zero
    real(kind=kind_noahmp)                :: TMPCM           ! temporary calculation for ExchCoeffMomAbvCan
    real(kind=kind_noahmp)                :: TMPCH           ! temporary calculation for CH
    real(kind=kind_noahmp)                :: FMNEW           ! stability correction factor, momentum, for current moz
    real(kind=kind_noahmp)                :: FHNEW           ! stability correction factor, sen heat, for current moz
    real(kind=kind_noahmp)                :: MOZOLD          ! Monin-Obukhov stability parameter from prior iteration
    real(kind=kind_noahmp)                :: TMP1,TMP2,TMP3,TMP4,TMP5 ! temporary calculation
    real(kind=kind_noahmp)                :: TVIR            ! temporary virtual temperature [K]
    real(kind=kind_noahmp)                :: TMPCM2          ! temporary calculation for CM2
    real(kind=kind_noahmp)                :: TMPCH2          ! temporary calculation for CH2
    real(kind=kind_noahmp)                :: FM2NEW          ! stability correction factor, momentum, for current moz
    real(kind=kind_noahmp)                :: FH2NEW          ! stability correction factor, sen heat, for current moz
    real(kind=kind_noahmp)                :: TMP12,TMP22,TMP32 ! temporary calculation
    real(kind=kind_noahmp)                :: CMFM, CHFH, CM2FM2, CH2FH2 ! temporary calculation

! --------------------------------------------------------------------
        associate(                                                                    &
                  TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight ,& ! in,    air temperature [K] at reference height
                  SpecHumidityRefHeight   => noahmp%forcing%SpecHumidityRefHeight ,& ! in,    specific humidity [kg/kg] at reference height
                  RefHeightAboveGrd       => noahmp%energy%state%RefHeightAboveGrd ,& ! in,    reference height [m] above ground
                  DensityAirRefHeight     => noahmp%energy%state%DensityAirRefHeight ,& ! in,    density air [kg/m3]
                  WindSpdRefHeight        => noahmp%energy%state%WindSpdRefHeight ,& ! in,    wind speed [m/s] at reference height
                  ZeroPlaneDispSfc        => noahmp%energy%state%ZeroPlaneDispSfc ,& ! in,    zero plane displacement [m]
                  RoughLenShCanopy        => noahmp%energy%state%RoughLenShCanopy ,& ! in,    roughness length [m], sensible heat, vegetated
                  RoughLenMomSfc          => noahmp%energy%state%RoughLenMomSfc   ,& ! in,    roughness length [m], momentum, surface
                  MoStabCorrMomAbvCan     => noahmp%energy%state%MoStabCorrMomAbvCan ,& ! inout, M-O momentum stability correction, above ZeroPlaneDispSfc, vegetated
                  MoStabCorrShAbvCan      => noahmp%energy%state%MoStabCorrShAbvCan ,& ! inout, M-O sen heat stability correction, above ZeroPlaneDispSfc, vegetated
                  MoStabCorrMomVeg2m      => noahmp%energy%state%MoStabCorrMomVeg2m ,& ! inout, M-O momentum stability correction, 2m, vegetated
                  MoStabCorrShVeg2m       => noahmp%energy%state%MoStabCorrShVeg2m  ,& ! inout, M-O sen heat stability correction, 2m, vegetated
                  MoStabParaAbvCan        => noahmp%energy%state%MoStabParaAbvCan   ,& ! inout, Monin-Obukhov stability (z/L), above ZeroPlaneDispSfc, vegetated
                  FrictionVelVeg          => noahmp%energy%state%FrictionVelVeg     ,& ! inout, friction velocity [m/s], vegetated
                  MoStabParaVeg2m         => noahmp%energy%state%MoStabParaVeg2m    ,& ! out,   Monin-Obukhov stability (z/L), 2m, vegetated
                  MoLengthAbvCan          => noahmp%energy%state%MoLengthAbvCan     ,& ! out,   Monin-Obukhov length [m], above ZeroPlaneDispSfc, vegetated
                  ExchCoeffMomAbvCan      => noahmp%energy%state%ExchCoeffMomAbvCan ,& ! out,   drag coefficient for momentum, above ZeroPlaneDispSfc, vegetated
                  ExchCoeffShAbvCan       => noahmp%energy%state%ExchCoeffShAbvCan  ,& ! out,   exchange coefficient for heat, above ZeroPlaneDispSfc, vegetated
                  ExchCoeffSh2mVegMo      => noahmp%energy%state%ExchCoeffSh2mVegMo ,& ! out,   exchange coefficient for heat, 2m, vegetated
                  ResistanceMomAbvCan     => noahmp%energy%state%ResistanceMomAbvCan ,& ! out,   aerodynamic resistance for momentum [s/m], above canopy
                  ResistanceShAbvCan      => noahmp%energy%state%ResistanceShAbvCan ,& ! out,   aerodynamic resistance for sensible heat [s/m], above canopy
                  ResistanceLhAbvCan      => noahmp%energy%state%ResistanceLhAbvCan  & ! out,   aerodynamic resistance for water vapor [s/m], above canopy
                 )

    !$acc loop collapse(2) gang vector private(MPE, TMPCM, TMPCH, FMNEW, FHNEW, MOZOLD, &
    !$acc TMP1, TMP2, TMP3, TMP4, TMP5, TVIR, TMPCM2, TMPCH2, FM2NEW, FH2NEW, TMP12, TMP22, TMP32, CMFM, CHFH, CM2FM2, &
    !$acc CH2FH2)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if ( .not. ((noahmp%energy%state%VegAreaIndEff(I,J) > 0.0 ) .and. (noahmp%energy%state%VegFrac(I,J) > 0)) ) cycle ! skip non-vegetated surface


        ! initialization
        MPE    = 1.0e-6
        MOZOLD = MoStabParaAbvCan(I,J)  ! M-O stability parameter for next iteration
        if ( RefHeightAboveGrd(I,J) <= ZeroPlaneDispSfc(I,J) ) then
#ifndef _OPENACC        
         write(*,*) "WARNING: critical problem: RefHeightAboveGrd(I,J) <= ZeroPlaneDispSfc(I,J); model stops"
#endif
         stop "Error in ResistanceAboveCanopyMostMod.F90"
        endif

        ! temporary drag coefficients
        TMPCM  = log((RefHeightAboveGrd(I,J) - ZeroPlaneDispSfc(I,J)) / RoughLenMomSfc(I,J))
        TMPCH  = log((RefHeightAboveGrd(I,J) - ZeroPlaneDispSfc(I,J)) / RoughLenShCanopy(I,J))
        TMPCM2 = log((2.0 + RoughLenMomSfc(I,J)) / RoughLenMomSfc(I,J))
        TMPCH2 = log((2.0 + RoughLenShCanopy(I,J)) / RoughLenShCanopy(I,J))

        ! compute M-O stability parameter
        if ( IterationInd == 1 ) then
           FrictionVelVeg(I,J)   = 0.0
           MoStabParaAbvCan(I,J) = 0.0
           MoLengthAbvCan(I,J)   = 0.0
           MoStabParaVeg2m(I,J)  = 0.0
        else
           TVIR = (1.0 + 0.61*SpecHumidityRefHeight(I,J)) * TemperatureAirRefHeight(I,J)
           TMP1 = ConstVonKarman * (ConstGravityAcc/TVIR) * HeatSensibleTmp(I,J) / (DensityAirRefHeight(I,J)*ConstHeatCapacAir)
           if ( abs(TMP1) <= MPE ) TMP1 = MPE
           MoLengthAbvCan(I,J)   = -1.0 * FrictionVelVeg(I,J)**3 / TMP1
           MoStabParaAbvCan(I,J) = min((RefHeightAboveGrd(I,J) - ZeroPlaneDispSfc(I,J)) / MoLengthAbvCan(I,J), 1.0)
           MoStabParaVeg2m(I,J)  = min((2.0 + RoughLenShCanopy(I,J)) / MoLengthAbvCan(I,J), 1.0)
        endif

        ! accumulate number of times moz changes sign.
        if ( MOZOLD*MoStabParaAbvCan(I,J) < 0.0 ) MoStabParaSgn(I,J) = MoStabParaSgn(I,J) + 1
        if ( MoStabParaSgn(I,J) >= 2 ) then
           MoStabParaAbvCan(I,J)    = 0.0
           MoStabCorrMomAbvCan(I,J) = 0.0
           MoStabCorrShAbvCan(I,J)  = 0.0
           MoStabParaVeg2m(I,J)     = 0.0
           MoStabCorrMomVeg2m(I,J)  = 0.0
           MoStabCorrShVeg2m(I,J)   = 0.0
        endif

        ! evaluate stability-dependent variables using moz from prior iteration
        if ( MoStabParaAbvCan(I,J) < 0.0 ) then
           TMP1   = (1.0 - 16.0 * MoStabParaAbvCan(I,J))**0.25
           TMP2   = log((1.0 + TMP1*TMP1) / 2.0)
           TMP3   = log((1.0 + TMP1) / 2.0)
           FMNEW  = 2.0 * TMP3 + TMP2 - 2.0 * atan(TMP1) + 1.5707963
           FHNEW  = 2 * TMP2
           ! 2-meter quantities
           TMP12  = (1.0 - 16.0 * MoStabParaVeg2m(I,J))**0.25
           TMP22  = log((1.0 + TMP12*TMP12) / 2.0)
           TMP32  = log((1.0 + TMP12) / 2.0)
           FM2NEW = 2.0 * TMP32 + TMP22 - 2.0 * atan(TMP12) + 1.5707963
           FH2NEW = 2 * TMP22
        else
           FMNEW  = -5.0 * MoStabParaAbvCan(I,J)
           FHNEW  = FMNEW
           FM2NEW = -5.0 * MoStabParaVeg2m(I,J)
           FH2NEW = FM2NEW
        endif

        ! except for first iteration, weight stability factors for previous
        ! iteration to help avoid flip-flops from one iteration to the next
        if ( IterationInd == 1 ) then
           MoStabCorrMomAbvCan(I,J) = FMNEW
           MoStabCorrShAbvCan(I,J)  = FHNEW
           MoStabCorrMomVeg2m(I,J)  = FM2NEW
           MoStabCorrShVeg2m(I,J)   = FH2NEW
        else
           MoStabCorrMomAbvCan(I,J) = 0.5 * (MoStabCorrMomAbvCan(I,J) + FMNEW)
           MoStabCorrShAbvCan(I,J)  = 0.5 * (MoStabCorrShAbvCan(I,J)  + FHNEW)
           MoStabCorrMomVeg2m(I,J)  = 0.5 * (MoStabCorrMomVeg2m(I,J)  + FM2NEW)
           MoStabCorrShVeg2m(I,J)   = 0.5 * (MoStabCorrShVeg2m(I,J)   + FH2NEW)
        endif

        ! exchange coefficients
        MoStabCorrShAbvCan(I,J)  = min(MoStabCorrShAbvCan(I,J) , 0.9*TMPCH)
        MoStabCorrMomAbvCan(I,J) = min(MoStabCorrMomAbvCan(I,J), 0.9*TMPCM)
        MoStabCorrShVeg2m(I,J)   = min(MoStabCorrShVeg2m(I,J)  , 0.9*TMPCH2)
        MoStabCorrMomVeg2m(I,J)  = min(MoStabCorrMomVeg2m(I,J) , 0.9*TMPCM2)
        CMFM   = TMPCM  - MoStabCorrMomAbvCan(I,J)
        CHFH   = TMPCH  - MoStabCorrShAbvCan(I,J)
        CM2FM2 = TMPCM2 - MoStabCorrMomVeg2m(I,J)
        CH2FH2 = TMPCH2 - MoStabCorrShVeg2m(I,J)
        if ( abs(CMFM) <= MPE )   CMFM   = MPE
        if ( abs(CHFH) <= MPE )   CHFH   = MPE
        if ( abs(CM2FM2) <= MPE ) CM2FM2 = MPE
        if ( abs(CH2FH2) <= MPE ) CH2FH2 = MPE
        ExchCoeffMomAbvCan(I,J) = ConstVonKarman * ConstVonKarman / (CMFM * CMFM)
        ExchCoeffShAbvCan(I,J)  = ConstVonKarman * ConstVonKarman / (CMFM * CHFH)
        !ExchCoeffSh2mVegMo = ConstVonKarman * ConstVonKarman / (CM2FM2 * CH2FH2)

        ! friction velocity
        FrictionVelVeg(I,J)     = WindSpdRefHeight(I,J) * sqrt(ExchCoeffMomAbvCan(I,J))
        ExchCoeffSh2mVegMo(I,J) = ConstVonKarman * FrictionVelVeg(I,J) / CH2FH2

        ! aerodynamic resistance
        ResistanceMomAbvCan(I,J) = max(1.0, 1.0/(ExchCoeffMomAbvCan(I,J)*WindSpdRefHeight(I,J)))
        ResistanceShAbvCan(I,J)  = max(1.0, 1.0/(ExchCoeffShAbvCan(I,J)*WindSpdRefHeight(I,J)))
        ResistanceLhAbvCan(I,J)  = ResistanceShAbvCan(I,J)


      enddo
    enddo


        end associate

  end subroutine ResistanceAboveCanopyMOST

end module ResistanceAboveCanopyMostMod
