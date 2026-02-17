module ResistanceCanopyStomataBallBerryMod

!!! Compute canopy stomatal resistance and foliage photosynthesis based on Ball-Berry scheme

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine ResistanceCanopyStomataBallBerry(noahmp, IndexShade)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: STOMATA
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    integer          , intent(in   ) :: IndexShade            ! index for sunlit/shaded (0=sunlit;1=shaded)
    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, JJ                  ! grid indices
    integer                          :: IndIter               ! iteration index
    integer, parameter               :: NumIter = 3           ! number of iterations
    real(kind=kind_noahmp)           :: RadPhotoActAbsTmp     ! temporary absorbed par for leaves [W/m2]
    real(kind=kind_noahmp)           :: ResistanceStomataTmp  ! temporary leaf stomatal resistance [s/m]
    real(kind=kind_noahmp)           :: PhotosynLeafTmp       ! temporary leaf photosynthesis [umol co2/m2/s]
    real(kind=kind_noahmp)           :: NitrogenFoliageFac    ! foliage nitrogen adjustment factor (0 to 1)
    real(kind=kind_noahmp)           :: CarboxylRateMax       ! maximum rate of carbonylation [umol co2/m2/s]
    real(kind=kind_noahmp)           :: MPE                   ! prevents overflow for division by zero
    real(kind=kind_noahmp)           :: RLB                   ! boundary layer resistance [s m2 / umol]
    real(kind=kind_noahmp)           :: TC                    ! foliage temperature [deg C]
    real(kind=kind_noahmp)           :: CS                    ! co2 concentration at leaf surface [Pa]
    real(kind=kind_noahmp)           :: KC                    ! co2 Michaelis-Menten constant [Pa]
    real(kind=kind_noahmp)           :: KO                    ! o2 Michaelis-Menten constant [Pa]
    real(kind=kind_noahmp)           :: A,B,C,Q               ! intermediate calculations for RS
    real(kind=kind_noahmp)           :: R1,R2                 ! roots for RS
    real(kind=kind_noahmp)           :: PPF                   ! absorb photosynthetic photon flux [umol photons/m2/s]
    real(kind=kind_noahmp)           :: WC                    ! Rubisco limited photosynthesis [umol co2/m2/s]
    real(kind=kind_noahmp)           :: WJ                    ! light limited photosynthesis [umol co2/m2/s]
    real(kind=kind_noahmp)           :: WE                    ! export limited photosynthesis [umol co2/m2/s]
    real(kind=kind_noahmp)           :: CP                    ! co2 compensation point [Pa]
    real(kind=kind_noahmp)           :: CI                    ! internal co2 [Pa]
    real(kind=kind_noahmp)           :: AWC                   ! intermediate calculation for wc
    real(kind=kind_noahmp)           :: J                     ! electron transport [umol co2/m2/s]
    real(kind=kind_noahmp)           :: CEA                   ! constrain ea or else model blows up
    real(kind=kind_noahmp)           :: CF                    ! [s m2/umol] -> [s/m]
    real(kind=kind_noahmp)           :: T                     ! temporary var
! local statement functions
    real(kind=kind_noahmp)           :: F1                    ! generic temperature response (statement function)
    real(kind=kind_noahmp)           :: F2                    ! generic temperature inhibition (statement function)
    real(kind=kind_noahmp)           :: AB                    ! used in statement functions
    real(kind=kind_noahmp)           :: BC                    ! used in statement functions
    F1(AB, BC) = AB**( (BC - 25.0) / 10.0 )
    F2(AB)     = 1.0 + exp( (-2.2e05 + 710.0 * (AB + 273.16)) / (8.314 * (AB + 273.16)) )

! --------------------------------------------------------------------
    associate(                                                                             &
              PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight              ,& ! in,  air pressure [Pa] at reference height
              TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight           ,& ! in,  air temperature [K] at reference height
              SoilTranspFacAcc        => noahmp%water%state%SoilTranspFacAcc              ,& ! in,  accumulated soil water transpiration factor (0 to 1)
              IndexGrowSeason         => noahmp%biochem%state%IndexGrowSeason             ,& ! in,  growing season index (0=off, 1=on)
              NitrogenConcFoliage     => noahmp%biochem%state%NitrogenConcFoliage         ,& ! in,  foliage nitrogen concentration [%]
              NitrogenConcFoliageMax  => noahmp%biochem%param%NitrogenConcFoliageMax      ,& ! in,  foliage nitrogen concentration when f(n)=1 [%]
              QuantumEfficiency25C    => noahmp%biochem%param%QuantumEfficiency25C         ,& ! in,  quantum efficiency at 25c [umol co2 / umol photon]
              CarboxylRateMax25C      => noahmp%biochem%param%CarboxylRateMax25C           ,& ! in,  maximum rate of carboxylation at 25c [umol co2/m**2/s]
              CarboxylRateMaxQ10      => noahmp%biochem%param%CarboxylRateMaxQ10           ,& ! in,  change in maximum rate of carboxylation for each 10C temp change
              PhotosynPathC3          => noahmp%biochem%param%PhotosynPathC3               ,& ! in,  C3 photosynthetic pathway indicator: 0. = c4, 1. = c3
              SlopeConductToPhotosyn  => noahmp%biochem%param%SlopeConductToPhotosyn       ,& ! in,  slope of conductance-to-photosynthesis relationship
              Co2MmConst25C           => noahmp%energy%param%Co2MmConst25C                ,& ! in,  co2 michaelis-menten constant at 25c [Pa]
              O2MmConst25C            => noahmp%energy%param%O2MmConst25C                 ,& ! in,  o2 michaelis-menten constant at 25c [Pa]
              Co2MmConstQ10           => noahmp%energy%param%Co2MmConstQ10                ,& ! in,  q10 for Co2MmConst25C
              O2MmConstQ10            => noahmp%energy%param%O2MmConstQ10                 ,& ! in,  q10 for ko25
              ConductanceLeafMin      => noahmp%energy%param%ConductanceLeafMin           ,& ! in,  minimum leaf conductance [umol/m**2/s]
              TemperatureCanopy       => noahmp%energy%state%TemperatureCanopy             ,& ! in,  vegetation temperature [K]
              VapPresSatCanopy        => noahmp%energy%state%VapPresSatCanopy              ,& ! in,  canopy saturation vapor pressure at TV [Pa]
              PressureVaporCanAir     => noahmp%energy%state%PressureVaporCanAir           ,& ! in,  canopy air vapor pressure [Pa]
              PressureAtmosO2         => noahmp%energy%state%PressureAtmosO2              ,& ! in,  atmospheric o2 pressure [Pa]
              PressureAtmosCO2        => noahmp%energy%state%PressureAtmosCO2             ,& ! in,  atmospheric co2 pressure [Pa]
              ResistanceLeafBoundary  => noahmp%energy%state%ResistanceLeafBoundary       ,& ! in,  leaf boundary layer resistance [s/m]
              VegFrac                 => noahmp%energy%state%VegFrac                      ,& ! in,  greeness vegetation fraction
              RadPhotoActAbsSunlit    => noahmp%energy%flux%RadPhotoActAbsSunlit           ,& ! in,  average absorbed par for sunlit leaves [W/m2]
              RadPhotoActAbsShade     => noahmp%energy%flux%RadPhotoActAbsShade            ,& ! in,  average absorbed par for shaded leaves [W/m2]
              ResistanceStomataSunlit => noahmp%energy%state%ResistanceStomataSunlit       ,& ! out, sunlit leaf stomatal resistance [s/m]
              ResistanceStomataShade  => noahmp%energy%state%ResistanceStomataShade        ,& ! out, shaded leaf stomatal resistance [s/m]
              PhotosynLeafSunlit      => noahmp%biochem%flux%PhotosynLeafSunlit            ,& ! out, sunlit leaf photosynthesis [umol co2/m2/s]
              PhotosynLeafShade       => noahmp%biochem%flux%PhotosynLeafShade              & ! out, shaded leaf photosynthesis [umol co2/m2/s]
             )

   !$acc parallel loop collapse(2) gang vector default(present) &
   !$acc private(IndIter,RadPhotoActAbsTmp,ResistanceStomataTmp,PhotosynLeafTmp,NitrogenFoliageFac) &
   !$acc private(CarboxylRateMax,MPE,RLB,TC,CS,KC,KO,A,B,C,Q,R1,R2,PPF,WC,WJ,WE,CP,CI,AWC,J,CEA,CF,T) &
   !$acc firstprivate(IndexShade)
    do JJ = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if ( .not. ((noahmp%energy%state%VegAreaIndEff(I,JJ) > 0.0 ) .and. (noahmp%energy%state%VegFrac(I,JJ) > 0)) ) cycle ! skip non-vegetated surface


    ! initialization
    MPE = 1.0e-6

    ! initialize ResistanceStomata=maximum value and photosynthesis=0 because will only do calculations
    ! for RadPhotoActAbs  > 0, in which case ResistanceStomata <= maximum value and photosynthesis >= 0
    CF = PressureAirRefHeight(I,JJ) / (8.314 * TemperatureAirRefHeight(I,JJ)) * 1.0e06  ! unit conversion factor
    ResistanceStomataTmp = 1.0 / ConductanceLeafMin(I,JJ) * CF
    PhotosynLeafTmp      = 0.0
    if ( IndexShade == 0 ) RadPhotoActAbsTmp = RadPhotoActAbsSunlit(I,JJ) / max(VegFrac(I,JJ),1.0e-6)  ! Sunlit case
    if ( IndexShade == 1 ) RadPhotoActAbsTmp = RadPhotoActAbsShade(I,JJ)  / max(VegFrac(I,JJ),1.0e-6)  ! Shaded case

    ! only compute when there is radiation absorption
    if ( RadPhotoActAbsTmp > 0.0 ) then

       NitrogenFoliageFac = min(NitrogenConcFoliage(I,JJ)/max(MPE, NitrogenConcFoliageMax(I,JJ)), 1.0)
       TC                 = TemperatureCanopy(I,JJ) - ConstFreezePoint
       PPF                = 4.6 * RadPhotoActAbsTmp
       J                  = PPF * QuantumEfficiency25C(I,JJ)
       KC                 = Co2MmConst25C(I,JJ) * F1(Co2MmConstQ10(I,JJ), TC)
       KO                 = O2MmConst25C(I,JJ) * F1(O2MmConstQ10(I,JJ), TC)
       AWC                = KC * ( 1.0 + PressureAtmosO2(I,JJ) / KO )
       CP                 = 0.5 * KC / KO * PressureAtmosO2(I,JJ) * 0.21
       CarboxylRateMax    = CarboxylRateMax25C(I,JJ) / F2(TC) * NitrogenFoliageFac * &
                            SoilTranspFacAcc(I,JJ) * F1(CarboxylRateMaxQ10(I,JJ), TC)
       ! first guess ci
       CI  = 0.7 * PressureAtmosCO2(I,JJ) * PhotosynPathC3(I,JJ) + 0.4 * PressureAtmosCO2(I,JJ) * (1.0 - PhotosynPathC3(I,JJ))
       ! ResistanceLeafBoundary: s/m -> s m**2 / umol
       RLB = ResistanceLeafBoundary(I,JJ) / CF
       ! constrain PressureVaporCanAir
       CEA = max(0.25*VapPresSatCanopy(I,JJ)*PhotosynPathC3(I,JJ) + 0.40*VapPresSatCanopy(I,JJ)*(1.0-PhotosynPathC3(I,JJ)), &
                 min(PressureVaporCanAir(I,JJ),VapPresSatCanopy(I,JJ)))

       ! ci iteration
       !$acc loop seq
       do IndIter = 1, NumIter
          WJ = max(CI-CP, 0.0) * J / (CI + 2.0*CP) * PhotosynPathC3(I,JJ) + J * (1.0 - PhotosynPathC3(I,JJ))
          WC = max(CI-CP, 0.0) * CarboxylRateMax / (CI + AWC) * PhotosynPathC3(I,JJ) + &
               CarboxylRateMax * (1.0 - PhotosynPathC3(I,JJ))
          WE = 0.5 * CarboxylRateMax * PhotosynPathC3(I,JJ) + &
               4000.0 * CarboxylRateMax * CI / PressureAirRefHeight(I,JJ) * (1.0 - PhotosynPathC3(I,JJ))
          PhotosynLeafTmp = min(WJ, WC, WE) * IndexGrowSeason(I,JJ)
          CS = max(PressureAtmosCO2(I,JJ)-1.37*RLB*PressureAirRefHeight(I,JJ)*PhotosynLeafTmp, MPE)
          A  = SlopeConductToPhotosyn(I,JJ) * PhotosynLeafTmp * PressureAirRefHeight(I,JJ) * CEA / &
               (CS * VapPresSatCanopy(I,JJ)) + ConductanceLeafMin(I,JJ)
          B  = (SlopeConductToPhotosyn(I,JJ) * PhotosynLeafTmp * PressureAirRefHeight(I,JJ) / CS + ConductanceLeafMin(I,JJ)) * &
               RLB - 1.0
          C  = -RLB
          if ( B >= 0.0 ) then
             Q = -0.5 * (B + sqrt(B*B-4.0*A*C))
          else
             Q = -0.5 * (B - sqrt(B*B-4.0*A*C))
          endif
          R1   = Q / A
          R2   = C / Q
          ResistanceStomataTmp = max(R1, R2)
          CI   = max(CS-PhotosynLeafTmp*PressureAirRefHeight(I,JJ)*1.65*ResistanceStomataTmp, 0.0)
       enddo

       ! ResistanceStomata:  s m**2 / umol -> s/m
       ResistanceStomataTmp = ResistanceStomataTmp * CF

    endif ! RadPhotoActAbsTmp > 0.0

    ! assign updated values
    ! Sunlit case
    if ( IndexShade == 0 ) then
       ResistanceStomataSunlit(I,JJ) = ResistanceStomataTmp
       PhotosynLeafSunlit(I,JJ)      = PhotosynLeafTmp
    endif
    ! Shaded case
    if ( IndexShade == 1 ) then
       ResistanceStomataShade(I,JJ)  = ResistanceStomataTmp
       PhotosynLeafShade(I,JJ)       = PhotosynLeafTmp
    endif


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine ResistanceCanopyStomataBallBerry

end module ResistanceCanopyStomataBallBerryMod
