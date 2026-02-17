module CarbonFluxNatureVegMod

!!! Main Carbon assimilation for natural/generic vegetation
!!! based on RE Dickinson et al.(1998), modifed by Guo-Yue Niu, 2004

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine CarbonFluxNatureVeg(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: CO2FLUX
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: I, J                 ! grid indices
    real(kind=kind_noahmp)           :: DeathCoeffTemp       ! temperature stress death coefficient
    real(kind=kind_noahmp)           :: DeathCoeffWater      ! water stress death coefficient
    real(kind=kind_noahmp)           :: NetPriProdLeafAdd    ! leaf assimil after resp. losses removed [gC/m2/s]
    real(kind=kind_noahmp)           :: NetPriProdStemAdd    ! stem assimil after resp. losses removed [gC/m2/s]
    real(kind=kind_noahmp)           :: RespTmp, Temp0       ! temperary vars for function below
    RespTmp(Temp0) = exp(0.08 * (Temp0 - 298.16))            ! Respiration as a function of temperature

    associate(                                                                           &
              VegType                  => noahmp%config%domain%VegType              ,& ! in,    vegetation type
              MainTimeStep             => noahmp%config%domain%MainTimeStep         ,& ! in,    main noahmp timestep [s]
              IndexEBLForest           => noahmp%config%domain%IndexEBLForest       ,& ! in,    flag for Evergreen Broadleaf Forest
              WoodToRootRatio          => noahmp%biochem%param%WoodToRootRatio      ,& ! in,    wood to root ratio
              TurnoverCoeffLeafVeg     => noahmp%biochem%param%TurnoverCoeffLeafVeg ,& ! in,    leaf turnover coefficient [1/s] for generic vegetation
              TemperaureLeafFreeze     => noahmp%biochem%param%TemperaureLeafFreeze ,& ! in,    characteristic temperature for leaf freezing [K]
              LeafDeathWaterCoeffVeg   => noahmp%biochem%param%LeafDeathWaterCoeffVeg,& ! in,    coeficient for leaf water stress death [1/s] for generic veg
              LeafDeathTempCoeffVeg    => noahmp%biochem%param%LeafDeathTempCoeffVeg,& ! in,    coeficient for leaf temp. stress death [1/s] for generic veg
              GrowthRespFrac           => noahmp%biochem%param%GrowthRespFrac       ,& ! in,    fraction of growth respiration
              TemperatureMinPhotosyn   => noahmp%biochem%param%TemperatureMinPhotosyn,& ! in,    minimum temperature for photosynthesis [K]
              MicroRespCoeff           => noahmp%biochem%param%MicroRespCoeff       ,& ! in,    microbial respiration parameter [umol CO2/kgC/s]
              NitrogenConcFoliageMax   => noahmp%biochem%param%NitrogenConcFoliageMax,& ! in,    foliage nitrogen concentration when f(n)=1 (%)
              RespMaintQ10             => noahmp%biochem%param%RespMaintQ10         ,& ! in,    q10 for maintenance respiration
              RespMaintLeaf25C         => noahmp%biochem%param%RespMaintLeaf25C     ,& ! in,    leaf maintenance respiration at 25c [umol CO2/m2/s]
              RespMaintRoot25C         => noahmp%biochem%param%RespMaintRoot25C     ,& ! in,    root maintenance respiration at 25c [umol CO2/kgC/s]
              RespMaintStem25C         => noahmp%biochem%param%RespMaintStem25C     ,& ! in,    stem maintenance respiration at 25c [umol CO2/kgC/s]
              WoodPoolIndex            => noahmp%biochem%param%WoodPoolIndex        ,& ! in,    wood pool index (0~1) depending on woody or not
              TurnoverCoeffRootVeg     => noahmp%biochem%param%TurnoverCoeffRootVeg ,& ! in,    root turnover coefficient [1/s] for generic vegetation
              WoodRespCoeff            => noahmp%biochem%param%WoodRespCoeff        ,& ! in,    wood respiration coeficient [1/s]
              WoodAllocFac             => noahmp%biochem%param%WoodAllocFac         ,& ! in,    parameter for present wood allocation
              WaterStressCoeff         => noahmp%biochem%param%WaterStressCoeff     ,& ! in,    water stress coeficient
              LeafAreaIndexMin         => noahmp%biochem%param%LeafAreaIndexMin     ,& ! in,    minimum leaf area index [m2/m2]
              StemAreaIndexMin         => noahmp%biochem%param%StemAreaIndexMin     ,& ! in,    minimum stem area index [m2/m2]
              IndexGrowSeason          => noahmp%biochem%state%IndexGrowSeason      ,& ! in,    growing season index (0=off, 1=on)
              NitrogenConcFoliage      => noahmp%biochem%state%NitrogenConcFoliage  ,& ! in,    foliage nitrogen concentration [%]
              LeafAreaPerMass          => noahmp%biochem%state%LeafAreaPerMass      ,& ! in,    leaf area per unit mass [m2/g]
              PhotosynTotal            => noahmp%biochem%flux%PhotosynTotal         ,& ! in,    total leaf photosynthesis [umolCO2/m2/s]
              SoilWaterRootZone        => noahmp%water%state%SoilWaterRootZone      ,& ! in,    root zone soil water
              SoilWaterStress          => noahmp%water%state%SoilWaterStress        ,& ! in,    water stress coeficient (1.0 for wilting)
              TemperatureSoilSnow      => noahmp%energy%state%TemperatureSoilSnow   ,& ! in,    snow and soil layer temperature [K]
              TemperatureCanopy        => noahmp%energy%state%TemperatureCanopy     ,& ! in,    vegetation temperature [K]
              LeafAreaIndex            => noahmp%energy%state%LeafAreaIndex         ,& ! inout, leaf area index
              StemAreaIndex            => noahmp%energy%state%StemAreaIndex         ,& ! inout, stem area index
              LeafMass                 => noahmp%biochem%state%LeafMass             ,& ! inout, leaf mass [gC/m2]
              RootMass                 => noahmp%biochem%state%RootMass             ,& ! inout, mass of fine roots [gC/m2]
              StemMass                 => noahmp%biochem%state%StemMass             ,& ! inout, stem mass [gC/m2]
              WoodMass                 => noahmp%biochem%state%WoodMass             ,& ! inout, mass of wood (incl. woody roots) [gC/m2]
              CarbonMassDeepSoil       => noahmp%biochem%state%CarbonMassDeepSoil   ,& ! inout, stable carbon in deep soil [gC/m2]
              CarbonMassShallowSoil    => noahmp%biochem%state%CarbonMassShallowSoil,& ! inout, short-lived carbon in shallow soil [gC/m2]
              CarbonMassSoilTot        => noahmp%biochem%state%CarbonMassSoilTot    ,& ! out,   total soil carbon [gC/m2]
              CarbonMassLiveTot        => noahmp%biochem%state%CarbonMassLiveTot    ,& ! out,   total living carbon ([gC/m2]
              LeafMassMin              => noahmp%biochem%state%LeafMassMin          ,& ! out,   minimum leaf mass [gC/m2]
              CarbonFracToLeaf         => noahmp%biochem%state%CarbonFracToLeaf     ,& ! out,   fraction of carbon allocated to leaves
              WoodCarbonFrac           => noahmp%biochem%state%WoodCarbonFrac       ,& ! out,   calculated wood to root ratio
              CarbonFracToWoodRoot     => noahmp%biochem%state%CarbonFracToWoodRoot ,& ! out,   fraction of carbon to root and wood
              CarbonFracToRoot         => noahmp%biochem%state%CarbonFracToRoot     ,& ! out,   fraction of carbon flux to roots
              CarbonFracToWood         => noahmp%biochem%state%CarbonFracToWood     ,& ! out,   fraction of carbon flux to wood
              CarbonFracToStem         => noahmp%biochem%state%CarbonFracToStem     ,& ! out,   fraction of carbon flux to stem
              MicroRespFactorSoilWater => noahmp%biochem%state%MicroRespFactorSoilWater,& ! out,   soil water factor for microbial respiration
              MicroRespFactorSoilTemp  => noahmp%biochem%state%MicroRespFactorSoilTemp,& ! out,   soil temperature factor for microbial respiration
              RespFacNitrogenFoliage   => noahmp%biochem%state%RespFacNitrogenFoliage,& ! out,   foliage nitrogen adjustemt to respiration (<= 1)
              RespFacTemperature       => noahmp%biochem%state%RespFacTemperature   ,& ! out,   temperature factor
              RespReductionFac         => noahmp%biochem%state%RespReductionFac     ,& ! out,   respiration reduction factor (<= 1)
              StemMassMin              => noahmp%biochem%state%StemMassMin          ,& ! out,   minimum stem mass [gC/m2]
              StemAreaPerMass          => noahmp%biochem%state%StemAreaPerMass      ,& ! out,   stem area per unit mass [m2/g]
              CarbonAssim              => noahmp%biochem%flux%CarbonAssim           ,& ! out,   carbon assimilated rate [gC/m2/s]
              GrossPriProduction       => noahmp%biochem%flux%GrossPriProduction    ,& ! out,   gross primary production [gC/m2/s]
              NetPriProductionTot      => noahmp%biochem%flux%NetPriProductionTot   ,& ! out,   total net primary productivity [gC/m2/s]
              NetEcoExchange           => noahmp%biochem%flux%NetEcoExchange        ,& ! out,   net ecosystem exchange [gCO2/m2/s]
              RespirationPlantTot      => noahmp%biochem%flux%RespirationPlantTot   ,& ! out,   total plant respiration [gC/m2/s]
              RespirationSoilOrg       => noahmp%biochem%flux%RespirationSoilOrg    ,& ! out,   soil organic respiration [gC/m2/s]
              CarbonToAtmos            => noahmp%biochem%flux%CarbonToAtmos         ,& ! out,   carbon flux to atmosphere [gC/m2/s]
              NetPriProductionLeaf     => noahmp%biochem%flux%NetPriProductionLeaf  ,& ! out,   leaf net primary productivity [gC/m2/s]
              NetPriProductionRoot     => noahmp%biochem%flux%NetPriProductionRoot  ,& ! out,   root net primary productivity [gC/m2/s]
              NetPriProductionWood     => noahmp%biochem%flux%NetPriProductionWood  ,& ! out,   wood net primary productivity [gC/m2/s]
              NetPriProductionStem     => noahmp%biochem%flux%NetPriProductionStem  ,& ! out,   stem net primary productivity [gC/m2/s]
              GrowthRespLeaf           => noahmp%biochem%flux%GrowthRespLeaf        ,& ! out,   growth respiration rate for leaf [gC/m2/s]
              GrowthRespRoot           => noahmp%biochem%flux%GrowthRespRoot        ,& ! out,   growth respiration rate for root [gC/m2/s]
              GrowthRespWood           => noahmp%biochem%flux%GrowthRespWood        ,& ! out,   growth respiration rate for wood [gC/m2/s]
              GrowthRespStem           => noahmp%biochem%flux%GrowthRespStem        ,& ! out,   growth respiration rate for stem [gC/m2/s]
              LeafMassMaxChg           => noahmp%biochem%flux%LeafMassMaxChg        ,& ! out,   maximum leaf mass available to change [gC/m2/s]
              CarbonDecayToStable      => noahmp%biochem%flux%CarbonDecayToStable   ,& ! out,   decay rate of fast carbon to slow carbon [gC/m2/s]
              RespirationLeaf          => noahmp%biochem%flux%RespirationLeaf       ,& ! out,   leaf respiration rate [umol CO2/m2/s]
              RespirationStem          => noahmp%biochem%flux%RespirationStem       ,& ! out,   stem respiration rate [gC/m2/s]
              RespirationWood          => noahmp%biochem%flux%RespirationWood       ,& ! out,   wood respiration rate [gC/m2/s]
              RespirationLeafMaint     => noahmp%biochem%flux%RespirationLeafMaint  ,& ! out,   leaf maintenance respiration rate [gC/m2/s]
              RespirationRoot          => noahmp%biochem%flux%RespirationRoot       ,& ! out,   fine root respiration rate [gC/m2/s]
              RespirationSoil          => noahmp%biochem%flux%RespirationSoil       ,& ! out,   soil respiration rate [gC/m2/s]
              DeathLeaf                => noahmp%biochem%flux%DeathLeaf             ,& ! out,   death rate of leaf mass [gC/m2/s]
              DeathStem                => noahmp%biochem%flux%DeathStem             ,& ! out,   death rate of stem mass [gC/m2/s]
              TurnoverLeaf             => noahmp%biochem%flux%TurnoverLeaf          ,& ! out,   leaf turnover rate [gC/m2/s]
              TurnoverStem             => noahmp%biochem%flux%TurnoverStem          ,& ! out,   stem turnover rate [gC/m2/s]
              TurnoverWood             => noahmp%biochem%flux%TurnoverWood          ,& ! out,   wood turnover rate [gC/m2/s]
              TurnoverRoot             => noahmp%biochem%flux%TurnoverRoot          ,& ! out,   root turnover rate [gC/m2/s]
              StemMassMaxChg           => noahmp%biochem%flux%StemMassMaxChg         & ! out,   maximum steam mass available to change [gC/m2/s]
             )

   !$acc parallel loop collapse(2) gang vector default(present) &
   !$acc private(DeathCoeffTemp, DeathCoeffWater, NetPriProdLeafAdd, NetPriProdStemAdd, Temp0)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! condition to cycle moved from NoahmpMainMod to here
    if ( noahmp%config%domain%FlagDynamicVeg(I,J) .eqv. .false. ) cycle

!-----------------------------------------------------------------------

    ! initialization
    StemAreaPerMass(I,J) = 3.0 * 0.001      ! m2/kg -->m2/g
    LeafMassMin(I,J)     = LeafAreaIndexMin(I,J) / LeafAreaPerMass(I,J)   ! gC/m2
    StemMassMin(I,J)     = StemAreaIndexMin(I,J) / StemAreaPerMass(I,J)   ! gC/m2

    ! respiration
    if ( IndexGrowSeason(I,J) == 0.0 ) then
       RespReductionFac(I,J) = 0.5
    else
       RespReductionFac(I,J) = 1.0
    endif
    RespFacNitrogenFoliage(I,J) = min(NitrogenConcFoliage(I,J) / max(1.0e-06,NitrogenConcFoliageMax(I,J)), 1.0)
    RespFacTemperature(I,J)     = RespMaintQ10(I,J)**((TemperatureCanopy(I,J) - 298.16) / 10.0)
    RespirationLeaf(I,J)        = RespMaintLeaf25C(I,J) * RespFacTemperature(I,J) * RespFacNitrogenFoliage(I,J) * &
                             LeafAreaIndex(I,J) * RespReductionFac(I,J) * (1.0 - SoilWaterStress(I,J))                               ! umol CO2/m2/s
    RespirationLeafMaint(I,J)   = min((LeafMass(I,J)-LeafMassMin(I,J))/MainTimeStep, RespirationLeaf(I,J)*12.0e-6)                        ! gC/m2/s
    RespirationRoot(I,J)        = RespMaintRoot25C(I,J) * (RootMass(I,J)*1.0e-3) * RespFacTemperature(I,J) * RespReductionFac(I,J) * 12.0e-6   ! gC/m2/s
    RespirationStem(I,J)        = RespMaintStem25C(I,J) * ((StemMass(I,J)-StemMassMin(I,J)) * 1.0e-3) * &
                             RespFacTemperature(I,J) * RespReductionFac(I,J) * 12.0e-6                                          ! gC/m2/s
    RespirationWood(I,J)        = WoodRespCoeff(I,J) * RespTmp(TemperatureCanopy(I,J)) * WoodMass(I,J) * WoodPoolIndex(I,J)                    ! gC/m2/s

    !!! carbon assimilation start
    ! 1 mole -> 12 g carbon or 44 g CO2; 1 umol -> 12.e-6 g carbon;
    CarbonAssim(I,J) = PhotosynTotal(I,J) * 12.0e-6      ! umol CO2/m2/s -> gC/m2/s

    ! fraction of carbon into leaf versus nonleaf
    CarbonFracToLeaf(I,J)     = exp(0.01 * (1.0 - exp(0.75*LeafAreaIndex(I,J))) * LeafAreaIndex(I,J))
    if ( VegType(I,J) == IndexEBLForest ) CarbonFracToLeaf(I,J) = exp(0.01 * (1.0 - exp(0.50*LeafAreaIndex(I,J))) * LeafAreaIndex(I,J))
    CarbonFracToWoodRoot(I,J) = 1.0 - CarbonFracToLeaf(I,J)
    CarbonFracToStem(I,J)     = LeafAreaIndex(I,J) / 10.0 * CarbonFracToLeaf(I,J)
    CarbonFracToLeaf(I,J)     = CarbonFracToLeaf(I,J) - CarbonFracToStem(I,J)

    !  fraction of carbon into wood versus root
    if ( WoodMass(I,J) > 1.0e-6 ) then
       WoodCarbonFrac(I,J) = (1.0 - exp(-WoodAllocFac(I,J) * (WoodToRootRatio(I,J)*RootMass(I,J)/WoodMass(I,J))) / WoodAllocFac(I,J)) * WoodPoolIndex(I,J)
    else
       WoodCarbonFrac(I,J) = WoodPoolIndex(I,J)
    endif
    CarbonFracToRoot(I,J)  = CarbonFracToWoodRoot(I,J) * (1.0 - WoodCarbonFrac(I,J))
    CarbonFracToWood(I,J)  = CarbonFracToWoodRoot(I,J) * WoodCarbonFrac(I,J)

    ! leaf and root turnover per time step
    TurnoverLeaf(I,J) = TurnoverCoeffLeafVeg(I,J) * 5.0e-7 * LeafMass(I,J)   ! gC/m2/s
    TurnoverStem(I,J) = TurnoverCoeffLeafVeg(I,J) * 5.0e-7 * StemMass(I,J)   ! gC/m2/s
    TurnoverRoot(I,J) = TurnoverCoeffRootVeg(I,J) * RootMass(I,J)            ! gC/m2/s
    TurnoverWood(I,J) = 9.5e-10 * WoodMass(I,J)                         ! gC/m2/s

    ! seasonal leaf die rate dependent on temp and water stress
    ! water stress is set to 1 at permanent wilting point
    DeathCoeffTemp  = exp(-0.3 * max(0.0, TemperatureCanopy(I,J)-TemperaureLeafFreeze(I,J))) * (LeafMass(I,J) / 120.0)
    DeathCoeffWater = exp((SoilWaterStress(I,J) - 1.0) * WaterStressCoeff(I,J))
    DeathLeaf(I,J)       = LeafMass(I,J) * 1.0e-6 * (LeafDeathWaterCoeffVeg(I,J) * DeathCoeffWater + LeafDeathTempCoeffVeg(I,J) * DeathCoeffTemp)  ! gC/m2/s
    DeathStem(I,J)       = StemMass(I,J) * 1.0e-6 * (LeafDeathWaterCoeffVeg(I,J) * DeathCoeffWater + LeafDeathTempCoeffVeg(I,J) * DeathCoeffTemp)  ! gC/m2/s

    ! calculate growth respiration for leaf, root and wood
    GrowthRespLeaf(I,J) = max(0.0, GrowthRespFrac(I,J) * (CarbonFracToLeaf(I,J)*CarbonAssim(I,J) - RespirationLeafMaint(I,J)))  ! gC/m2/s
    GrowthRespStem(I,J) = max(0.0, GrowthRespFrac(I,J) * (CarbonFracToStem(I,J)*CarbonAssim(I,J) - RespirationStem(I,J)))       ! gC/m2/s
    GrowthRespRoot(I,J) = max(0.0, GrowthRespFrac(I,J) * (CarbonFracToRoot(I,J)*CarbonAssim(I,J) - RespirationRoot(I,J)))       ! gC/m2/s
    GrowthRespWood(I,J) = max(0.0, GrowthRespFrac(I,J) * (CarbonFracToWood(I,J)*CarbonAssim(I,J) - RespirationWood(I,J)))       ! gC/m2/s

    ! Impose lower T limit for photosynthesis
    NetPriProdLeafAdd = max(0.0, CarbonFracToLeaf(I,J)*CarbonAssim(I,J) - GrowthRespLeaf(I,J) - RespirationLeafMaint(I,J)) ! gC/m2/s
    NetPriProdStemAdd = max(0.0, CarbonFracToStem(I,J)*CarbonAssim(I,J) - GrowthRespStem(I,J) - RespirationStem(I,J))      ! gC/m2/s
   !NetPriProdLeafAdd = CarbonFracToLeaf(I,J)*CarbonAssim(I,J) - GrowthRespLeaf(I,J) - RespirationLeafMaint(I,J)  ! MB: test Kjetil
   !NetPriProdStemAdd = CarbonFracToStem(I,J)*CarbonAssim(I,J) - GrowthRespStem(I,J) - RespirationStem(I,J)       ! MB: test Kjetil
    if ( TemperatureCanopy(I,J) < TemperatureMinPhotosyn(I,J) ) NetPriProdLeafAdd = 0.0
    if ( TemperatureCanopy(I,J) < TemperatureMinPhotosyn(I,J) ) NetPriProdStemAdd = 0.0

    ! update leaf, root, and wood carbon
    ! avoid reducing leaf mass below its minimum value but conserve mass
    LeafMassMaxChg(I,J) = (LeafMass(I,J) - LeafMassMin(I,J)) / MainTimeStep                        ! gC/m2/s
    StemMassMaxChg(I,J) = (StemMass(I,J) - StemMassMin(I,J)) / MainTimeStep                        ! gC/m2/s
    DeathLeaf(I,J)      = min(DeathLeaf(I,J), LeafMassMaxChg(I,J)+NetPriProdLeafAdd-TurnoverLeaf(I,J))  ! gC/m2/s
    DeathStem(I,J)      = min(DeathStem(I,J), StemMassMaxChg(I,J)+NetPriProdStemAdd-TurnoverStem(I,J))  ! gC/m2/s

    ! net primary productivities
    NetPriProductionLeaf(I,J) = max(NetPriProdLeafAdd, -LeafMassMaxChg(I,J))                             ! gC/m2/s
    NetPriProductionStem(I,J) = max(NetPriProdStemAdd, -StemMassMaxChg(I,J))                             ! gC/m2/s
    NetPriProductionRoot(I,J) = CarbonFracToRoot(I,J) * CarbonAssim(I,J) - RespirationRoot(I,J) - GrowthRespRoot(I,J)   ! gC/m2/s
    NetPriProductionWood(I,J) = CarbonFracToWood(I,J) * CarbonAssim(I,J) - RespirationWood(I,J) - GrowthRespWood(I,J)   ! gC/m2/s

    ! masses of plant components
    LeafMass(I,J) = LeafMass(I,J) + (NetPriProductionLeaf(I,J) - TurnoverLeaf(I,J) - DeathLeaf(I,J)) * MainTimeStep   ! gC/m2
    StemMass(I,J) = StemMass(I,J) + (NetPriProductionStem(I,J) - TurnoverStem(I,J) - DeathStem(I,J)) * MainTimeStep   ! gC/m2
    RootMass(I,J) = RootMass(I,J) + (NetPriProductionRoot(I,J) - TurnoverRoot(I,J)) * MainTimeStep               ! gC/m2
    if ( RootMass(I,J) < 0.0 ) then
       TurnoverRoot(I,J) = NetPriProductionRoot(I,J)
       RootMass(I,J)     = 0.0
    endif
    WoodMass(I,J) = (WoodMass(I,J) + (NetPriProductionWood(I,J) - TurnoverWood(I,J)) * MainTimeStep ) * WoodPoolIndex(I,J)  ! gC/m2

    ! soil carbon budgets
    CarbonMassShallowSoil(I,J)    = CarbonMassShallowSoil(I,J) + &
                               (TurnoverRoot(I,J)+TurnoverLeaf(I,J)+TurnoverStem(I,J)+TurnoverWood(I,J)+DeathLeaf(I,J)+DeathStem(I,J)) * MainTimeStep  ! gC/m2, MB: add DeathStem v3.7
    MicroRespFactorSoilTemp(I,J)  = 2.0**( (TemperatureSoilSnow(I,1,J) - 283.16) / 10.0 )
    MicroRespFactorSoilWater(I,J) = SoilWaterRootZone(I,J) / (0.20 + SoilWaterRootZone(I,J)) * 0.23 / (0.23 + SoilWaterRootZone(I,J))
    RespirationSoil(I,J)          = MicroRespFactorSoilWater(I,J) * MicroRespFactorSoilTemp(I,J) * &
                               MicroRespCoeff(I,J) * max(0.0, CarbonMassShallowSoil(I,J)*1.0e-3) * 12.0e-6              ! gC/m2/s
    CarbonDecayToStable(I,J)      = 0.1 * RespirationSoil(I,J)                                                          ! gC/m2/s
    CarbonMassShallowSoil(I,J)    = CarbonMassShallowSoil(I,J) - (RespirationSoil(I,J) + CarbonDecayToStable(I,J)) * MainTimeStep ! gC/m2
    CarbonMassDeepSoil(I,J)       = CarbonMassDeepSoil(I,J) + CarbonDecayToStable(I,J) * MainTimeStep                        ! gC/m2

    !  total carbon flux ! MB: add RespirationStem,GrowthRespStem,0.9*RespirationSoil v3.7
    CarbonToAtmos(I,J)       = - CarbonAssim(I,J) + RespirationLeafMaint(I,J) + RespirationRoot(I,J) + RespirationWood(I,J) + RespirationStem(I,J) + &
                          0.9*RespirationSoil(I,J) + GrowthRespLeaf(I,J) + GrowthRespRoot(I,J) + GrowthRespWood(I,J) + GrowthRespStem(I,J)        ! gC/m2/s

    ! for outputs ! MB: add RespirationStem, GrowthRespStem in RespirationPlantTot v3.7
    GrossPriProduction(I,J)  = CarbonAssim(I,J)                                                                                 ! gC/m2/s
    NetPriProductionTot(I,J) = NetPriProductionLeaf(I,J) + NetPriProductionWood(I,J) + NetPriProductionRoot(I,J) + NetPriProductionStem(I,J)   ! gC/m2/s
    RespirationPlantTot(I,J) = RespirationRoot(I,J) + RespirationWood(I,J) + RespirationLeafMaint(I,J) + RespirationStem(I,J) + &
                          GrowthRespLeaf(I,J) + GrowthRespRoot(I,J) + GrowthRespWood(I,J) + GrowthRespStem(I,J)                           ! gC/m2/s
    RespirationSoilOrg(I,J)  = 0.9 * RespirationSoil(I,J)                                                                       ! gC/m2/s MB: add 0.9* v3.7
    NetEcoExchange(I,J)      = (RespirationPlantTot(I,J) + RespirationSoilOrg(I,J) - GrossPriProduction(I,J)) * 44.0 / 12.0               ! gCO2/m2/s
    CarbonMassSoilTot(I,J)   = CarbonMassShallowSoil(I,J) + CarbonMassDeepSoil(I,J)                                                  ! gC/m2
    CarbonMassLiveTot(I,J)   = LeafMass(I,J) + RootMass(I,J) + StemMass(I,J) + WoodMass(I,J)                                                   ! gC/m2   MB: add StemMass v3.7

    ! leaf area index and stem area index
    LeafAreaIndex(I,J)       = max(LeafMass(I,J)*LeafAreaPerMass(I,J), LeafAreaIndexMin(I,J))
    StemAreaIndex(I,J)       = max(StemMass(I,J)*StemAreaPerMass(I,J), StemAreaIndexMin(I,J))


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine CarbonFluxNatureVeg

end module CarbonFluxNatureVegMod
