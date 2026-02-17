module CarbonFluxCropMod

!!! Main Carbon assimilation for crops
!!! based on RE Dickinson et al.(1998), modifed by Guo-Yue Niu, 2004
!!! Modified by Xing Liu, 2014
        
  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none
        
contains

  subroutine CarbonFluxCrop(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: CO2FLUX_CROP
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------
 
    implicit none
        
    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                 ! grid indices
    real(kind=kind_noahmp)           :: DeathCoeffTemp       ! temperature stress death coefficient
    real(kind=kind_noahmp)           :: DeathCoeffWater      ! water stress death coefficient
    real(kind=kind_noahmp)           :: NetPriProdLeafAdd    ! leaf assimil after resp. losses removed [gCH2O/m2/s]
    real(kind=kind_noahmp)           :: NetPriProdStemAdd    ! stem assimil after resp. losses removed [gCH2O/m2/s]
   !real(kind=kind_noahmp)           :: RespTmp, Temp0       ! temperary vars for function below
   !RespTmp(Temp0) = exp(0.08 * (Temp0 - 298.16))            ! Respiration as a function of temperature

!------------------------------------------------------------------------
    associate(                                                                           &
              MainTimeStep             => noahmp%config%domain%MainTimeStep             ,& ! in,    main noahmp timestep [s]
              WaterStressCoeff         => noahmp%biochem%param%WaterStressCoeff         ,& ! in,    water stress coeficient
              LeafAreaIndexMin         => noahmp%biochem%param%LeafAreaIndexMin         ,& ! in,    minimum leaf area index [m2/m2]
              StemAreaIndexMin         => noahmp%biochem%param%StemAreaIndexMin         ,& ! in,    minimum stem area index [m2/m2]
              NitrogenConcFoliageMax   => noahmp%biochem%param%NitrogenConcFoliageMax   ,& ! in,    foliage nitrogen concentration when f(n)=1 [%]
              RespMaintQ10             => noahmp%biochem%param%RespMaintQ10             ,& ! in,    change in maintenance respiration for each 10C temp. change
              RespMaintLeaf25C         => noahmp%biochem%param%RespMaintLeaf25C         ,& ! in,    leaf maintenance respiration at 25C [umol CO2/m2/s]
              RespMaintRoot25C         => noahmp%biochem%param%RespMaintRoot25C         ,& ! in,    root maintenance respiration at 25C [umol CO2/kgCH2O/s]
              RespMaintStem25C         => noahmp%biochem%param%RespMaintStem25C         ,& ! in,    stem maintenance respiration at 25C [umol CO2/kgCH2O/s]
              RespMaintGrain25C        => noahmp%biochem%param%RespMaintGrain25C        ,& ! in,    grain maintenance respiration at 25C [umol CO2/kgCH2O/s]
              GrowthRespFrac           => noahmp%biochem%param%GrowthRespFrac           ,& ! in,    fraction of growth respiration
              CarbohydrFracToLeaf      => noahmp%biochem%param%CarbohydrFracToLeaf      ,& ! in,    fraction of carbohydrate flux to leaf
              CarbohydrFracToStem      => noahmp%biochem%param%CarbohydrFracToStem      ,& ! in,    fraction of carbohydrate flux to stem
              CarbohydrFracToRoot      => noahmp%biochem%param%CarbohydrFracToRoot      ,& ! in,    fraction of carbohydrate flux to root
              CarbohydrFracToGrain     => noahmp%biochem%param%CarbohydrFracToGrain     ,& ! in,    fraction of carbohydrate flux to grain
              TurnoverCoeffLeafCrop    => noahmp%biochem%param%TurnoverCoeffLeafCrop    ,& ! in,    leaf turnover coefficient [1/s] for crop
              TurnoverCoeffRootCrop    => noahmp%biochem%param%TurnoverCoeffRootCrop    ,& ! in,    root tunrover coefficient [1/s] for crop
              TurnoverCoeffStemCrop    => noahmp%biochem%param%TurnoverCoeffStemCrop    ,& ! in,    stem turnover coefficient [1/s] for crop
              TemperaureLeafFreeze     => noahmp%biochem%param%TemperaureLeafFreeze     ,& ! in,    characteristic temperature for leaf freezing [K]
              LeafDeathWaterCoeffCrop  => noahmp%biochem%param%LeafDeathWaterCoeffCrop  ,& ! in,    coeficient for water leaf stress death [1/s] for crop
              LeafDeathTempCoeffCrop   => noahmp%biochem%param%LeafDeathTempCoeffCrop   ,& ! in,    coeficient for temperature leaf stress death [1/s] for crop
              CarbohydrLeafToGrain     => noahmp%biochem%param%CarbohydrLeafToGrain     ,& ! in,    fraction of carbohydrate translocation from leaf to grain
              CarbohydrStemToGrain     => noahmp%biochem%param%CarbohydrStemToGrain     ,& ! in,    fraction of carbohydrate translocation from stem to grain
              CarbohydrRootToGrain     => noahmp%biochem%param%CarbohydrRootToGrain     ,& ! in,    fraction of carbohydrate translocation from root to grain
              MicroRespCoeff           => noahmp%biochem%param%MicroRespCoeff           ,& ! in,    microbial respiration parameter [umol CO2/kgC/s]
              LeafAreaPerBiomass       => noahmp%biochem%param%LeafAreaPerBiomass       ,& ! in,    leaf area per living leaf biomass [m2/g]
              SoilWaterRootZone        => noahmp%water%state%SoilWaterRootZone          ,& ! in,    root zone soil water
              SoilWaterStress          => noahmp%water%state%SoilWaterStress            ,& ! in,    water stress coeficient (1.0 for wilting)
              PhotosynTotal            => noahmp%biochem%flux%PhotosynTotal             ,& ! in,    total leaf photosynthesis [umol CO2/m2/s]
              NitrogenConcFoliage      => noahmp%biochem%state%NitrogenConcFoliage      ,& ! in,    foliage nitrogen concentration [%]
              IndexPlanting            => noahmp%biochem%state%IndexPlanting            ,& ! in,    Planting index
              PlantGrowStage           => noahmp%biochem%state%PlantGrowStage           ,& ! in,    plant growing stage
              TemperatureSoilSnow      => noahmp%energy%state%TemperatureSoilSnow       ,& ! in,    snow and soil layer temperature [K]
              TemperatureCanopy        => noahmp%energy%state%TemperatureCanopy         ,& ! in,    vegetation temperature [K]
              LeafAreaIndex            => noahmp%energy%state%LeafAreaIndex             ,& ! inout, leaf area index
              StemAreaIndex            => noahmp%energy%state%StemAreaIndex             ,& ! inout, stem area index
              LeafMass                 => noahmp%biochem%state%LeafMass                 ,& ! inout, leaf mass [gCH2O/m2]
              RootMass                 => noahmp%biochem%state%RootMass                 ,& ! inout, mass of fine roots [gCH2O/m2]
              StemMass                 => noahmp%biochem%state%StemMass                 ,& ! inout, stem mass [gCH2O/m2]
              CarbonMassDeepSoil       => noahmp%biochem%state%CarbonMassDeepSoil       ,& ! inout, stable carbon in deep soil [gC/m2]
              CarbonMassShallowSoil    => noahmp%biochem%state%CarbonMassShallowSoil    ,& ! inout, short-lived carbon in shallow soil [gC/m2]
              GrainMass                => noahmp%biochem%state%GrainMass                ,& ! inout, mass of grain [gCH2O/m2]
              RespFacNitrogenFoliage   => noahmp%biochem%state%RespFacNitrogenFoliage   ,& ! out,   foliage nitrogen adjustemt to respiration (<= 1)
              MicroRespFactorSoilWater => noahmp%biochem%state%MicroRespFactorSoilWater ,& ! out,   soil water factor for microbial respiration
              MicroRespFactorSoilTemp  => noahmp%biochem%state%MicroRespFactorSoilTemp  ,& ! out,   soil temperature factor for microbial respiration
              LeafMassMin              => noahmp%biochem%state%LeafMassMin              ,& ! out,   minimum leaf mass [gCH2O/m2]
              StemMassMin              => noahmp%biochem%state%StemMassMin              ,& ! out,   minimum stem mass [gCH2O/m2]
              StemAreaPerMass          => noahmp%biochem%state%StemAreaPerMass          ,& ! out,   stem area per unit mass [m2/g]
              RespFacTemperature       => noahmp%biochem%state%RespFacTemperature       ,& ! out,   temperature factor
              CarbonMassSoilTot        => noahmp%biochem%state%CarbonMassSoilTot        ,& ! out,   total soil carbon [gC/m2]
              CarbonMassLiveTot        => noahmp%biochem%state%CarbonMassLiveTot        ,& ! out,   total living carbon [gC/m2]
              CarbonAssim              => noahmp%biochem%flux%CarbonAssim               ,& ! out,   carbon assimilated rate [gC/m2/s]
              CarbohydrAssim           => noahmp%biochem%flux%CarbohydrAssim            ,& ! out,   carbohydrate assimilated rate [gCH2O/m2/s]
              TurnoverLeaf             => noahmp%biochem%flux%TurnoverLeaf              ,& ! out,   leaf turnover rate [gCH2O/m2/s]
              TurnoverStem             => noahmp%biochem%flux%TurnoverStem              ,& ! out,   stem turnover rate [gCH2O/m2/s]
              TurnoverRoot             => noahmp%biochem%flux%TurnoverRoot              ,& ! out,   root carbon loss rate by turnover [gCH2O/m2/s]
              ConvLeafToGrain          => noahmp%biochem%flux%ConvLeafToGrain           ,& ! out,   leaf to grain conversion [gCH2O/m2]
              ConvRootToGrain          => noahmp%biochem%flux%ConvRootToGrain           ,& ! out,   root to grain conversion [gCH2O/m2]
              ConvStemToGrain          => noahmp%biochem%flux%ConvStemToGrain           ,& ! out,   stem to grain conversion [gCH2O/m2]
              RespirationPlantTot      => noahmp%biochem%flux%RespirationPlantTot       ,& ! out,   total plant respiration [gC/m2/s C]
              CarbonToAtmos            => noahmp%biochem%flux%CarbonToAtmos             ,& ! out,   carbon flux to atmosphere [gC/m2/s]
              GrossPriProduction       => noahmp%biochem%flux%GrossPriProduction        ,& ! out,   gross primary production [gC/m2/s]
              NetPriProductionTot      => noahmp%biochem%flux%NetPriProductionTot       ,& ! out,   total net primary productivity [gC/m2/s]
              NetPriProductionLeaf     => noahmp%biochem%flux%NetPriProductionLeaf      ,& ! out,   leaf net primary productivity [gCH2O/m2/s]
              NetPriProductionRoot     => noahmp%biochem%flux%NetPriProductionRoot      ,& ! out,   root net primary productivity [gCH2O/m2/s]
              NetPriProductionStem     => noahmp%biochem%flux%NetPriProductionStem      ,& ! out,   stem net primary productivity [gCH2O/m2/s]
              NetPriProductionGrain    => noahmp%biochem%flux%NetPriProductionGrain     ,& ! out,   grain net primary productivity [gCH2O/m2/s]
              NetEcoExchange           => noahmp%biochem%flux%NetEcoExchange            ,& ! out,   net ecosystem exchange [gCO2/m2/s]
              GrowthRespGrain          => noahmp%biochem%flux%GrowthRespGrain           ,& ! out,   growth respiration rate for grain [gCH2O/m2/s]
              GrowthRespLeaf           => noahmp%biochem%flux%GrowthRespLeaf            ,& ! out,   growth respiration rate for leaf [gCH2O/m2/s]
              GrowthRespRoot           => noahmp%biochem%flux%GrowthRespRoot            ,& ! out,   growth respiration rate for root [gCH2O/m2/s]
              GrowthRespStem           => noahmp%biochem%flux%GrowthRespStem            ,& ! out,   growth respiration rate for stem [gCH2O/m2/s]
              RespirationSoilOrg       => noahmp%biochem%flux%RespirationSoilOrg        ,& ! out,   soil organic respiration rate [gC/m2/s]
              LeafMassMaxChg           => noahmp%biochem%flux%LeafMassMaxChg            ,& ! out,   maximum leaf mass available to change [gCH2O/m2/s]
              StemMassMaxChg           => noahmp%biochem%flux%StemMassMaxChg            ,& ! out,   maximum steam  mass available to change [gCH2O/m2/s]
              RespirationLeaf          => noahmp%biochem%flux%RespirationLeaf           ,& ! out,   leaf respiration rate [umol CO2/m2/s]
              RespirationStem          => noahmp%biochem%flux%RespirationStem           ,& ! out,   stem respiration rate [gCH2O/m2/s]
              RespirationLeafMaint     => noahmp%biochem%flux%RespirationLeafMaint      ,& ! out,   leaf maintenance respiration rate [gCH2O/m2/s]
              RespirationRoot          => noahmp%biochem%flux%RespirationRoot           ,& ! out,   fine root respiration rate [gCH2O/m2/s]
              RespirationSoil          => noahmp%biochem%flux%RespirationSoil           ,& ! out,   soil respiration rate [gCH2O/m2/s]
              RespirationGrain         => noahmp%biochem%flux%RespirationGrain          ,& ! out,   grain respiration rate [gCH2O/m2/s]
              DeathLeaf                => noahmp%biochem%flux%DeathLeaf                 ,& ! out,   death rate of leaf mass [gCH2O/m2/s]
              CarbonDecayToStable      => noahmp%biochem%flux%CarbonDecayToStable        & ! out,   decay rate of fast carbon to slow carbon [gCH2O/m2/s]
             )

!$acc parallel loop collapse(2) gang vector default(present) &
!$acc private(DeathCoeffTemp, DeathCoeffWater, NetPriProdLeafAdd, NetPriProdStemAdd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! condition to cycle moved from NoahmpMainMod to here
    if ( .not. (noahmp%config%domain%FlagDynamicCrop(I,J) .and. (noahmp%config%nmlist%OptCropModel == 1) ) ) cycle


    ! initialization
    StemAreaPerMass(I,J) = 3.0 * 0.001         ! m2/kg -->m2/g
    LeafMassMin(I,J)     = LeafAreaIndexMin(I,J) / 0.035
    StemMassMin(I,J)     = StemAreaIndexMin(I,J) / StemAreaPerMass(I,J)

    !!! carbon assimilation starts
    ! 1 mole -> 12 g carbon or 44 g CO2 or 30 g CH20
    CarbonAssim(I,J)     = PhotosynTotal(I,J) * 12.0e-6   !*IndexPlanting(I,J)   !umol co2 /m2/ s -> g/m2/s C
    CarbohydrAssim(I,J)  = PhotosynTotal(I,J) * 30.0e-6   !*IndexPlanting(I,J)   !umol co2 /m2/ s -> g/m2/s CH2O

    ! mainteinance respiration
    RespFacNitrogenFoliage(I,J) = min(NitrogenConcFoliage(I,J) / max(1.0e-06, NitrogenConcFoliageMax(I,J)), 1.0)
    RespFacTemperature(I,J)     = RespMaintQ10(I,J)**((TemperatureCanopy(I,J) - 298.16) / 10.0)
    RespirationLeaf(I,J)        = RespMaintLeaf25C(I,J) * RespFacTemperature(I,J) * RespFacNitrogenFoliage(I,J) * &
                             LeafAreaIndex(I,J) * (1.0 - SoilWaterStress(I,J))                                       ! umolCO2/m2/s
    RespirationLeafMaint(I,J)   = min((LeafMass(I,J) - LeafMassMin(I,J)) / MainTimeStep, RespirationLeaf(I,J)*30.0e-6)         ! gCH2O/m2/s
    RespirationRoot(I,J)        = RespMaintRoot25C(I,J) * (RootMass(I,J) * 1.0e-3) * RespFacTemperature(I,J) * 30.0e-6         ! gCH2O/m2/s
    RespirationStem(I,J)        = RespMaintStem25C(I,J) * (StemMass(I,J) * 1.0e-3) * RespFacTemperature(I,J) * 30.0e-6         ! gCH2O/m2/s
    RespirationGrain(I,J)       = RespMaintGrain25C(I,J) * (GrainMass(I,J) * 1.0e-3) * RespFacTemperature(I,J) * 30.0e-6       ! gCH2O/m2/s

    ! calculate growth respiration for leaf, root and grain
    GrowthRespLeaf(I,J)  = max(0.0, GrowthRespFrac(I,J) * (CarbohydrFracToLeaf(I,PlantGrowStage(I,J),J)*CarbohydrAssim(I,J) - RespirationLeafMaint(I,J)))  ! gCH2O/m2/s
    GrowthRespStem(I,J)  = max(0.0, GrowthRespFrac(I,J) * (CarbohydrFracToStem(I,PlantGrowStage(I,J),J)*CarbohydrAssim(I,J) - RespirationStem(I,J)))       ! gCH2O/m2/s
    GrowthRespRoot(I,J)  = max(0.0, GrowthRespFrac(I,J) * (CarbohydrFracToRoot(I,PlantGrowStage(I,J),J)*CarbohydrAssim(I,J) - RespirationRoot(I,J)))       ! gCH2O/m2/s
    GrowthRespGrain(I,J) = max(0.0, GrowthRespFrac(I,J) * (CarbohydrFracToGrain(I,PlantGrowStage(I,J),J)*CarbohydrAssim(I,J) - RespirationGrain(I,J)))     ! gCH2O/m2/s

    ! leaf turnover, stem turnover, root turnover and leaf death caused by soil water and soil temperature stress
    TurnoverLeaf(I,J)    = TurnoverCoeffLeafCrop(I,PlantGrowStage(I,J),J) * 1.0e-6 * LeafMass(I,J)     ! gCH2O/m2/s
    TurnoverRoot(I,J)    = TurnoverCoeffRootCrop(I,PlantGrowStage(I,J),J) * 1.0e-6 * RootMass(I,J)     ! gCH2O/m2/s
    TurnoverStem(I,J)    = TurnoverCoeffStemCrop(I,PlantGrowStage(I,J),J) * 1.0e-6 * StemMass(I,J)     ! gCH2O/m2/s
    DeathCoeffTemp  = exp(-0.3 * max(0.0, TemperatureCanopy(I,J)-TemperaureLeafFreeze(I,J))) * (LeafMass(I,J)/120.0)
    DeathCoeffWater = exp((SoilWaterStress(I,J) - 1.0) * WaterStressCoeff(I,J))
    DeathLeaf(I,J)       = LeafMass(I,J) * 1.0e-6 * (LeafDeathWaterCoeffCrop(I,PlantGrowStage(I,J),J) * DeathCoeffWater + &
                                           LeafDeathTempCoeffCrop(I,PlantGrowStage(I,J),J) * DeathCoeffTemp)      ! gCH2O/m2/s

    ! Allocation of CarbohydrAssim to leaf, stem, root and grain at each growth stage
    !NetPriProdLeafAdd = max(0.0, CarbohydrFracToLeaf(I,PlantGrowStage,J)*CarbohydrAssim - GrowthRespLeaf - RespirationLeafMaint) ! gCH2O/m2/s
    NetPriProdLeafAdd = CarbohydrFracToLeaf(I,PlantGrowStage(I,J),J)*CarbohydrAssim(I,J) - GrowthRespLeaf(I,J) - RespirationLeafMaint(I,J)            ! gCH2O/m2/s
    !NetPriProdStemAdd = max(0.0, CarbohydrFracToStem(I,PlantGrowStage,J)*CarbohydrAssim - GrowthRespStem - RespirationStem)      ! gCH2O/m2/s
    NetPriProdStemAdd = CarbohydrFracToStem(I,PlantGrowStage(I,J),J)*CarbohydrAssim(I,J) - GrowthRespStem(I,J) - RespirationStem(I,J)                 ! gCH2O/m2/s
    
    ! avoid reducing leaf mass below its minimum value but conserve mass
    LeafMassMaxChg(I,J) = (LeafMass(I,J) - LeafMassMin(I,J)) / MainTimeStep                         ! gCH2O/m2/s
    StemMassMaxChg(I,J) = (StemMass(I,J) - StemMassMin(I,J)) / MainTimeStep                         ! gCH2O/m2/s
    TurnoverLeaf(I,J)   = min(TurnoverLeaf(I,J), LeafMassMaxChg(I,J)+NetPriProdLeafAdd)             ! gCH2O/m2/s
    TurnoverStem(I,J)   = min(TurnoverStem(I,J), StemMassMaxChg(I,J)+NetPriProdStemAdd)             ! gCH2O/m2/s
    DeathLeaf(I,J)      = min(DeathLeaf(I,J), LeafMassMaxChg(I,J)+NetPriProdLeafAdd-TurnoverLeaf(I,J))   ! gCH2O/m2/s

    ! net primary productivities
    !NetPriProductionLeaf  = max(NetPriProdLeafAdd, -LeafMassMaxChg)    ! gCH2O/m2/s
    NetPriProductionLeaf(I,J)  = NetPriProdLeafAdd                           ! gCH2O/m2/s
    !NetPriProductionStem  = max(NetPriProdStemAdd, -StemMassMaxChg)    ! gCH2O/m2/s
    NetPriProductionStem(I,J)  = NetPriProdStemAdd                           ! gCH2O/m2/s
    NetPriProductionRoot(I,J)  = CarbohydrFracToRoot(I,PlantGrowStage(I,J),J) * CarbohydrAssim(I,J) - RespirationRoot(I,J) - GrowthRespRoot(I,J)     ! gCH2O/m2/s
    NetPriProductionGrain(I,J) = CarbohydrFracToGrain(I,PlantGrowStage(I,J),J) * CarbohydrAssim(I,J) - RespirationGrain(I,J) - GrowthRespGrain(I,J)  ! gCH2O/m2/s

    ! masses of plant components
    LeafMass(I,J)           = LeafMass(I,J) + (NetPriProductionLeaf(I,J) - TurnoverLeaf(I,J) - DeathLeaf(I,J)) * MainTimeStep ! gCH2O/m2
    StemMass(I,J)           = StemMass(I,J) + (NetPriProductionStem(I,J) - TurnoverStem(I,J)) * MainTimeStep             ! gCH2O/m2
    RootMass(I,J)           = RootMass(I,J) + (NetPriProductionRoot(I,J) - TurnoverRoot(I,J)) * MainTimeStep             ! gCH2O/m2
    GrainMass(I,J)          = GrainMass(I,J) + NetPriProductionGrain(I,J) * MainTimeStep                            ! gCH2O/m2
    GrossPriProduction(I,J) = CarbohydrAssim(I,J) * 0.4    ! gC/m2/s  0.4=12/30, CH20 to C

    ! carbon convert to grain ! Zhe Zhang 2020-07-13
    ConvLeafToGrain(I,J) = 0.0
    ConvStemToGrain(I,J) = 0.0
    ConvRootToGrain(I,J) = 0.0
    ConvLeafToGrain(I,J) = LeafMass(I,J) * (CarbohydrLeafToGrain(I,PlantGrowStage(I,J),J) * MainTimeStep / 3600.0)      ! gCH2O/m2
    ConvStemToGrain(I,J) = StemMass(I,J) * (CarbohydrStemToGrain(I,PlantGrowStage(I,J),J) * MainTimeStep / 3600.0)      ! gCH2O/m2
    ConvRootToGrain(I,J) = RootMass(I,J) * (CarbohydrRootToGrain(I,PlantGrowStage(I,J),J) * MainTimeStep / 3600.0)      ! gCH2O/m2
    LeafMass(I,J)        = LeafMass(I,J) - ConvLeafToGrain(I,J)                                                     ! gCH2O/m2
    StemMass(I,J)        = StemMass(I,J) - ConvStemToGrain(I,J)                                                     ! gCH2O/m2
    RootMass(I,J)        = RootMass(I,J) - ConvRootToGrain(I,J)                                                     ! gCH2O/m2
    GrainMass(I,J)       = GrainMass(I,J) + ConvStemToGrain(I,J) + ConvRootToGrain(I,J) + ConvLeafToGrain(I,J)                ! gCH2O/m2
    !if ( PlantGrowStage==6 ) then
    !   ConvStemToGrain = StemMass * (0.00005 * MainTimeStep / 3600.0)    ! gCH2O/m2
    !   StemMass        = StemMass - ConvStemToGrain                      ! gCH2O/m2
    !   ConvRootToGrain = RootMass * (0.0005 * MainTimeStep / 3600.0)     ! gCH2O/m2
    !   RootMass        = RootMass - ConvRootToGrain                      ! gCH2O/m2
    !   GrainMass       = GrainMass + ConvStemToGrain + ConvRootToGrain   ! gCH2O/m2
    !endif
    
    if ( RootMass(I,J) < 0.0 ) then
       TurnoverRoot(I,J) = NetPriProductionRoot(I,J)
       RootMass(I,J)     = 0.0
    endif
    if ( GrainMass(I,J) < 0.0 ) then
       GrainMass(I,J)    = 0.0
    endif

    ! soil carbon budgets
    !if ( (PlantGrowStage == 1) .or. (PlantGrowStage == 2) .or. (PlantGrowStage == 8) ) then
    !   CarbonMassShallowSoil = 1000
    !else
    CarbonMassShallowSoil(I,J) = CarbonMassShallowSoil(I,J) + &
                            (TurnoverRoot(I,J)+TurnoverLeaf(I,J)+TurnoverStem(I,J)+DeathLeaf(I,J)) * MainTimeStep * 0.4  ! 0.4: gCH2O/m2 -> gC/m2 
    !endif
    MicroRespFactorSoilTemp(I,J)  = 2.0**((TemperatureSoilSnow(I,1,J) - 283.16) / 10.0)
    MicroRespFactorSoilWater(I,J) = SoilWaterRootZone(I,J) / (0.20 + SoilWaterRootZone(I,J)) * 0.23 / (0.23 + SoilWaterRootZone(I,J))
    RespirationSoil(I,J)          = MicroRespFactorSoilWater(I,J) * MicroRespFactorSoilTemp(I,J) * &
                               MicroRespCoeff(I,J) * max(0.0, CarbonMassShallowSoil(I,J)*1.0e-3) * 30.0e-6     ! gCH2O/m2/s
    CarbonDecayToStable(I,J)      = 0.1 * RespirationSoil(I,J)                                                 ! gCH2O/m2/s
    CarbonMassShallowSoil(I,J)    = CarbonMassShallowSoil(I,J) - (RespirationSoil(I,J) + CarbonDecayToStable(I,J)) * MainTimeStep * 0.4    ! 0.4: gCH2O/m2 -> gC/m2
    CarbonMassDeepSoil(I,J)       = CarbonMassDeepSoil(I,J) + CarbonDecayToStable(I,J) * MainTimeStep * 0.4                           ! 0.4: gCH2O/m2 -> gC/m2
 
    !  total carbon flux
    CarbonToAtmos(I,J)  = - CarbonAssim(I,J) + (RespirationLeafMaint(I,J) + RespirationRoot(I,J) + RespirationStem(I,J) + RespirationGrain(I,J) + &
                       0.9*RespirationSoil(I,J) + GrowthRespLeaf(I,J) + GrowthRespRoot(I,J) + GrowthRespStem(I,J) + GrowthRespGrain(I,J)) * 0.4 ! gC/m2/s 0.4=12/30, CH20 to C

    ! for outputs
    NetPriProductionTot(I,J) = (NetPriProductionLeaf(I,J) + NetPriProductionStem(I,J) + &
                           NetPriProductionRoot(I,J) + NetPriProductionGrain(I,J)) * 0.4                              ! gC/m2/s  0.4=12/30, CH20 to C 
    RespirationPlantTot(I,J) = (RespirationRoot(I,J) + RespirationGrain(I,J) + RespirationLeafMaint(I,J) + RespirationStem(I,J) + &
                           GrowthRespLeaf(I,J) + GrowthRespRoot(I,J) + GrowthRespGrain(I,J) + GrowthRespStem(I,J)) * 0.4        ! gC/m2/s  0.4=12/30, CH20 to C
    RespirationSoilOrg(I,J)  = 0.9 * RespirationSoil(I,J) * 0.4                                                       ! gC/m2/s  0.4=12/30, CH20 to C
    NetEcoExchange(I,J)      = (RespirationPlantTot(I,J) + RespirationSoilOrg(I,J) - GrossPriProduction(I,J)) * 44.0 / 12.0     ! gCO2/m2/s
    CarbonMassSoilTot(I,J)   = CarbonMassShallowSoil(I,J) + CarbonMassDeepSoil(I,J)                                        ! gC/m2
    CarbonMassLiveTot(I,J)   = (LeafMass(I,J) + RootMass(I,J) + StemMass(I,J) + GrainMass(I,J)) * 0.4                                ! gC/m2 0.4=12/30, CH20 to C
 
    ! leaf area index and stem area index
    LeafAreaIndex(I,J) = 0.9!max(LeafMass(I,J)*LeafAreaPerBiomass(I,J), LeafAreaIndexMin(I,J))
    StemAreaIndex(I,J) = max(StemMass(I,J)*StemAreaPerMass(I,J), StemAreaIndexMin(I,J))
   
    ! After harversting
    !if ( PlantGrowStage == 8 ) then
    !   LeafMass  = 0.62
    !   StemMass  = 0.0
    !   GrainMass = 0.0
    !endif

    !if ( (PlantGrowStage == 1) .or. (PlantGrowStage == 2) .or. (PlantGrowStage == 8) ) then
    if ( (PlantGrowStage(I,J) == 8) .and. &
         ((GrainMass(I,J) > 0) .or. (LeafMass(I,J) > 0) .or. (StemMass(I,J) > 0) .or. (RootMass(I,J) > 0)) ) then
       LeafAreaIndex(I,J) = 0.7!0.05
       StemAreaIndex(I,J) = 0.05
       LeafMass(I,J)      = LeafMassMin(I,J)
       StemMass(I,J)      = StemMassMin(I,J)
       RootMass(I,J)      = 0.0
       GrainMass(I,J)     = 0.0
    endif


      end do
    end do
!$acc end parallel loop


    end associate

  end subroutine CarbonFluxCrop

end module CarbonFluxCropMod
