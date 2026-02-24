module BiochemVarInitMod

!!! Initialize column (1-D) Noah-MP biochemistry (carbon,nitrogen,etc) variables
!!! Biochemistry variables should be first defined in BiochemVarType.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpVarType

  implicit none

contains

!=== initialize with default values

  subroutine BiochemVarInitDefault(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer :: I, J      ! grid indices
    integer :: LoopInd   ! loop index

    ! Allocate 3D crop parameter arrays and transfer to GPU
    associate(                                                                      &
              NumCropGrowStage => noahmp%config%domain%NumCropGrowStage ,&
               ITS => noahmp%config%domain%ITS, ITE => noahmp%config%domain%ITE ,&
               JTS => noahmp%config%domain%JTS, JTE => noahmp%config%domain%JTE  &
             )


    if ( .not. allocated(noahmp%biochem%param%LeafDeathTempCoeffCrop) ) then
       allocate( noahmp%biochem%param%LeafDeathTempCoeffCrop(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafDeathWaterCoeffCrop) ) then
       allocate( noahmp%biochem%param%LeafDeathWaterCoeffCrop(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrLeafToGrain) ) then
       allocate( noahmp%biochem%param%CarbohydrLeafToGrain(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrStemToGrain) ) then
       allocate( noahmp%biochem%param%CarbohydrStemToGrain(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrRootToGrain) ) then
       allocate( noahmp%biochem%param%CarbohydrRootToGrain(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrFracToLeaf) ) then
       allocate( noahmp%biochem%param%CarbohydrFracToLeaf(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrFracToStem) ) then
       allocate( noahmp%biochem%param%CarbohydrFracToStem(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrFracToRoot) ) then
       allocate( noahmp%biochem%param%CarbohydrFracToRoot(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrFracToGrain) ) then
       allocate( noahmp%biochem%param%CarbohydrFracToGrain(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TurnoverCoeffLeafCrop) ) then
       allocate( noahmp%biochem%param%TurnoverCoeffLeafCrop(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TurnoverCoeffStemCrop) ) then
       allocate( noahmp%biochem%param%TurnoverCoeffStemCrop(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TurnoverCoeffRootCrop) ) then
       allocate( noahmp%biochem%param%TurnoverCoeffRootCrop(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
    endif

    ! Allocate 2D biochem state arrays and transfer to GPU
    if ( .not. allocated(noahmp%biochem%state%PlantGrowStage) ) then
       allocate( noahmp%biochem%state%PlantGrowStage(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%IndexPlanting) ) then
       allocate( noahmp%biochem%state%IndexPlanting(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%IndexHarvest) ) then
       allocate( noahmp%biochem%state%IndexHarvest(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%IndexGrowSeason) ) then
       allocate( noahmp%biochem%state%IndexGrowSeason(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%NitrogenConcFoliage) ) then
       allocate( noahmp%biochem%state%NitrogenConcFoliage(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%LeafMass) ) then
       allocate( noahmp%biochem%state%LeafMass(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%RootMass) ) then
       allocate( noahmp%biochem%state%RootMass(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%StemMass) ) then
       allocate( noahmp%biochem%state%StemMass(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%WoodMass) ) then
       allocate( noahmp%biochem%state%WoodMass(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%GrainMass) ) then
       allocate( noahmp%biochem%state%GrainMass(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonMassDeepSoil) ) then
       allocate( noahmp%biochem%state%CarbonMassDeepSoil(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonMassShallowSoil) ) then
       allocate( noahmp%biochem%state%CarbonMassShallowSoil(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonMassSoilTot) ) then
       allocate( noahmp%biochem%state%CarbonMassSoilTot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonMassLiveTot) ) then
       allocate( noahmp%biochem%state%CarbonMassLiveTot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%LeafAreaPerMass) ) then
       allocate( noahmp%biochem%state%LeafAreaPerMass(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%StemAreaPerMass) ) then
       allocate( noahmp%biochem%state%StemAreaPerMass(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%LeafMassMin) ) then
       allocate( noahmp%biochem%state%LeafMassMin(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%StemMassMin) ) then
       allocate( noahmp%biochem%state%StemMassMin(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonFracToLeaf) ) then
       allocate( noahmp%biochem%state%CarbonFracToLeaf(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonFracToRoot) ) then
       allocate( noahmp%biochem%state%CarbonFracToRoot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonFracToWood) ) then
       allocate( noahmp%biochem%state%CarbonFracToWood(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonFracToStem) ) then
       allocate( noahmp%biochem%state%CarbonFracToStem(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%WoodCarbonFrac) ) then
       allocate( noahmp%biochem%state%WoodCarbonFrac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonFracToWoodRoot) ) then
       allocate( noahmp%biochem%state%CarbonFracToWoodRoot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%MicroRespFactorSoilWater) ) then
       allocate( noahmp%biochem%state%MicroRespFactorSoilWater(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%MicroRespFactorSoilTemp) ) then
       allocate( noahmp%biochem%state%MicroRespFactorSoilTemp(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%RespFacNitrogenFoliage) ) then
       allocate( noahmp%biochem%state%RespFacNitrogenFoliage(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%RespFacTemperature) ) then
       allocate( noahmp%biochem%state%RespFacTemperature(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%RespReductionFac) ) then
       allocate( noahmp%biochem%state%RespReductionFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%state%GrowDegreeDay) ) then
       allocate( noahmp%biochem%state%GrowDegreeDay(ITS:ITE,JTS:JTE) )
    endif

    ! Allocate 2D biochem flux arrays and transfer to GPU
    if ( .not. allocated(noahmp%biochem%flux%PhotosynTotal) ) then
       allocate( noahmp%biochem%flux%PhotosynTotal(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%PhotosynLeafSunlit) ) then
       allocate( noahmp%biochem%flux%PhotosynLeafSunlit(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%PhotosynLeafShade) ) then
       allocate( noahmp%biochem%flux%PhotosynLeafShade(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%PhotosynCrop) ) then
       allocate( noahmp%biochem%flux%PhotosynCrop(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrossPriProduction) ) then
       allocate( noahmp%biochem%flux%GrossPriProduction(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetEcoExchange) ) then
       allocate( noahmp%biochem%flux%NetEcoExchange(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionTot) ) then
       allocate( noahmp%biochem%flux%NetPriProductionTot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionLeaf) ) then
       allocate( noahmp%biochem%flux%NetPriProductionLeaf(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionRoot) ) then
       allocate( noahmp%biochem%flux%NetPriProductionRoot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionWood) ) then
       allocate( noahmp%biochem%flux%NetPriProductionWood(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionStem) ) then
       allocate( noahmp%biochem%flux%NetPriProductionStem(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionGrain) ) then
       allocate( noahmp%biochem%flux%NetPriProductionGrain(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationPlantTot) ) then
       allocate( noahmp%biochem%flux%RespirationPlantTot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationSoilOrg) ) then
       allocate( noahmp%biochem%flux%RespirationSoilOrg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%CarbonToAtmos) ) then
       allocate( noahmp%biochem%flux%CarbonToAtmos(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrowthRespLeaf) ) then
       allocate( noahmp%biochem%flux%GrowthRespLeaf(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrowthRespRoot) ) then
       allocate( noahmp%biochem%flux%GrowthRespRoot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrowthRespWood) ) then
       allocate( noahmp%biochem%flux%GrowthRespWood(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrowthRespStem) ) then
       allocate( noahmp%biochem%flux%GrowthRespStem(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrowthRespGrain) ) then
       allocate( noahmp%biochem%flux%GrowthRespGrain(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%LeafMassMaxChg) ) then
       allocate( noahmp%biochem%flux%LeafMassMaxChg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%StemMassMaxChg) ) then
       allocate( noahmp%biochem%flux%StemMassMaxChg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%CarbonDecayToStable) ) then
       allocate( noahmp%biochem%flux%CarbonDecayToStable(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationLeaf) ) then
       allocate( noahmp%biochem%flux%RespirationLeaf(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationStem) ) then
       allocate( noahmp%biochem%flux%RespirationStem(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationWood) ) then
       allocate( noahmp%biochem%flux%RespirationWood(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationLeafMaint) ) then
       allocate( noahmp%biochem%flux%RespirationLeafMaint(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationRoot) ) then
       allocate( noahmp%biochem%flux%RespirationRoot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationSoil) ) then
       allocate( noahmp%biochem%flux%RespirationSoil(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationGrain) ) then
       allocate( noahmp%biochem%flux%RespirationGrain(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%ConvRootToGrain) ) then
       allocate( noahmp%biochem%flux%ConvRootToGrain(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%ConvStemToGrain) ) then
       allocate( noahmp%biochem%flux%ConvStemToGrain(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%ConvLeafToGrain) ) then
       allocate( noahmp%biochem%flux%ConvLeafToGrain(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%TurnoverLeaf) ) then
       allocate( noahmp%biochem%flux%TurnoverLeaf(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%TurnoverStem) ) then
       allocate( noahmp%biochem%flux%TurnoverStem(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%TurnoverWood) ) then
       allocate( noahmp%biochem%flux%TurnoverWood(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%TurnoverRoot) ) then
       allocate( noahmp%biochem%flux%TurnoverRoot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%TurnoverGrain) ) then
       allocate( noahmp%biochem%flux%TurnoverGrain(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%DeathLeaf) ) then
       allocate( noahmp%biochem%flux%DeathLeaf(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%DeathStem) ) then
       allocate( noahmp%biochem%flux%DeathStem(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%CarbonAssim) ) then
       allocate( noahmp%biochem%flux%CarbonAssim(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%flux%CarbohydrAssim) ) then
       allocate( noahmp%biochem%flux%CarbohydrAssim(ITS:ITE,JTS:JTE) )
    endif

    ! Allocate 2D biochem param arrays and transfer to GPU
    if ( .not. allocated(noahmp%biochem%param%DatePlanting) ) then
       allocate( noahmp%biochem%param%DatePlanting(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%DateHarvest) ) then
       allocate( noahmp%biochem%param%DateHarvest(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%QuantumEfficiency25C) ) then
       allocate( noahmp%biochem%param%QuantumEfficiency25C(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%CarboxylRateMax25C) ) then
       allocate( noahmp%biochem%param%CarboxylRateMax25C(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%CarboxylRateMaxQ10) ) then
       allocate( noahmp%biochem%param%CarboxylRateMaxQ10(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%PhotosynPathC3) ) then
       allocate( noahmp%biochem%param%PhotosynPathC3(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%SlopeConductToPhotosyn) ) then
       allocate( noahmp%biochem%param%SlopeConductToPhotosyn(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TemperatureMinPhotosyn) ) then
       allocate( noahmp%biochem%param%TemperatureMinPhotosyn(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafAreaPerMass1side) ) then
       allocate( noahmp%biochem%param%LeafAreaPerMass1side(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%NitrogenConcFoliageMax) ) then
       allocate( noahmp%biochem%param%NitrogenConcFoliageMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%WoodToRootRatio) ) then
       allocate( noahmp%biochem%param%WoodToRootRatio(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%WoodPoolIndex) ) then
       allocate( noahmp%biochem%param%WoodPoolIndex(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TurnoverCoeffLeafVeg) ) then
       allocate( noahmp%biochem%param%TurnoverCoeffLeafVeg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafDeathWaterCoeffVeg) ) then
       allocate( noahmp%biochem%param%LeafDeathWaterCoeffVeg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafDeathTempCoeffVeg) ) then
       allocate( noahmp%biochem%param%LeafDeathTempCoeffVeg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%MicroRespCoeff) ) then
       allocate( noahmp%biochem%param%MicroRespCoeff(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%RespMaintQ10) ) then
       allocate( noahmp%biochem%param%RespMaintQ10(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%RespMaintLeaf25C) ) then
       allocate( noahmp%biochem%param%RespMaintLeaf25C(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%RespMaintStem25C) ) then
       allocate( noahmp%biochem%param%RespMaintStem25C(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%RespMaintRoot25C) ) then
       allocate( noahmp%biochem%param%RespMaintRoot25C(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%RespMaintGrain25C) ) then
       allocate( noahmp%biochem%param%RespMaintGrain25C(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowthRespFrac) ) then
       allocate( noahmp%biochem%param%GrowthRespFrac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TemperaureLeafFreeze) ) then
       allocate( noahmp%biochem%param%TemperaureLeafFreeze(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafAreaPerBiomass) ) then
       allocate( noahmp%biochem%param%LeafAreaPerBiomass(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TempBaseGrowDegDay) ) then
       allocate( noahmp%biochem%param%TempBaseGrowDegDay(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TempMaxGrowDegDay) ) then
       allocate( noahmp%biochem%param%TempMaxGrowDegDay(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowDegDayEmerg) ) then
       allocate( noahmp%biochem%param%GrowDegDayEmerg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowDegDayInitVeg) ) then
       allocate( noahmp%biochem%param%GrowDegDayInitVeg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowDegDayPostVeg) ) then
       allocate( noahmp%biochem%param%GrowDegDayPostVeg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowDegDayInitReprod) ) then
       allocate( noahmp%biochem%param%GrowDegDayInitReprod(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowDegDayMature) ) then
       allocate( noahmp%biochem%param%GrowDegDayMature(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%PhotosynRadFrac) ) then
       allocate( noahmp%biochem%param%PhotosynRadFrac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TempMinCarbonAssim) ) then
       allocate( noahmp%biochem%param%TempMinCarbonAssim(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TempMaxCarbonAssim) ) then
       allocate( noahmp%biochem%param%TempMaxCarbonAssim(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TempMaxCarbonAssimMax) ) then
       allocate( noahmp%biochem%param%TempMaxCarbonAssimMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbonAssimRefMax) ) then
       allocate( noahmp%biochem%param%CarbonAssimRefMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%LightExtCoeff) ) then
       allocate( noahmp%biochem%param%LightExtCoeff(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%LightUseEfficiency) ) then
       allocate( noahmp%biochem%param%LightUseEfficiency(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbonAssimReducFac) ) then
       allocate( noahmp%biochem%param%CarbonAssimReducFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%StemAreaIndexMin) ) then
       allocate( noahmp%biochem%param%StemAreaIndexMin(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%WoodAllocFac) ) then
       allocate( noahmp%biochem%param%WoodAllocFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%WaterStressCoeff) ) then
       allocate( noahmp%biochem%param%WaterStressCoeff(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafAreaIndexMin) ) then
       allocate( noahmp%biochem%param%LeafAreaIndexMin(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%TurnoverCoeffRootVeg) ) then
       allocate( noahmp%biochem%param%TurnoverCoeffRootVeg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%biochem%param%WoodRespCoeff) ) then
       allocate( noahmp%biochem%param%WoodRespCoeff(ITS:ITE,JTS:JTE) )
    endif
    end associate

    ! Allocate 3D crop parameter arrays and transfer to GPU
    associate(                                                                      &
              NumCropGrowStage => noahmp%config%domain%NumCropGrowStage ,&
                  PlantGrowStage           => noahmp%biochem%state%PlantGrowStage ,&
                  IndexPlanting            => noahmp%biochem%state%IndexPlanting ,&
                  IndexHarvest             => noahmp%biochem%state%IndexHarvest ,&
                  IndexGrowSeason          => noahmp%biochem%state%IndexGrowSeason ,&
                  NitrogenConcFoliage      => noahmp%biochem%state%NitrogenConcFoliage ,&
                  LeafMass                 => noahmp%biochem%state%LeafMass ,&
                  RootMass                 => noahmp%biochem%state%RootMass ,&
                  StemMass                 => noahmp%biochem%state%StemMass ,&
                  WoodMass                 => noahmp%biochem%state%WoodMass ,&
                  CarbonMassDeepSoil       => noahmp%biochem%state%CarbonMassDeepSoil ,&
                  CarbonMassShallowSoil    => noahmp%biochem%state%CarbonMassShallowSoil ,&
                  CarbonMassSoilTot        => noahmp%biochem%state%CarbonMassSoilTot ,&
                  CarbonMassLiveTot        => noahmp%biochem%state%CarbonMassLiveTot ,&
                  LeafAreaPerMass          => noahmp%biochem%state%LeafAreaPerMass ,&
                  StemAreaPerMass          => noahmp%biochem%state%StemAreaPerMass ,&
                  LeafMassMin              => noahmp%biochem%state%LeafMassMin ,&
                  StemMassMin              => noahmp%biochem%state%StemMassMin ,&
                  CarbonFracToLeaf         => noahmp%biochem%state%CarbonFracToLeaf ,&
                  CarbonFracToRoot         => noahmp%biochem%state%CarbonFracToRoot ,&
                  CarbonFracToWood         => noahmp%biochem%state%CarbonFracToWood ,&
                  CarbonFracToStem         => noahmp%biochem%state%CarbonFracToStem ,&
                  WoodCarbonFrac           => noahmp%biochem%state%WoodCarbonFrac ,&
                  CarbonFracToWoodRoot     => noahmp%biochem%state%CarbonFracToWoodRoot ,&
                  MicroRespFactorSoilWater => noahmp%biochem%state%MicroRespFactorSoilWater ,&
                  MicroRespFactorSoilTemp  => noahmp%biochem%state%MicroRespFactorSoilTemp ,&
                  RespFacNitrogenFoliage   => noahmp%biochem%state%RespFacNitrogenFoliage ,&
                  RespFacTemperature       => noahmp%biochem%state%RespFacTemperature ,&
                  RespReductionFac         => noahmp%biochem%state%RespReductionFac ,&
                  GrainMass                => noahmp%biochem%state%GrainMass ,&
                  GrowDegreeDay            => noahmp%biochem%state%GrowDegreeDay ,&
                  PhotosynLeafSunlit       => noahmp%biochem%flux%PhotosynLeafSunlit ,&
                  PhotosynLeafShade        => noahmp%biochem%flux%PhotosynLeafShade ,&
                  PhotosynCrop             => noahmp%biochem%flux%PhotosynCrop ,&
                  PhotosynTotal            => noahmp%biochem%flux%PhotosynTotal ,&
                  GrossPriProduction       => noahmp%biochem%flux%GrossPriProduction ,&
                  NetPriProductionTot      => noahmp%biochem%flux%NetPriProductionTot ,&
                  NetEcoExchange           => noahmp%biochem%flux%NetEcoExchange ,&
                  RespirationPlantTot      => noahmp%biochem%flux%RespirationPlantTot ,&
                  RespirationSoilOrg       => noahmp%biochem%flux%RespirationSoilOrg ,&
                  CarbonToAtmos            => noahmp%biochem%flux%CarbonToAtmos ,&
                  NetPriProductionLeaf     => noahmp%biochem%flux%NetPriProductionLeaf ,&
                  NetPriProductionRoot     => noahmp%biochem%flux%NetPriProductionRoot ,&
                  NetPriProductionWood     => noahmp%biochem%flux%NetPriProductionWood ,&
                  NetPriProductionStem     => noahmp%biochem%flux%NetPriProductionStem ,&
                  GrowthRespLeaf           => noahmp%biochem%flux%GrowthRespLeaf ,&
                  GrowthRespRoot           => noahmp%biochem%flux%GrowthRespRoot ,&
                  GrowthRespWood           => noahmp%biochem%flux%GrowthRespWood ,&
                  GrowthRespStem           => noahmp%biochem%flux%GrowthRespStem ,&
                  LeafMassMaxChg           => noahmp%biochem%flux%LeafMassMaxChg ,&
                  StemMassMaxChg           => noahmp%biochem%flux%StemMassMaxChg ,&
                  CarbonDecayToStable      => noahmp%biochem%flux%CarbonDecayToStable ,&
                  RespirationLeaf          => noahmp%biochem%flux%RespirationLeaf ,&
                  RespirationStem          => noahmp%biochem%flux%RespirationStem ,&
                  GrowthRespGrain          => noahmp%biochem%flux%GrowthRespGrain ,&
                  NetPriProductionGrain    => noahmp%biochem%flux%NetPriProductionGrain ,&
                  ConvRootToGrain          => noahmp%biochem%flux%ConvRootToGrain ,&
                  ConvStemToGrain          => noahmp%biochem%flux%ConvStemToGrain ,&
                  RespirationWood          => noahmp%biochem%flux%RespirationWood ,&
                  RespirationLeafMaint     => noahmp%biochem%flux%RespirationLeafMaint ,&
                  RespirationRoot          => noahmp%biochem%flux%RespirationRoot ,&
                  DeathLeaf                => noahmp%biochem%flux%DeathLeaf ,&
                  DeathStem                => noahmp%biochem%flux%DeathStem ,&
                  CarbonAssim              => noahmp%biochem%flux%CarbonAssim ,&
                  TurnoverLeaf             => noahmp%biochem%flux%TurnoverLeaf ,&
                  TurnoverStem             => noahmp%biochem%flux%TurnoverStem ,&
                  TurnoverWood             => noahmp%biochem%flux%TurnoverWood ,&
                  RespirationSoil          => noahmp%biochem%flux%RespirationSoil ,&
                  TurnoverRoot             => noahmp%biochem%flux%TurnoverRoot ,&
                  CarbohydrAssim           => noahmp%biochem%flux%CarbohydrAssim ,&
                  TurnoverGrain            => noahmp%biochem%flux%TurnoverGrain ,&
                  ConvLeafToGrain          => noahmp%biochem%flux%ConvLeafToGrain ,&
                  RespirationGrain         => noahmp%biochem%flux%RespirationGrain ,&
                  DatePlanting             => noahmp%biochem%param%DatePlanting ,&
                  DateHarvest              => noahmp%biochem%param%DateHarvest ,&
                  QuantumEfficiency25C     => noahmp%biochem%param%QuantumEfficiency25C ,&
                  CarboxylRateMax25C       => noahmp%biochem%param%CarboxylRateMax25C ,&
                  CarboxylRateMaxQ10       => noahmp%biochem%param%CarboxylRateMaxQ10 ,&
                  PhotosynPathC3           => noahmp%biochem%param%PhotosynPathC3 ,&
                  SlopeConductToPhotosyn   => noahmp%biochem%param%SlopeConductToPhotosyn ,&
                  TemperatureMinPhotosyn   => noahmp%biochem%param%TemperatureMinPhotosyn ,&
                  LeafAreaPerMass1side     => noahmp%biochem%param%LeafAreaPerMass1side ,&
                  NitrogenConcFoliageMax   => noahmp%biochem%param%NitrogenConcFoliageMax ,&
                  WoodToRootRatio          => noahmp%biochem%param%WoodToRootRatio ,&
                  WoodPoolIndex            => noahmp%biochem%param%WoodPoolIndex ,&
                  TurnoverCoeffLeafVeg     => noahmp%biochem%param%TurnoverCoeffLeafVeg ,&
                  LeafDeathWaterCoeffVeg   => noahmp%biochem%param%LeafDeathWaterCoeffVeg ,&
                  LeafDeathTempCoeffVeg    => noahmp%biochem%param%LeafDeathTempCoeffVeg ,&
                  MicroRespCoeff           => noahmp%biochem%param%MicroRespCoeff ,&
                  RespMaintQ10             => noahmp%biochem%param%RespMaintQ10 ,&
                  RespMaintLeaf25C         => noahmp%biochem%param%RespMaintLeaf25C ,&
                  RespMaintStem25C         => noahmp%biochem%param%RespMaintStem25C ,&
                  RespMaintRoot25C         => noahmp%biochem%param%RespMaintRoot25C ,&
                  RespMaintGrain25C        => noahmp%biochem%param%RespMaintGrain25C ,&
                  GrowthRespFrac           => noahmp%biochem%param%GrowthRespFrac ,&
                  TemperaureLeafFreeze     => noahmp%biochem%param%TemperaureLeafFreeze ,&
                  LeafAreaPerBiomass       => noahmp%biochem%param%LeafAreaPerBiomass ,&
                  TempBaseGrowDegDay       => noahmp%biochem%param%TempBaseGrowDegDay ,&
                  TempMaxGrowDegDay        => noahmp%biochem%param%TempMaxGrowDegDay ,&
                  GrowDegDayEmerg          => noahmp%biochem%param%GrowDegDayEmerg ,&
                  GrowDegDayInitVeg        => noahmp%biochem%param%GrowDegDayInitVeg ,&
                  GrowDegDayPostVeg        => noahmp%biochem%param%GrowDegDayPostVeg ,&
                  GrowDegDayInitReprod     => noahmp%biochem%param%GrowDegDayInitReprod ,&
                  GrowDegDayMature         => noahmp%biochem%param%GrowDegDayMature ,&
                  PhotosynRadFrac          => noahmp%biochem%param%PhotosynRadFrac ,&
                  TempMinCarbonAssim       => noahmp%biochem%param%TempMinCarbonAssim ,&
                  TempMaxCarbonAssim       => noahmp%biochem%param%TempMaxCarbonAssim ,&
                  TempMaxCarbonAssimMax    => noahmp%biochem%param%TempMaxCarbonAssimMax ,&
                  CarbonAssimRefMax        => noahmp%biochem%param%CarbonAssimRefMax ,&
                  LightExtCoeff            => noahmp%biochem%param%LightExtCoeff ,&
                  LightUseEfficiency       => noahmp%biochem%param%LightUseEfficiency ,&
                  CarbonAssimReducFac      => noahmp%biochem%param%CarbonAssimReducFac ,&
                  StemAreaIndexMin         => noahmp%biochem%param%StemAreaIndexMin ,&
                  WoodAllocFac             => noahmp%biochem%param%WoodAllocFac ,&
                  WaterStressCoeff         => noahmp%biochem%param%WaterStressCoeff ,&
                  LeafAreaIndexMin         => noahmp%biochem%param%LeafAreaIndexMin ,&
                  TurnoverCoeffRootVeg     => noahmp%biochem%param%TurnoverCoeffRootVeg ,&
                  WoodRespCoeff            => noahmp%biochem%param%WoodRespCoeff  &
             )


    ! Initialize 2D state, flux, and parameter arrays (including 3D crop arrays)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%biochem%param%LeafDeathTempCoeffCrop,   &
    !$acc   noahmp%biochem%param%LeafDeathWaterCoeffCrop,   &
    !$acc   noahmp%biochem%param%CarbohydrLeafToGrain,   &
    !$acc   noahmp%biochem%param%CarbohydrStemToGrain,   &
    !$acc   noahmp%biochem%param%CarbohydrRootToGrain,   &
    !$acc   noahmp%biochem%param%CarbohydrFracToLeaf,   &
    !$acc   noahmp%biochem%param%CarbohydrFracToStem,   &
    !$acc   noahmp%biochem%param%CarbohydrFracToRoot,   &
    !$acc   noahmp%biochem%param%CarbohydrFracToGrain,   &
    !$acc   noahmp%biochem%param%TurnoverCoeffLeafCrop,   &
    !$acc   noahmp%biochem%param%TurnoverCoeffStemCrop,   &
    !$acc   noahmp%biochem%param%TurnoverCoeffRootCrop,   &
    !$acc   noahmp%biochem%state%PlantGrowStage,   &
    !$acc   noahmp%biochem%state%IndexPlanting,   &
    !$acc   noahmp%biochem%state%IndexHarvest,   &
    !$acc   noahmp%biochem%state%IndexGrowSeason,   &
    !$acc   noahmp%biochem%state%NitrogenConcFoliage,   &
    !$acc   noahmp%biochem%state%LeafMass,   &
    !$acc   noahmp%biochem%state%RootMass,   &
    !$acc   noahmp%biochem%state%StemMass,   &
    !$acc   noahmp%biochem%state%WoodMass,   &
    !$acc   noahmp%biochem%state%GrainMass,   &
    !$acc   noahmp%biochem%state%CarbonMassDeepSoil,   &
    !$acc   noahmp%biochem%state%CarbonMassShallowSoil,   &
    !$acc   noahmp%biochem%state%CarbonMassSoilTot,   &
    !$acc   noahmp%biochem%state%CarbonMassLiveTot,   &
    !$acc   noahmp%biochem%state%LeafAreaPerMass,   &
    !$acc   noahmp%biochem%state%StemAreaPerMass,   &
    !$acc   noahmp%biochem%state%LeafMassMin,   &
    !$acc   noahmp%biochem%state%StemMassMin    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%biochem%state%CarbonFracToLeaf,   &
    !$acc   noahmp%biochem%state%CarbonFracToRoot,   &
    !$acc   noahmp%biochem%state%CarbonFracToWood,   &
    !$acc   noahmp%biochem%state%CarbonFracToStem,   &
    !$acc   noahmp%biochem%state%WoodCarbonFrac,   &
    !$acc   noahmp%biochem%state%CarbonFracToWoodRoot,   &
    !$acc   noahmp%biochem%state%MicroRespFactorSoilWater,   &
    !$acc   noahmp%biochem%state%MicroRespFactorSoilTemp,   &
    !$acc   noahmp%biochem%state%RespFacNitrogenFoliage,   &
    !$acc   noahmp%biochem%state%RespFacTemperature,   &
    !$acc   noahmp%biochem%state%RespReductionFac,   &
    !$acc   noahmp%biochem%state%GrowDegreeDay,   &
    !$acc   noahmp%biochem%flux%PhotosynTotal,   &
    !$acc   noahmp%biochem%flux%PhotosynLeafSunlit,   &
    !$acc   noahmp%biochem%flux%PhotosynLeafShade,   &
    !$acc   noahmp%biochem%flux%PhotosynCrop,   &
    !$acc   noahmp%biochem%flux%GrossPriProduction,   &
    !$acc   noahmp%biochem%flux%NetEcoExchange,   &
    !$acc   noahmp%biochem%flux%NetPriProductionTot,   &
    !$acc   noahmp%biochem%flux%NetPriProductionLeaf,   &
    !$acc   noahmp%biochem%flux%NetPriProductionRoot,   &
    !$acc   noahmp%biochem%flux%NetPriProductionWood,   &
    !$acc   noahmp%biochem%flux%NetPriProductionStem,   &
    !$acc   noahmp%biochem%flux%NetPriProductionGrain,   &
    !$acc   noahmp%biochem%flux%RespirationPlantTot,   &
    !$acc   noahmp%biochem%flux%RespirationSoilOrg,   &
    !$acc   noahmp%biochem%flux%CarbonToAtmos,   &
    !$acc   noahmp%biochem%flux%GrowthRespLeaf,   &
    !$acc   noahmp%biochem%flux%GrowthRespRoot,   &
    !$acc   noahmp%biochem%flux%GrowthRespWood    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%biochem%flux%GrowthRespStem,   &
    !$acc   noahmp%biochem%flux%GrowthRespGrain,   &
    !$acc   noahmp%biochem%flux%LeafMassMaxChg,   &
    !$acc   noahmp%biochem%flux%StemMassMaxChg,   &
    !$acc   noahmp%biochem%flux%CarbonDecayToStable,   &
    !$acc   noahmp%biochem%flux%RespirationLeaf,   &
    !$acc   noahmp%biochem%flux%RespirationStem,   &
    !$acc   noahmp%biochem%flux%RespirationWood,   &
    !$acc   noahmp%biochem%flux%RespirationLeafMaint,   &
    !$acc   noahmp%biochem%flux%RespirationRoot,   &
    !$acc   noahmp%biochem%flux%RespirationSoil,   &
    !$acc   noahmp%biochem%flux%RespirationGrain,   &
    !$acc   noahmp%biochem%flux%ConvRootToGrain,   &
    !$acc   noahmp%biochem%flux%ConvStemToGrain,   &
    !$acc   noahmp%biochem%flux%ConvLeafToGrain,   &
    !$acc   noahmp%biochem%flux%TurnoverLeaf,   &
    !$acc   noahmp%biochem%flux%TurnoverStem,   &
    !$acc   noahmp%biochem%flux%TurnoverWood,   &
    !$acc   noahmp%biochem%flux%TurnoverRoot,   &
    !$acc   noahmp%biochem%flux%TurnoverGrain,   &
    !$acc   noahmp%biochem%flux%DeathLeaf,   &
    !$acc   noahmp%biochem%flux%DeathStem,   &
    !$acc   noahmp%biochem%flux%CarbonAssim,   &
    !$acc   noahmp%biochem%flux%CarbohydrAssim,   &
    !$acc   noahmp%biochem%param%DatePlanting,   &
    !$acc   noahmp%biochem%param%DateHarvest,   &
    !$acc   noahmp%biochem%param%QuantumEfficiency25C,   &
    !$acc   noahmp%biochem%param%CarboxylRateMax25C,   &
    !$acc   noahmp%biochem%param%CarboxylRateMaxQ10,   &
    !$acc   noahmp%biochem%param%PhotosynPathC3    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%biochem%param%SlopeConductToPhotosyn,   &
    !$acc   noahmp%biochem%param%TemperatureMinPhotosyn,   &
    !$acc   noahmp%biochem%param%LeafAreaPerMass1side,   &
    !$acc   noahmp%biochem%param%NitrogenConcFoliageMax,   &
    !$acc   noahmp%biochem%param%WoodToRootRatio,   &
    !$acc   noahmp%biochem%param%WoodPoolIndex,   &
    !$acc   noahmp%biochem%param%TurnoverCoeffLeafVeg,   &
    !$acc   noahmp%biochem%param%LeafDeathWaterCoeffVeg,   &
    !$acc   noahmp%biochem%param%LeafDeathTempCoeffVeg,   &
    !$acc   noahmp%biochem%param%MicroRespCoeff,   &
    !$acc   noahmp%biochem%param%RespMaintQ10,   &
    !$acc   noahmp%biochem%param%RespMaintLeaf25C,   &
    !$acc   noahmp%biochem%param%RespMaintStem25C,   &
    !$acc   noahmp%biochem%param%RespMaintRoot25C,   &
    !$acc   noahmp%biochem%param%RespMaintGrain25C,   &
    !$acc   noahmp%biochem%param%GrowthRespFrac,   &
    !$acc   noahmp%biochem%param%TemperaureLeafFreeze,   &
    !$acc   noahmp%biochem%param%LeafAreaPerBiomass,   &
    !$acc   noahmp%biochem%param%TempBaseGrowDegDay,   &
    !$acc   noahmp%biochem%param%TempMaxGrowDegDay,   &
    !$acc   noahmp%biochem%param%GrowDegDayEmerg,   &
    !$acc   noahmp%biochem%param%GrowDegDayInitVeg,   &
    !$acc   noahmp%biochem%param%GrowDegDayPostVeg,   &
    !$acc   noahmp%biochem%param%GrowDegDayInitReprod,   &
    !$acc   noahmp%biochem%param%GrowDegDayMature,   &
    !$acc   noahmp%biochem%param%PhotosynRadFrac,   &
    !$acc   noahmp%biochem%param%TempMinCarbonAssim,   &
    !$acc   noahmp%biochem%param%TempMaxCarbonAssim,   &
    !$acc   noahmp%biochem%param%TempMaxCarbonAssimMax,   &
    !$acc   noahmp%biochem%param%CarbonAssimRefMax    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%biochem%param%LightExtCoeff,   &
    !$acc   noahmp%biochem%param%LightUseEfficiency,   &
    !$acc   noahmp%biochem%param%CarbonAssimReducFac,   &
    !$acc   noahmp%biochem%param%StemAreaIndexMin,   &
    !$acc   noahmp%biochem%param%WoodAllocFac,   &
    !$acc   noahmp%biochem%param%WaterStressCoeff,   &
    !$acc   noahmp%biochem%param%LeafAreaIndexMin,   &
    !$acc   noahmp%biochem%param%TurnoverCoeffRootVeg,   &
    !$acc   noahmp%biochem%param%WoodRespCoeff    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    !$acc wait(NOAHMP_ACC_QUEUE)

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


        ! biochem state variables
        PlantGrowStage(I,J)             = undefined_int
        IndexPlanting(I,J)              = undefined_int
        IndexHarvest(I,J)               = undefined_int
        IndexGrowSeason(I,J)            = undefined_real
        NitrogenConcFoliage(I,J)        = undefined_real
        LeafMass(I,J)                   = undefined_real
        RootMass(I,J)                   = undefined_real
        StemMass(I,J)                   = undefined_real
        WoodMass(I,J)                   = undefined_real
        CarbonMassDeepSoil(I,J)         = undefined_real
        CarbonMassShallowSoil(I,J)      = undefined_real
        CarbonMassSoilTot(I,J)          = undefined_real
        CarbonMassLiveTot(I,J)          = undefined_real
        LeafAreaPerMass(I,J)            = undefined_real
        StemAreaPerMass(I,J)            = undefined_real
        LeafMassMin(I,J)                = undefined_real
        StemMassMin(I,J)                = undefined_real
        CarbonFracToLeaf(I,J)           = undefined_real
        CarbonFracToRoot(I,J)           = undefined_real
        CarbonFracToWood(I,J)           = undefined_real
        CarbonFracToStem(I,J)           = undefined_real
        WoodCarbonFrac(I,J)             = undefined_real
        CarbonFracToWoodRoot(I,J)       = undefined_real
        MicroRespFactorSoilWater(I,J)   = undefined_real
        MicroRespFactorSoilTemp(I,J)    = undefined_real
        RespFacNitrogenFoliage(I,J)     = undefined_real
        RespFacTemperature(I,J)         = undefined_real
        RespReductionFac(I,J)           = undefined_real
        GrainMass(I,J)                  = undefined_real
        GrowDegreeDay(I,J)              = undefined_real

        ! biochem flux variables
        PhotosynLeafSunlit(I,J)         = undefined_real
        PhotosynLeafShade(I,J)          = undefined_real
        PhotosynCrop(I,J)               = undefined_real
        PhotosynTotal(I,J)              = undefined_real
        GrossPriProduction(I,J)         = undefined_real
        NetPriProductionTot(I,J)        = undefined_real
        NetEcoExchange(I,J)             = undefined_real
        RespirationPlantTot(I,J)        = undefined_real
        RespirationSoilOrg(I,J)         = undefined_real
        CarbonToAtmos(I,J)              = undefined_real
        NetPriProductionLeaf(I,J)       = undefined_real
        NetPriProductionRoot(I,J)       = undefined_real
        NetPriProductionWood(I,J)       = undefined_real
        NetPriProductionStem(I,J)       = undefined_real
        GrowthRespLeaf(I,J)             = undefined_real
        GrowthRespRoot(I,J)             = undefined_real
        GrowthRespWood(I,J)             = undefined_real
        GrowthRespStem(I,J)             = undefined_real
        LeafMassMaxChg(I,J)             = undefined_real
        StemMassMaxChg(I,J)             = undefined_real
        CarbonDecayToStable(I,J)        = undefined_real
        RespirationLeaf(I,J)            = undefined_real
        RespirationStem(I,J)            = undefined_real
        GrowthRespGrain(I,J)            = undefined_real
        NetPriProductionGrain(I,J)      = undefined_real
        ConvRootToGrain(I,J)            = undefined_real
        ConvStemToGrain(I,J)            = undefined_real
        RespirationWood(I,J)            = undefined_real
        RespirationLeafMaint(I,J)       = undefined_real
        RespirationRoot(I,J)            = undefined_real
        DeathLeaf(I,J)                  = undefined_real
        DeathStem(I,J)                  = undefined_real
        CarbonAssim(I,J)                = undefined_real
        TurnoverLeaf(I,J)               = undefined_real
        TurnoverStem(I,J)               = undefined_real
        TurnoverWood(I,J)               = undefined_real
        RespirationSoil(I,J)            = undefined_real
        TurnoverRoot(I,J)               = undefined_real
        CarbohydrAssim(I,J)             = undefined_real
        TurnoverGrain(I,J)              = undefined_real
        ConvLeafToGrain(I,J)            = undefined_real
        RespirationGrain(I,J)           = undefined_real

        ! biochem parameter variables
        DatePlanting(I,J)               = undefined_int
        DateHarvest(I,J)                = undefined_int
        QuantumEfficiency25C(I,J)       = undefined_real
        CarboxylRateMax25C(I,J)         = undefined_real
        CarboxylRateMaxQ10(I,J)         = undefined_real
        PhotosynPathC3(I,J)             = undefined_real
        SlopeConductToPhotosyn(I,J)     = undefined_real
        TemperatureMinPhotosyn(I,J)     = undefined_real
        LeafAreaPerMass1side(I,J)       = undefined_real
        NitrogenConcFoliageMax(I,J)     = undefined_real
        WoodToRootRatio(I,J)            = undefined_real
        WoodPoolIndex(I,J)              = undefined_real
        TurnoverCoeffLeafVeg(I,J)       = undefined_real
        LeafDeathWaterCoeffVeg(I,J)     = undefined_real
        LeafDeathTempCoeffVeg(I,J)      = undefined_real
        MicroRespCoeff(I,J)             = undefined_real
        RespMaintQ10(I,J)               = undefined_real
        RespMaintLeaf25C(I,J)           = undefined_real
        RespMaintStem25C(I,J)           = undefined_real
        RespMaintRoot25C(I,J)           = undefined_real
        RespMaintGrain25C(I,J)          = undefined_real
        GrowthRespFrac(I,J)             = undefined_real
        TemperaureLeafFreeze(I,J)       = undefined_real
        LeafAreaPerBiomass(I,J)         = undefined_real
        TempBaseGrowDegDay(I,J)         = undefined_real
        TempMaxGrowDegDay(I,J)          = undefined_real
        GrowDegDayEmerg(I,J)            = undefined_real
        GrowDegDayInitVeg(I,J)          = undefined_real
        GrowDegDayPostVeg(I,J)          = undefined_real
        GrowDegDayInitReprod(I,J)       = undefined_real
        GrowDegDayMature(I,J)           = undefined_real
        PhotosynRadFrac(I,J)            = undefined_real
        TempMinCarbonAssim(I,J)         = undefined_real
        TempMaxCarbonAssim(I,J)         = undefined_real
        TempMaxCarbonAssimMax(I,J)      = undefined_real
        CarbonAssimRefMax(I,J)          = undefined_real
        LightExtCoeff(I,J)              = undefined_real
        LightUseEfficiency(I,J)         = undefined_real
        CarbonAssimReducFac(I,J)        = undefined_real
        StemAreaIndexMin(I,J)           = undefined_real
        WoodAllocFac(I,J)               = undefined_real
        WaterStressCoeff(I,J)           = undefined_real
        LeafAreaIndexMin(I,J)           = undefined_real
        TurnoverCoeffRootVeg(I,J)       = undefined_real
        WoodRespCoeff(I,J)              = undefined_real

        ! Initialize 3D crop parameter arrays (I,LoopInd,J)
        !$acc loop seq
        do LoopInd = 1, noahmp%config%domain%NumCropGrowStage
           noahmp%biochem%param%LeafDeathTempCoeffCrop (I,LoopInd,J) = undefined_real
           noahmp%biochem%param%LeafDeathWaterCoeffCrop(I,LoopInd,J) = undefined_real
           noahmp%biochem%param%CarbohydrLeafToGrain   (I,LoopInd,J) = undefined_real
           noahmp%biochem%param%CarbohydrStemToGrain   (I,LoopInd,J) = undefined_real
           noahmp%biochem%param%CarbohydrRootToGrain   (I,LoopInd,J) = undefined_real
           noahmp%biochem%param%CarbohydrFracToLeaf    (I,LoopInd,J) = undefined_real
           noahmp%biochem%param%CarbohydrFracToStem    (I,LoopInd,J) = undefined_real
           noahmp%biochem%param%CarbohydrFracToRoot    (I,LoopInd,J) = undefined_real
           noahmp%biochem%param%CarbohydrFracToGrain   (I,LoopInd,J) = undefined_real
           noahmp%biochem%param%TurnoverCoeffLeafCrop  (I,LoopInd,J) = undefined_real
           noahmp%biochem%param%TurnoverCoeffStemCrop  (I,LoopInd,J) = undefined_real
           noahmp%biochem%param%TurnoverCoeffRootCrop  (I,LoopInd,J) = undefined_real
        enddo


      end do
    end do
    !$acc end parallel loop


    end associate

  end subroutine BiochemVarInitDefault


  subroutine BiochemVarExitDevice(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    !$acc exit data delete(               &
    !$acc   noahmp%biochem%param%LeafDeathTempCoeffCrop,   &
    !$acc   noahmp%biochem%param%LeafDeathWaterCoeffCrop,   &
    !$acc   noahmp%biochem%param%CarbohydrLeafToGrain,   &
    !$acc   noahmp%biochem%param%CarbohydrStemToGrain,   &
    !$acc   noahmp%biochem%param%CarbohydrRootToGrain,   &
    !$acc   noahmp%biochem%param%CarbohydrFracToLeaf,   &
    !$acc   noahmp%biochem%param%CarbohydrFracToStem,   &
    !$acc   noahmp%biochem%param%CarbohydrFracToRoot,   &
    !$acc   noahmp%biochem%param%CarbohydrFracToGrain,   &
    !$acc   noahmp%biochem%param%TurnoverCoeffLeafCrop,   &
    !$acc   noahmp%biochem%param%TurnoverCoeffStemCrop,   &
    !$acc   noahmp%biochem%param%TurnoverCoeffRootCrop,   &
    !$acc   noahmp%biochem%state%PlantGrowStage,   &
    !$acc   noahmp%biochem%state%IndexPlanting,   &
    !$acc   noahmp%biochem%state%IndexHarvest,   &
    !$acc   noahmp%biochem%state%IndexGrowSeason,   &
    !$acc   noahmp%biochem%state%NitrogenConcFoliage,   &
    !$acc   noahmp%biochem%state%LeafMass,   &
    !$acc   noahmp%biochem%state%RootMass,   &
    !$acc   noahmp%biochem%state%StemMass,   &
    !$acc   noahmp%biochem%state%WoodMass,   &
    !$acc   noahmp%biochem%state%GrainMass,   &
    !$acc   noahmp%biochem%state%CarbonMassDeepSoil,   &
    !$acc   noahmp%biochem%state%CarbonMassShallowSoil,   &
    !$acc   noahmp%biochem%state%CarbonMassSoilTot,   &
    !$acc   noahmp%biochem%state%CarbonMassLiveTot,   &
    !$acc   noahmp%biochem%state%LeafAreaPerMass,   &
    !$acc   noahmp%biochem%state%StemAreaPerMass,   &
    !$acc   noahmp%biochem%state%LeafMassMin,   &
    !$acc   noahmp%biochem%state%StemMassMin    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%biochem%state%CarbonFracToLeaf,   &
    !$acc   noahmp%biochem%state%CarbonFracToRoot,   &
    !$acc   noahmp%biochem%state%CarbonFracToWood,   &
    !$acc   noahmp%biochem%state%CarbonFracToStem,   &
    !$acc   noahmp%biochem%state%WoodCarbonFrac,   &
    !$acc   noahmp%biochem%state%CarbonFracToWoodRoot,   &
    !$acc   noahmp%biochem%state%MicroRespFactorSoilWater,   &
    !$acc   noahmp%biochem%state%MicroRespFactorSoilTemp,   &
    !$acc   noahmp%biochem%state%RespFacNitrogenFoliage,   &
    !$acc   noahmp%biochem%state%RespFacTemperature,   &
    !$acc   noahmp%biochem%state%RespReductionFac,   &
    !$acc   noahmp%biochem%state%GrowDegreeDay,   &
    !$acc   noahmp%biochem%flux%PhotosynTotal,   &
    !$acc   noahmp%biochem%flux%PhotosynLeafSunlit,   &
    !$acc   noahmp%biochem%flux%PhotosynLeafShade,   &
    !$acc   noahmp%biochem%flux%PhotosynCrop,   &
    !$acc   noahmp%biochem%flux%GrossPriProduction,   &
    !$acc   noahmp%biochem%flux%NetEcoExchange,   &
    !$acc   noahmp%biochem%flux%NetPriProductionTot,   &
    !$acc   noahmp%biochem%flux%NetPriProductionLeaf,   &
    !$acc   noahmp%biochem%flux%NetPriProductionRoot,   &
    !$acc   noahmp%biochem%flux%NetPriProductionWood,   &
    !$acc   noahmp%biochem%flux%NetPriProductionStem,   &
    !$acc   noahmp%biochem%flux%NetPriProductionGrain,   &
    !$acc   noahmp%biochem%flux%RespirationPlantTot,   &
    !$acc   noahmp%biochem%flux%RespirationSoilOrg,   &
    !$acc   noahmp%biochem%flux%CarbonToAtmos,   &
    !$acc   noahmp%biochem%flux%GrowthRespLeaf,   &
    !$acc   noahmp%biochem%flux%GrowthRespRoot,   &
    !$acc   noahmp%biochem%flux%GrowthRespWood    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%biochem%flux%GrowthRespStem,   &
    !$acc   noahmp%biochem%flux%GrowthRespGrain,   &
    !$acc   noahmp%biochem%flux%LeafMassMaxChg,   &
    !$acc   noahmp%biochem%flux%StemMassMaxChg,   &
    !$acc   noahmp%biochem%flux%CarbonDecayToStable,   &
    !$acc   noahmp%biochem%flux%RespirationLeaf,   &
    !$acc   noahmp%biochem%flux%RespirationStem,   &
    !$acc   noahmp%biochem%flux%RespirationWood,   &
    !$acc   noahmp%biochem%flux%RespirationLeafMaint,   &
    !$acc   noahmp%biochem%flux%RespirationRoot,   &
    !$acc   noahmp%biochem%flux%RespirationSoil,   &
    !$acc   noahmp%biochem%flux%RespirationGrain,   &
    !$acc   noahmp%biochem%flux%ConvRootToGrain,   &
    !$acc   noahmp%biochem%flux%ConvStemToGrain,   &
    !$acc   noahmp%biochem%flux%ConvLeafToGrain,   &
    !$acc   noahmp%biochem%flux%TurnoverLeaf,   &
    !$acc   noahmp%biochem%flux%TurnoverStem,   &
    !$acc   noahmp%biochem%flux%TurnoverWood,   &
    !$acc   noahmp%biochem%flux%TurnoverRoot,   &
    !$acc   noahmp%biochem%flux%TurnoverGrain,   &
    !$acc   noahmp%biochem%flux%DeathLeaf,   &
    !$acc   noahmp%biochem%flux%DeathStem,   &
    !$acc   noahmp%biochem%flux%CarbonAssim,   &
    !$acc   noahmp%biochem%flux%CarbohydrAssim,   &
    !$acc   noahmp%biochem%param%DatePlanting,   &
    !$acc   noahmp%biochem%param%DateHarvest,   &
    !$acc   noahmp%biochem%param%QuantumEfficiency25C,   &
    !$acc   noahmp%biochem%param%CarboxylRateMax25C,   &
    !$acc   noahmp%biochem%param%CarboxylRateMaxQ10,   &
    !$acc   noahmp%biochem%param%PhotosynPathC3    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%biochem%param%SlopeConductToPhotosyn,   &
    !$acc   noahmp%biochem%param%TemperatureMinPhotosyn,   &
    !$acc   noahmp%biochem%param%LeafAreaPerMass1side,   &
    !$acc   noahmp%biochem%param%NitrogenConcFoliageMax,   &
    !$acc   noahmp%biochem%param%WoodToRootRatio,   &
    !$acc   noahmp%biochem%param%WoodPoolIndex,   &
    !$acc   noahmp%biochem%param%TurnoverCoeffLeafVeg,   &
    !$acc   noahmp%biochem%param%LeafDeathWaterCoeffVeg,   &
    !$acc   noahmp%biochem%param%LeafDeathTempCoeffVeg,   &
    !$acc   noahmp%biochem%param%MicroRespCoeff,   &
    !$acc   noahmp%biochem%param%RespMaintQ10,   &
    !$acc   noahmp%biochem%param%RespMaintLeaf25C,   &
    !$acc   noahmp%biochem%param%RespMaintStem25C,   &
    !$acc   noahmp%biochem%param%RespMaintRoot25C,   &
    !$acc   noahmp%biochem%param%RespMaintGrain25C,   &
    !$acc   noahmp%biochem%param%GrowthRespFrac,   &
    !$acc   noahmp%biochem%param%TemperaureLeafFreeze,   &
    !$acc   noahmp%biochem%param%LeafAreaPerBiomass,   &
    !$acc   noahmp%biochem%param%TempBaseGrowDegDay,   &
    !$acc   noahmp%biochem%param%TempMaxGrowDegDay,   &
    !$acc   noahmp%biochem%param%GrowDegDayEmerg,   &
    !$acc   noahmp%biochem%param%GrowDegDayInitVeg,   &
    !$acc   noahmp%biochem%param%GrowDegDayPostVeg,   &
    !$acc   noahmp%biochem%param%GrowDegDayInitReprod,   &
    !$acc   noahmp%biochem%param%GrowDegDayMature,   &
    !$acc   noahmp%biochem%param%PhotosynRadFrac,   &
    !$acc   noahmp%biochem%param%TempMinCarbonAssim,   &
    !$acc   noahmp%biochem%param%TempMaxCarbonAssim,   &
    !$acc   noahmp%biochem%param%TempMaxCarbonAssimMax,   &
    !$acc   noahmp%biochem%param%CarbonAssimRefMax    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%biochem%param%LightExtCoeff,   &
    !$acc   noahmp%biochem%param%LightUseEfficiency,   &
    !$acc   noahmp%biochem%param%CarbonAssimReducFac,   &
    !$acc   noahmp%biochem%param%StemAreaIndexMin,   &
    !$acc   noahmp%biochem%param%WoodAllocFac,   &
    !$acc   noahmp%biochem%param%WaterStressCoeff,   &
    !$acc   noahmp%biochem%param%LeafAreaIndexMin,   &
    !$acc   noahmp%biochem%param%TurnoverCoeffRootVeg,   &
    !$acc   noahmp%biochem%param%WoodRespCoeff    &
    !$acc   )

  end subroutine BiochemVarExitDevice
end module BiochemVarInitMod
