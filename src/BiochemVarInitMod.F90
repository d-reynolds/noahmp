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
    associate( NumCropGrowStage => noahmp%config%domain%NumCropGrowStage ,&
               ITS => noahmp%config%domain%ITS, ITE => noahmp%config%domain%ITE ,&
               JTS => noahmp%config%domain%JTS, JTE => noahmp%config%domain%JTE )

    if ( .not. allocated(noahmp%biochem%param%LeafDeathTempCoeffCrop) ) then
       allocate( noahmp%biochem%param%LeafDeathTempCoeffCrop(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%LeafDeathTempCoeffCrop)
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafDeathWaterCoeffCrop) ) then
       allocate( noahmp%biochem%param%LeafDeathWaterCoeffCrop(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%LeafDeathWaterCoeffCrop)
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrLeafToGrain) ) then
       allocate( noahmp%biochem%param%CarbohydrLeafToGrain(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%CarbohydrLeafToGrain)
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrStemToGrain) ) then
       allocate( noahmp%biochem%param%CarbohydrStemToGrain(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%CarbohydrStemToGrain)
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrRootToGrain) ) then
       allocate( noahmp%biochem%param%CarbohydrRootToGrain(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%CarbohydrRootToGrain)
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrFracToLeaf) ) then
       allocate( noahmp%biochem%param%CarbohydrFracToLeaf(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%CarbohydrFracToLeaf)
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrFracToStem) ) then
       allocate( noahmp%biochem%param%CarbohydrFracToStem(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%CarbohydrFracToStem)
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrFracToRoot) ) then
       allocate( noahmp%biochem%param%CarbohydrFracToRoot(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%CarbohydrFracToRoot)
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbohydrFracToGrain) ) then
       allocate( noahmp%biochem%param%CarbohydrFracToGrain(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%CarbohydrFracToGrain)
    endif
    if ( .not. allocated(noahmp%biochem%param%TurnoverCoeffLeafCrop) ) then
       allocate( noahmp%biochem%param%TurnoverCoeffLeafCrop(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TurnoverCoeffLeafCrop)
    endif
    if ( .not. allocated(noahmp%biochem%param%TurnoverCoeffStemCrop) ) then
       allocate( noahmp%biochem%param%TurnoverCoeffStemCrop(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TurnoverCoeffStemCrop)
    endif
    if ( .not. allocated(noahmp%biochem%param%TurnoverCoeffRootCrop) ) then
       allocate( noahmp%biochem%param%TurnoverCoeffRootCrop(ITS:ITE,1:NumCropGrowStage,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TurnoverCoeffRootCrop)
    endif

    ! Allocate 2D biochem state arrays and transfer to GPU
    if ( .not. allocated(noahmp%biochem%state%PlantGrowStage) ) then
       allocate( noahmp%biochem%state%PlantGrowStage(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%PlantGrowStage)
    endif
    if ( .not. allocated(noahmp%biochem%state%IndexPlanting) ) then
       allocate( noahmp%biochem%state%IndexPlanting(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%IndexPlanting)
    endif
    if ( .not. allocated(noahmp%biochem%state%IndexHarvest) ) then
       allocate( noahmp%biochem%state%IndexHarvest(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%IndexHarvest)
    endif
    if ( .not. allocated(noahmp%biochem%state%IndexGrowSeason) ) then
       allocate( noahmp%biochem%state%IndexGrowSeason(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%IndexGrowSeason)
    endif
    if ( .not. allocated(noahmp%biochem%state%NitrogenConcFoliage) ) then
       allocate( noahmp%biochem%state%NitrogenConcFoliage(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%NitrogenConcFoliage)
    endif
    if ( .not. allocated(noahmp%biochem%state%LeafMass) ) then
       allocate( noahmp%biochem%state%LeafMass(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%LeafMass)
    endif
    if ( .not. allocated(noahmp%biochem%state%RootMass) ) then
       allocate( noahmp%biochem%state%RootMass(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%RootMass)
    endif
    if ( .not. allocated(noahmp%biochem%state%StemMass) ) then
       allocate( noahmp%biochem%state%StemMass(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%StemMass)
    endif
    if ( .not. allocated(noahmp%biochem%state%WoodMass) ) then
       allocate( noahmp%biochem%state%WoodMass(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%WoodMass)
    endif
    if ( .not. allocated(noahmp%biochem%state%GrainMass) ) then
       allocate( noahmp%biochem%state%GrainMass(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%GrainMass)
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonMassDeepSoil) ) then
       allocate( noahmp%biochem%state%CarbonMassDeepSoil(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%CarbonMassDeepSoil)
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonMassShallowSoil) ) then
       allocate( noahmp%biochem%state%CarbonMassShallowSoil(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%CarbonMassShallowSoil)
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonMassSoilTot) ) then
       allocate( noahmp%biochem%state%CarbonMassSoilTot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%CarbonMassSoilTot)
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonMassLiveTot) ) then
       allocate( noahmp%biochem%state%CarbonMassLiveTot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%CarbonMassLiveTot)
    endif
    if ( .not. allocated(noahmp%biochem%state%LeafAreaPerMass) ) then
       allocate( noahmp%biochem%state%LeafAreaPerMass(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%LeafAreaPerMass)
    endif
    if ( .not. allocated(noahmp%biochem%state%StemAreaPerMass) ) then
       allocate( noahmp%biochem%state%StemAreaPerMass(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%StemAreaPerMass)
    endif
    if ( .not. allocated(noahmp%biochem%state%LeafMassMin) ) then
       allocate( noahmp%biochem%state%LeafMassMin(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%LeafMassMin)
    endif
    if ( .not. allocated(noahmp%biochem%state%StemMassMin) ) then
       allocate( noahmp%biochem%state%StemMassMin(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%StemMassMin)
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonFracToLeaf) ) then
       allocate( noahmp%biochem%state%CarbonFracToLeaf(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%CarbonFracToLeaf)
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonFracToRoot) ) then
       allocate( noahmp%biochem%state%CarbonFracToRoot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%CarbonFracToRoot)
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonFracToWood) ) then
       allocate( noahmp%biochem%state%CarbonFracToWood(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%CarbonFracToWood)
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonFracToStem) ) then
       allocate( noahmp%biochem%state%CarbonFracToStem(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%CarbonFracToStem)
    endif
    if ( .not. allocated(noahmp%biochem%state%WoodCarbonFrac) ) then
       allocate( noahmp%biochem%state%WoodCarbonFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%WoodCarbonFrac)
    endif
    if ( .not. allocated(noahmp%biochem%state%CarbonFracToWoodRoot) ) then
       allocate( noahmp%biochem%state%CarbonFracToWoodRoot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%CarbonFracToWoodRoot)
    endif
    if ( .not. allocated(noahmp%biochem%state%MicroRespFactorSoilWater) ) then
       allocate( noahmp%biochem%state%MicroRespFactorSoilWater(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%MicroRespFactorSoilWater)
    endif
    if ( .not. allocated(noahmp%biochem%state%MicroRespFactorSoilTemp) ) then
       allocate( noahmp%biochem%state%MicroRespFactorSoilTemp(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%MicroRespFactorSoilTemp)
    endif
    if ( .not. allocated(noahmp%biochem%state%RespFacNitrogenFoliage) ) then
       allocate( noahmp%biochem%state%RespFacNitrogenFoliage(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%RespFacNitrogenFoliage)
    endif
    if ( .not. allocated(noahmp%biochem%state%RespFacTemperature) ) then
       allocate( noahmp%biochem%state%RespFacTemperature(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%RespFacTemperature)
    endif
    if ( .not. allocated(noahmp%biochem%state%RespReductionFac) ) then
       allocate( noahmp%biochem%state%RespReductionFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%RespReductionFac)
    endif
    if ( .not. allocated(noahmp%biochem%state%GrowDegreeDay) ) then
       allocate( noahmp%biochem%state%GrowDegreeDay(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%state%GrowDegreeDay)
    endif

    ! Allocate 2D biochem flux arrays and transfer to GPU
    if ( .not. allocated(noahmp%biochem%flux%PhotosynTotal) ) then
       allocate( noahmp%biochem%flux%PhotosynTotal(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%PhotosynTotal)
    endif
    if ( .not. allocated(noahmp%biochem%flux%PhotosynLeafSunlit) ) then
       allocate( noahmp%biochem%flux%PhotosynLeafSunlit(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%PhotosynLeafSunlit)
    endif
    if ( .not. allocated(noahmp%biochem%flux%PhotosynLeafShade) ) then
       allocate( noahmp%biochem%flux%PhotosynLeafShade(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%PhotosynLeafShade)
    endif
    if ( .not. allocated(noahmp%biochem%flux%PhotosynCrop) ) then
       allocate( noahmp%biochem%flux%PhotosynCrop(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%PhotosynCrop)
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrossPriProduction) ) then
       allocate( noahmp%biochem%flux%GrossPriProduction(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%GrossPriProduction)
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetEcoExchange) ) then
       allocate( noahmp%biochem%flux%NetEcoExchange(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%NetEcoExchange)
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionTot) ) then
       allocate( noahmp%biochem%flux%NetPriProductionTot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%NetPriProductionTot)
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionLeaf) ) then
       allocate( noahmp%biochem%flux%NetPriProductionLeaf(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%NetPriProductionLeaf)
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionRoot) ) then
       allocate( noahmp%biochem%flux%NetPriProductionRoot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%NetPriProductionRoot)
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionWood) ) then
       allocate( noahmp%biochem%flux%NetPriProductionWood(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%NetPriProductionWood)
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionStem) ) then
       allocate( noahmp%biochem%flux%NetPriProductionStem(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%NetPriProductionStem)
    endif
    if ( .not. allocated(noahmp%biochem%flux%NetPriProductionGrain) ) then
       allocate( noahmp%biochem%flux%NetPriProductionGrain(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%NetPriProductionGrain)
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationPlantTot) ) then
       allocate( noahmp%biochem%flux%RespirationPlantTot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%RespirationPlantTot)
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationSoilOrg) ) then
       allocate( noahmp%biochem%flux%RespirationSoilOrg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%RespirationSoilOrg)
    endif
    if ( .not. allocated(noahmp%biochem%flux%CarbonToAtmos) ) then
       allocate( noahmp%biochem%flux%CarbonToAtmos(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%CarbonToAtmos)
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrowthRespLeaf) ) then
       allocate( noahmp%biochem%flux%GrowthRespLeaf(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%GrowthRespLeaf)
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrowthRespRoot) ) then
       allocate( noahmp%biochem%flux%GrowthRespRoot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%GrowthRespRoot)
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrowthRespWood) ) then
       allocate( noahmp%biochem%flux%GrowthRespWood(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%GrowthRespWood)
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrowthRespStem) ) then
       allocate( noahmp%biochem%flux%GrowthRespStem(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%GrowthRespStem)
    endif
    if ( .not. allocated(noahmp%biochem%flux%GrowthRespGrain) ) then
       allocate( noahmp%biochem%flux%GrowthRespGrain(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%GrowthRespGrain)
    endif
    if ( .not. allocated(noahmp%biochem%flux%LeafMassMaxChg) ) then
       allocate( noahmp%biochem%flux%LeafMassMaxChg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%LeafMassMaxChg)
    endif
    if ( .not. allocated(noahmp%biochem%flux%StemMassMaxChg) ) then
       allocate( noahmp%biochem%flux%StemMassMaxChg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%StemMassMaxChg)
    endif
    if ( .not. allocated(noahmp%biochem%flux%CarbonDecayToStable) ) then
       allocate( noahmp%biochem%flux%CarbonDecayToStable(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%CarbonDecayToStable)
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationLeaf) ) then
       allocate( noahmp%biochem%flux%RespirationLeaf(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%RespirationLeaf)
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationStem) ) then
       allocate( noahmp%biochem%flux%RespirationStem(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%RespirationStem)
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationWood) ) then
       allocate( noahmp%biochem%flux%RespirationWood(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%RespirationWood)
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationLeafMaint) ) then
       allocate( noahmp%biochem%flux%RespirationLeafMaint(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%RespirationLeafMaint)
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationRoot) ) then
       allocate( noahmp%biochem%flux%RespirationRoot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%RespirationRoot)
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationSoil) ) then
       allocate( noahmp%biochem%flux%RespirationSoil(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%RespirationSoil)
    endif
    if ( .not. allocated(noahmp%biochem%flux%RespirationGrain) ) then
       allocate( noahmp%biochem%flux%RespirationGrain(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%RespirationGrain)
    endif
    if ( .not. allocated(noahmp%biochem%flux%ConvRootToGrain) ) then
       allocate( noahmp%biochem%flux%ConvRootToGrain(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%ConvRootToGrain)
    endif
    if ( .not. allocated(noahmp%biochem%flux%ConvStemToGrain) ) then
       allocate( noahmp%biochem%flux%ConvStemToGrain(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%ConvStemToGrain)
    endif
    if ( .not. allocated(noahmp%biochem%flux%ConvLeafToGrain) ) then
       allocate( noahmp%biochem%flux%ConvLeafToGrain(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%ConvLeafToGrain)
    endif
    if ( .not. allocated(noahmp%biochem%flux%TurnoverLeaf) ) then
       allocate( noahmp%biochem%flux%TurnoverLeaf(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%TurnoverLeaf)
    endif
    if ( .not. allocated(noahmp%biochem%flux%TurnoverStem) ) then
       allocate( noahmp%biochem%flux%TurnoverStem(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%TurnoverStem)
    endif
    if ( .not. allocated(noahmp%biochem%flux%TurnoverWood) ) then
       allocate( noahmp%biochem%flux%TurnoverWood(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%TurnoverWood)
    endif
    if ( .not. allocated(noahmp%biochem%flux%TurnoverRoot) ) then
       allocate( noahmp%biochem%flux%TurnoverRoot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%TurnoverRoot)
    endif
    if ( .not. allocated(noahmp%biochem%flux%TurnoverGrain) ) then
       allocate( noahmp%biochem%flux%TurnoverGrain(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%TurnoverGrain)
    endif
    if ( .not. allocated(noahmp%biochem%flux%DeathLeaf) ) then
       allocate( noahmp%biochem%flux%DeathLeaf(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%DeathLeaf)
    endif
    if ( .not. allocated(noahmp%biochem%flux%DeathStem) ) then
       allocate( noahmp%biochem%flux%DeathStem(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%DeathStem)
    endif
    if ( .not. allocated(noahmp%biochem%flux%CarbonAssim) ) then
       allocate( noahmp%biochem%flux%CarbonAssim(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%CarbonAssim)
    endif
    if ( .not. allocated(noahmp%biochem%flux%CarbohydrAssim) ) then
       allocate( noahmp%biochem%flux%CarbohydrAssim(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%flux%CarbohydrAssim)
    endif

    ! Allocate 2D biochem param arrays and transfer to GPU
    if ( .not. allocated(noahmp%biochem%param%DatePlanting) ) then
       allocate( noahmp%biochem%param%DatePlanting(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%DatePlanting)
    endif
    if ( .not. allocated(noahmp%biochem%param%DateHarvest) ) then
       allocate( noahmp%biochem%param%DateHarvest(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%DateHarvest)
    endif
    if ( .not. allocated(noahmp%biochem%param%QuantumEfficiency25C) ) then
       allocate( noahmp%biochem%param%QuantumEfficiency25C(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%QuantumEfficiency25C)
    endif
    if ( .not. allocated(noahmp%biochem%param%CarboxylRateMax25C) ) then
       allocate( noahmp%biochem%param%CarboxylRateMax25C(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%CarboxylRateMax25C)
    endif
    if ( .not. allocated(noahmp%biochem%param%CarboxylRateMaxQ10) ) then
       allocate( noahmp%biochem%param%CarboxylRateMaxQ10(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%CarboxylRateMaxQ10)
    endif
    if ( .not. allocated(noahmp%biochem%param%PhotosynPathC3) ) then
       allocate( noahmp%biochem%param%PhotosynPathC3(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%PhotosynPathC3)
    endif
    if ( .not. allocated(noahmp%biochem%param%SlopeConductToPhotosyn) ) then
       allocate( noahmp%biochem%param%SlopeConductToPhotosyn(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%SlopeConductToPhotosyn)
    endif
    if ( .not. allocated(noahmp%biochem%param%TemperatureMinPhotosyn) ) then
       allocate( noahmp%biochem%param%TemperatureMinPhotosyn(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TemperatureMinPhotosyn)
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafAreaPerMass1side) ) then
       allocate( noahmp%biochem%param%LeafAreaPerMass1side(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%LeafAreaPerMass1side)
    endif
    if ( .not. allocated(noahmp%biochem%param%NitrogenConcFoliageMax) ) then
       allocate( noahmp%biochem%param%NitrogenConcFoliageMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%NitrogenConcFoliageMax)
    endif
    if ( .not. allocated(noahmp%biochem%param%WoodToRootRatio) ) then
       allocate( noahmp%biochem%param%WoodToRootRatio(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%WoodToRootRatio)
    endif
    if ( .not. allocated(noahmp%biochem%param%WoodPoolIndex) ) then
       allocate( noahmp%biochem%param%WoodPoolIndex(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%WoodPoolIndex)
    endif
    if ( .not. allocated(noahmp%biochem%param%TurnoverCoeffLeafVeg) ) then
       allocate( noahmp%biochem%param%TurnoverCoeffLeafVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TurnoverCoeffLeafVeg)
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafDeathWaterCoeffVeg) ) then
       allocate( noahmp%biochem%param%LeafDeathWaterCoeffVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%LeafDeathWaterCoeffVeg)
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafDeathTempCoeffVeg) ) then
       allocate( noahmp%biochem%param%LeafDeathTempCoeffVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%LeafDeathTempCoeffVeg)
    endif
    if ( .not. allocated(noahmp%biochem%param%MicroRespCoeff) ) then
       allocate( noahmp%biochem%param%MicroRespCoeff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%MicroRespCoeff)
    endif
    if ( .not. allocated(noahmp%biochem%param%RespMaintQ10) ) then
       allocate( noahmp%biochem%param%RespMaintQ10(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%RespMaintQ10)
    endif
    if ( .not. allocated(noahmp%biochem%param%RespMaintLeaf25C) ) then
       allocate( noahmp%biochem%param%RespMaintLeaf25C(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%RespMaintLeaf25C)
    endif
    if ( .not. allocated(noahmp%biochem%param%RespMaintStem25C) ) then
       allocate( noahmp%biochem%param%RespMaintStem25C(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%RespMaintStem25C)
    endif
    if ( .not. allocated(noahmp%biochem%param%RespMaintRoot25C) ) then
       allocate( noahmp%biochem%param%RespMaintRoot25C(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%RespMaintRoot25C)
    endif
    if ( .not. allocated(noahmp%biochem%param%RespMaintGrain25C) ) then
       allocate( noahmp%biochem%param%RespMaintGrain25C(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%RespMaintGrain25C)
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowthRespFrac) ) then
       allocate( noahmp%biochem%param%GrowthRespFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%GrowthRespFrac)
    endif
    if ( .not. allocated(noahmp%biochem%param%TemperaureLeafFreeze) ) then
       allocate( noahmp%biochem%param%TemperaureLeafFreeze(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TemperaureLeafFreeze)
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafAreaPerBiomass) ) then
       allocate( noahmp%biochem%param%LeafAreaPerBiomass(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%LeafAreaPerBiomass)
    endif
    if ( .not. allocated(noahmp%biochem%param%TempBaseGrowDegDay) ) then
       allocate( noahmp%biochem%param%TempBaseGrowDegDay(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TempBaseGrowDegDay)
    endif
    if ( .not. allocated(noahmp%biochem%param%TempMaxGrowDegDay) ) then
       allocate( noahmp%biochem%param%TempMaxGrowDegDay(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TempMaxGrowDegDay)
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowDegDayEmerg) ) then
       allocate( noahmp%biochem%param%GrowDegDayEmerg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%GrowDegDayEmerg)
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowDegDayInitVeg) ) then
       allocate( noahmp%biochem%param%GrowDegDayInitVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%GrowDegDayInitVeg)
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowDegDayPostVeg) ) then
       allocate( noahmp%biochem%param%GrowDegDayPostVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%GrowDegDayPostVeg)
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowDegDayInitReprod) ) then
       allocate( noahmp%biochem%param%GrowDegDayInitReprod(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%GrowDegDayInitReprod)
    endif
    if ( .not. allocated(noahmp%biochem%param%GrowDegDayMature) ) then
       allocate( noahmp%biochem%param%GrowDegDayMature(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%GrowDegDayMature)
    endif
    if ( .not. allocated(noahmp%biochem%param%PhotosynRadFrac) ) then
       allocate( noahmp%biochem%param%PhotosynRadFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%PhotosynRadFrac)
    endif
    if ( .not. allocated(noahmp%biochem%param%TempMinCarbonAssim) ) then
       allocate( noahmp%biochem%param%TempMinCarbonAssim(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TempMinCarbonAssim)
    endif
    if ( .not. allocated(noahmp%biochem%param%TempMaxCarbonAssim) ) then
       allocate( noahmp%biochem%param%TempMaxCarbonAssim(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TempMaxCarbonAssim)
    endif
    if ( .not. allocated(noahmp%biochem%param%TempMaxCarbonAssimMax) ) then
       allocate( noahmp%biochem%param%TempMaxCarbonAssimMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TempMaxCarbonAssimMax)
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbonAssimRefMax) ) then
       allocate( noahmp%biochem%param%CarbonAssimRefMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%CarbonAssimRefMax)
    endif
    if ( .not. allocated(noahmp%biochem%param%LightExtCoeff) ) then
       allocate( noahmp%biochem%param%LightExtCoeff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%LightExtCoeff)
    endif
    if ( .not. allocated(noahmp%biochem%param%LightUseEfficiency) ) then
       allocate( noahmp%biochem%param%LightUseEfficiency(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%LightUseEfficiency)
    endif
    if ( .not. allocated(noahmp%biochem%param%CarbonAssimReducFac) ) then
       allocate( noahmp%biochem%param%CarbonAssimReducFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%CarbonAssimReducFac)
    endif
    if ( .not. allocated(noahmp%biochem%param%StemAreaIndexMin) ) then
       allocate( noahmp%biochem%param%StemAreaIndexMin(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%StemAreaIndexMin)
    endif
    if ( .not. allocated(noahmp%biochem%param%WoodAllocFac) ) then
       allocate( noahmp%biochem%param%WoodAllocFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%WoodAllocFac)
    endif
    if ( .not. allocated(noahmp%biochem%param%WaterStressCoeff) ) then
       allocate( noahmp%biochem%param%WaterStressCoeff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%WaterStressCoeff)
    endif
    if ( .not. allocated(noahmp%biochem%param%LeafAreaIndexMin) ) then
       allocate( noahmp%biochem%param%LeafAreaIndexMin(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%LeafAreaIndexMin)
    endif
    if ( .not. allocated(noahmp%biochem%param%TurnoverCoeffRootVeg) ) then
       allocate( noahmp%biochem%param%TurnoverCoeffRootVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%TurnoverCoeffRootVeg)
    endif
    if ( .not. allocated(noahmp%biochem%param%WoodRespCoeff) ) then
       allocate( noahmp%biochem%param%WoodRespCoeff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%biochem%param%WoodRespCoeff)
    endif

    end associate

    ! Initialize 2D state, flux, and parameter arrays (including 3D crop arrays)
    !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                          &
                  PlantGrowStage           => noahmp%biochem%state%PlantGrowStage(I,J)     ,&
                  IndexPlanting            => noahmp%biochem%state%IndexPlanting(I,J)      ,&
                  IndexHarvest             => noahmp%biochem%state%IndexHarvest(I,J)       ,&
                  IndexGrowSeason          => noahmp%biochem%state%IndexGrowSeason(I,J)    ,&
                  NitrogenConcFoliage      => noahmp%biochem%state%NitrogenConcFoliage(I,J),&
                  LeafMass                 => noahmp%biochem%state%LeafMass(I,J)           ,&
                  RootMass                 => noahmp%biochem%state%RootMass(I,J)           ,&
                  StemMass                 => noahmp%biochem%state%StemMass(I,J)           ,&
                  WoodMass                 => noahmp%biochem%state%WoodMass(I,J)           ,&
                  CarbonMassDeepSoil       => noahmp%biochem%state%CarbonMassDeepSoil(I,J) ,&
                  CarbonMassShallowSoil    => noahmp%biochem%state%CarbonMassShallowSoil(I,J),&
                  CarbonMassSoilTot        => noahmp%biochem%state%CarbonMassSoilTot(I,J)  ,&
                  CarbonMassLiveTot        => noahmp%biochem%state%CarbonMassLiveTot(I,J)  ,&
                  LeafAreaPerMass          => noahmp%biochem%state%LeafAreaPerMass(I,J)    ,&
                  StemAreaPerMass          => noahmp%biochem%state%StemAreaPerMass(I,J)    ,&
                  LeafMassMin              => noahmp%biochem%state%LeafMassMin(I,J)        ,&
                  StemMassMin              => noahmp%biochem%state%StemMassMin(I,J)        ,&
                  CarbonFracToLeaf         => noahmp%biochem%state%CarbonFracToLeaf(I,J)   ,&
                  CarbonFracToRoot         => noahmp%biochem%state%CarbonFracToRoot(I,J)   ,&
                  CarbonFracToWood         => noahmp%biochem%state%CarbonFracToWood(I,J)   ,&
                  CarbonFracToStem         => noahmp%biochem%state%CarbonFracToStem(I,J)   ,&
                  WoodCarbonFrac           => noahmp%biochem%state%WoodCarbonFrac(I,J)     ,&
                  CarbonFracToWoodRoot     => noahmp%biochem%state%CarbonFracToWoodRoot(I,J),&
                  MicroRespFactorSoilWater => noahmp%biochem%state%MicroRespFactorSoilWater(I,J),&
                  MicroRespFactorSoilTemp  => noahmp%biochem%state%MicroRespFactorSoilTemp(I,J),&
                  RespFacNitrogenFoliage   => noahmp%biochem%state%RespFacNitrogenFoliage(I,J),&
                  RespFacTemperature       => noahmp%biochem%state%RespFacTemperature(I,J),&
                  RespReductionFac         => noahmp%biochem%state%RespReductionFac(I,J)  ,&
                  GrainMass                => noahmp%biochem%state%GrainMass(I,J)         ,&
                  GrowDegreeDay            => noahmp%biochem%state%GrowDegreeDay(I,J)     ,&
                  PhotosynLeafSunlit       => noahmp%biochem%flux%PhotosynLeafSunlit(I,J) ,&
                  PhotosynLeafShade        => noahmp%biochem%flux%PhotosynLeafShade(I,J)  ,&
                  PhotosynCrop             => noahmp%biochem%flux%PhotosynCrop(I,J)       ,&
                  PhotosynTotal            => noahmp%biochem%flux%PhotosynTotal(I,J)      ,&
                  GrossPriProduction       => noahmp%biochem%flux%GrossPriProduction(I,J) ,&
                  NetPriProductionTot      => noahmp%biochem%flux%NetPriProductionTot(I,J),&
                  NetEcoExchange           => noahmp%biochem%flux%NetEcoExchange(I,J)     ,&
                  RespirationPlantTot      => noahmp%biochem%flux%RespirationPlantTot(I,J),&
                  RespirationSoilOrg       => noahmp%biochem%flux%RespirationSoilOrg(I,J) ,&
                  CarbonToAtmos            => noahmp%biochem%flux%CarbonToAtmos(I,J)      ,&
                  NetPriProductionLeaf     => noahmp%biochem%flux%NetPriProductionLeaf(I,J),&
                  NetPriProductionRoot     => noahmp%biochem%flux%NetPriProductionRoot(I,J),&
                  NetPriProductionWood     => noahmp%biochem%flux%NetPriProductionWood(I,J),&
                  NetPriProductionStem     => noahmp%biochem%flux%NetPriProductionStem(I,J),&
                  GrowthRespLeaf           => noahmp%biochem%flux%GrowthRespLeaf(I,J)     ,&
                  GrowthRespRoot           => noahmp%biochem%flux%GrowthRespRoot(I,J)     ,&
                  GrowthRespWood           => noahmp%biochem%flux%GrowthRespWood(I,J)     ,&
                  GrowthRespStem           => noahmp%biochem%flux%GrowthRespStem(I,J)     ,&
                  LeafMassMaxChg           => noahmp%biochem%flux%LeafMassMaxChg(I,J)     ,&
                  StemMassMaxChg           => noahmp%biochem%flux%StemMassMaxChg(I,J)     ,&
                  CarbonDecayToStable      => noahmp%biochem%flux%CarbonDecayToStable(I,J),&
                  RespirationLeaf          => noahmp%biochem%flux%RespirationLeaf(I,J)    ,&
                  RespirationStem          => noahmp%biochem%flux%RespirationStem(I,J)    ,&
                  GrowthRespGrain          => noahmp%biochem%flux%GrowthRespGrain(I,J)    ,&
                  NetPriProductionGrain    => noahmp%biochem%flux%NetPriProductionGrain(I,J),&
                  ConvRootToGrain          => noahmp%biochem%flux%ConvRootToGrain(I,J)    ,&
                  ConvStemToGrain          => noahmp%biochem%flux%ConvStemToGrain(I,J)    ,&
                  RespirationWood          => noahmp%biochem%flux%RespirationWood(I,J)    ,&
                  RespirationLeafMaint     => noahmp%biochem%flux%RespirationLeafMaint(I,J),&
                  RespirationRoot          => noahmp%biochem%flux%RespirationRoot(I,J)    ,&
                  DeathLeaf                => noahmp%biochem%flux%DeathLeaf(I,J)          ,&
                  DeathStem                => noahmp%biochem%flux%DeathStem(I,J)          ,&
                  CarbonAssim              => noahmp%biochem%flux%CarbonAssim(I,J)        ,&
                  TurnoverLeaf             => noahmp%biochem%flux%TurnoverLeaf(I,J)       ,&
                  TurnoverStem             => noahmp%biochem%flux%TurnoverStem(I,J)       ,&
                  TurnoverWood             => noahmp%biochem%flux%TurnoverWood(I,J)       ,&
                  RespirationSoil          => noahmp%biochem%flux%RespirationSoil(I,J)    ,&
                  TurnoverRoot             => noahmp%biochem%flux%TurnoverRoot(I,J)       ,&
                  CarbohydrAssim           => noahmp%biochem%flux%CarbohydrAssim(I,J)     ,&
                  TurnoverGrain            => noahmp%biochem%flux%TurnoverGrain(I,J)      ,&
                  ConvLeafToGrain          => noahmp%biochem%flux%ConvLeafToGrain(I,J)    ,&
                  RespirationGrain         => noahmp%biochem%flux%RespirationGrain(I,J)   ,&
                  DatePlanting             => noahmp%biochem%param%DatePlanting(I,J)      ,&
                  DateHarvest              => noahmp%biochem%param%DateHarvest(I,J)       ,&
                  QuantumEfficiency25C     => noahmp%biochem%param%QuantumEfficiency25C(I,J),&
                  CarboxylRateMax25C       => noahmp%biochem%param%CarboxylRateMax25C(I,J),&
                  CarboxylRateMaxQ10       => noahmp%biochem%param%CarboxylRateMaxQ10(I,J),&
                  PhotosynPathC3           => noahmp%biochem%param%PhotosynPathC3(I,J)    ,&
                  SlopeConductToPhotosyn   => noahmp%biochem%param%SlopeConductToPhotosyn(I,J),&
                  TemperatureMinPhotosyn   => noahmp%biochem%param%TemperatureMinPhotosyn(I,J),&
                  LeafAreaPerMass1side     => noahmp%biochem%param%LeafAreaPerMass1side(I,J),&
                  NitrogenConcFoliageMax   => noahmp%biochem%param%NitrogenConcFoliageMax(I,J),&
                  WoodToRootRatio          => noahmp%biochem%param%WoodToRootRatio(I,J)   ,&
                  WoodPoolIndex            => noahmp%biochem%param%WoodPoolIndex(I,J)     ,&
                  TurnoverCoeffLeafVeg     => noahmp%biochem%param%TurnoverCoeffLeafVeg(I,J),&
                  LeafDeathWaterCoeffVeg   => noahmp%biochem%param%LeafDeathWaterCoeffVeg(I,J),&
                  LeafDeathTempCoeffVeg    => noahmp%biochem%param%LeafDeathTempCoeffVeg(I,J),&
                  MicroRespCoeff           => noahmp%biochem%param%MicroRespCoeff(I,J)    ,&
                  RespMaintQ10             => noahmp%biochem%param%RespMaintQ10(I,J)      ,&
                  RespMaintLeaf25C         => noahmp%biochem%param%RespMaintLeaf25C(I,J)  ,&
                  RespMaintStem25C         => noahmp%biochem%param%RespMaintStem25C(I,J)  ,&
                  RespMaintRoot25C         => noahmp%biochem%param%RespMaintRoot25C(I,J)  ,&
                  RespMaintGrain25C        => noahmp%biochem%param%RespMaintGrain25C(I,J) ,&
                  GrowthRespFrac           => noahmp%biochem%param%GrowthRespFrac(I,J)    ,&
                  TemperaureLeafFreeze     => noahmp%biochem%param%TemperaureLeafFreeze(I,J),&
                  LeafAreaPerBiomass       => noahmp%biochem%param%LeafAreaPerBiomass(I,J),&
                  TempBaseGrowDegDay       => noahmp%biochem%param%TempBaseGrowDegDay(I,J),&
                  TempMaxGrowDegDay        => noahmp%biochem%param%TempMaxGrowDegDay(I,J) ,&
                  GrowDegDayEmerg          => noahmp%biochem%param%GrowDegDayEmerg(I,J)   ,&
                  GrowDegDayInitVeg        => noahmp%biochem%param%GrowDegDayInitVeg(I,J) ,&
                  GrowDegDayPostVeg        => noahmp%biochem%param%GrowDegDayPostVeg(I,J) ,&
                  GrowDegDayInitReprod     => noahmp%biochem%param%GrowDegDayInitReprod(I,J),&
                  GrowDegDayMature         => noahmp%biochem%param%GrowDegDayMature(I,J)  ,&
                  PhotosynRadFrac          => noahmp%biochem%param%PhotosynRadFrac(I,J)   ,&
                  TempMinCarbonAssim       => noahmp%biochem%param%TempMinCarbonAssim(I,J),&
                  TempMaxCarbonAssim       => noahmp%biochem%param%TempMaxCarbonAssim(I,J),&
                  TempMaxCarbonAssimMax    => noahmp%biochem%param%TempMaxCarbonAssimMax(I,J),&
                  CarbonAssimRefMax        => noahmp%biochem%param%CarbonAssimRefMax(I,J) ,&
                  LightExtCoeff            => noahmp%biochem%param%LightExtCoeff(I,J)     ,&
                  LightUseEfficiency       => noahmp%biochem%param%LightUseEfficiency(I,J),&
                  CarbonAssimReducFac      => noahmp%biochem%param%CarbonAssimReducFac(I,J),&
                  StemAreaIndexMin         => noahmp%biochem%param%StemAreaIndexMin(I,J)  ,&
                  WoodAllocFac             => noahmp%biochem%param%WoodAllocFac(I,J)      ,&
                  WaterStressCoeff         => noahmp%biochem%param%WaterStressCoeff(I,J)  ,&
                  LeafAreaIndexMin         => noahmp%biochem%param%LeafAreaIndexMin(I,J)  ,&
                  TurnoverCoeffRootVeg     => noahmp%biochem%param%TurnoverCoeffRootVeg(I,J),&
                  WoodRespCoeff            => noahmp%biochem%param%WoodRespCoeff(I,J)      &
                 )

        ! biochem state variables
        PlantGrowStage             = undefined_int
        IndexPlanting              = undefined_int
        IndexHarvest               = undefined_int
        IndexGrowSeason            = undefined_real
        NitrogenConcFoliage        = undefined_real
        LeafMass                   = undefined_real
        RootMass                   = undefined_real
        StemMass                   = undefined_real
        WoodMass                   = undefined_real
        CarbonMassDeepSoil         = undefined_real
        CarbonMassShallowSoil      = undefined_real
        CarbonMassSoilTot          = undefined_real
        CarbonMassLiveTot          = undefined_real
        LeafAreaPerMass            = undefined_real
        StemAreaPerMass            = undefined_real
        LeafMassMin                = undefined_real
        StemMassMin                = undefined_real
        CarbonFracToLeaf           = undefined_real
        CarbonFracToRoot           = undefined_real
        CarbonFracToWood           = undefined_real
        CarbonFracToStem           = undefined_real
        WoodCarbonFrac             = undefined_real
        CarbonFracToWoodRoot       = undefined_real
        MicroRespFactorSoilWater   = undefined_real
        MicroRespFactorSoilTemp    = undefined_real
        RespFacNitrogenFoliage     = undefined_real
        RespFacTemperature         = undefined_real
        RespReductionFac           = undefined_real
        GrainMass                  = undefined_real
        GrowDegreeDay              = undefined_real

        ! biochem flux variables
        PhotosynLeafSunlit         = undefined_real
        PhotosynLeafShade          = undefined_real
        PhotosynCrop               = undefined_real
        PhotosynTotal              = undefined_real
        GrossPriProduction         = undefined_real
        NetPriProductionTot        = undefined_real
        NetEcoExchange             = undefined_real
        RespirationPlantTot        = undefined_real
        RespirationSoilOrg         = undefined_real
        CarbonToAtmos              = undefined_real
        NetPriProductionLeaf       = undefined_real
        NetPriProductionRoot       = undefined_real
        NetPriProductionWood       = undefined_real
        NetPriProductionStem       = undefined_real
        GrowthRespLeaf             = undefined_real
        GrowthRespRoot             = undefined_real
        GrowthRespWood             = undefined_real
        GrowthRespStem             = undefined_real
        LeafMassMaxChg             = undefined_real
        StemMassMaxChg             = undefined_real
        CarbonDecayToStable        = undefined_real
        RespirationLeaf            = undefined_real
        RespirationStem            = undefined_real
        GrowthRespGrain            = undefined_real
        NetPriProductionGrain      = undefined_real
        ConvRootToGrain            = undefined_real
        ConvStemToGrain            = undefined_real
        RespirationWood            = undefined_real
        RespirationLeafMaint       = undefined_real
        RespirationRoot            = undefined_real
        DeathLeaf                  = undefined_real
        DeathStem                  = undefined_real
        CarbonAssim                = undefined_real
        TurnoverLeaf               = undefined_real
        TurnoverStem               = undefined_real
        TurnoverWood               = undefined_real
        RespirationSoil            = undefined_real
        TurnoverRoot               = undefined_real
        CarbohydrAssim             = undefined_real
        TurnoverGrain              = undefined_real
        ConvLeafToGrain            = undefined_real
        RespirationGrain           = undefined_real

        ! biochem parameter variables
        DatePlanting               = undefined_int
        DateHarvest                = undefined_int
        QuantumEfficiency25C       = undefined_real
        CarboxylRateMax25C         = undefined_real
        CarboxylRateMaxQ10         = undefined_real
        PhotosynPathC3             = undefined_real
        SlopeConductToPhotosyn     = undefined_real
        TemperatureMinPhotosyn     = undefined_real
        LeafAreaPerMass1side       = undefined_real
        NitrogenConcFoliageMax     = undefined_real
        WoodToRootRatio            = undefined_real
        WoodPoolIndex              = undefined_real
        TurnoverCoeffLeafVeg       = undefined_real
        LeafDeathWaterCoeffVeg     = undefined_real
        LeafDeathTempCoeffVeg      = undefined_real
        MicroRespCoeff             = undefined_real
        RespMaintQ10               = undefined_real
        RespMaintLeaf25C           = undefined_real
        RespMaintStem25C           = undefined_real
        RespMaintRoot25C           = undefined_real
        RespMaintGrain25C          = undefined_real
        GrowthRespFrac             = undefined_real
        TemperaureLeafFreeze       = undefined_real
        LeafAreaPerBiomass         = undefined_real
        TempBaseGrowDegDay         = undefined_real
        TempMaxGrowDegDay          = undefined_real
        GrowDegDayEmerg            = undefined_real
        GrowDegDayInitVeg          = undefined_real
        GrowDegDayPostVeg          = undefined_real
        GrowDegDayInitReprod       = undefined_real
        GrowDegDayMature           = undefined_real
        PhotosynRadFrac            = undefined_real
        TempMinCarbonAssim         = undefined_real
        TempMaxCarbonAssim         = undefined_real
        TempMaxCarbonAssimMax      = undefined_real
        CarbonAssimRefMax          = undefined_real
        LightExtCoeff              = undefined_real
        LightUseEfficiency         = undefined_real
        CarbonAssimReducFac        = undefined_real
        StemAreaIndexMin           = undefined_real
        WoodAllocFac               = undefined_real
        WaterStressCoeff           = undefined_real
        LeafAreaIndexMin           = undefined_real
        TurnoverCoeffRootVeg       = undefined_real
        WoodRespCoeff              = undefined_real

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

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine BiochemVarInitDefault

end module BiochemVarInitMod
