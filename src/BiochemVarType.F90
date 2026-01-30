module BiochemVarType

!!! Define 2D Noah-MP Biochemistry (carbon,nitrogen,etc) variables
!!! Biochemistry variable initialization is done in BiochemVarInitMod.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

  use Machine

  implicit none
  save
  private

!=== define "flux" sub-type of biochem (biochem%flux%variable)
  type :: flux_type

    ! All flux fields are now 2D arrays (I,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PhotosynTotal              ! total leaf photosynthesis [umol co2/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PhotosynLeafSunlit         ! sunlit leaf photosynthesis [umol co2/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PhotosynLeafShade          ! shaded leaf photosynthesis [umol co2/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PhotosynCrop               ! crop photosynthesis rate [umol co2/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrossPriProduction         ! gross primary production [g/m2/s C]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: NetEcoExchange             ! net ecosystem exchange [g/m2/s CO2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: NetPriProductionTot        ! total net primary production [g/m2/s C]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: NetPriProductionLeaf       ! leaf net primary production [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: NetPriProductionRoot       ! root net primary production [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: NetPriProductionWood       ! wood net primary production [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: NetPriProductionStem       ! stem net primary production [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: NetPriProductionGrain      ! grain net primary production [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespirationPlantTot        ! total plant respiration (leaf,stem,root,wood,grain) [g/m2/s C]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespirationSoilOrg         ! soil heterotrophic (organic) respiration [g/m2/s C]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonToAtmos              ! carbon flux to atmosphere [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowthRespLeaf             ! growth respiration rate for leaf [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowthRespRoot             ! growth respiration rate for root [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowthRespWood             ! growth respiration rate for wood [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowthRespStem             ! growth respiration rate for stem [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowthRespGrain            ! growth respiration rate for grain [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafMassMaxChg             ! maximum leaf mass available to change [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: StemMassMaxChg             ! maximum stem mass available to change [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonDecayToStable        ! decay rate of fast carbon to slow carbon [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespirationLeaf            ! leaf respiration [umol CO2/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespirationStem            ! stem respiration [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespirationWood            ! wood respiration rate [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespirationLeafMaint       ! leaf maintenance respiration rate [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespirationRoot            ! fine root respiration rate [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespirationSoil            ! soil respiration rate [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespirationGrain           ! grain respiration rate [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ConvRootToGrain            ! root to grain conversion [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ConvStemToGrain            ! stem to grain conversion [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ConvLeafToGrain            ! leaf to grain conversion [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TurnoverLeaf               ! leaf turnover rate [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TurnoverStem               ! stem turnover rate [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TurnoverWood               ! wood turnover rate [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TurnoverRoot               ! root turnover rate [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TurnoverGrain              ! grain turnover rate [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DeathLeaf                  ! death rate of leaf mass [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DeathStem                  ! death rate of stem mass [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonAssim                ! carbon assimilated rate [g/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbohydrAssim             ! carbohydrate assimilated rate [g/m2/s]

  end type flux_type


!=== define "state" sub-type of biochem (biochem%state%variable)
  type :: state_type

    ! All state fields are now 2D arrays (I,J)
    integer, allocatable, dimension(:,:) :: PlantGrowStage             ! plant growing stage
    integer, allocatable, dimension(:,:) :: IndexPlanting              ! Planting index (0=off, 1=on)
    integer, allocatable, dimension(:,:) :: IndexHarvest               ! Harvest index (0=on,1=off)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IndexGrowSeason            ! growing season index (0=off, 1=on)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: NitrogenConcFoliage        ! foliage nitrogen concentration [%]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafMass                   ! leaf mass [g/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RootMass                   ! mass of fine roots [g/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: StemMass                   ! stem mass [g/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WoodMass                   ! mass of wood (include woody roots) [g/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrainMass                  ! mass of grain [g/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonMassDeepSoil         ! stable carbon in deep soil [g/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonMassShallowSoil      ! short-lived carbon in shallow soil [g/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonMassSoilTot          ! total soil carbon mass [g/m2 C]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonMassLiveTot          ! total living carbon mass ([g/m2 C]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafAreaPerMass            ! leaf area per unit mass [m2/g]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: StemAreaPerMass            ! stem area per unit mass (m2/g)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafMassMin                ! minimum leaf mass [g/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: StemMassMin                ! minimum stem mass [g/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonFracToLeaf           ! fraction of carbon flux allocated to leaves
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonFracToRoot           ! fraction of carbon flux allocated to roots
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonFracToWood           ! fraction of carbon flux allocated to wood
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonFracToStem           ! fraction of carbon flux allocated to stem
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WoodCarbonFrac             ! wood carbon fraction in (root + wood) carbon
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonFracToWoodRoot       ! fraction of carbon to root and wood
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MicroRespFactorSoilWater   ! soil water factor for microbial respiration
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MicroRespFactorSoilTemp    ! soil temperature factor for microbial respiration
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespFacNitrogenFoliage     ! foliage nitrogen adjustemt factor to respiration (<= 1)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespFacTemperature         ! temperature factor for respiration
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespReductionFac           ! respiration reduction factor (<= 1)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowDegreeDay              ! growing degree days

  end type state_type


!=== define "parameter" sub-type of biochem (biochem%param%variable)
  type :: parameter_type

    ! All parameter fields are now 2D arrays (I,J) for spatial flexibility
    integer, allocatable, dimension(:,:) :: DatePlanting               ! planting date
    integer, allocatable, dimension(:,:) :: DateHarvest                ! harvest date
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: QuantumEfficiency25C       ! quantum efficiency at 25c [umol CO2/umol photon]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarboxylRateMax25C         ! maximum rate of carboxylation at 25c [umol CO2/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarboxylRateMaxQ10         ! change in maximum rate of carboxylation for every 10-deg C temperature change
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PhotosynPathC3             ! C3 photosynthetic pathway indicator: 0.0 = c4, 1.0 = c3
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SlopeConductToPhotosyn     ! slope of conductance-to-photosynthesis relationship
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureMinPhotosyn     ! minimum temperature for photosynthesis [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafAreaPerMass1side       ! single-side leaf area per mass [m2/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: NitrogenConcFoliageMax     ! foliage nitrogen concentration when f(n)=1 (%)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WoodToRootRatio            ! wood to root ratio
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WoodPoolIndex              ! wood pool index (0~1) depending on woody or not
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TurnoverCoeffLeafVeg       ! leaf turnover coefficient [1/s] for generic vegetation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafDeathWaterCoeffVeg     ! coeficient for leaf water stress death [1/s] for generic vegetation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafDeathTempCoeffVeg      ! coeficient for leaf temperature stress death [1/s] for generic vegetation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MicroRespCoeff             ! microbial respiration coefficient [umol co2 /kg c/ s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespMaintQ10               ! change in maintenance respiration for every 10-deg C temperature change
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespMaintLeaf25C           ! leaf maintenance respiration at 25C [umol CO2/m2  /s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespMaintStem25C           ! stem maintenance respiration at 25C [umol CO2/kg bio/s], bio: C or CH2O
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespMaintRoot25C           ! root maintenance respiration at 25C [umol CO2/kg bio/s], bio: C or CH2O
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RespMaintGrain25C          ! grain maintenance respiration at 25C [umol CO2/kg bio/s], bio: C or CH2O
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowthRespFrac             ! fraction of growth respiration
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperaureLeafFreeze       ! characteristic temperature for leaf freezing [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafAreaPerBiomass         ! leaf area per living leaf biomass [m2/g]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TempBaseGrowDegDay         ! Base temperature for growing degree day (GDD) accumulation [C]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TempMaxGrowDegDay          ! Maximum temperature for growing degree day (GDD) accumulation [C]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowDegDayEmerg            ! growing degree day (GDD) from seeding to emergence
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowDegDayInitVeg          ! growing degree day (GDD) from seeding to initial vegetative
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowDegDayPostVeg          ! growing degree day (GDD) from seeding to post vegetative
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowDegDayInitReprod       ! growing degree day (GDD) from seeding to intial reproductive
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GrowDegDayMature           ! growing degree day (GDD) from seeding to pysical maturity
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PhotosynRadFrac            ! Fraction of incoming solar radiation to photosynthetically active radiation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TempMinCarbonAssim         ! Minimum temperature for CO2 assimulation [C]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TempMaxCarbonAssim         ! CO2 assimulation linearly increasing until reaching this temperature [C]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TempMaxCarbonAssimMax      ! CO2 assmilation rate remain at CarbonAssimRefMax until reaching this temperature [C]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonAssimRefMax          ! reference maximum CO2 assimilation rate [g co2/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LightExtCoeff              ! light extinction coefficient
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LightUseEfficiency         ! initial light use efficiency
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CarbonAssimReducFac        ! CO2 assimilation reduction factor(0-1) (caused by non-modeling part,e.g.pest,weeds)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: StemAreaIndexMin           ! minimum stem area index [m2/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WoodAllocFac               ! present wood allocation factor
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterStressCoeff           ! water stress coeficient
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafAreaIndexMin           ! minimum leaf area index [m2/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TurnoverCoeffRootVeg       ! root turnover coefficient [1/s] for generic vegetation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WoodRespCoeff              ! wood respiration coeficient [1/s]

    ! Crop-specific growth stage parameters (I,stage,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: LeafDeathTempCoeffCrop      ! coeficient for leaf temperature stress death [1/s] for crop
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: LeafDeathWaterCoeffCrop     ! coeficient for leaf water stress death [1/s] for crop
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: CarbohydrLeafToGrain        ! fraction of carbohydrate flux transallocate from leaf to grain
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: CarbohydrStemToGrain        ! fraction of carbohydrate flux transallocate from stem to grain
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: CarbohydrRootToGrain        ! fraction of carbohydrate flux transallocate from root to grain
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: CarbohydrFracToLeaf         ! fraction of carbohydrate flux to leaf for crop
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: CarbohydrFracToStem         ! fraction of carbohydrate flux to stem for crop
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: CarbohydrFracToRoot         ! fraction of carbohydrate flux to root for crop
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: CarbohydrFracToGrain        ! fraction of carbohydrate flux to grain for crop
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: TurnoverCoeffLeafCrop       ! leaf turnover coefficient [1/s] for crop
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: TurnoverCoeffStemCrop       ! stem turnover coefficient [1/s] for crop
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: TurnoverCoeffRootCrop       ! root tunrover coefficient [1/s] for crop

  end type parameter_type


!=== define biochem type that includes 3 subtypes (flux,state,parameter)
  type, public :: biochem_type

    type(flux_type)      :: flux
    type(state_type)     :: state
    type(parameter_type) :: param

  end type biochem_type

end module BiochemVarType
