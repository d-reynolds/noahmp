module WaterVarType

!!! Define 2D Noah-MP Water variables
!!! Water variable initialization is done in WaterVarInitMod.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

  use Machine

  implicit none
  save
  private

!=== define "flux" sub-type of water (water%flux%variable)
  type :: flux_type

    ! All flux fields are now 2D arrays (I,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RainfallRefHeight          ! liquid rainfall rate [mm/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowfallRefHeight          ! snowfall rate [mm/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PrecipTotRefHeight         ! total precipitation [mm/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PrecipConvTotRefHeight     ! total convective precipitation [mm/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PrecipLargeSclRefHeight    ! large-scale precipitation [mm/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EvapCanopyNet              ! net evaporation of canopy intercepted total water [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: Transpiration              ! transpiration rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EvapCanopyLiq              ! canopy liquid water evaporation rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DewCanopyLiq               ! canopy water dew rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FrostCanopyIce             ! canopy ice frost rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SublimCanopyIce            ! canopy ice sublimation rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MeltCanopyIce              ! canopy ice melting rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FreezeCanopyLiq            ! canopy water freezing rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowfallGround             ! snowfall on the ground (below canopy) [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowDepthIncr              ! snow depth increasing rate [m/s] due to snowfall
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FrostSnowSfcIce            ! snow surface ice frost rate[mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SublimSnowSfcIce           ! snow surface ice sublimation rate[mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RainfallGround             ! ground surface rain rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowBotOutflow             ! total water (snowmelt + rain through pack) out of snowpack bottom [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GlacierExcessFlow          ! glacier excess flow [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrrigationRateFlood        ! flood irrigation water rate [m/timestep]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrrigationRateMicro        ! micro irrigation water rate [m/timestep]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrrigationRateSprinkler    ! sprinkler irrigation water rate [m/timestep]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrriEvapLossSprinkler      ! loss of irrigation water to evaporation,sprinkler [m/timestep]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilSfcInflow              ! water input on soil surface [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RunoffSurface              ! surface runoff [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RunoffSubsurface           ! subsurface runoff [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: InfilRateSfc               ! infiltration rate at surface [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EvapSoilSfcLiq             ! soil surface water evaporation [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DrainSoilBot               ! soil bottom drainage [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TileDrain                  ! tile drainage [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RechargeGw                 ! groundwater recharge rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DischargeGw                ! groundwater discharge rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VaporizeGrd                ! ground vaporize rate total (evap+sublim) [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CondenseVapGrd             ! ground vapor condense rate total (dew+frost) [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DewSoilSfcLiq              ! soil surface water dew rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EvapIrriSprinkler          ! evaporation of irrigation water, sprinkler [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: InterceptCanopyRain        ! interception rate for rain [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DripCanopyRain             ! drip rate for intercepted rain [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ThroughfallRain            ! throughfall for rain [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: InterceptCanopySnow        ! interception (loading) rate for snowfall [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DripCanopySnow             ! drip (unloading) rate for intercepted snow [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ThroughfallSnow            ! throughfall of snowfall [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EvapGroundNet              ! net ground (soil/snow) evaporation [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MeltGroundSnow             ! ground snow melting rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterToAtmosTotal          ! total surface water vapor flux to atmosphere [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EvapSoilSfcLiqAcc          ! accumulated soil surface water evaporation per soil timestep [m/s * dt_soil/dt_main]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilSfcInflowAcc           ! accumulated water input on soil surface per soil timestep [m/s * dt_soil/dt_main]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SfcWaterTotChgAcc          ! accumulated snow,soil,canopy water change per soil timestep [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PrecipTotAcc               ! accumulated precipitation per soil timestep [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EvapCanopyNetAcc           ! accumulated net evaporation of canopy intercepted water per soil timestep [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TranspirationAcc           ! accumulated transpiration per soil timestep [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EvapGroundNetAcc           ! accumulated net ground (soil/snow) evaporation per soil timestep [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GlacierExcessFlowAcc       ! accumulated glacier excessive flow [mm] per soil timestep
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EvapSoilSfcLiqMean         ! mean soil surface water evaporation during soil timestep [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilSfcInflowMean          ! mean water input on soil surface during soil timestep [m/s]

    ! Layer-specific arrays (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: TranspWatLossSoil     ! transpiration water loss from soil layers [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: TranspWatLossSoilAcc  ! accumulated transpiration water loss from soil per soil timestep [m/s * dt_soil/dt_main]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: TranspWatLossSoilMean ! mean transpiration water loss from soil during soil timestep [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: CompactionSnowAging   ! rate of snow compaction due to destructive metamorphism/aging [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: CompactionSnowBurden  ! rate of snow compaction due to overburden [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: CompactionSnowMelt    ! rate of snow compaction due to melt [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: CompactionSnowTot     ! rate of total snow compaction [fraction/timestep]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: OutflowSnowLayer      ! water flow out of each snow layer [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowFreezeRate        ! rate of snow freezing [mm/s]

  end type flux_type


!=== define "state" sub-type of water (water%state%variable)
  type :: state_type

    ! All state fields are now 2D arrays (I,J)
    integer, allocatable, dimension(:,:) :: IrrigationCntSprinkler     ! irrigation event number, Sprinkler
    integer, allocatable, dimension(:,:) :: IrrigationCntMicro         ! irrigation event number, Micro
    integer, allocatable, dimension(:,:) :: IrrigationCntFlood         ! irrigation event number, Flood
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyTotalWater           ! total (liquid+ice) canopy intercepted water [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyWetFrac              ! wetted or snowed fraction of the canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowfallDensity            ! bulk density of snowfall (kg/m3)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyLiqWater             ! intercepted canopy liquid water [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyIce                  ! intercepted canopy ice [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyIceMax               ! canopy capacity for snow interception [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyLiqWaterMax          ! canopy capacity for rain interception [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowDepth                  ! snow depth [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowWaterEquiv             ! snow water equivalent (ice+liquid) [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowWaterEquivPrev         ! snow water equivalent at previous time step (mm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PondSfcThinSnwMelt         ! surface ponding [mm] from snowmelt when snow has no layer
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PondSfcThinSnwComb         ! surface ponding [mm] from liquid in thin snow layer combination
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PondSfcThinSnwTrans        ! surface ponding [mm] from thin snow liquid during transition from multilayer to no layer
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrrigationFracFlood        ! fraction of grid under flood irrigation (0 to 1)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrrigationAmtFlood         ! flood irrigation water amount [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrrigationFracMicro        ! fraction of grid under micro irrigation (0 to 1)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrrigationAmtMicro         ! micro irrigation water amount [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrrigationFracSprinkler    ! fraction of grid under sprinkler irrigation (0 to 1)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrrigationAmtSprinkler     ! sprinkler irrigation water amount [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterTableDepth            ! water table depth [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilIceMax                 ! maximum soil ice content [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilLiqWaterMin            ! minimum soil liquid water content [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilSaturateFrac           ! fractional saturated area for soil moisture
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilImpervFracMax          ! maximum soil imperviousness fraction
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilMoistureToWT           ! soil moisture between bottom of the soil and the water table
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RechargeGwDeepWT           ! groundwater recharge to or from the water table when deep [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RechargeGwShallowWT        ! groundwater recharge to or from shallow water table [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilSaturationExcess       ! saturation excess of the total soil [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterTableHydro            ! water table depth estimated in WRF-Hydro fine grids [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TileDrainFrac              ! tile drainage fraction
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterStorageAquifer        ! water storage in aquifer [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterStorageSoilAqf        ! water storage in aquifer + saturated soil [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterStorageLake           ! water storage in lake (can be negative) [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterStorageWetland        ! water storage in wetland [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterHeadSfc               ! surface water head [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrrigationFracGrid         ! total irrigation fraction from input for a grid
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PrecipAreaFrac             ! fraction of the gridcell that receives precipitation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCoverFrac              ! snow cover fraction
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilTranspFacAcc           ! accumulated soil water transpiration factor (0 to 1)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FrozenPrecipFrac           ! fraction of frozen precip in total precipitation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilWaterRootZone          ! root zone soil water
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilWaterStress            ! soil water stress
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterStorageTotBeg         ! total water storage [mm] at the begining before NoahMP process
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterBalanceError          ! water balance error [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterStorageTotEnd         ! total water storage [mm] at the end of NoahMP process
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowRadiusFresh            ! fresh snow radius [microns]

    ! Layer-specific arrays (I,layer,J)
    integer               , allocatable, dimension(:,:,:) :: IndexPhaseChange      ! phase change index (0-none;1-melt;2-refreeze)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowIce               ! snow layer ice [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowLiqWater          ! snow layer liquid water [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowIceFracPrev       ! ice fraction in snow layers at previous timestep
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowIceFrac           ! ice fraction in snow layers at current timestep
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilIceFrac           ! ice fraction in soil layers at current timestep
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowEffPorosity       ! snow effective porosity [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilLiqWater          ! soil liquid moisture [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilIce               ! soil ice moisture [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilMoisture          ! total soil moisture [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilImpervFrac        ! fraction of imperviousness due to frozen soil
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilWatConductivity   ! soil hydraulic/water conductivity [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilWatDiffusivity    ! soil water diffusivity [m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilEffPorosity       ! soil effective porosity [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilMoistureEqui      ! equilibrium soil water  content [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilTranspFac         ! soil water transpiration factor (0 to 1)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowIceVol            ! partial volume of snow ice [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowLiqWaterVol       ! partial volume of snow liquid water [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilSupercoolWater    ! supercooled water in soil [kg/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilMatPotential      ! soil matric potential [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowRadius            ! snow effective grain radius [microns, m-6]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassBChydropho        ! mass of hydrophobic Black Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassBChydrophi        ! mass of hydrophillic Black Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassOChydropho        ! mass of hydrophobic Organic Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassOChydrophi        ! mass of hydrophillic Organic Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassDust1             ! mass of dust species 1 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassDust2             ! mass of dust species 2 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassDust3             ! mass of dust species 3 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassDust4             ! mass of dust species 4 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassDust5             ! mass of dust species 5 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassConcBChydropho    ! mass concentration of hydrophobic Black Carbon in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassConcBChydrophi    ! mass concentration of hydrophillic Black Carbon in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassConcOChydropho    ! mass concentration of hydrophobic Organic Carbon in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassConcOChydrophi    ! mass concentration of hydrophillic Organic Carbon in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassConcDust1         ! mass concentration of dust species 1 in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassConcDust2         ! mass concentration of dust species 2 in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassConcDust3         ! mass concentration of dust species 3 in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassConcDust4         ! mass concentration of dust species 4 in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassConcDust5         ! mass concentration of dust species 5 in snow [kg/kg]

  end type state_type


!=== define "parameter" sub-type of water (water%param%variable)
  type :: parameter_type

    ! All parameter fields are now 2D arrays (I,J) for spatial flexibility
    integer, allocatable, dimension(:,:) :: DrainSoilLayerInd          ! starting soil layer for drainage
    integer, allocatable, dimension(:,:) :: TileDrainTubeDepth         ! depth [m] of drain tube from the soil surface for simple scheme
    integer, allocatable, dimension(:,:) :: NumSoilLayerRoot           ! number of soil layers with root present
    integer, allocatable, dimension(:,:) :: IrriStopDayBfHarvest       ! number of days before harvest date to stop irrigation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyLiqHoldCap           ! maximum canopy intercepted liquid water per unit veg area index [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCompactBurdenFac       ! overburden snow compaction parameter [m3/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCompactAgingFac1       ! snow desctructive metamorphism compaction parameter1 [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCompactAgingFac2       ! snow desctructive metamorphism compaction parameter2 [1/k]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCompactAgingFac3       ! snow desctructive metamorphism compaction parameter3
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCompactAgingMax        ! upper Limit on destructive metamorphism compaction [kg/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowViscosityCoeff         ! snow viscosity coefficient [kg-s/m2], Anderson1979: 0.52e6~1.38e6
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCompactmAR24           ! snow compaction m parameter for linear sfc temp fitting from AR24
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCompactbAR24           ! snow compaction b parameter for linear sfc temp fitting from AR24
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCompactP1AR24          ! lower constraint for SnowCompactBurdenFac for high pressure bin from AR24
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCompactP2AR24          ! lower constraint for SnowCompactBurdenFac for mid pressure bin from AR24
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCompactP3AR24          ! lower constraint for SnowCompactBurdenFac for low pressure bin from AR24
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: BurdenFacUpAR24            ! upper constraint on SnowCompactBurdenFac from AR24
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCoverM1AR25            ! SCFm1 ground SCF parameter from AR2025
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCoverM2AR25            ! SCFm2 ground SCF parameter from AR2025
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCoverFac1AR25          ! SCfac1 ground SCF parameter from AR2025
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCoverFac2AR25          ! SCfac2 ground SCF parameter from AR2025
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowLiqFracMax             ! maximum liquid water fraction in snow
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowLiqHoldCap             ! liquid water holding capacity for snowpack [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowLiqReleaseFac          ! snowpack water release timescale factor [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrriFloodRateFac           ! flood irrigation application rate factor
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrriMicroRate              ! micro irrigation rate [mm/hr]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilInfilMaxCoeff          ! parameter to calculate maximum soil infiltration rate
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilImpervFracCoeff        ! parameter to calculate frozen soil impermeable fraction
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: InfilFacVic                ! VIC model infiltration parameter
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TensionWatDistrInfl        ! Tension water distribution inflection parameter
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TensionWatDistrShp         ! Tension water distribution shape parameter
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FreeWatDistrShp            ! Free water distribution shape parameter
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: InfilHeteroDynVic          ! DVIC heterogeniety parameter for infiltration
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: InfilCapillaryDynVic       ! DVIC Mean Capillary Drive (m) for infiltration models
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: InfilFacDynVic             ! DVIC model infiltration parameter
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilDrainSlope             ! slope index for soil drainage
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TileDrainCoeffSp           ! drainage coefficient [mm d^-1] for simple scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DrainFacSoilWat            ! drainage factor for soil moisture
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TileDrainCoeff             ! drainage coefficent [m d^-1] for Hooghoudt scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DrainDepthToImperv         ! actual depth of tile drainage to impermeable layer form surface
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LateralWatCondFac          ! multiplication factor to determine lateral hydraulic conductivity
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TileDrainDepth             ! Depth of drain [m] for Hooghoudt scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DrainTubeDist              ! distance between two drain tubes or tiles [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DrainTubeRadius            ! effective radius of drain tubes [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DrainWatDepToImperv        ! depth to impervious layer from drain water level [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RunoffDecayFac             ! runoff decay factor [m^-1]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: BaseflowCoeff              ! baseflow coefficient [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GridTopoIndex              ! gridcell mean topgraphic index (global mean)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilSfcSatFracMax          ! maximum surface soil saturated fraction (global mean)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SpecYieldGw                ! specific yield [-] for Niu et al. 2007 groundwater scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MicroPoreContent           ! microprore content (0.0-1.0), 0.0: close to free drainage
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WaterStorageLakeMax        ! maximum lake water storage [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnoWatEqvMaxGlacier        ! Maximum SWE allowed at glaciers [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilConductivityRef        ! reference Soil Conductivity parameter (used in runoff formulation)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilInfilFacRef            ! reference Soil Infiltration Parameter (used in runoff formulation)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GroundFrzCoeff             ! frozen ground parameter to compute frozen soil impervious fraction
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrriTriggerLaiMin          ! minimum lai to trigger irrigation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilWatDeficitAllow        ! management allowable deficit (0-1)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrriFloodLossFrac          ! factor of flood irrigation loss
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrriSprinklerRate          ! sprinkler irrigation rate [mm/h]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrriFracThreshold          ! irrigation Fraction threshold in a grid
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: IrriStopPrecipThr          ! precipitation threshold [mm/hr] to stop irrigation trigger
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowfallDensityMax         ! maximum fresh snowfall density [kg/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowMassFullCoverOld       ! new snow mass to fully cover old snow [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilMatPotentialWilt       ! soil metric potential for wilting point [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowMeltFac                ! snowmelt m parameter in snow cover fraction calculation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowCoverFac               ! snow cover factor [m] (originally hard-coded 2.5*z0 in SCF formulation)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WetlandCapMax              ! maximum wetland capacity [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowRadiusMin              ! minimum allowed snow effective radius for SNICAR (also cold "fresh snow" value) [microns]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FreshSnowRadiusMax         ! maximum warm fresh snow effective radius [microns]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowRadiusRefrz            ! Effective radius of re-frozen snow [microns]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ScavEffMeltScale           ! Scaling factor modifying scavenging factors for aerosol in meltwater (-)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ScavEffMeltBCphi           ! scavenging factor for hydrophillic BC inclusion in meltwater [frc]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ScavEffMeltBCpho           ! scavenging factor for hydrophobic BC inclusion in meltwater  [frc]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ScavEffMeltOCphi           ! scavenging factor for hydrophillic OC inclusion in meltwater [frc]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ScavEffMeltOCpho           ! scavenging factor for hydrophobic OC inclusion in meltwater  [frc]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ScavEffMeltDust1           ! scavenging factor for dust species 1 inclusion in meltwater  [frc]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ScavEffMeltDust2           ! scavenging factor for dust species 2 inclusion in meltwater  [frc]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ScavEffMeltDust3           ! scavenging factor for dust species 3 inclusion in meltwater  [frc]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ScavEffMeltDust4           ! scavenging factor for dust species 4 inclusion in meltwater  [frc]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ScavEffMeltDust5           ! scavenging factor for dust species 5 inclusion in meltwater  [frc]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowRadiusMax              ! maximum allowed snow effective radius [microns]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowWetAgeC1Brun89         ! constant for liquid water grain growth [m3 s-1], from Brun89
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowWetAgeC2Brun89         ! constant for liquid water grain growth [m3 s-1], from Brun89: corrected for LWC
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowAgeScaleFac            ! arbitrary tuning/scaling factor applied to snow aging rate (-)

    ! Soil layer parameters (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilMoistureSat        ! saturated value of soil moisture [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilMoistureWilt       ! wilting point soil moisture [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilMoistureFieldCap   ! reference soil moisture (field capacity) [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilMoistureDry        ! dry soil moisture threshold [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilWatDiffusivitySat  ! saturated soil hydraulic diffusivity [m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilWatConductivitySat ! saturated soil hydraulic conductivity [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilExpCoeffB          ! soil exponent B paramete
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilMatPotentialSat    ! saturated soil matric potential [m]

    ! Snow aging lookup tables (I,K1,K2,J) - Note: these may remain 3D depending on implementation
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: snowage_tau            ! snowage tau from table [hours]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: snowage_kappa          ! snowage kappa from table [unitless]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: snowage_drdt0          ! snowage dr/dt_0 from table [m2/kg/hr]

  end type parameter_type


!=== define water type that includes 3 subtypes (flux,state,parameter)
  type, public :: water_type

    type(flux_type)      :: flux
    type(state_type)     :: state
    type(parameter_type) :: param

  end type water_type

end module WaterVarType
