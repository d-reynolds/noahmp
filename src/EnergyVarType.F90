module EnergyVarType

!!! Define 2D Noah-MP Energy variables
!!! Energy variable initialization is done in EnergyVarInitMod.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

  use Machine

  implicit none
  save
  private

!=== define "flux" sub-type of energy (energy%flux%variable)
  type :: flux_type

    ! All scalar fields now 2D arrays (I,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatLatentCanopy            ! canopy latent heat flux [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatLatentTransp            ! latent heat flux from transpiration [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatLatentGrd               ! total ground latent heat [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatLatentIrriEvap          ! latent heating due to sprinkler irrigation evaporation [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatPrecipAdvCanopy         ! precipitation advected heat - canopy net [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatPrecipAdvVegGrd         ! precipitation advected heat - vegetated ground net [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatPrecipAdvBareGrd        ! precipitation advected heat - bare ground net [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatPrecipAdvSfc            ! precipitation advected heat - total [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatSensibleCanopy          ! canopy sensible heat flux [W/m2]     (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatLatentCanEvap           ! canopy evaporation heat flux [W/m2]  (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatSensibleVegGrd          ! vegetated ground sensible heat flux [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatSensibleSfc             ! total sensible heat [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatLatentVegGrd            ! vegetated ground latent heat flux [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatLatentCanTransp         ! canopy transpiration latent heat flux [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatGroundVegGrd            ! vegetated ground heat flux [W/m2] (+ to soil/snow)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatSensibleBareGrd         ! bare ground sensible heat flux [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatLatentBareGrd           ! bare ground latent heat flux [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatGroundBareGrd           ! bare ground heat flux [W/m2] (+ to soil/snow)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatGroundTot               ! total ground heat flux [W/m2] (+ to soil/snow)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatGroundTotMean           ! total ground heat flux [W/m2] averaged over soil timestep
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatFromSoilBot             ! energy influx from soil bottom [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatCanStorageChg           ! canopy heat storage change [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatGroundTotAcc            ! accumulated total ground heat flux per soil timestep [W/m2 * dt_soil/dt_main] (+ to soil/snow)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadPhotoActAbsSunlit        ! absorbed photosyn. active radiation for sunlit leaves [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadPhotoActAbsShade         ! absorbed photosyn. active radiation  for shaded leaves [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadSwAbsVeg                 ! solar radiation absorbed by vegetation [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadSwAbsGrd                 ! solar radiation absorbed by ground [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadSwAbsSfc                 ! total absorbed solar radiation [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadSwReflSfc                ! total reflected solar radiation [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadSwReflVeg                ! reflected solar radiation by vegetation [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadSwReflGrd                ! reflected solar radiation by ground [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadLwNetCanopy              ! canopy net longwave radiation [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadLwNetSfc                 ! total net longwave radiation [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadPhotoActAbsCan           ! total photosyn. active energy [W/m2] absorbed by canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadLwEmitSfc                ! emitted outgoing longwave radiation [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadLwNetVegGrd              ! vegetated ground net longwave radiation [W/m2] (+ to atm)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadLwNetBareGrd             ! bare ground net longwave rad [W/m2] (+ to atm)

    ! Spectral/band arrays now 3D: (I,K,J) where K is band dimension
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwAbsVegDir        ! solar flux absorbed by veg per unit direct flux (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwAbsVegDif        ! solar flux absorbed by veg per unit diffuse flux (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwDirTranGrdDir    ! transmitted direct flux below veg per unit direct flux (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwDirTranGrdDif    ! transmitted direct flux below veg per unit diffuse flux (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwDifTranGrdDir    ! transmitted diffuse flux below veg per unit direct flux (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwDifTranGrdDif    ! transmitted diffuse flux below veg per unit diffuse flux (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwReflVegDir       ! solar flux reflected by veg layer per unit direct flux (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwReflVegDif       ! solar flux reflected by veg layer per unit diffuse flux (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwReflGrdDir       ! solar flux reflected by ground per unit direct flux (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwReflGrdDif       ! solar flux reflected by ground per unit diffuse flux (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwDownDir          ! incoming direct solar radiation [W/m2] (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwDownDif          ! incoming diffuse solar radiation [W/m2] (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwPenetrateGrd     ! light penetrating through soil/snow water [W/m2] (I,band,J)
    ! Snow layer spectral absorption: (I, K_layer, K_band, J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:,:) :: FracRadSwAbsSnowDir   ! direct spectral solar flux absorbed by snow layer [frc] (I,layer,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:,:) :: FracRadSwAbsSnowDif   ! diffuse spectral solar flux absorbed by snow layer [frc] (I,layer,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:)   :: RadSwAbsSnowSoilLayer ! solar flux absorbed by snow/soil for each layer [W/m2] (I,layer,J)


  end type flux_type


!=== define "state" sub-type of energy (energy%state%variable)
  type :: state_type

    ! All scalar fields now 2D arrays (I,J)
    logical, allocatable, dimension(:,:) :: FlagFrozenCanopy            ! frozen canopy flag used to define latent heat pathway
    logical, allocatable, dimension(:,:) :: FlagFrozenGround            ! frozen ground flag used to define latent heat pathway
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafAreaIndEff              ! effective leaf area index, after burying by snow
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: StemAreaIndEff              ! effective stem area index, after burying by snow
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafAreaIndex               ! leaf area index
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: StemAreaIndex               ! stem area index
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VegAreaIndEff               ! one-sided leaf+stem area index [m2/m2], after burying by snow
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VegFrac                     ! greeness vegetation fraction
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureGrd              ! ground temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureCanopy           ! vegetation/canopy temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureSfc              ! surface temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureRootZone         ! root-zone averaged temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PressureVaporRefHeight      ! vapor pressure air [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowAgeFac                  ! snow age factor
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowAgeNondim               ! non-dimensional snow age
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: AlbedoSnowPrev              ! snow albedo at last time step
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VegAreaProjDir              ! projected leaf+stem area in solar direction
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GapBtwCanopy                ! between canopy gap fraction for beam
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GapInCanopy                 ! within canopy gap fraction for beam
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GapCanopyDif                ! gap fraction for diffue light
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: GapCanopyDir                ! total gap fraction for beam (<=1-shafac)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopySunlitFrac            ! sunlit fraction of canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyShadeFrac             ! shaded fraction of canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafAreaIndSunlit           ! sunlit leaf area
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafAreaIndShade            ! shaded leaf area
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VapPresSatCanopy            ! canopy saturation vapor pressure at veg temperature [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VapPresSatGrdVeg            ! below-canopy saturation vapor pressure at ground temperature [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VapPresSatGrdBare           ! bare ground saturation vapor pressure at ground temperature [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VapPresSatCanTempD          ! canopy saturation vapor pressure derivative with temperature at veg temp. [Pa/K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VapPresSatGrdVegTempD       ! below-canopy saturation vapor pressure derivative with temperature at ground temp. [Pa/K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VapPresSatGrdBareTempD      ! bare ground saturation vapor pressure derivative with temperature at ground temp. [Pa/K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PressureVaporCanAir         ! canopy air vapor pressure [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PressureAtmosCO2            ! atmospheric co2 partial pressure [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PressureAtmosO2             ! atmospheric o2 partial pressure [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceStomataSunlit     ! sunlit leaf stomatal resistance [s/m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceStomataShade      ! shaded leaf stomatal resistance [s/m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DensityAirRefHeight         ! density air [kg/m3] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureCanopyAir        ! canopy air temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ZeroPlaneDispSfc            ! surface zero plane displacement [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ZeroPlaneDispGrd            ! ground zero plane displacement [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RoughLenMomGrd              ! roughness length, momentum, ground [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RoughLenMomSfc              ! roughness length, momentum, surface [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RoughLenShCanopy            ! roughness length, sensible heat, canopy [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RoughLenShVegGrd            ! roughness length, sensible heat, ground, below canopy [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RoughLenShBareGrd           ! roughness length, sensible heat, bare ground [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyHeight                ! canopy height [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WindSpdCanopyTop            ! wind speed at top of canopy [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FrictionVelVeg              ! friction velocity [m/s], vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FrictionVelBare             ! friction velocity [m/s], bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WindExtCoeffCanopy          ! canopy wind extinction coefficient
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabParaUndCan            ! M-O stability parameter ground, below canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabParaAbvCan            ! M-O stability parameter (z/L), above ZeroPlaneDisp, vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabParaBare              ! M-O stability parameter (z/L), above ZeroPlaneDisp, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabParaVeg2m             ! M-O stability parameter (2/L), 2m, vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabParaBare2m            ! M-O stability parameter (2/L), 2m, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoLengthUndCan              ! M-O length [m], ground, below canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoLengthAbvCan              ! M-O length [m], above ZeroPlaneDisp, vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoLengthBare                ! M-O length [m], above ZeroPlaneDisp, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabCorrShUndCan          ! M-O stability correction ground, below canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabCorrMomAbvCan         ! M-O momentum stability correction, above ZeroPlaneDisp, vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabCorrShAbvCan          ! M-O sensible heat stability correction, above ZeroPlaneDisp, vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabCorrMomVeg2m          ! M-O momentum stability correction, 2m, vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabCorrShVeg2m           ! M-O sensible heat stability correction, 2m, vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabCorrShBare            ! M-O sensible heat stability correction, above ZeroPlaneDisp, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabCorrMomBare           ! M-O momentum stability correction, above ZeroPlaneDisp, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabCorrMomBare2m         ! M-O momentum stability correction, 2m, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MoStabCorrShBare2m          ! M-O sensible heat stability correction, 2m, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffMomSfc             ! exchange coefficient [m/s] for momentum, surface, grid mean
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffMomAbvCan          ! exchange coefficient [m/s] for momentum, above ZeroPlaneDisp, vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffMomBare            ! exchange coefficient [m/s] for momentum, above ZeroPlaneDisp, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffShSfc              ! exchange coefficient [m/s] for sensible heat, surface, grid mean
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffShAbvCan           ! exchange coefficient [m/s] for sensible heat, above ZeroPlaneDisp, vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffShBare             ! exchange coefficient [m/s] for sensible heat, above ZeroPlaneDisp, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffSh2mVegMo          ! exchange coefficient [m/s] for sensible heat, 2m, vegetated (M-O)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffSh2mBareMo         ! exchange coefficient [m/s] for sensible heat, 2m, bare ground (M-O)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffSh2mVeg            ! exchange coefficient [m/s] for sensible heat, 2m, vegetated (diagnostic)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffLhAbvCan           ! exchange coefficient [m/s] for latent heat, canopy air to ref height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffLhTransp           ! exchange coefficient [m/s] for transpiration, leaf to canopy air
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffLhEvap             ! exchange coefficient [m/s] for leaf evaporation, leaf to canopy air
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffLhUndCan           ! exchange coefficient [m/s] for latent heat, ground to canopy air
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceMomUndCan         ! aerodynamic resistance [s/m] for momentum, ground, below canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceShUndCan          ! aerodynamic resistance [s/m] for sensible heat, ground, below canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceLhUndCan          ! aerodynamic resistance [s/m] for water vapor, ground, below canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceMomAbvCan         ! aerodynamic resistance [s/m] for momentum, above canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceShAbvCan          ! aerodynamic resistance [s/m] for sensible heat, above canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceLhAbvCan          ! aerodynamic resistance [s/m] for water vapor, above canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceMomBareGrd        ! aerodynamic resistance [s/m] for momentum, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceShBareGrd         ! aerodynamic resistance [s/m] for sensible heat, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceLhBareGrd         ! aerodynamic resistance [s/m] for water vapor, bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceLeafBoundary      ! bulk leaf boundary layer resistance [s/m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperaturePotRefHeight     ! potential temp at reference height [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WindSpdRefHeight            ! wind speed [m/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FrictionVelVertVeg          ! friction velocity in vertical direction [m/s], vegetated (only for Chen97)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FrictionVelVertBare         ! friction velocity in vertical direction [m/s], bare ground (only for Chen97)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EmissivityVeg               ! vegetation emissivity
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EmissivityGrd               ! ground emissivity
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceGrdEvap           ! ground surface resistance [s/m] to evaporation/sublimation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PsychConstCanopy            ! psychrometric constant [Pa/K], canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LatHeatVapCanopy            ! latent heat of vaporization/subli [J/kg], canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PsychConstGrd               ! psychrometric constant [Pa/K], ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LatHeatVapGrd               ! latent heat of vaporization/subli [J/kg], ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RelHumidityGrd              ! raltive humidity in surface soil/snow air space (-)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SpecHumiditySfc             ! specific humidity at surface (bare or vegetated or urban)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SpecHumiditySfcMean         ! specific humidity at surface grid mean
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SpecHumidity2mVeg           ! specific humidity at 2m vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SpecHumidity2mBare          ! specific humidity at 2m bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SpecHumidity2m              ! specific humidity at 2m grid mean
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureGrdVeg           ! vegetated ground (below-canopy) temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureGrdBare          ! bare ground temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WindStressEwVeg             ! wind stress [N/m2]: east-west above canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WindStressNsVeg             ! wind stress [N/m2]: north-south above canopy
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WindStressEwBare            ! wind stress [N/m2]: east-west bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WindStressNsBare            ! wind stress [N/m2]: north-south bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WindStressEwSfc             ! wind stress [N/m2]: east-west grid mean
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WindStressNsSfc             ! wind stress [N/m2]: north-south grid mean
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureAir2mVeg         ! 2 m height air temperature [K], vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureAir2mBare        ! 2 m height air temperature [K], bare ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureAir2m            ! 2 m height air temperature [K], grid mean
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffShLeaf             ! leaf sensible heat exchange coefficient [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffShUndCan           ! under canopy sensible heat exchange coefficient [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ExchCoeffSh2mBare           ! bare ground 2-m sensible heat exchange coefficient [m/s] (diagnostic)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RefHeightAboveGrd           ! reference height [m] above ground
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyFracSnowBury          ! fraction of canopy buried by snow
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DepthSoilTempBotToSno       ! depth of soil temperature lower boundary condition from snow surface [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RoughLenMomSfcToAtm         ! roughness length, momentum, surface, sent to coupled atmos model
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureRadSfc           ! radiative temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EmissivitySfc               ! surface emissivity
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: AlbedoSfc                   ! total surface albedo
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EnergyBalanceError          ! error in surface energy balance [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadSwBalanceError           ! error in shortwave radiation balance [W/m2]

    ! Layer arrays now 3D: (I,K,J) where K is layer dimension
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: TemperatureSoilSnow   ! snow and soil layer temperature [K] (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: HeatCapacVolSnow      ! snow layer volumetric specific heat capacity [J/m3/K] (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ThermConductSnow      ! snow layer thermal conductivity [W/m/K] (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: HeatCapacVolSoil      ! soil layer volumetric specific heat capacity [J/m3/K] (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ThermConductSoil      ! soil layer thermal conductivity [W/m/K] (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: HeatCapacGlaIce       ! glacier ice layer volumetric specific heat [J/m3/K] (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ThermConductGlaIce    ! glacier ice thermal conductivity [W/m/K] (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ThermConductSoilSnow  ! thermal conductivity for all soil and snow layers [W/m/K] (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: HeatCapacSoilSnow     ! heat capacity for all snow and soil layers [J/m3/K] (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: PhaseChgFacSoilSnow   ! energy factor for soil and snow phase change (I,layer,J)
    ! Band/spectral arrays now 3D: (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoSnowDir         ! snow albedo for direct(1=vis, 2=nir) (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoSnowDif         ! snow albedo for diffuse(1=vis, 2=nir) (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoSoilDir         ! soil albedo (direct) (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoSoilDif         ! soil albedo (diffuse) (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoGrdDir          ! ground albedo (direct beam: vis, nir) (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoGrdDif          ! ground albedo (diffuse: vis, nir) (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ReflectanceVeg        ! leaf/stem reflectance weighted by LeafAreaIndex and StemAreaIndex (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: TransmittanceVeg      ! leaf/stem transmittance weighted by LeafAreaIndex and StemAreaIndex (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoSfcDir          ! surface albedo (direct) (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoSfcDif          ! surface albedo (diffuse) (I,band,J)

  end type state_type


!=== define "parameter" sub-type of energy (energy%param%variable)
  type :: parameter_type

    ! Parameters are typically spatially uniform, but keeping as 2D for flexibility
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TreeCrownRadius             ! tree crown radius [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeightCanopyTop             ! height of canopy top [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeightCanopyBot             ! height of canopy bottom [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RoughLenMomVeg              ! momentum roughness length [m] vegetated
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TreeDensity                 ! tree density [no. of trunks per m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyOrientIndex           ! leaf/stem orientation index
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: UpscatterCoeffSnowDir       ! Upscattering parameters for snow for direct radiation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: UpscatterCoeffSnowDif       ! Upscattering parameters for snow for diffuse radiation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SoilHeatCapacity            ! volumetric soil heat capacity [j/m3/K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowAgeFacBats              ! snow aging parameter for BATS snow albedo
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowGrowVapFacBats          ! vapor diffusion snow growth factor for BATS snow albedo
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowSootFacBats             ! dirt and soot effect factor for BATS snow albedo
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowGrowFrzFacBats          ! extra snow growth factor near freezing for BATS snow albedo
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SolarZenithAdjBats          ! zenith angle snow albedo adjustment
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FreshSnoAlbVisBats          ! new snow visible albedo for BATS
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FreshSnoAlbNirBats          ! new snow NIR albedo for BATS
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnoAgeFacDifVisBats         ! age factor for diffuse visible snow albedo for BATS
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnoAgeFacDifNirBats         ! age factor for diffuse NIR snow albedo for BATS
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SzaFacDirVisBats            ! cosz factor for direct visible snow albedo for BATS
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SzaFacDirNirBats            ! cosz factor for direct NIR snow albedo for BATS
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowAlbRefClass             ! reference snow albedo in CLASS scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowAgeFacClass             ! snow aging e-folding time [s] in CLASS albedo scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SnowAlbFreshClass           ! fresh snow albedo in CLASS albedo scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ConductanceLeafMin          ! minimum leaf conductance [umol/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: Co2MmConst25C               ! co2 michaelis-menten constant at 25c [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: O2MmConst25C                ! o2 michaelis-menten constant at 25c [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: Co2MmConstQ10               ! change in co2 Michaelis-Menten constant for every 10-deg C temperature change
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: O2MmConstQ10                ! change in o2 michaelis-menten constant for every 10-deg C temperature change
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadiationStressFac          ! Parameter used in radiation stress function in Jarvis scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceStomataMin        ! Minimum stomatal resistance [s/m] in Jarvis scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceStomataMax        ! Maximal stomatal resistance [s/m] in Jarvis scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: AirTempOptimTransp          ! Optimum transpiration air temperature [K] in Jarvis scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VaporPresDeficitFac         ! Parameter used in vapor pressure deficit function in Jarvis scheme
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LeafDimLength               ! characteristic leaf dimension [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ZilitinkevichCoeff          ! Zilitinkevich coefficient for heat exchange coefficient calculation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EmissivitySnow              ! snow emissivity
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: CanopyWindExtFac            ! empirical canopy wind extinction parameter
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RoughLenMomSnow             ! snow surface roughness length [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RoughLenMomSoil             ! Bare-soil roughness length [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RoughLenMomLake             ! lake surface roughness length [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: EmissivityIceSfc            ! ice surface emissivity
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceSoilExp           ! exponent in the shape parameter for soil resistance option 1
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ResistanceSnowSfc           ! surface resistance for snow [s/m]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VegFracGreen                ! green vegetation fraction
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: VegFracAnnMax               ! annual maximum vegetation fraction
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: HeatCapacCanFac             ! canopy biomass heat capacity parameter [m]

    ! Parameter arrays: 1D arrays (e.g., monthly, spectral) become 3D: (I,K,J) where K is parameter dimension
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: LeafAreaIndexMon      ! monthly leaf area index, one-sided (I,month,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: StemAreaIndexMon      ! monthly stem area index, one-sided (I,month,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilQuartzFrac        ! soil quartz content (I,layer,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoSoilSat         ! saturated soil albedos: 1=vis, 2=nir (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoSoilDry         ! dry soil albedos: 1=vis, 2=nir (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoLakeFrz         ! albedo frozen lakes: 1=vis, 2=nir (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ScatterCoeffSnow      ! Scattering coefficient for snow (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ReflectanceLeaf       ! leaf reflectance: 1=vis, 2=nir (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ReflectanceStem       ! stem reflectance: 1=vis, 2=nir (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: TransmittanceLeaf     ! leaf transmittance: 1=vis, 2=nir (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: TransmittanceStem     ! stem transmittance: 1=vis, 2=nir (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: EmissivitySoilLake    ! emissivity soil surface: 1=soil, 2=lake (I,type,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AlbedoLandIce         ! land/glacier ice albedo: 1=vis, 2=nir (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwWgtDir           ! downward solar radiation spectral weights (direct) (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: RadSwWgtDif           ! downward solar radiation spectral weights (diffuse) (I,band,J)
    ! 2D parameter lookup tables become 4D: (I, K1, K2, J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:,:) :: SsAlbSnwRadDir        ! Mie single scatter albedos for direct-beam ice (I,band,radius,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:,:) :: AsyPrmSnwRadDir       ! asymmetry parameter of direct-beam ice (I,band,radius,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:,:) :: ExtCffMassSnwRadDir   ! mass extinction coefficient for direct-beam ice [m2/kg] (I,band,radius,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:,:) :: SsAlbSnwRadDif        ! Mie single scatter albedos for diffuse ice (I,band,radius,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:,:) :: AsyPrmSnwRadDif       ! asymmetry parameter of diffuse ice (I,band,radius,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:,:) :: ExtCffMassSnwRadDif   ! mass extinction coefficient for diffuse ice [m2/kg] (I,band,radius,J)
    ! Aerosol optical properties: 1D spectral arrays become 3D
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SsAlbBCphi            ! Mie single scatter albedos for hydrophillic BC (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AsyPrmBCphi           ! asymmetry parameter for hydrophillic BC (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ExtCffMassBCphi       ! mass extinction coefficient for hydrophillic BC [m2/kg] (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SsAlbBCpho            ! Mie single scatter albedos for hydrophobic BC (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AsyPrmBCpho           ! asymmetry parameter for hydrophobic BC (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ExtCffMassBCpho       ! mass extinction coefficient for hydrophobic BC [m2/kg] (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SsAlbOCphi            ! Mie single scatter albedos for hydrophillic OC (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AsyPrmOCphi           ! asymmetry parameter for hydrophillic OC (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ExtCffMassOCphi       ! mass extinction coefficient for hydrophillic OC [m2/kg] (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SsAlbOCpho            ! Mie single scatter albedos for hydrophobic OC (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AsyPrmOCpho           ! asymmetry parameter for hydrophobic OC (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ExtCffMassOCpho       ! mass extinction coefficient for hydrophobic OC [m2/kg] (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SsAlbDustB1           ! Mie single scatter albedos for dust species 1 (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AsyPrmDustB1          ! asymmetry parameter for dust species 1 (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ExtCffMassDustB1      ! mass extinction coefficient for dust species 1 [m2/kg] (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SsAlbDustB2           ! Mie single scatter albedos for dust species 2 (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AsyPrmDustB2          ! asymmetry parameter for dust species 2 (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ExtCffMassDustB2      ! mass extinction coefficient for dust species 2 [m2/kg] (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SsAlbDustB3           ! Mie single scatter albedos for dust species 3 (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AsyPrmDustB3          ! asymmetry parameter for dust species 3 (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ExtCffMassDustB3      ! mass extinction coefficient for dust species 3 [m2/kg] (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SsAlbDustB4           ! Mie single scatter albedos for dust species 4 (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AsyPrmDustB4          ! asymmetry parameter for dust species 4 (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ExtCffMassDustB4      ! mass extinction coefficient for dust species 4 [m2/kg] (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SsAlbDustB5           ! Mie single scatter albedos for dust species 5 (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: AsyPrmDustB5          ! asymmetry parameter for dust species 5 (I,band,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ExtCffMassDustB5      ! mass extinction coefficient for dust species 5 [m2/kg] (I,band,J)

  end type parameter_type


!=== define energy type that includes 3 subtypes (flux,state,parameter)
  type, public :: energy_type

    type(flux_type)      :: flux
    type(state_type)     :: state
    type(parameter_type) :: param

  end type energy_type

end module EnergyVarType
