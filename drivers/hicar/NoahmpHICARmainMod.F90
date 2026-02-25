module NoahmpHICARmainMod 

! -------------------------------------------------------------
! this is the interface for NoahMP and HICAR variable remapping
! and calling the main NoahMP driver: NoahmpDriverMain(NoahmpIO)
! adapted from original module_sf_noahmpdrv.F file
!
! Coder: Cenlin He (NCAR), December 2025
! Adapted for HICAR code: Dylan Reynolds (EPFL), January 2026
! -------------------------------------------------------------

contains

  subroutine NoahmpHICARmain(NoahmpIO, ITIMESTEP, YR, JULIAN, COSZIN, XLAT, XLONG, & ! IN : Time/Space-related
                   DZ8W,          DT,        DZS,    NSOIL,                      & ! IN : Model configuration 
                   IVGTYP,    ISLTYP,     VEGFRA,   VEGMAX,      TMN,            & ! IN : Vegetation/Soil characteristics
                   XLAND,       XICE,              CROPCAT,                      & ! IN : Vegetation/Soil characteristics
                   PLANTING, HARVEST, SEASON_GDD,                                & ! IN : agriculture input
                   SOILCOMP, SOILCL1,    SOILCL2,  SOILCL3,  SOILCL4,            & ! IN : User options
                   T3D, QV3D,  U_PHY,      V_PHY,   SWDOWN,   SWDDIR,            & ! IN : forcing
                   SWDDIF,       GLW,                                            & ! IN : Forcing
                   P8W3D,PRECIP_IN,           SR,                                & ! IN : Forcing
                   IRFRACT,  SIFRACT,    MIFRACT,  FIFRACT,                      & ! IN : irrigation
                   TSK,          HFX,        QFX,       LH,   GRDFLX,   SMSTAV,  & ! IN/OUT LSM eqv
                   SMSTOT, SFCRUNOFF,   UDRUNOFF,   ALBEDO,    SNOWC,    SMOIS,  & ! IN/OUT LSM eqv
                   SH2O,        TSLB,       SNOW,    SNOWH,   CANWAT,   ACSNOM,  & ! IN/OUT LSM eqv
                   ACSNOW,     EMISS,       QSFC,                                & ! IN/OUT LSM eqv
                   Z0,                                                           & ! IN/OUT LSM eqv
                   IRNUMSI,  IRNUMMI,    IRNUMFI,  IRWATSI,  IRWATMI,  IRWATFI,  & ! IN/OUT irrigation
                   IRELOSS,  IRSIVOL,    IRMIVOL,  IRFIVOL,  IRRSPLH, LLANDUSE,  & ! IN/OUT irrigation
                   ISNOWXY,     TVXY,       TGXY, CANICEXY, CANLIQXY,    EAHXY,  & ! IN/OUT Noah MP only
                   TAHXY,       CMXY,       CHXY,   FWETXY, SNEQVOXY, ALBOLDXY,  & ! IN/OUT Noah MP only
                   QSNOWXY,  QRAINXY,   WSLAKEXY,    ZWTXY,WAXY,WTXY,   TSNOXY,  & ! IN/OUT Noah MP only
                   ZSNSOXY,  SNICEXY,    SNLIQXY, LFMASSXY, RTMASSXY, STMASSXY,  & ! IN/OUT Noah MP only
                   WOODXY,  STBLCPXY,   FASTCPXY,   XLAIXY,   XSAIXY,  TAUSSXY,  & ! IN/OUT Noah MP only
                   SMOISEQ, SMCWTDXY, DEEPRECHXY,   RECHXY,GRAINXY,GDDXY,PGSXY,  & ! IN/OUT Noah MP only
                !    QTDRAIN,   TD_FRACTION,                                       & ! IN/OUT tile drainage
                   T2MVXY,    T2MBXY,     Q2MVXY,   Q2MBXY,                      & ! OUT Noah MP only
                   TRADXY,     NEEXY,      GPPXY,    NPPXY,   FVEGXY,  RUNSFXY,  & ! OUT Noah MP only
                   RUNSBXY,   ECANXY,     EDIRXY,  ETRANXY,    FSAXY,   FIRAXY,  & ! OUT Noah MP only
                   APARXY,     PSNXY,      SAVXY,    SAGXY,  RSSUNXY,  RSSHAXY,  & ! OUT Noah MP only
                   BGAPXY,    WGAPXY,      TGVXY,    TGBXY,    CHVXY,    CHBXY,  & ! OUT Noah MP only
                   SHGXY,      SHCXY,      SHBXY,    EVGXY,    EVBXY,    GHVXY,  & ! OUT Noah MP only
                   GHBXY,      IRGXY,      IRCXY,    IRBXY,     TRXY,    EVCXY,  & ! OUT Noah MP only
                   CHLEAFXY,  CHUCXY,     CHV2XY,   CHB2XY,       RS,            & ! OUT Noah MP only
                !    QINTSXY,  QINTRXY,   QDRIPSXY,                                & ! OUT Noah MP only
                !    QDRIPRXY,QTHROSXY,   QTHRORXY,                                & ! OUT Noah MP only
                !    QSNSUBXY,QSNFROXY,    QSUBCXY,                                & ! OUT Noah MP only
                !    QFROCXY,  QEVACXY,    QDEWCXY,  QFRZCXY, QMELTCXY,            & ! OUT Noah MP only
                !    QSNBOTXY, QMELTXY,  PONDINGXY,  PAHXY,PAHGXY,PAHVXY, PAHBXY,  & ! OUT Noah MP only
                !    FPICEXY,RAINLSM,SNOWLSM,FORCTLSM,FORCQLSM,FORCPLSM,FORCZLSM,  & ! OUT Noah MP only
                !    FORCWLSM,ACC_SSOILXY,ACC_QINSURXY,ACC_QSEVAXY, ACC_ETRANIXY,  & ! IN/OUT Noah MP
                !    EFLXBXY, SOILENERGY, SNOWENERGY, CANHSXY,                     & ! OUT Noah MP only
                !    ACC_DWATERXY, ACC_PRCPXY, ACC_ECANXY,ACC_ETRANXY,ACC_EDIRXY,  & ! IN/OUT Noah MP
                   FSATXY, WSURFXY,                                              & ! IN/OUT Noah MP
                   SNRDSXY, SNFRXY, BCPHIXY, BCPHOXY,  & ! SNICAR variable
                   OCPHIXY, OCPHOXY, DUST1XY, DUST2XY, DUST3XY, DUST4XY, DUST5XY,& ! SNICAR variable
                   MassConcBCPHIXY, MassConcBCPHOXY, MassConcOCPHIXY,            & ! SNICAR variable
                   MassConcOCPHOXY, MassConcDUST1XY, MassConcDUST2XY,            & ! SNICAR variable
                   MassConcDUST3XY, MassConcDUST4XY, MassConcDUST5XY,            & ! SNICAR variable
                   ALBSOILDIRXY, ALBSOILDIFXY,                                   & ! SNICAR variable
                 ! BEXP_3D,SMCDRY_3D,SMCWLT_3D,SMCREF_3D,SMCMAX_3D,              & ! placeholders to activate 3D soil
                 ! DKSAT_3D,DWSAT_3D,PSISAT_3D,QUARTZ_3D,                        & ! placeholders to activate 3D soil
                 ! REFDK_2D,REFKDT_2D,                                           & ! placeholders to activate 3D soil
                 ! IRR_FRAC_2D,IRR_HAR_2D,IRR_LAI_2D,IRR_MAD_2D,FILOSS_2D,       & ! placeholders to activate 3D soil
                 ! SPRIR_RATE_2D,MICIR_RATE_2D,FIRTFAC_2D,IR_RAIN_2D,            & ! placeholders to activate 3D soil
                 ! BVIC_2D,AXAJ_2D,BXAJ_2D,XXAJ_2D,BDVIC_2D,GDVIC_2D,BBVIC_2D,   & ! placeholders to activate 3D soil
                 ! KLAT_FAC,TDSMC_FAC,TD_DC,TD_DCOEF,TD_DDRAIN,TD_RADI,TD_SPAC,  & ! placeholders to activate 3D soil
! #ifdef WRF_HYDRO
!                    sfcheadrt,INFXSRT,soldrain,qtiledrain,ZWATBLE2D,              & ! OUT WRF-Hydro only
! #endif
                   ids,ide,  jds,jde,  kds,kde,                                  & ! IN: HICAR dimension
                   ims,ime,  jms,jme,  kms,kme,                                  & ! IN: HICAR dimension
                   its,ite,  jts,jte,  kts,kte,                                  & ! IN: HICAR dimension
                   MP_RAINC,MP_RAINNC,MP_SHCV,MP_SNOW,MP_GRAUP,MP_HAIL           ) ! IN: HICAR forcing

!----------------------------------------------------------------

    use NoahmpIOVarType, only : NoahmpIO_type
    use NoahmpIOVarInitMod
    use NoahmpReadTableMod
   !  use SnowInputSnicarMod
    use NoahmpDriverMainMod
    use LanduseConvertMod
    ! use module_sf_urban,    only: IRI_SCHEME

    implicit none

    type(NoahmpIO_type),                             intent(inout) ::  NoahmpIO

    ! input
    INTEGER,                                         INTENT(IN   ) ::  ids,ide, jds,jde, kds,kde,  &  ! d -> domain
                                                                       ims,ime, jms,jme, kms,kme,  &  ! m -> memory
                                                                       its,ite, jts,jte, kts,kte      ! t -> tile
    INTEGER,                                         INTENT(IN   ) ::  ITIMESTEP    ! timestep number
    INTEGER,                                         INTENT(IN   ) ::  YR           ! 4-digit year
    REAL,                                            INTENT(IN   ) ::  JULIAN       ! Julian day
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  COSZIN       ! cosine zenith angle
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLAT         ! latitude [rad]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLONG        ! latitude [rad]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  DZ8W         ! thickness of atmo layers [m]
    REAL,                                            INTENT(IN   ) ::  DT           ! timestep [s]
    REAL,    DIMENSION(1:nsoil),                     INTENT(IN   ) ::  DZS          ! thickness of soil layers [m]
    INTEGER,                                         INTENT(IN   ) ::  NSOIL        ! number of soil layers
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  IVGTYP       ! vegetation type
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  ISLTYP       ! soil type
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  VEGFRA       ! vegetation fraction []
    REAL,    DIMENSION( ims:ime ,         jms:jme ), INTENT(IN   ) ::  VEGMAX       ! annual max vegetation fraction []
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  TMN          ! deep soil temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLAND        ! =2 ocean; =1 land/seaice
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XICE         ! fraction of grid that is seaice
    REAL,    DIMENSION( ims:ime, nsoil*2, jms:jme ), INTENT(IN   ) ::  SOILCOMP     ! soil sand and clay percentage
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL1      ! soil texture in layer 1
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL2      ! soil texture in layer 2
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL3      ! soil texture in layer 3
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL4      ! soil texture in layer 4
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  T3D          ! 3D atmospheric temperature valid at mid-levels [K]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  QV3D         ! 3D water vapor mixing ratio [kg/kg_dry]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  U_PHY        ! 3D U wind component [m/s]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  V_PHY        ! 3D V wind component [m/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDOWN       ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDDIF       ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDDIR       ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  GLW          ! longwave down at surface [W m-2]
    REAL,    DIMENSION( ims:ime, kms:kme+1, jms:jme ), INTENT(IN   ) ::  P8W3D        ! 3D pressure, valid at interface [Pa]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  PRECIP_IN    ! total input precipitation [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SR           ! frozen precipitation ratio [-]
    ! Optional Detailed Precipitation Partitioning Inputs
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_RAINC  ! convective precipitation entering land model [mm] ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_RAINNC ! large-scale precipitation entering land model [mm]! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_SHCV   ! shallow conv precip entering land model [mm]      ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_SNOW   ! snow precipitation entering land model [mm]       ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_GRAUP  ! graupel precipitation entering land model [mm]    ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_HAIL   ! hail precipitation entering land model [mm]       ! MB/AN : v3.7
    ! Crop Model
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  CROPCAT     ! crop catagory
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  PLANTING    ! planting date
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  HARVEST     ! harvest date
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SEASON_GDD  ! growing season GDD
    ! Tile drain variables    
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  TD_FRACTION
    !2D inout irrigation variables 
    CHARACTER(LEN=256),                              INTENT(IN   ) ::  LLANDUSE    ! landuse data name (USGS or MODIS_IGBP)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  IRFRACT     ! irrigation fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SIFRACT     ! sprinkler irrigation fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  MIFRACT     ! micro irrigation fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  FIFRACT     ! flood irrigation fraction   

    ! placeholders for 2D/3D soil
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  BEXP_3D       ! C-H B exponent
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCDRY_3D     ! Soil Moisture Limit: Dry
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCWLT_3D     ! Soil Moisture Limit: Wilt
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCREF_3D     ! Soil Moisture Limit: Reference
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCMAX_3D     ! Soil Moisture Limit: Max
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  DKSAT_3D      ! Saturated Soil Conductivity
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  DWSAT_3D      ! Saturated Soil Diffusivity
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  PSISAT_3D     ! Saturated Matric Potential
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  QUARTZ_3D     ! Soil quartz content
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  REFDK_2D      ! Reference Soil Conductivity
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  REFKDT_2D     ! Soil Infiltration Parameter
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BVIC_2D       ! VIC model infiltration parameter [-] for opt_run=6
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  AXAJ_2D       ! Xinanjiang: Tension water distribution inflection parameter [-] for opt_run=7
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BXAJ_2D       ! Xinanjiang: Tension water distribution shape parameter [-] for opt_run=7
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  XXAJ_2D       ! Xinanjiang: Free water distribution shape parameter [-] for opt_run=7
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BDVIC_2D      ! VIC model infiltration parameter [-]  
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  GDVIC_2D      ! Mean Capillary Drive (m) for infiltration models
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BBVIC_2D      ! DVIC heterogeniety paramater [-]

    ! placeholders for 2D irrigation parameters
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_FRAC_2D   ! irrigation Fraction
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_HAR_2D    ! number of days before harvest date to stop irrigation 
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_LAI_2D    ! Minimum lai to trigger irrigation
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_MAD_2D    ! management allowable deficit (0-1)
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  FILOSS_2D     ! fraction of flood irrigation loss (0-1) 
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  SPRIR_RATE_2D ! mm/h, sprinkler irrigation rate
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  MICIR_RATE_2D ! mm/h, micro irrigation rate
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  FIRTFAC_2D    ! flood application rate factor
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IR_RAIN_2D    ! maximum precipitation to stop irrigation trigger

    ! placeholders for 2D tile drainage parameters
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  KLAT_FAC      ! factor multiplier to hydraulic conductivity
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TDSMC_FAC     ! factor multiplier to field capacity
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_DC         ! drainage coefficient for simple
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_DCOEF      ! drainge coefficient for Hooghoudt 
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_DDRAIN     ! depth of drain
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_RADI       ! tile radius
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_SPAC       ! tile spacing

    ! INOUT (with generic LSM equivalent)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TSK          ! surface radiative temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  HFX          ! sensible heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QFX          ! latent heat flux [kg s-1 m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  LH           ! latent heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GRDFLX       ! ground/snow heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMSTAV       ! soil moisture avail. [not used]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMSTOT       ! total soil water [mm][not used]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SFCRUNOFF    ! accumulated surface runoff [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  UDRUNOFF     ! accumulated sub-surface runoff [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ALBEDO       ! total grid albedo []
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOWC        ! snow cover fraction []
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SMOIS        ! volumetric soil moisture [m3/m3]
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SH2O         ! volumetric liquid soil moisture [m3/m3]
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  TSLB         ! soil temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOW         ! snow water equivalent [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOWH        ! physical snow depth [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANWAT       ! total canopy water + ice [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACSNOM       ! accumulated snow melt (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACSNOW       ! accumulated snow on grid
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  EMISS        ! surface bulk emissivity
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QSFC         ! bulk surface specific humidity
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  Z0           ! combined z0 sent to coupled model
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RS           ! Total stomatal resistance (s/m)
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ISNOWXY      ! actual no. of snow layers
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TVXY         ! vegetation leaf temperature
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TGXY         ! bulk ground surface temperature
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANICEXY     ! canopy-intercepted ice (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANLIQXY     ! canopy-intercepted liquid water (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  EAHXY        ! canopy air vapor pressure (pa)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TAHXY        ! canopy air temperature (k)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CMXY         ! bulk momentum drag coefficient
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CHXY         ! bulk sensible heat exchange coefficient
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  FWETXY       ! wetted or snowed fraction of the canopy (-)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNEQVOXY     ! snow mass at last time step(mm h2o)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ALBOLDXY     ! snow albedo at last time step (-)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QSNOWXY      ! snowfall on the ground [mm/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QRAINXY      ! rainfall on the ground [mm/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WSLAKEXY     ! lake water storage [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ZWTXY        ! water table depth [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WAXY         ! water in the "aquifer" [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WTXY         ! groundwater storage [mm]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  TSNOXY       ! snow temperature [K]
    REAL,    DIMENSION( ims:ime,-2:NSOIL, jms:jme ), INTENT(INOUT) ::  ZSNSOXY      ! snow layer depth [m]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNICEXY      ! snow layer ice [mm]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNLIQXY      ! snow layer liquid water [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  LFMASSXY     ! leaf mass [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RTMASSXY     ! mass of fine roots [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  STMASSXY     ! stem mass [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WOODXY       ! mass of wood (incl. woody roots) [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  STBLCPXY     ! stable carbon in deep soil [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  FASTCPXY     ! short-lived carbon, shallow soil [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  XLAIXY       ! leaf area index
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  XSAIXY       ! stem area index
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TAUSSXY      ! snow age factor
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SMOISEQ      ! eq volumetric soil moisture [m3/m3]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMCWTDXY     ! soil moisture content in the layer to the water table when deep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  DEEPRECHXY   ! recharge to the water table when deep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RECHXY       ! recharge to the water table (diagnostic) 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GRAINXY      ! mass of grain XING [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GDDXY        ! growing degree days XING (based on 10C) 
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  PGSXY        ! growing stage
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_SSOILXY  ! m/s * soil_dt/main_dt
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_QINSURXY ! m/s * soil_dt/main_dt
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_QSEVAXY  ! m/s * soil_dt/main_dt
    ! REAL,    DIMENSION( ims:ime, 1:NSOIL, jms:jme ), INTENT(INOUT) ::  ACC_ETRANIXY ! m/s * soil_dt/main_dt
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_DWATERXY ! m/s * soil_dt/main_dt
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_PRCPXY   ! m/s * soil_dt/main_dt
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_ECANXY   ! m/s * soil_dt/main_dt
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_ETRANXY  ! m/s * soil_dt/main_dt
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_EDIRXY   ! m/s * soil_dt/main_dt
    !2D inout irrigation variables 
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRNUMSI      ! irrigation event number, Sprinkler
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRNUMMI      ! irrigation event number, Micro
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRNUMFI      ! irrigation event number, Flood 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRWATSI      ! irrigation water amount [m] to be applied, Sprinkler
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRWATMI      ! irrigation water amount [m] to be applied, Micro
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRWATFI      ! irrigation water amount [m] to be applied, Flood
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRELOSS      ! loss of irrigation water to evaporation,sprinkler [m/timestep]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRSIVOL      ! amount of irrigation by sprinkler (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRMIVOL      ! amount of irrigation by micro (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRFIVOL      ! amount of irrigation by micro (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRRSPLH      ! latent heating from sprinkler evaporation (w/m2)    
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QTDRAIN      ! Tile drain
    ! wetland state varible
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  FSATXY       ! saturated fraction of the grid (-)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WSURFXY      ! wetland water storage [mm]
    ! SNICAR state variable
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNRDSXY         ! snow layer effective grain radius [microns, m-6]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNFRXY          ! snow layer rate of snow freezing [mm/s]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  BCPHIXY         ! mass of hydrophillic Black Carbon in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  BCPHOXY         ! mass of hydrophobic Black Carbon in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  OCPHIXY         ! mass of hydrophillic Organic Carbon in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  OCPHOXY         ! mass of hydrophobic Organic Carbon in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  DUST1XY         ! mass of dust species 1 in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  DUST2XY         ! mass of dust species 2 in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  DUST3XY         ! mass of dust species 3 in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  DUST4XY         ! mass of dust species 4 in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  DUST5XY         ! mass of dust species 5 in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcBCPHIXY ! mass concentration of hydrophillic Black Carbon in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcBCPHOXY ! mass concentration of hydrophobic Black Carbon in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcOCPHIXY ! mass concentration of hydrophillic Organic Carbon in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcOCPHOXY ! mass concentration of hydrophobic Organic Carbon in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcDUST1XY ! mass concentration of dust species 1 in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcDUST2XY ! mass concentration of dust species 2 in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcDUST3XY ! mass concentration of dust species 3 in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcDUST4XY ! mass concentration of dust species 4 in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcDUST5XY ! mass concentration of dust species 5 in snow [kg/kg]
    REAL,    DIMENSION(ims:ime,1:2,       jms:jme),  INTENT(INOUT) ::  ALBSOILDIRXY    ! soil albedo direct
    REAL,    DIMENSION(ims:ime,1:2,       jms:jme),  INTENT(INOUT) ::  ALBSOILDIFXY    ! soil albedo diffuse

    ! OUT (with no Noah LSM equivalent)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  T2MVXY       ! 2m temperature of vegetation part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  T2MBXY       ! 2m temperature of bare ground part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  Q2MVXY       ! 2m mixing ratio of vegetation part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  Q2MBXY       ! 2m mixing ratio of bare ground part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TRADXY       ! surface radiative temperature (k)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  NEEXY        ! net ecosys exchange (g/m2/s CO2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GPPXY        ! gross primary assimilation [g/m2/s C]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  NPPXY        ! net primary productivity [g/m2/s C]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FVEGXY       ! Noah-MP vegetation fraction [-]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RUNSFXY      ! surface runoff [mm] per soil timestep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RUNSBXY      ! subsurface runoff [mm] per soil timestep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  ECANXY       ! evaporation of intercepted water (mm/s)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EDIRXY       ! soil surface evaporation rate (mm/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  ETRANXY      ! transpiration rate (mm/s)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FSAXY        ! total absorbed solar radiation (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FIRAXY       ! total net longwave rad (w/m2) [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  APARXY       ! photosyn active energy by canopy (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PSNXY        ! total photosynthesis (umol co2/m2/s) [+]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SAVXY        ! solar rad absorbed by veg. (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SAGXY        ! solar rad absorbed by ground (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RSSUNXY      ! sunlit leaf stomatal resistance (s/m)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RSSHAXY      ! shaded leaf stomatal resistance (s/m)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  BGAPXY       ! between gap fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  WGAPXY       ! within gap fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TGVXY        ! under canopy ground temperature[K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TGBXY        ! bare ground temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHVXY        ! sensible heat exchange coefficient vegetated
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHBXY        ! sensible heat exchange coefficient bare-ground
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHGXY        ! veg ground sen. heat [w/m2]   [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHCXY        ! canopy sen. heat [w/m2]   [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHBXY        ! bare sensible heat [w/m2]     [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVGXY        ! veg ground evap. heat [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVBXY        ! bare soil evaporation [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GHVXY        ! veg ground heat flux [w/m2]  [+ to soil]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GHBXY        ! bare ground heat flux [w/m2] [+ to soil]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRGXY        ! veg ground net LW rad. [w/m2] [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRCXY        ! canopy net LW rad. [w/m2] [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRBXY        ! bare net longwave rad. [w/m2] [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TRXY         ! transpiration [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVCXY        ! canopy evaporation heat [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHLEAFXY     ! leaf exchange coefficient 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHUCXY       ! under canopy exchange coefficient 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHV2XY       ! veg 2m exchange coefficient 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHB2XY       ! bare 2m exchange coefficient
    ! additional output variables
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHXY        ! precipitation advected heat 
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHGXY       ! precipitation advected heat 
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHBXY       ! precipitation advected heat 
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHVXY       ! precipitation advected heat 
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QINTSXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QINTRXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDRIPSXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDRIPRXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QTHROSXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QTHRORXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNSUBXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNFROXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSUBCXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QFROCXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QEVACXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDEWCXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QFRZCXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QMELTCXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNBOTXY     ! total liquid water (snowmelt + rain through pack)out of snowpack bottom [mm/s]
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QMELTXY      ! snowmelt due to phase change (mm/s)
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PONDINGXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FPICEXY      ! fraction of ice in precip
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RAINLSM      ! rain rate                   (mm/s)  AJN
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SNOWLSM      ! liquid equivalent snow rate (mm/s)  AJN
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCTLSM
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCQLSM
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCPLSM
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCZLSM
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCWLSM
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EFLXBXY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SOILENERGY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SNOWENERGY
    ! REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CANHSXY
    ! local
    INTEGER :: I, J
    INTEGER :: LoopInd    ! loop index for array section expansion

! ----------------------------------------------------------------------

    !--------- Input variable mapping start ---------

    NoahmpIO%YR                 = YR
    NoahmpIO%JULIAN             = JULIAN
    NoahmpIO%DTBL               = DT
    NoahmpIO%DZS                = DZS
    ! NoahmpIO%IRI_URBAN          = IRI_SCHEME
    NoahmpIO%ITIMESTEP          = ITIMESTEP

    !$acc update device(NoahmpIO%YR, NoahmpIO%JULIAN, NoahmpIO%DTBL, NoahmpIO%DZS, NoahmpIO%ITIMESTEP)

    ! input HICAR variables mapped to NoahmpIO variables
    ! 2D/3D variables
    !$acc parallel loop gang vector collapse(2) default(present) private(LoopInd) firstprivate(kte, kts, NSOIL)
    do J = jts, jte
    do I = its, ite

    if ( NoahmpIO%IOPT_SOIL > 1 ) then
       !$acc loop seq
       do LoopInd = 1, nsoil*2
          NoahmpIO%SOILCOMP(I,LoopInd,J)        = SOILCOMP(I,LoopInd,J)
       enddo
       NoahmpIO%SOILCL1(I,J)           = SOILCL1(I,J)
       NoahmpIO%SOILCL2(I,J)           = SOILCL2(I,J)
       NoahmpIO%SOILCL3(I,J)           = SOILCL3(I,J)
       NoahmpIO%SOILCL4(I,J)           = SOILCL4(I,J)
    endif
    NoahmpIO%IVGTYP(I,J)               = IVGTYP(I,J)
    NoahmpIO%ISLTYP(I,J)               = ISLTYP(I,J)
    NoahmpIO%VEGFRA(I,J)               = VEGFRA(I,J)
    NoahmpIO%GVFMAX(I,J)               = VEGMAX(I,J)
    NoahmpIO%TMN(I,J)                  = TMN(I,J)
    NoahmpIO%XLAND(I,J)                = XLAND(I,J)
    NoahmpIO%XICE(I,J)                 = XICE(I,J)
    NoahmpIO%CROPCAT(I,J)              = CROPCAT(I,J)
    NoahmpIO%PLANTING(I,J)             = PLANTING(I,J)
    NoahmpIO%HARVEST(I,J)              = HARVEST(I,J)
    NoahmpIO%SEASON_GDD(I,J)           = SEASON_GDD(I,J)
    NoahmpIO%XLAT(I,J)                 = XLAT(I,J)
    NoahmpIO%XLONG(I,J)                = XLONG(I,J)
    NoahmpIO%COSZEN(I,J)               = COSZIN(I,J)
    !$acc loop seq
    do LoopInd = kts, kte
       NoahmpIO%DZ8W(I,LoopInd,J)               = DZ8W(I,LoopInd,J)
       NoahmpIO%T_PHY(I,LoopInd,J)              = T3D(I,LoopInd,J)
       NoahmpIO%QV_CURR(I,LoopInd,J)            = QV3D(I,LoopInd,J)
       NoahmpIO%U_PHY(I,LoopInd,J)              = U_PHY(I,LoopInd,J)
       NoahmpIO%V_PHY(I,LoopInd,J)              = V_PHY(I,LoopInd,J)
       NoahmpIO%P8W(I,LoopInd,J) = P8W3D(I,LoopInd,J)
    enddo
    NoahmpIO%SWDOWN(I,J)               = SWDOWN(I,J)
    NoahmpIO%SWDDIR(I,J)               = SWDDIR(I,J)
    NoahmpIO%SWDDIF(I,J)               = SWDDIF(I,J)
    NoahmpIO%RadSwDirFrac(I,J)         = min(SWDDIR(I,J),SWDOWN(I,J))/max(SWDOWN(I,J), 0.001)
    NoahmpIO%GLW(I,J)                  = GLW(I,J)
    NoahmpIO%RAINBL(I,J)               = PRECIP_IN(I,J)
    NoahmpIO%SR(I,J)                   = SR(I,J)
    NoahmpIO%IRFRACT(I,J)              = IRFRACT(I,J)
    NoahmpIO%SIFRACT(I,J)              = SIFRACT(I,J)
    NoahmpIO%MIFRACT(I,J)              = MIFRACT(I,J)
    NoahmpIO%FIFRACT(I,J)              = FIFRACT(I,J)
    ! NoahmpIO%TD_FRACTION(I,J)          = TD_FRACTION(I,J)
    if (present(MP_RAINC) .and. present(MP_RAINNC) .and. &
        present(MP_SHCV)  .and. present(MP_SNOW)   .and. &
        present(MP_GRAUP) .and. present(MP_HAIL) ) then
       NoahmpIO%MP_RAINC(I,J)          = MP_RAINC(I,J)
       NoahmpIO%MP_RAINNC(I,J)         = MP_RAINNC(I,J)
       NoahmpIO%MP_SHCV(I,J)           = MP_SHCV(I,J)
       NoahmpIO%MP_SNOW(I,J)           = MP_SNOW(I,J)
       NoahmpIO%MP_GRAUP(I,J)          = MP_GRAUP(I,J)
       NoahmpIO%MP_HAIL(I,J)           = MP_HAIL(I,J)
    endif
    ! NoahmpIO%BEXP_3D(I,:,J)          = BEXP_3D(I,:,J)
    ! NoahmpIO%SMCDRY_3D(I,:,J)        = SMCDRY_3D(I,:,J)
    ! NoahmpIO%SMCWLT_3D(I,:,J)        = SMCWLT_3D(I,:,J)
    ! NoahmpIO%SMCREF_3D(I,:,J)        = SMCREF_3D(I,:,J)
    ! NoahmpIO%SMCMAX_3D(I,:,J)        = SMCMAX_3D(I,:,J)
    ! NoahmpIO%DKSAT_3D(I,:,J)         = DKSAT_3D(I,:,J)
    ! NoahmpIO%DWSAT_3D(I,:,J)         = DWSAT_3D(I,:,J)
    ! NoahmpIO%PSISAT_3D(I,:,J)        = PSISAT_3D(I,:,J)
    ! NoahmpIO%QUARTZ_3D(I,:,J)        = QUARTZ_3D(I,:,J)
    ! NoahmpIO%REFDK_2D(I,J)           = REFDK_2D(I,J)
    ! NoahmpIO%REFKDT_2D(I,J)          = REFKDT_2D(I,J)
    ! NoahmpIO%IRR_FRAC_2D(I,J)        = IRR_FRAC_2D(I,J)
    ! NoahmpIO%IRR_HAR_2D(I,J)         = IRR_HAR_2D(I,J)
    ! NoahmpIO%IRR_LAI_2D(I,J)         = IRR_LAI_2D(I,J)
    ! NoahmpIO%IRR_MAD_2D(I,J)         = IRR_MAD_2D(I,J)
    ! NoahmpIO%FILOSS_2D(I,J)          = FILOSS_2D(I,J)
    ! NoahmpIO%SPRIR_RATE_2D(I,J)      = SPRIR_RATE_2D(I,J)
    ! NoahmpIO%MICIR_RATE_2D(I,J)      = MICIR_RATE_2D(I,J)
    ! NoahmpIO%FIRTFAC_2D(I,J)         = FIRTFAC_2D(I,J)
    ! NoahmpIO%IR_RAIN_2D(I,J)         = IR_RAIN_2D(I,J)
    ! NoahmpIO%BVIC_2D(I,J)            = BVIC_2D(I,J)
    ! NoahmpIO%AXAJ_2D(I,J)            = AXAJ_2D(I,J)
    ! NoahmpIO%BXAJ_2D(I,J)            = BXAJ_2D(I,J)
    ! NoahmpIO%XXAJ_2D(I,J)            = XXAJ_2D(I,J)
    ! NoahmpIO%BDVIC_2D(I,J)           = BDVIC_2D(I,J)
    ! NoahmpIO%GDVIC_2D(I,J)           = GDVIC_2D(I,J)
    ! NoahmpIO%BBVIC_2D(I,J)           = BBVIC_2D(I,J)
    ! NoahmpIO%KLAT_FAC(I,J)           = KLAT_FAC(I,J)
    ! NoahmpIO%TDSMC_FAC(I,J)          = TDSMC_FAC(I,J)
    ! NoahmpIO%TD_DC(I,J)              = TD_DC(I,J)
    ! NoahmpIO%TD_DCOEF(I,J)           = TD_DCOEF(I,J)
    ! NoahmpIO%TD_DDRAIN(I,J)          = TD_DDRAIN(I,J)
    ! NoahmpIO%TD_RADI(I,J)            = TD_RADI(I,J)
    ! NoahmpIO%TD_SPAC(I,J)            = TD_SPAC(I,J)
    
    ! in/out HICAR variables mapped to NoahmpIO variables
    NoahmpIO%TSK(I,J)                  = TSK(I,J)
    NoahmpIO%HFX(I,J)                  = HFX(I,J)
    NoahmpIO%QFX(I,J)                  = QFX(I,J)
    NoahmpIO%LH(I,J)                   = LH(I,J)
    NoahmpIO%GRDFLX(I,J)               = GRDFLX(I,J)
    NoahmpIO%SMSTAV(I,J)               = SMSTAV(I,J)
    NoahmpIO%SMSTOT(I,J)               = SMSTOT(I,J)
    NoahmpIO%SFCRUNOFF(I,J)            = SFCRUNOFF(I,J)
    NoahmpIO%UDRUNOFF(I,J)             = UDRUNOFF(I,J)
    NoahmpIO%ALBEDO(I,J)               = ALBEDO(I,J)
    NoahmpIO%SNOWC(I,J)                = SNOWC(I,J)
    !$acc loop seq
    do LoopInd = 1, NSOIL
       NoahmpIO%SMOIS(I,LoopInd,J)              = SMOIS(I,LoopInd,J)
       NoahmpIO%SH2O(I,LoopInd,J)               = SH2O(I,LoopInd,J)
       NoahmpIO%TSLB(I,LoopInd,J)               = TSLB(I,LoopInd,J)
    enddo
    NoahmpIO%SNOW(I,J)                 = SNOW(I,J)
    NoahmpIO%SNOWH(I,J)                = SNOWH(I,J)
    NoahmpIO%CANWAT(I,J)               = CANWAT(I,J)
    NoahmpIO%CANICEXY(I,J)             = CANICEXY(I,J)
    NoahmpIO%CANLIQXY(I,J)             = CANLIQXY(I,J)
    NoahmpIO%ACSNOM(I,J)               = ACSNOM(I,J)
    NoahmpIO%ACSNOW(I,J)               = ACSNOW(I,J)
    NoahmpIO%EMISS(I,J)                = EMISS(I,J)
    NoahmpIO%QSFC(I,J)                 = QSFC(I,J)
    NoahmpIO%Z0(I,J)                   = Z0(I,J)
    NoahmpIO%IRNUMSI(I,J)              = IRNUMSI(I,J)
    NoahmpIO%IRNUMMI(I,J)              = IRNUMMI(I,J)
    NoahmpIO%IRNUMFI(I,J)              = IRNUMFI(I,J)
    NoahmpIO%IRWATSI(I,J)              = IRWATSI(I,J)
    NoahmpIO%IRWATMI(I,J)              = IRWATMI(I,J)
    NoahmpIO%IRWATFI(I,J)              = IRWATFI(I,J)
    NoahmpIO%IRELOSS(I,J)              = IRELOSS(I,J)
    NoahmpIO%IRSIVOL(I,J)              = IRSIVOL(I,J)
    NoahmpIO%IRMIVOL(I,J)              = IRMIVOL(I,J)
    NoahmpIO%IRFIVOL(I,J)              = IRFIVOL(I,J)
    NoahmpIO%IRRSPLH(I,J)              = IRRSPLH(I,J)
    NoahmpIO%ISNOWXY(I,J)              = ISNOWXY(I,J)
    NoahmpIO%TVXY(I,J)                 = TVXY(I,J)
    NoahmpIO%TGXY(I,J)                 = TGXY(I,J)
    NoahmpIO%EAHXY(I,J)                = EAHXY(I,J)
    NoahmpIO%TAHXY(I,J)                = TAHXY(I,J)
    NoahmpIO%CMXY(I,J)                 = CMXY(I,J)
    NoahmpIO%CHXY(I,J)                 = CHXY(I,J)
    NoahmpIO%FWETXY(I,J)               = FWETXY(I,J)
    NoahmpIO%SNEQVOXY(I,J)             = SNEQVOXY(I,J)
    NoahmpIO%ALBOLDXY(I,J)             = ALBOLDXY(I,J)
    NoahmpIO%QSNOWXY(I,J)              = QSNOWXY(I,J)
    NoahmpIO%QRAINXY(I,J)              = QRAINXY(I,J)
    NoahmpIO%WSLAKEXY(I,J)             = WSLAKEXY(I,J)
    NoahmpIO%ZWTXY(I,J)                = ZWTXY(I,J)
    NoahmpIO%WAXY(I,J)                 = WAXY(I,J)
    NoahmpIO%WTXY(I,J)                 = WTXY(I,J)
    !$acc loop seq
    do LoopInd = -2, 0
       NoahmpIO%TSNOXY(I,LoopInd,J)             = TSNOXY(I,LoopInd,J)
       NoahmpIO%ZSNSOXY(I,LoopInd,J)            = ZSNSOXY(I,LoopInd,J)
       NoahmpIO%SNICEXY(I,LoopInd,J)            = SNICEXY(I,LoopInd,J)
       NoahmpIO%SNLIQXY(I,LoopInd,J)            = SNLIQXY(I,LoopInd,J)
    enddo
    NoahmpIO%LFMASSXY(I,J)             = LFMASSXY(I,J)
    NoahmpIO%RTMASSXY(I,J)             = RTMASSXY(I,J)
    NoahmpIO%STMASSXY(I,J)             = STMASSXY(I,J)
    NoahmpIO%WOODXY(I,J)               = WOODXY(I,J)
    NoahmpIO%STBLCPXY(I,J)             = STBLCPXY(I,J)
    NoahmpIO%FASTCPXY(I,J)             = FASTCPXY(I,J)
    NoahmpIO%LAI(I,J)                  = XLAIXY(I,J)
    NoahmpIO%XSAIXY(I,J)               = XSAIXY(I,J)
    NoahmpIO%TAUSSXY(I,J)              = TAUSSXY(I,J)
    !$acc loop seq
    do LoopInd = 1, NSOIL
       NoahmpIO%SMOISEQ(I,LoopInd,J)            = SMOISEQ(I,LoopInd,J)
    enddo
    NoahmpIO%SMCWTDXY(I,J)             = SMCWTDXY(I,J)
    NoahmpIO%DEEPRECHXY(I,J)           = DEEPRECHXY(I,J)
    NoahmpIO%RECHXY(I,J)               = RECHXY(I,J)
    NoahmpIO%GRAINXY(I,J)              = GRAINXY(I,J)
    NoahmpIO%GDDXY(I,J)                = GDDXY(I,J)
    NoahmpIO%PGSXY(I,J)                = PGSXY(I,J)
    ! NoahmpIO%QTDRAIN(I,J)              = QTDRAIN(I,J)
    NoahmpIO%RS(I,J)                   = RS(I,J)
    ! NoahmpIO%ACC_SSOILXY(I,J)          = ACC_SSOILXY(I,J)
    ! NoahmpIO%ACC_QINSURXY(I,J)         = ACC_QINSURXY(I,J)
    ! NoahmpIO%ACC_QSEVAXY(I,J)          = ACC_QSEVAXY(I,J)
    ! NoahmpIO%ACC_ETRANIXY(I,:,J)       = ACC_ETRANIXY(I,:,J)
    ! NoahmpIO%ACC_DWATERXY(I,J)         = ACC_DWATERXY(I,J)
    ! NoahmpIO%ACC_PRCPXY(I,J)           = ACC_PRCPXY(I,J)
    ! NoahmpIO%ACC_ECANXY(I,J)           = ACC_ECANXY(I,J)
    ! NoahmpIO%ACC_ETRANXY(I,J)          = ACC_ETRANXY(I,J)
    ! NoahmpIO%ACC_EDIRXY(I,J)           = ACC_EDIRXY(I,J)
    !$acc loop seq
    do LoopInd = 1, 2
       NoahmpIO%ALBSOILDIRXY(I,LoopInd,J)       = ALBSOILDIRXY(I,LoopInd,J)
       NoahmpIO%ALBSOILDIFXY(I,LoopInd,J)       = ALBSOILDIFXY(I,LoopInd,J)
    enddo
    if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       NoahmpIO%FSATXY(I,J)            = FSATXY(I,J)
       NoahmpIO%WSURFXY(I,J)           = WSURFXY(I,J)
    endif
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       !$acc loop seq
       do LoopInd = -2, 0
          NoahmpIO%SNRDSXY(I,LoopInd,J)         = SNRDSXY(I,LoopInd,J)
          NoahmpIO%SNFRXY(I,LoopInd,J)          = SNFRXY(I,LoopInd,J)
          NoahmpIO%BCPHIXY(I,LoopInd,J)         = BCPHIXY(I,LoopInd,J)
          NoahmpIO%BCPHOXY(I,LoopInd,J)         = BCPHOXY(I,LoopInd,J)
          NoahmpIO%OCPHIXY(I,LoopInd,J)         = OCPHIXY(I,LoopInd,J)
          NoahmpIO%OCPHOXY(I,LoopInd,J)         = OCPHOXY(I,LoopInd,J)
          NoahmpIO%DUST1XY(I,LoopInd,J)         = DUST1XY(I,LoopInd,J)
          NoahmpIO%DUST2XY(I,LoopInd,J)         = DUST2XY(I,LoopInd,J)
          NoahmpIO%DUST3XY(I,LoopInd,J)         = DUST3XY(I,LoopInd,J)
          NoahmpIO%DUST4XY(I,LoopInd,J)         = DUST4XY(I,LoopInd,J)
          NoahmpIO%DUST5XY(I,LoopInd,J)         = DUST5XY(I,LoopInd,J)
          NoahmpIO%MassConcBCPHIXY(I,LoopInd,J) = MassConcBCPHIXY(I,LoopInd,J)
          NoahmpIO%MassConcBCPHOXY(I,LoopInd,J) = MassConcBCPHOXY(I,LoopInd,J)
          NoahmpIO%MassConcOCPHIXY(I,LoopInd,J) = MassConcOCPHIXY(I,LoopInd,J)
          NoahmpIO%MassConcOCPHOXY(I,LoopInd,J) = MassConcOCPHOXY(I,LoopInd,J)
          NoahmpIO%MassConcDUST1XY(I,LoopInd,J) = MassConcDUST1XY(I,LoopInd,J)
          NoahmpIO%MassConcDUST2XY(I,LoopInd,J) = MassConcDUST2XY(I,LoopInd,J)
          NoahmpIO%MassConcDUST3XY(I,LoopInd,J) = MassConcDUST3XY(I,LoopInd,J)
          NoahmpIO%MassConcDUST4XY(I,LoopInd,J) = MassConcDUST4XY(I,LoopInd,J)
          NoahmpIO%MassConcDUST5XY(I,LoopInd,J) = MassConcDUST5XY(I,LoopInd,J)
       enddo
    endif
! #ifdef WRF_HYDRO
!     NoahmpIO%sfcheadrt(I,J)            = sfcheadrt(I,J)
!     NoahmpIO%INFXSRT(I,J)              = INFXSRT(I,J)
!     NoahmpIO%soldrain(I,J)             = soldrain(I,J)
!     NoahmpIO%qtiledrain(I,J)           = qtiledrain(I,J)
!     NoahmpIO%ZWATBLE2D(I,J)            = ZWATBLE2D(I,J)              
! #endif

    enddo ! I
    enddo ! J
    !--------- Input variable mapping end ---------


    !-------- call main Noah-MP driver ------------
    call NoahmpDriverMain(NoahmpIO)


    !--------- Output variable mapping start ---------
    !$acc parallel loop gang vector collapse(2) default(present) private(LoopInd)
    do J = jts, jte
    do I = its, ite

    ! in/out NoahmpIO variables mapped to HICAR variables
    TSK(I,J)            = NoahmpIO%TSK(I,J)
    HFX(I,J)            = NoahmpIO%HFX(I,J)
    QFX(I,J)            = NoahmpIO%QFX(I,J)
    LH(I,J)             = NoahmpIO%LH(I,J)
    GRDFLX(I,J)         = NoahmpIO%GRDFLX(I,J)
    SMSTAV(I,J)         = NoahmpIO%SMSTAV(I,J)
    SMSTOT(I,J)         = NoahmpIO%SMSTOT(I,J)
    SFCRUNOFF(I,J)      = NoahmpIO%SFCRUNOFF(I,J)
    UDRUNOFF(I,J)       = NoahmpIO%UDRUNOFF(I,J)
    if (NoahmpIO%ALBEDO(I,J) > 0.0) ALBEDO(I,J)         = NoahmpIO%ALBEDO(I,J)
    SNOWC(I,J)          = NoahmpIO%SNOWC(I,J)
    SMOIS(I,:,J)        = NoahmpIO%SMOIS(I,:,J)
    SH2O(I,:,J)         = NoahmpIO%SH2O(I,:,J)
    TSLB(I,:,J)         = NoahmpIO%TSLB(I,:,J)
    SNOW(I,J)           = NoahmpIO%SNOW(I,J)
    SNOWH(I,J)          = NoahmpIO%SNOWH(I,J)
    CANWAT(I,J)         = NoahmpIO%CANWAT(I,J)
    CANICEXY(I,J)       = NoahmpIO%CANICEXY(I,J)
    CANLIQXY(I,J)       = NoahmpIO%CANLIQXY(I,J)
    ACSNOM(I,J)         = NoahmpIO%ACSNOM(I,J)
    ACSNOW(I,J)         = NoahmpIO%ACSNOW(I,J)
    EMISS(I,J)          = NoahmpIO%EMISS(I,J)
    QSFC(I,J)           = NoahmpIO%QSFC(I,J)
    Z0(I,J)             = NoahmpIO%Z0(I,J)
    IRNUMSI(I,J)        = NoahmpIO%IRNUMSI(I,J)
    IRNUMMI(I,J)        = NoahmpIO%IRNUMMI(I,J)
    IRNUMFI(I,J)        = NoahmpIO%IRNUMFI(I,J)
    IRWATSI(I,J)        = NoahmpIO%IRWATSI(I,J)
    IRWATMI(I,J)        = NoahmpIO%IRWATMI(I,J)
    IRWATFI(I,J)        = NoahmpIO%IRWATFI(I,J)
    IRELOSS(I,J)        = NoahmpIO%IRELOSS(I,J)
    IRSIVOL(I,J)        = NoahmpIO%IRSIVOL(I,J)
    IRMIVOL(I,J)        = NoahmpIO%IRMIVOL(I,J)
    IRFIVOL(I,J)        = NoahmpIO%IRFIVOL(I,J)
    IRRSPLH(I,J)        = NoahmpIO%IRRSPLH(I,J)
    ISNOWXY(I,J)        = NoahmpIO%ISNOWXY(I,J)
    TVXY(I,J)           = NoahmpIO%TVXY(I,J)
    TGXY(I,J)           = NoahmpIO%TGXY(I,J)
    EAHXY(I,J)          = NoahmpIO%EAHXY(I,J)
    TAHXY(I,J)          = NoahmpIO%TAHXY(I,J)
    CMXY(I,J)           = NoahmpIO%CMXY(I,J)
    CHXY(I,J)           = NoahmpIO%CHXY(I,J)
    FWETXY(I,J)         = NoahmpIO%FWETXY(I,J)
    SNEQVOXY(I,J)       = NoahmpIO%SNEQVOXY(I,J)
    ALBOLDXY(I,J)       = NoahmpIO%ALBOLDXY(I,J)
    QSNOWXY(I,J)        = NoahmpIO%QSNOWXY(I,J)
    QRAINXY(I,J)        = NoahmpIO%QRAINXY(I,J)
    WSLAKEXY(I,J)       = NoahmpIO%WSLAKEXY(I,J)
    ZWTXY(I,J)          = NoahmpIO%ZWTXY(I,J)
    WAXY(I,J)           = NoahmpIO%WAXY(I,J)
    WTXY(I,J)           = NoahmpIO%WTXY(I,J)
    TSNOXY(I,:,J)       = NoahmpIO%TSNOXY(I,:,J)
    ZSNSOXY(I,:,J)      = NoahmpIO%ZSNSOXY(I,:,J)
    SNICEXY(I,:,J)      = NoahmpIO%SNICEXY(I,:,J)
    SNLIQXY(I,:,J)      = NoahmpIO%SNLIQXY(I,:,J)
    LFMASSXY(I,J)       = NoahmpIO%LFMASSXY(I,J)
    RTMASSXY(I,J)       = NoahmpIO%RTMASSXY(I,J)
    STMASSXY(I,J)       = NoahmpIO%STMASSXY(I,J)
    WOODXY(I,J)         = NoahmpIO%WOODXY(I,J)
    STBLCPXY(I,J)       = NoahmpIO%STBLCPXY(I,J)
    FASTCPXY(I,J)       = NoahmpIO%FASTCPXY(I,J)
    XLAIXY(I,J)         = NoahmpIO%LAI(I,J)
    XSAIXY(I,J)         = NoahmpIO%XSAIXY(I,J)
    TAUSSXY(I,J)        = NoahmpIO%TAUSSXY(I,J)
    SMOISEQ(I,:,J)      = NoahmpIO%SMOISEQ(I,:,J)
    SMCWTDXY(I,J)       = NoahmpIO%SMCWTDXY(I,J)
    DEEPRECHXY(I,J)     = NoahmpIO%DEEPRECHXY(I,J)
    RECHXY(I,J)         = NoahmpIO%RECHXY(I,J)
    GRAINXY(I,J)        = NoahmpIO%GRAINXY(I,J)
    GDDXY(I,J)          = NoahmpIO%GDDXY(I,J)
    PGSXY(I,J)          = NoahmpIO%PGSXY(I,J)
    ! QTDRAIN(I,J)        = NoahmpIO%QTDRAIN(I,J)
    RS(I,J)             = NoahmpIO%RS(I,J)
    ! ACC_SSOILXY(I,J)    = NoahmpIO%ACC_SSOILXY(I,J)
    ! ACC_QINSURXY(I,J)   = NoahmpIO%ACC_QINSURXY(I,J)
    ! ACC_QSEVAXY(I,J)    = NoahmpIO%ACC_QSEVAXY(I,J)
    ! ACC_ETRANIXY(I,:,J) = NoahmpIO%ACC_ETRANIXY(I,:,J)
    ! ACC_DWATERXY(I,J)   = NoahmpIO%ACC_DWATERXY(I,J)
    ! ACC_PRCPXY(I,J)     = NoahmpIO%ACC_PRCPXY(I,J)
    ! ACC_ECANXY(I,J)     = NoahmpIO%ACC_ECANXY(I,J)
    ! ACC_ETRANXY(I,J)    = NoahmpIO%ACC_ETRANXY(I,J)
    ! ACC_EDIRXY(I,J)     = NoahmpIO%ACC_EDIRXY(I,J)
    ALBSOILDIRXY(I,:,J) = NoahmpIO%ALBSOILDIRXY(I,:,J)
    ALBSOILDIFXY(I,:,J) = NoahmpIO%ALBSOILDIFXY(I,:,J)
    if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       FSATXY(I,J)      = NoahmpIO%FSATXY(I,J)
       WSURFXY(I,J)     = NoahmpIO%WSURFXY(I,J)
    endif
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       SNRDSXY(I,:,J)         = NoahmpIO%SNRDSXY(I,:,J)
       SNFRXY(I,:,J)          = NoahmpIO%SNFRXY(I,:,J)
       BCPHIXY(I,:,J)         = NoahmpIO%BCPHIXY(I,:,J)
       BCPHOXY(I,:,J)         = NoahmpIO%BCPHOXY(I,:,J)
       OCPHIXY(I,:,J)         = NoahmpIO%OCPHIXY(I,:,J)
       OCPHOXY(I,:,J)         = NoahmpIO%OCPHOXY(I,:,J)
       DUST1XY(I,:,J)         = NoahmpIO%DUST1XY(I,:,J)
       DUST2XY(I,:,J)         = NoahmpIO%DUST2XY(I,:,J)
       DUST3XY(I,:,J)         = NoahmpIO%DUST3XY(I,:,J)
       DUST4XY(I,:,J)         = NoahmpIO%DUST4XY(I,:,J)
       DUST5XY(I,:,J)         = NoahmpIO%DUST5XY(I,:,J)
       MassConcBCPHIXY(I,:,J) = NoahmpIO%MassConcBCPHIXY(I,:,J)
       MassConcBCPHOXY(I,:,J) = NoahmpIO%MassConcBCPHOXY(I,:,J)
       MassConcOCPHIXY(I,:,J) = NoahmpIO%MassConcOCPHIXY(I,:,J)
       MassConcOCPHOXY(I,:,J) = NoahmpIO%MassConcOCPHOXY(I,:,J)
       MassConcDUST1XY(I,:,J) = NoahmpIO%MassConcDUST1XY(I,:,J)
       MassConcDUST2XY(I,:,J) = NoahmpIO%MassConcDUST2XY(I,:,J)
       MassConcDUST3XY(I,:,J) = NoahmpIO%MassConcDUST3XY(I,:,J)
       MassConcDUST4XY(I,:,J) = NoahmpIO%MassConcDUST4XY(I,:,J)
       MassConcDUST5XY(I,:,J) = NoahmpIO%MassConcDUST5XY(I,:,J)
    endif
! #ifdef WRF_HYDRO
!     sfcheadrt(I,J)    = NoahmpIO%sfcheadrt(I,J)
!     INFXSRT(I,J)      = NoahmpIO%INFXSRT(I,J)
!     soldrain(I,J)     = NoahmpIO%soldrain(I,J)
!     qtiledrain(I,J)   = NoahmpIO%qtiledrain(I,J)
!     ZWATBLE2D(I,J)    = NoahmpIO%ZWATBLE2D(I,J)
! #endif

    ! output NoahmpIO variables mapped to HICAR variables
    T2MVXY(I,J)       = NoahmpIO%T2MVXY(I,J)
    T2MBXY(I,J)       = NoahmpIO%T2MBXY(I,J)
    Q2MVXY(I,J)       = NoahmpIO%Q2MVXY(I,J)
    Q2MBXY(I,J)       = NoahmpIO%Q2MBXY(I,J)
    TRADXY(I,J)       = NoahmpIO%TRADXY(I,J)
    NEEXY(I,J)        = NoahmpIO%NEEXY(I,J)
    GPPXY(I,J)        = NoahmpIO%GPPXY(I,J)
    NPPXY(I,J)        = NoahmpIO%NPPXY(I,J)
    FVEGXY(I,J)       = NoahmpIO%FVEGXY(I,J)
    RUNSFXY(I,J)      = NoahmpIO%RUNSFXY(I,J)
    RUNSBXY(I,J)      = NoahmpIO%RUNSBXY(I,J)
    ECANXY(I,J)       = NoahmpIO%ECANXY(I,J)
    EDIRXY(I,J)       = NoahmpIO%EDIRXY(I,J)
    ETRANXY(I,J)      = NoahmpIO%ETRANXY(I,J)
    FSAXY(I,J)        = NoahmpIO%FSAXY(I,J)
    FIRAXY(I,J)       = NoahmpIO%FIRAXY(I,J)
    APARXY(I,J)       = NoahmpIO%APARXY(I,J)
    PSNXY(I,J)        = NoahmpIO%PSNXY(I,J)
    SAVXY(I,J)        = NoahmpIO%SAVXY(I,J)
    SAGXY(I,J)        = NoahmpIO%SAGXY(I,J)
    RSSUNXY(I,J)      = NoahmpIO%RSSUNXY(I,J)
    RSSHAXY(I,J)      = NoahmpIO%RSSHAXY(I,J)
    BGAPXY(I,J)       = NoahmpIO%BGAPXY(I,J)
    WGAPXY(I,J)       = NoahmpIO%WGAPXY(I,J)
    TGVXY(I,J)        = NoahmpIO%TGVXY(I,J)
    TGBXY(I,J)        = NoahmpIO%TGBXY(I,J)
    CHVXY(I,J)        = NoahmpIO%CHVXY(I,J)
    CHBXY(I,J)        = NoahmpIO%CHBXY(I,J)
    SHGXY(I,J)        = NoahmpIO%SHGXY(I,J)
    SHCXY(I,J)        = NoahmpIO%SHCXY(I,J)
    SHBXY(I,J)        = NoahmpIO%SHBXY(I,J)
    EVGXY(I,J)        = NoahmpIO%EVGXY(I,J)
    EVBXY(I,J)        = NoahmpIO%EVBXY(I,J)
    GHVXY(I,J)        = NoahmpIO%GHVXY(I,J)
    GHBXY(I,J)        = NoahmpIO%GHBXY(I,J)
    IRGXY(I,J)        = NoahmpIO%IRGXY(I,J)
    IRCXY(I,J)        = NoahmpIO%IRCXY(I,J)
    IRBXY(I,J)        = NoahmpIO%IRBXY(I,J)
    TRXY(I,J)         = NoahmpIO%TRXY(I,J)
    EVCXY(I,J)        = NoahmpIO%EVCXY(I,J)
    CHLEAFXY(I,J)     = NoahmpIO%CHLEAFXY(I,J)
    CHUCXY(I,J)       = NoahmpIO%CHUCXY(I,J)
    CHV2XY(I,J)       = NoahmpIO%CHV2XY(I,J)
    CHB2XY(I,J)       = NoahmpIO%CHB2XY(I,J)
    ! QINTSXY(I,J)      = NoahmpIO%QINTSXY(I,J)
    ! QINTRXY(I,J)      = NoahmpIO%QINTRXY(I,J)
    ! QDRIPSXY(I,J)     = NoahmpIO%QDRIPSXY(I,J)
    ! QDRIPRXY(I,J)     = NoahmpIO%QDRIPRXY(I,J)
    ! QTHROSXY(I,J)     = NoahmpIO%QTHROSXY(I,J)
    ! QTHRORXY(I,J)     = NoahmpIO%QTHRORXY(I,J)
    ! QSNSUBXY(I,J)     = NoahmpIO%QSNSUBXY(I,J)
    ! QSNFROXY(I,J)     = NoahmpIO%QSNFROXY(I,J)
    ! QSUBCXY(I,J)      = NoahmpIO%QSUBCXY(I,J)
    ! QFROCXY(I,J)      = NoahmpIO%QFROCXY(I,J)
    ! QEVACXY(I,J)      = NoahmpIO%QEVACXY(I,J)
    ! QDEWCXY(I,J)      = NoahmpIO%QDEWCXY(I,J)
    ! QFRZCXY(I,J)      = NoahmpIO%QFRZCXY(I,J)
    ! QMELTCXY(I,J)     = NoahmpIO%QMELTCXY(I,J)
    ! QSNBOTXY(I,J)     = NoahmpIO%QSNBOTXY(I,J)
    ! QMELTXY(I,J)      = NoahmpIO%QMELTXY(I,J)
    ! PONDINGXY(I,J)    = NoahmpIO%PONDINGXY(I,J)
    ! PAHXY(I,J)        = NoahmpIO%PAHXY(I,J)
    ! PAHGXY(I,J)       = NoahmpIO%PAHGXY(I,J)
    ! PAHVXY(I,J)       = NoahmpIO%PAHVXY(I,J)
    ! PAHBXY(I,J)       = NoahmpIO%PAHBXY(I,J)
    ! FPICEXY(I,J)      = NoahmpIO%FPICEXY(I,J)
    ! RAINLSM(I,J)      = NoahmpIO%RAINLSM(I,J)
    ! SNOWLSM(I,J)      = NoahmpIO%SNOWLSM(I,J)
    ! FORCTLSM(I,J)     = NoahmpIO%FORCTLSM(I,J)
    ! FORCQLSM(I,J)     = NoahmpIO%FORCQLSM(I,J)
    ! FORCPLSM(I,J)     = NoahmpIO%FORCPLSM(I,J)
    ! FORCZLSM(I,J)     = NoahmpIO%FORCZLSM(I,J)
    ! FORCWLSM(I,J)     = NoahmpIO%FORCWLSM(I,J)
    ! EFLXBXY(I,J)      = NoahmpIO%EFLXBXY(I,J)
    ! SOILENERGY(I,J)   = NoahmpIO%SOILENERGY(I,J)
    ! SNOWENERGY(I,J)   = NoahmpIO%SNOWENERGY(I,J)
    ! CANHSXY(I,J)      = NoahmpIO%CANHSXY(I,J)

    enddo ! I
    enddo ! J

    !--------- Output variable mapping end ---------

  end subroutine NoahmpHICARmain

end module NoahmpHICARmainMod
