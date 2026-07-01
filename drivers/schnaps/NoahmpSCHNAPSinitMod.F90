module NoahmpSCHNAPSinitMod

! --------------------------------------------------------------------------
! this is for NoahmpIO variable mapping & initialization in WRF physics_init
! and calling the main NoahMP initialization module: NoahmpInitMain(NoahmpIO) 
! adapted from NOAHMP_INIT in original module_sf_noahmpdrv.F file
!
! Coder: Cenlin He (NCAR), December 2025
! ---------------------------------------------------------------------------

contains

  subroutine NoahmpSCHNAPSinit(NoahmpIO, MMINLU, SNOW, SNOWH, CANWAT, ISLTYP, IVGTYP, XLAT, &
                   TSLB,  SMOIS, SH2O,   DZS, FNDSOILW, FNDSNOWH,                       &
                   TSK, isnowxy, tvxy,  tgxy, canicexy,      TMN,   XICE,               &
                   canliqxy,    eahxy, tahxy,     cmxy,     chxy,                       &
                   fwetxy, sneqvoxy, alboldxy, qsnowxy, qrainxy, wslakexy, zwtxy, waxy, &
                   wtxy, tsnoxy, zsnsoxy, snicexy, snliqxy, lfmassxy, rtmassxy,         &
                   stmassxy, woodxy, stblcpxy, fastcpxy, xsaixy, lai,                   &
                   grainxy,   gddxy,                                                    &
                   croptype, cropcat,                                                   &
                   irnumsi, irnummi, irnumfi, irwatsi,                                  &
                   irwatmi, irwatfi, ireloss, irsivol,                                  &
                   irmivol, irfivol, irrsplh,                                           &
                   t2mvxy,   t2mbxy, fsatxy, wsurfxy,                                   &
                   snrdsxy, snfrxy, bcphixy, bcphoxy, ocphixy, ocphoxy, dust1xy,        &
                   dust2xy, dust3xy, dust4xy, dust5xy, massconcbcphixy, massconcbcphoxy,&
                   massconcocphixy, massconcocphoxy, massconcdust1xy, massconcdust2xy,  &
                   massconcdust3xy, massconcdust4xy, massconcdust5xy,                   &
                   ALBSOILDIRXY, ALBSOILDIFXY,                                          &
                   NSOIL,   NSNOW, restart, allowed_to_read, XICE_THRES, DX,            &
                   IDVEG,   IOPT_CRS, IOPT_BTR, IOPT_RUNSUB,IOPT_SFC, IOPT_FRZ,  & ! IN : User options
                   IOPT_INF,IOPT_RAD,   IOPT_ALB, IOPT_SNF,IOPT_TBOT, IOPT_STC,  & ! IN : User options
                   IOPT_GLA,IOPT_RSF,  IOPT_SOIL,IOPT_PEDO,IOPT_CROP, IOPT_IRR,  & ! IN : User options
                   IOPT_IRRM,IOPT_INFDV,IOPT_TDRN, soiltstep,                    & ! IN : User options
                   IOPT_RUNSRF, IOPT_TKSNO, IOPT_COMPACT, IOPT_SCF, IOPT_WETLAND, & ! IN : User options                 
                   IZ0TLND, SF_URBAN_PHYSICS,                                    & ! IN : User options
                   SNICAR_BANDNUMBER_OPT, SNICAR_SOLARSPEC_OPT,                  & ! SNICAR variable
                   SNICAR_SNOWOPTICS_OPT, SNICAR_DUSTOPTICS_OPT,                 & ! SNICAR variable
                   SNICAR_RTSOLVER_OPT, SNICAR_SNOWSHAPE_OPT,                    & ! SNICAR variable
                   SNICAR_USE_AEROSOL, SNICAR_SNOWBC_INTMIX,                     & ! SNICAR variable
                   SNICAR_SNOWDUST_INTMIX, SNICAR_USE_OC,                        & ! SNICAR variable
                   SNICAR_AEROSOL_READTABLE,                                            &
                   ids,ide, jds,jde, kds,kde,                                           &
                   ims,ime, jms,jme, kms,kme,                                           &
                   its,ite, jts,jte, kts,kte,                                           &
                   smoiseq,smcwtdxy,rechxy,deeprechxy,qtdrain,areaxy,dy,msftx,msfty, & ! Optional groundwater
                   wtddt,   stepwtd, dt, qrfsxy, qspringsxy, qslatxy,                   & ! Optional groundwater
                   fdepthxy, ht, riverbedxy, eqzwt, rivercondxy, pexpxy, rechclim       ) ! Optional groundwater

! ---------------------------------------------------------------------------------------

    use NoahmpIOVarType, only : NoahmpIO_type
    use NoahmpIOVarInitMod
    use NoahmpReadTableMod
    use NoahmpInitMainMod
    use LanduseConvertMod
    use SnowInputSnicarMod
    use NoahmpDriverMainMod, only : NoahmpDriverInit, noahmp_initialized

    implicit none

    type(NoahmpIO_type), intent(inout)                         :: NoahmpIO

    ! input only
    INTEGER, INTENT(IN)                                        :: ids,ide, jds,jde, kds,kde,  &
                                                                  ims,ime, jms,jme, kms,kme,  &
                                                                  its,ite, jts,jte, kts,kte
    INTEGER, INTENT(IN)                                        :: NSOIL, NSNOW
    INTEGER, INTENT(IN)                                        ::  IOPT_CRS     ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
    INTEGER, INTENT(IN)                                        ::  IOPT_BTR     ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
    INTEGER, INTENT(IN)                                        ::  IOPT_RUNSUB  ! subsurface runoff and groundwater (currently keep the same as surface runoff option)
    INTEGER, INTENT(IN)                                        ::  IOPT_RUNSRF  ! surface runoff (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS; 5->MMF; 6->VIC; 7->XianAnJiang; 8->DynVIC)
    INTEGER, INTENT(IN)                                        ::  IOPT_COMPACT ! snowpack compaction (1->Anderson1976; 2->Abolafia-Rosenzweig2024)
    INTEGER, INTENT(IN)                                        ::  IOPT_TKSNO   ! snow thermal conductivity: 1 -> Stieglitz(yen,1965) scheme (default), 2 -> Anderson, 1976 scheme, 3 -> constant, 4 -> Verseghy (1991) scheme, 5 -> Douvill(Yen, 1981) scheme
    INTEGER, INTENT(IN)                                        ::  IOPT_SCF     ! snow cover fraction (1->NiuYang07; 2->Abolafia-Rosenzweig2025)
    INTEGER, INTENT(IN)                                        ::  IOPT_WETLAND ! wetland model option (0->off; 1->Zhang2022 fixed parameter; 2->Zhang2022 read in 2D parameter)
    INTEGER, INTENT(IN)                                        ::  IOPT_SFC     ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
    INTEGER, INTENT(IN)                                        ::  IOPT_FRZ     ! supercooled liquid water (1-> NY06; 2->Koren99)
    INTEGER, INTENT(IN)                                        ::  IOPT_INF     ! frozen soil permeability (1-> NY06; 2->Koren99)
    INTEGER, INTENT(IN)                                        ::  IOPT_RAD     ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
    INTEGER, INTENT(IN)                                        ::  IOPT_ALB     ! snow surface albedo (1->BATS; 2->CLASS; 3->SNICAR)
    INTEGER, INTENT(IN)                                        ::  IOPT_SNF     ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
    INTEGER, INTENT(IN)                                        ::  IOPT_TBOT    ! lower boundary of soil temperature (1->zero-flux; 2->Noah)
    INTEGER, INTENT(IN)                                        ::  IOPT_STC     ! snow/soil temperature time scheme
    INTEGER, INTENT(IN)                                        ::  IOPT_GLA     ! glacier option (1->phase change; 2->simple)
    INTEGER, INTENT(IN)                                        ::  IOPT_RSF     ! surface resistance (1->Sakaguchi/Zeng; 2->Seller; 3->mod Sellers; 4->1+snow)
    INTEGER, INTENT(IN)                                        ::  IOPT_SOIL    ! soil configuration option
    INTEGER, INTENT(IN)                                        ::  IOPT_PEDO    ! soil pedotransfer function option
    INTEGER, INTENT(IN)                                        ::  IOPT_CROP    ! crop model option (0->none; 1->Liu et al.; 2->Gecros)
    INTEGER, INTENT(IN)                                        ::  IOPT_IRR     ! irrigation scheme (0->none; >1 irrigation scheme ON)
    INTEGER, INTENT(IN)                                        ::  IOPT_IRRM    ! irrigation method
    INTEGER, INTENT(IN)                                        ::  IOPT_INFDV   ! infiltration options for dynamic VIC infiltration (1->Philip; 2-> Green-Ampt;3->Smith-Parlange)
    INTEGER, INTENT(IN)                                        ::  IOPT_TDRN    ! tile drainage (0-> no tile drainage; 1-> simple tile drainage;2->Hooghoudt's)
    REAL,                                            INTENT(IN   ) ::  XICE_THRES   ! fraction of grid determining seaice
    INTEGER,                                         INTENT(IN   ) ::  IDVEG        ! dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1      
    REAL,                                            INTENT(IN   ) ::  soiltstep    ! soil timestep (s), default:0->same as main model timestep
    INTEGER,                                         INTENT(IN   ) ::  IZ0TLND      ! option of Chen adjustment of Czil (not used)
    INTEGER,                                         INTENT(IN   ) ::  sf_urban_physics ! urban physics option
    REAL,                                            INTENT(IN   ) ::  DX           ! horizontal grid spacing [m]


    LOGICAL, INTENT(IN)                                        :: restart, allowed_to_read
    LOGICAL, INTENT(IN)                                        :: FNDSOILW, FNDSNOWH
    CHARACTER(LEN=*),                    INTENT(IN)            :: MMINLU
    REAL,    DIMENSION(NSOIL), INTENT(IN)                      :: DZS                 ! Thickness of the soil layers [m]
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: ISLTYP, IVGTYP
    REAL   , DIMENSION(ims:ime,5,jms:jme),INTENT(IN)           :: croptype            ! crop type fraction
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: XLAT                ! latitude
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: TSK                 ! skin temperature (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: XICE                ! sea ice fraction
    REAL,                                INTENT(IN), OPTIONAL  :: DT, WTDDT
    REAL,                                INTENT(IN), OPTIONAL  :: DY
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: FDEPTHXY            ! efolding depth for transmissivity (m)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: HT                  ! terrain height (m)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: MSFTX, MSFTY
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: rechclim

    ! SNICAR snow albedo variables
    INTEGER, INTENT(IN),                             OPTIONAL  :: SNICAR_BANDNUMBER_OPT, &
                                                                  SNICAR_SOLARSPEC_OPT, &
                                                                  SNICAR_SNOWOPTICS_OPT,  &
                                                                  SNICAR_DUSTOPTICS_OPT,  &
                                                                  SNICAR_RTSOLVER_OPT,    &
                                                                  SNICAR_SNOWSHAPE_OPT
    LOGICAL, INTENT(IN),                             OPTIONAL  :: SNICAR_USE_AEROSOL,     &
                                                                  SNICAR_SNOWBC_INTMIX,   &
                                                                  SNICAR_SNOWDUST_INTMIX, &
                                                                  SNICAR_USE_OC,          &
                                                                  SNICAR_AEROSOL_READTABLE

    ! in/out
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: TMN                 ! deep soil temperature (k)
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: isnowxy             ! actual no. of snow layers
    REAL,    DIMENSION(ims:ime,1:NSOIL,jms:jme), INTENT(INOUT) :: SMOIS, SH2O, TSLB
    REAL,    DIMENSION(ims:ime, jms:jme), INTENT(INOUT)        :: SNOW, SNOWH, CANWAT
    REAL,    DIMENSION(ims:ime,1:(NSNOW+NSOIL),jms:jme),INTENT(INOUT) :: zsnsoxy             ! snow layer depth [m]
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: tsnoxy              ! snow temperature [K]
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: snicexy             ! snow layer ice [mm]
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: snliqxy             ! snow layer liquid water [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: tvxy                ! vegetation canopy temperature
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: tgxy                ! ground surface temperature
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: canicexy            ! canopy-intercepted ice (mm)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: canliqxy            ! canopy-intercepted liquid water (mm)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: eahxy               ! canopy air vapor pressure (pa)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: tahxy               ! canopy air temperature (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: cmxy                ! momentum drag coefficient
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: chxy                ! sensible heat exchange coefficient
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: fwetxy              ! wetted or snowed fraction of the canopy (-)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: sneqvoxy            ! snow mass at last time step(mm h2o)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: alboldxy            ! snow albedo at last time step (-)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: qsnowxy             ! snowfall on the ground [mm/s]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: qrainxy             ! rainfall on the ground [mm/s]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wslakexy            ! lake water storage [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: zwtxy               ! water table depth [m]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: waxy                ! water in the "aquifer" [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wtxy                ! groundwater storage [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: lfmassxy            ! leaf mass [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: rtmassxy            ! mass of fine roots [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: stmassxy            ! stem mass [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: woodxy              ! mass of wood (incl. woody roots) [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: grainxy             ! mass of grain [g/m2] !XING
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: gddxy               ! growing degree days !XING
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: stblcpxy            ! stable carbon in deep soil [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: fastcpxy            ! short-lived carbon, shallow soil [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: xsaixy              ! stem area index
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: lai                 ! leaf area index
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irnumsi             ! irrigation number
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irnummi             ! irrigation number
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irnumfi             ! irrigation number
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irwatsi             ! irrigation water amount
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irwatmi             ! irrigation water amount
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irwatfi             ! irrigation water amount
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: ireloss             ! irrigation loss
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irsivol             ! irrigation water volume
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irmivol             ! irrigation water volume
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irfivol             ! irrigation water volume
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irrsplh             ! irrigation evaporation heat
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: t2mvxy              ! 2m temperature vegetation part (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: t2mbxy              ! 2m temperature bare ground part (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: fsatxy              ! saturation fraction of grid (-)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wsurfxy             ! wetland water storage (mm)
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: snrdsxy             ! SNICAR snow radius
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: snfrxy              ! SNICAR snow freezing rate
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: bcphixy             ! SNICAR BCPHI mass in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: bcphoxy             ! SNICAR BCPHO mass in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: ocphixy             ! SNICAR OCPHI mass in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: ocphoxy             ! SNICAR OCPHO mass in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: dust1xy             ! SNICAR DUST1 mass in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: dust2xy             ! SNICAR DUST2 mass in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: dust3xy             ! SNICAR DUST3 mass in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: dust4xy             ! SNICAR DUST4 mass in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: dust5xy             ! SNICAR DUST5 mass in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: massconcbcphixy     ! SNICAR BCPHI mass conc in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: massconcbcphoxy     ! SNICAR BCPHO mass conc in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: massconcocphixy     ! SNICAR OCPHI mass conc in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: massconcocphoxy     ! SNICAR OCPHO mass conc in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: massconcdust1xy     ! SNICAR DUST1 mass conc in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: massconcdust2xy     ! SNICAR DUST2 mass conc in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: massconcdust3xy     ! SNICAR DUST3 mass conc in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: massconcdust4xy     ! SNICAR DUST4 mass conc in snow
    REAL,    DIMENSION(ims:ime,1:NSNOW,jms:jme), INTENT(INOUT)    :: massconcdust5xy     ! SNICAR DUST5 mass conc in snow
    REAL,    DIMENSION(ims:ime,1:2,jms:jme),  INTENT(INOUT)    :: ALBSOILDIRXY        ! soil albedo direct
    REAL,    DIMENSION(ims:ime,1:2,jms:jme),  INTENT(INOUT)    :: ALBSOILDIFXY        ! soil albedo diffuse
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: qtdrain      ! tile drainage (mm)
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: smcwtdxy     ! deep soil moisture content [m3m-3]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: deeprechxy   ! deep recharge [m]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: rechxy       ! accumulated recharge [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: qrfsxy       ! accumulated flux from groundwater to rivers [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: qspringsxy   ! accumulated seeping water [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: qslatxy      ! accumulated lateral flow [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: areaxy       ! grid cell area [m2]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: RIVERBEDXY   ! riverbed depth (m)
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: EQZWT        ! equilibrium water table depth (m)
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: RIVERCONDXY  ! river conductance
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: PEXPXY       ! factor for river conductance
    REAL,    DIMENSION(ims:ime,1:NSOIL,jms:jme), INTENT(INOUT), OPTIONAL :: smoiseq      ! equilibrium soil moisture content [m3m-3]
    ! output only
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)           :: cropcat             ! crop type
    INTEGER,                             INTENT(OUT), OPTIONAL :: STEPWTD
    ! local
    integer :: itf, jtf, I, J
    integer :: LoopInd    ! loop index for array section expansion
! ----------------------------------------------------------------------------------

    ! initialize NoahmpIO dimension and key config variables
    NoahmpIO%ids                = ids
    NoahmpIO%ide                = ide
    NoahmpIO%jds                = jds
    NoahmpIO%jde                = jde
    NoahmpIO%kds                = kds
    NoahmpIO%kde                = kde
    NoahmpIO%ims                = ims
    NoahmpIO%ime                = ime
    NoahmpIO%jms                = jms
    NoahmpIO%jme                = jme
    NoahmpIO%kms                = kms
    NoahmpIO%kme                = kme
    NoahmpIO%its                = its
    NoahmpIO%ite                = ite
    NoahmpIO%jts                = jts
    NoahmpIO%jte                = jte
    NoahmpIO%kts                = kts
    NoahmpIO%kte                = kte
    NoahmpIO%xstart             = ims
    NoahmpIO%xend               = ime
    NoahmpIO%ystart             = jms
    NoahmpIO%yend               = jme    
    NoahmpIO%NSOIL              = NSOIL
    NoahmpIO%DX                 = DX
    NoahmpIO%DY                 = DX
    NoahmpIO%IOPT_DVEG          = IDVEG
    NoahmpIO%IOPT_CRS           = IOPT_CRS
    NoahmpIO%IOPT_BTR           = IOPT_BTR
    NoahmpIO%IOPT_SFC           = IOPT_SFC
    NoahmpIO%IOPT_FRZ           = IOPT_FRZ
    NoahmpIO%IOPT_INF           = IOPT_INF
    NoahmpIO%IOPT_RAD           = IOPT_RAD
    NoahmpIO%IOPT_ALB           = IOPT_ALB
    NoahmpIO%IOPT_SNF           = IOPT_SNF
    NoahmpIO%IOPT_TBOT          = IOPT_TBOT
    NoahmpIO%IOPT_STC           = IOPT_STC
    NoahmpIO%IOPT_GLA           = IOPT_GLA
    NoahmpIO%IOPT_RSF           = IOPT_RSF
    NoahmpIO%IOPT_SOIL          = IOPT_SOIL
    NoahmpIO%IOPT_PEDO          = IOPT_PEDO
    NoahmpIO%IOPT_CROP          = IOPT_CROP
    NoahmpIO%IOPT_IRR           = IOPT_IRR
    NoahmpIO%IOPT_IRRM          = IOPT_IRRM
    NoahmpIO%IOPT_INFDV         = IOPT_INFDV
    NoahmpIO%IOPT_TDRN          = IOPT_TDRN
    NoahmpIO%IOPT_RUNSRF        = IOPT_RUNSRF
    NoahmpIO%IOPT_RUNSUB        = IOPT_RUNSUB
    NoahmpIO%IOPT_TKSNO         = IOPT_TKSNO
    NoahmpIO%IOPT_COMPACT       = IOPT_COMPACT
    NoahmpIO%IOPT_SCF           = IOPT_SCF
    NoahmpIO%IOPT_WETLAND       = IOPT_WETLAND
    NoahmpIO%SOILTSTEP          = SOILTSTEP
    NoahmpIO%SF_URBAN_PHYSICS   = SF_URBAN_PHYSICS
    NoahmpIO%IZ0TLND            = IZ0TLND
    NoahmpIO%LLANDUSE           = LanduseConvert(MMINLU)
    NoahmpIO%XICE_THRESHOLD     = XICE_THRES
    NoahmpIO%DZS                = DZS
    ! NoahmpIO%IRI_URBAN          = IRI_SCHEME
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       NoahmpIO%SNICAR_BANDNUMBER_OPT    = SNICAR_BANDNUMBER_OPT
       NoahmpIO%SNICAR_SOLARSPEC_OPT     = SNICAR_SOLARSPEC_OPT
       NoahmpIO%SNICAR_SNOWOPTICS_OPT    = SNICAR_SNOWOPTICS_OPT
       NoahmpIO%SNICAR_DUSTOPTICS_OPT    = SNICAR_DUSTOPTICS_OPT
       NoahmpIO%SNICAR_RTSOLVER_OPT      = SNICAR_RTSOLVER_OPT
       NoahmpIO%SNICAR_SNOWSHAPE_OPT     = SNICAR_SNOWSHAPE_OPT
       NoahmpIO%SNICAR_USE_AEROSOL       = SNICAR_USE_AEROSOL
       NoahmpIO%SNICAR_SNOWBC_INTMIX     = SNICAR_SNOWBC_INTMIX
       NoahmpIO%SNICAR_SNOWDUST_INTMIX   = SNICAR_SNOWDUST_INTMIX
       NoahmpIO%SNICAR_USE_OC            = SNICAR_USE_OC
       NoahmpIO%SNICAR_AEROSOL_READTABLE = SNICAR_AEROSOL_READTABLE
    endif

    ! $acc update device(NoahmpIO%ids, NoahmpIO%ide, NoahmpIO%jds, NoahmpIO%jde, NoahmpIO%kds, NoahmpIO%kde) &
    ! $acc        device(NoahmpIO%ims, NoahmpIO%ime, NoahmpIO%jms, NoahmpIO%jme, NoahmpIO%kms, NoahmpIO%kme) &
    ! $acc        device(NoahmpIO%its, NoahmpIO%ite, NoahmpIO%jts, NoahmpIO%jte, NoahmpIO%kts, NoahmpIO%kte) &
    ! $acc        device(NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend) &
    ! $acc        device(NoahmpIO%YR, NoahmpIO%JULIAN, NoahmpIO%DTBL, NoahmpIO%NSOIL, NoahmpIO%DX, NoahmpIO%DY) &
    ! $acc        device(NoahmpIO%IOPT_DVEG, NoahmpIO%IOPT_CRS, NoahmpIO%IOPT_BTR, NoahmpIO%IOPT_SFC) &
    ! $acc        device(NoahmpIO%IOPT_FRZ, NoahmpIO%IOPT_INF, NoahmpIO%IOPT_RAD, NoahmpIO%IOPT_ALB) &
    ! $acc        device(NoahmpIO%IOPT_SNF, NoahmpIO%IOPT_TBOT, NoahmpIO%IOPT_STC, NoahmpIO%IOPT_GLA) &
    ! $acc        device(NoahmpIO%IOPT_RSF, NoahmpIO%IOPT_SOIL, NoahmpIO%IOPT_PEDO, NoahmpIO%IOPT_CROP) &
    ! $acc        device(NoahmpIO%IOPT_IRR, NoahmpIO%IOPT_IRRM, NoahmpIO%IOPT_INFDV, NoahmpIO%IOPT_TDRN) &
    ! $acc        device(NoahmpIO%IOPT_RUNSRF, NoahmpIO%IOPT_RUNSUB, NoahmpIO%IOPT_TKSNO, NoahmpIO%IOPT_COMPACT) &
    ! $acc        device(NoahmpIO%IOPT_SCF, NoahmpIO%IOPT_WETLAND, NoahmpIO%SF_URBAN_PHYSICS, NoahmpIO%IZ0TLND) &
    ! $acc        device(NoahmpIO%LLANDUSE, NoahmpIO%SOILTSTEP, NoahmpIO%XICE_THRESHOLD, NoahmpIO%DZS) &
    ! $acc        device(NoahmpIO%ITIMESTEP) &
    ! $acc        device(NoahmpIO%SNICAR_BANDNUMBER_OPT, NoahmpIO%SNICAR_SOLARSPEC_OPT) &
    ! $acc        device(NoahmpIO%SNICAR_SNOWOPTICS_OPT, NoahmpIO%SNICAR_DUSTOPTICS_OPT) &
    ! $acc        device(NoahmpIO%SNICAR_RTSOLVER_OPT, NoahmpIO%SNICAR_SNOWSHAPE_OPT) &
    ! $acc        device(NoahmpIO%SNICAR_USE_AEROSOL, NoahmpIO%SNICAR_SNOWBC_INTMIX) &
    ! $acc        device(NoahmpIO%SNICAR_SNOWDUST_INTMIX, NoahmpIO%SNICAR_USE_OC) &
    ! $acc        device(NoahmpIO%SNICAR_AEROSOL_READTABLE)

    ! initialze all NoahmpIO variables with default values
    call NoahmpIOVarInitDefault(NoahmpIO)

    ! read in Noahmp table parameters
    call NoahmpReadTable(NoahmpIO)

    ! read in SNICAR parameter netcdif file
    if ( NoahmpIO%IOPT_ALB == 3 ) then
      call SnowInputSnicar(NoahmpIO, snicar_optic_flnm="snicar_optics_5bnd_c013122.nc", snicar_age_flnm="snicar_drdt_bst_fit_60_c070416.nc")
    endif

    !--------- WRF variables mapped to NoahmpIO variables

    ! non-2D variables
    itf = min0(ite, ide-1)
    jtf = min0(jte, jde-1)
    NoahmpIO%DZS                = DZS
    NoahmpIO%FNDSNOWH           = FNDSNOWH
    NoahmpIO%restart_flag       = restart
    NoahmpIO%DX                 = DX
    if(present(DT)) NoahmpIO%DTBL               = DT
    if(present(WTDDT)) NoahmpIO%WTDDT              = WTDDT
    if(present(DY)) NoahmpIO%DY                 = DY

    ! If re-initializing, delete old noahmpIO device data first
    if (noahmp_initialized) then
       !$acc exit data delete(noahmpIO)
    endif
    !$acc enter data copyin(noahmpIO)

    ! 2D/3D variables
    !$acc parallel loop gang vector collapse(2) default(present) private(LoopInd) firstprivate(NSNOW, NSOIL)
    do J = jts, jtf
    do I = its, itf
    
    ! input only
    NoahmpIO%IVGTYP(I,J)               = IVGTYP(I,J)
    NoahmpIO%ISLTYP(I,J)               = ISLTYP(I,J)
    NoahmpIO%XLAT(I,J)                 = XLAT(I,J)
    NoahmpIO%TSK(I,J)                  = TSK(I,J)
    NoahmpIO%XICE(I,J)                 = XICE(I,J)
    !$acc loop seq
    do LoopInd = 1, 5
       NoahmpIO%CROPTYPE(I,LoopInd,J)      = CROPTYPE(I,LoopInd,J)
    enddo
    ! in/out variables
    !$acc loop seq
    do LoopInd = 1, NSOIL
       NoahmpIO%SMOIS(I,LoopInd,J) = SMOIS(I,LoopInd,J)
       NoahmpIO%SH2O(I,LoopInd,J) = SH2O(I,LoopInd,J)
       NoahmpIO%TSLB(I,LoopInd,J) = TSLB(I,LoopInd,J)
    enddo
    NoahmpIO%SNOW(I,J)                 = SNOW(I,J)    
    NoahmpIO%SNOWH(I,J)                = SNOWH(I,J)   
    NoahmpIO%CANWAT(I,J)               = CANWAT(I,J)  
    NoahmpIO%CANICEXY(I,J)             = CANICEXY(I,J)
    NoahmpIO%CANLIQXY(I,J)             = CANLIQXY(I,J)
    NoahmpIO%TMN(I,J)                  = TMN(I,J)
    NoahmpIO%ISNOWXY(I,J)              = ISNOWXY(I,J)
    !$acc loop seq
    do LoopInd = -NSNOW+1, NSOIL
       NoahmpIO%ZSNSOXY(I,LoopInd,J)      = ZSNSOXY(I,LoopInd+NSNOW,J)
    enddo 
    !$acc loop seq
    do LoopInd = -NSNOW+1, 0
       NoahmpIO%TSNOXY(I,LoopInd,J) = TSNOXY(I,LoopInd+NSNOW,J)
       NoahmpIO%SNICEXY(I,LoopInd,J) = SNICEXY(I,LoopInd+NSNOW,J)
       NoahmpIO%SNLIQXY(I,LoopInd,J) = SNLIQXY(I,LoopInd+NSNOW,J)
    enddo
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
    NoahmpIO%LFMASSXY(I,J)             = LFMASSXY(I,J)
    NoahmpIO%RTMASSXY(I,J)             = RTMASSXY(I,J)
    NoahmpIO%STMASSXY(I,J)             = STMASSXY(I,J)
    NoahmpIO%WOODXY(I,J)               = WOODXY(I,J)
    NoahmpIO%GRAINXY(I,J)              = GRAINXY(I,J)
    NoahmpIO%GDDXY(I,J)                = GDDXY(I,J)
    NoahmpIO%STBLCPXY(I,J)             = STBLCPXY(I,J)
    NoahmpIO%FASTCPXY(I,J)             = FASTCPXY(I,J)
    NoahmpIO%LAI(I,J)                  = LAI(I,J)
    NoahmpIO%XSAIXY(I,J)               = XSAIXY(I,J)
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
    NoahmpIO%T2MVXY(I,J)               = T2MVXY(I,J)
    NoahmpIO%T2MBXY(I,J)               = T2MBXY(I,J)
    !$acc loop seq
    do LoopInd = 1, 2
       NoahmpIO%ALBSOILDIRXY(I,LoopInd,J) = ALBSOILDIRXY(I,LoopInd,J)
       NoahmpIO%ALBSOILDIFXY(I,LoopInd,J) = ALBSOILDIFXY(I,LoopInd,J)
    enddo
    if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       NoahmpIO%FSATXY(I,J)            = FSATXY(I,J)
       NoahmpIO%WSURFXY(I,J)           = WSURFXY(I,J)
    endif
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       !$acc loop seq
       do LoopInd = -NSNOW+1, 0
          NoahmpIO%SNRDSXY(I,LoopInd,J) = SNRDSXY(I,LoopInd+NSNOW,J)
          NoahmpIO%SNFRXY(I,LoopInd,J) = SNFRXY(I,LoopInd+NSNOW,J)
          NoahmpIO%BCPHIXY(I,LoopInd,J) = BCPHIXY(I,LoopInd+NSNOW,J)
          NoahmpIO%BCPHOXY(I,LoopInd,J) = BCPHOXY(I,LoopInd+NSNOW,J)
          NoahmpIO%OCPHIXY(I,LoopInd,J) = OCPHIXY(I,LoopInd+NSNOW,J)
          NoahmpIO%OCPHOXY(I,LoopInd,J) = OCPHOXY(I,LoopInd+NSNOW,J)
          NoahmpIO%DUST1XY(I,LoopInd,J) = DUST1XY(I,LoopInd+NSNOW,J)
          NoahmpIO%DUST2XY(I,LoopInd,J) = DUST2XY(I,LoopInd+NSNOW,J)
          NoahmpIO%DUST3XY(I,LoopInd,J) = DUST3XY(I,LoopInd+NSNOW,J)
          NoahmpIO%DUST4XY(I,LoopInd,J) = DUST4XY(I,LoopInd+NSNOW,J)
          NoahmpIO%DUST5XY(I,LoopInd,J) = DUST5XY(I,LoopInd+NSNOW,J)
          NoahmpIO%MassConcBCPHIXY(I,LoopInd,J) = MassConcBCPHIXY(I,LoopInd+NSNOW,J)
          NoahmpIO%MassConcBCPHOXY(I,LoopInd,J) = MassConcBCPHOXY(I,LoopInd+NSNOW,J)
          NoahmpIO%MassConcOCPHIXY(I,LoopInd,J) = MassConcOCPHIXY(I,LoopInd+NSNOW,J)
          NoahmpIO%MassConcOCPHOXY(I,LoopInd,J) = MassConcOCPHOXY(I,LoopInd+NSNOW,J)
          NoahmpIO%MassConcDUST1XY(I,LoopInd,J) = MassConcDUST1XY(I,LoopInd+NSNOW,J)
          NoahmpIO%MassConcDUST2XY(I,LoopInd,J) = MassConcDUST2XY(I,LoopInd+NSNOW,J)
          NoahmpIO%MassConcDUST3XY(I,LoopInd,J) = MassConcDUST3XY(I,LoopInd+NSNOW,J)
          NoahmpIO%MassConcDUST4XY(I,LoopInd,J) = MassConcDUST4XY(I,LoopInd+NSNOW,J)
          NoahmpIO%MassConcDUST5XY(I,LoopInd,J) = MassConcDUST5XY(I,LoopInd+NSNOW,J)
       enddo
    endif

    NoahmpIO%CROPCAT(I,J) = CROPCAT(I,J)

    enddo ! I
    enddo ! J

    ! Optional argument copies - must be outside GPU parallel region
    ! because if(present(...)) cannot safely execute on the GPU
    if(present(FDEPTHXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%FDEPTHXY(I,J) = FDEPTHXY(I,J)
      enddo; enddo
    endif
    if(present(MSFTX)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%MSFTX(I,J) = MSFTX(I,J)
      enddo; enddo
    endif
    if(present(MSFTY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%MSFTY(I,J) = MSFTY(I,J)
      enddo; enddo
    endif
    if(present(HT)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%TERRAIN(I,J) = HT(I,J)
      enddo; enddo
    endif
    if(present(RECHCLIM)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%RECHCLIM(I,J) = RECHCLIM(I,J)
      enddo; enddo
    endif
    if(present(QTDRAIN)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%QTDRAIN(I,J) = QTDRAIN(I,J)
      enddo; enddo
    endif
    if(present(SMCWTDXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%SMCWTDXY(I,J) = SMCWTDXY(I,J)
      enddo; enddo
    endif
    if(present(DEEPRECHXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%DEEPRECHXY(I,J) = DEEPRECHXY(I,J)
      enddo; enddo
    endif
    if(present(RECHXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%RECHXY(I,J) = RECHXY(I,J)
      enddo; enddo
    endif
    if(present(QRFSXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%QRFSXY(I,J) = QRFSXY(I,J)
      enddo; enddo
    endif
    if(present(QSPRINGSXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%QSPRINGSXY(I,J) = QSPRINGSXY(I,J)
      enddo; enddo
    endif
    if(present(QSLATXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%QSLATXY(I,J) = QSLATXY(I,J)
      enddo; enddo
    endif
    if(present(AREAXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%AREAXY(I,J) = AREAXY(I,J)
      enddo; enddo
    endif
    if(present(RIVERBEDXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%RIVERBEDXY(I,J) = RIVERBEDXY(I,J)
      enddo; enddo
    endif
    if(present(EQZWT)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%EQZWT(I,J) = EQZWT(I,J)
      enddo; enddo
    endif
    if(present(RIVERCONDXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%RIVERCONDXY(I,J) = RIVERCONDXY(I,J)
      enddo; enddo
    endif
    if(present(PEXPXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        NoahmpIO%PEXPXY(I,J) = PEXPXY(I,J)
      enddo; enddo
    endif
    if(present(SMOISEQ)) then
      !$acc parallel loop gang vector collapse(2) default(present) firstprivate(LoopInd, NSOIL)
      do J = jts, jtf; do I = its, itf
        !$acc loop seq
        do LoopInd = 1, NSOIL
           NoahmpIO%SMOISEQ(I,LoopInd,J) = SMOISEQ(I,LoopInd,J)
        enddo
      enddo; enddo
    endif

    !--------- WRF -> NoahmpIO variables mapping ends

    !--------- initialize noahmp data types and create device arrays (batched async)
    call NoahmpDriverInit(NoahmpIO)
    !---------

    !--------- main Noahmp initialization module
    call NoahmpInitMain(NoahmpIO)
    !---------

    !--------- initialized NoahmpIO variable mapped to WRF variables
    !$acc parallel loop gang vector collapse(2) default(present) firstprivate(LoopInd, NSNOW, NSOIL)
    do J = jts, jtf
    do I = its, itf

    ! in/out variables
    !$acc loop seq
    do LoopInd = 1, NSOIL
       SMOIS(I,LoopInd,J) = NoahmpIO%SMOIS(I,LoopInd,J)
       SH2O(I,LoopInd,J) = NoahmpIO%SH2O(I,LoopInd,J)
       TSLB(I,LoopInd,J) = NoahmpIO%TSLB(I,LoopInd,J)
    enddo
    SNOW(I,J)           = NoahmpIO%SNOW(I,J)
    SNOWH(I,J)          = NoahmpIO%SNOWH(I,J) 
    CANWAT(I,J)         = NoahmpIO%CANWAT(I,J)
    CANICEXY(I,J)       = NoahmpIO%CANICEXY(I,J)
    CANLIQXY(I,J)       = NoahmpIO%CANLIQXY(I,J)
    TMN(I,J)            = NoahmpIO%TMN(I,J)
    ISNOWXY(I,J)        = NoahmpIO%ISNOWXY(I,J)
    !$acc loop seq
    do LoopInd = 1, NSNOW+NSOIL
       ZSNSOXY(I,LoopInd,J) = NoahmpIO%ZSNSOXY(I,LoopInd-NSNOW,J)
    enddo
    !$acc loop seq
    do LoopInd = 1, NSNOW
       TSNOXY(I,LoopInd,J) = NoahmpIO%TSNOXY(I,LoopInd-NSNOW,J)
       SNICEXY(I,LoopInd,J) = NoahmpIO%SNICEXY(I,LoopInd-NSNOW,J)
       SNLIQXY(I,LoopInd,J) = NoahmpIO%SNLIQXY(I,LoopInd-NSNOW,J)
    enddo
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
    LFMASSXY(I,J)       = NoahmpIO%LFMASSXY(I,J)
    RTMASSXY(I,J)       = NoahmpIO%RTMASSXY(I,J)
    STMASSXY(I,J)       = NoahmpIO%STMASSXY(I,J)
    WOODXY(I,J)         = NoahmpIO%WOODXY(I,J)
    GRAINXY(I,J)        = NoahmpIO%GRAINXY(I,J)
    GDDXY(I,J)          = NoahmpIO%GDDXY(I,J)
    STBLCPXY(I,J)       = NoahmpIO%STBLCPXY(I,J)
    FASTCPXY(I,J)       = NoahmpIO%FASTCPXY(I,J)
    LAI(I,J)            = NoahmpIO%LAI(I,J)
    XSAIXY(I,J)         = NoahmpIO%XSAIXY(I,J)
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
    T2MVXY(I,J)         = NoahmpIO%T2MVXY(I,J)
    T2MBXY(I,J)         = NoahmpIO%T2MBXY(I,J)
    !$acc loop seq
    do LoopInd = 1, 2
       ALBSOILDIRXY(I,LoopInd,J) = NoahmpIO%ALBSOILDIRXY(I,LoopInd,J)
       ALBSOILDIFXY(I,LoopInd,J) = NoahmpIO%ALBSOILDIFXY(I,LoopInd,J)
    enddo
    if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       FSATXY(I,J)      = NoahmpIO%FSATXY(I,J)
       WSURFXY(I,J)     = NoahmpIO%WSURFXY(I,J)
    endif
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       !$acc loop seq
       do LoopInd = 1, NSNOW
          SNRDSXY(I,LoopInd,J) = NoahmpIO%SNRDSXY(I,LoopInd-NSNOW,J)
          SNFRXY(I,LoopInd,J) = NoahmpIO%SNFRXY(I,LoopInd-NSNOW,J)
          BCPHIXY(I,LoopInd,J) = NoahmpIO%BCPHIXY(I,LoopInd-NSNOW,J)
          BCPHOXY(I,LoopInd,J) = NoahmpIO%BCPHOXY(I,LoopInd-NSNOW,J)
          OCPHIXY(I,LoopInd,J) = NoahmpIO%OCPHIXY(I,LoopInd-NSNOW,J)
          OCPHOXY(I,LoopInd,J) = NoahmpIO%OCPHOXY(I,LoopInd-NSNOW,J)
          DUST1XY(I,LoopInd,J) = NoahmpIO%DUST1XY(I,LoopInd-NSNOW,J)
          DUST2XY(I,LoopInd,J) = NoahmpIO%DUST2XY(I,LoopInd-NSNOW,J)
          DUST3XY(I,LoopInd,J) = NoahmpIO%DUST3XY(I,LoopInd-NSNOW,J)
          DUST4XY(I,LoopInd,J) = NoahmpIO%DUST4XY(I,LoopInd-NSNOW,J)
          DUST5XY(I,LoopInd,J) = NoahmpIO%DUST5XY(I,LoopInd-NSNOW,J)
          MassConcBCPHIXY(I,LoopInd,J) = NoahmpIO%MassConcBCPHIXY(I,LoopInd-NSNOW,J)
          MassConcBCPHOXY(I,LoopInd,J) = NoahmpIO%MassConcBCPHOXY(I,LoopInd-NSNOW,J)
          MassConcOCPHIXY(I,LoopInd,J) = NoahmpIO%MassConcOCPHIXY(I,LoopInd-NSNOW,J)
          MassConcOCPHOXY(I,LoopInd,J) = NoahmpIO%MassConcOCPHOXY(I,LoopInd-NSNOW,J)
          MassConcDUST1XY(I,LoopInd,J) = NoahmpIO%MassConcDUST1XY(I,LoopInd-NSNOW,J)
          MassConcDUST2XY(I,LoopInd,J) = NoahmpIO%MassConcDUST2XY(I,LoopInd-NSNOW,J)
          MassConcDUST3XY(I,LoopInd,J) = NoahmpIO%MassConcDUST3XY(I,LoopInd-NSNOW,J)
          MassConcDUST4XY(I,LoopInd,J) = NoahmpIO%MassConcDUST4XY(I,LoopInd-NSNOW,J)
          MassConcDUST5XY(I,LoopInd,J) = NoahmpIO%MassConcDUST5XY(I,LoopInd-NSNOW,J)
       enddo
    endif

    CROPCAT(I,J) = NoahmpIO%CROPCAT(I,J)

    enddo ! I
    enddo ! J

    ! Optional argument copies back - must be outside GPU parallel region
    if(present(QTDRAIN)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        QTDRAIN(I,J) = NoahmpIO%QTDRAIN(I,J)
      enddo; enddo
    endif
    if(present(SMCWTDXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        SMCWTDXY(I,J) = NoahmpIO%SMCWTDXY(I,J)
      enddo; enddo
    endif
    if(present(DEEPRECHXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        DEEPRECHXY(I,J) = NoahmpIO%DEEPRECHXY(I,J)
      enddo; enddo
    endif
    if(present(RECHXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        RECHXY(I,J) = NoahmpIO%RECHXY(I,J)
      enddo; enddo
    endif
    if(present(QRFSXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        QRFSXY(I,J) = NoahmpIO%QRFSXY(I,J)
      enddo; enddo
    endif
    if(present(QSPRINGSXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        QSPRINGSXY(I,J) = NoahmpIO%QSPRINGSXY(I,J)
      enddo; enddo
    endif
    if(present(QSLATXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        QSLATXY(I,J) = NoahmpIO%QSLATXY(I,J)
      enddo; enddo
    endif
    if(present(AREAXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        AREAXY(I,J) = NoahmpIO%AREAXY(I,J)
      enddo; enddo
    endif
    if(present(RIVERBEDXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        RIVERBEDXY(I,J) = NoahmpIO%RIVERBEDXY(I,J)
      enddo; enddo
    endif
    if(present(EQZWT)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        EQZWT(I,J) = NoahmpIO%EQZWT(I,J)
      enddo; enddo
    endif
    if(present(RIVERCONDXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        RIVERCONDXY(I,J) = NoahmpIO%RIVERCONDXY(I,J)
      enddo; enddo
    endif
    if(present(PEXPXY)) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do J = jts, jtf; do I = its, itf
        PEXPXY(I,J) = NoahmpIO%PEXPXY(I,J)
      enddo; enddo
    endif
    if(present(SMOISEQ)) then
      !$acc parallel loop gang vector collapse(2) default(present) firstprivate(LoopInd, NSOIL)
      do J = jts, jtf; do I = its, itf
        !$acc loop seq
        do LoopInd = 1, NSOIL
           SMOISEQ(I,LoopInd,J) = NoahmpIO%SMOISEQ(I,LoopInd,J)
        enddo
      enddo; enddo
    endif

    if(present(STEPWTD)) STEPWTD      = NoahmpIO%STEPWTD

    !--------- NoahmpIO -> WRF variables mapping ends

  end subroutine NoahmpSCHNAPSinit

end module NoahmpSCHNAPSinitMod
