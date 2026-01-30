module WaterVarInTransferMod

!!! Transfer input 2-D NoahmpIO Water variables to 1-D column variable
!!! 1-D variables should be first defined in /src/WaterVarType.F90
!!! 2-D variables should be first defined in NoahmpIOVarType.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type
  use NoahmpVarType
  use PedoTransferSR2006Mod

  implicit none

contains

!=== initialize with input data or table values

  subroutine WaterVarInTransfer(noahmp, NoahmpIO)

    implicit none

    type(noahmp_type),   intent(inout) :: noahmp
    type(NoahmpIO_type), intent(inout) :: NoahmpIO

    ! local variables 
    integer                            :: I, J
    integer                            :: IndexSoilLayer, LoopInd
    real(kind=kind_noahmp) :: SoilSand(1:noahmp%config%domain%NumSoilLayer)
    real(kind=kind_noahmp) :: SoilClay(1:noahmp%config%domain%NumSoilLayer)
    real(kind=kind_noahmp) :: SoilOrg(1:noahmp%config%domain%NumSoilLayer)

    !$acc parallel loop collapse(2) present(noahmp, NoahmpIO) private(IndexSoilLayer, SoilSand, SoilClay, SoilOrg)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
! -------------------------------------------------------------------------
    associate(                                                         &
              NumSnowLayerMax => noahmp%config%domain%NumSnowLayerMax ,&
              NumSoilLayer    => noahmp%config%domain%NumSoilLayer    ,&
              VegType         => noahmp%config%domain%VegType(I,J)         ,&
              SoilType        => noahmp%config%domain%SoilType             ,&
              FlagUrban       => noahmp%config%domain%FlagUrban(I,J)       ,&
              RunoffSlopeType => noahmp%config%domain%RunoffSlopeType      ,&
              NumSnowLayerNeg => noahmp%config%domain%NumSnowLayerNeg(I,J)  &
             )
! -------------------------------------------------------------------------

    ! water state variables
    noahmp%water%state%CanopyLiqWater(I,J)                     = NoahmpIO%CANLIQXY   (I,J)
    noahmp%water%state%CanopyIce(I,J)                          = NoahmpIO%CANICEXY   (I,J)
    noahmp%water%state%CanopyWetFrac(I,J)                      = NoahmpIO%FWETXY     (I,J)
    noahmp%water%state%SnowWaterEquiv(I,J)                     = NoahmpIO%SNOW       (I,J)
    noahmp%water%state%SnowWaterEquivPrev(I,J)                 = NoahmpIO%SNEQVOXY   (I,J) 
    noahmp%water%state%SnowDepth(I,J)                          = NoahmpIO%SNOWH      (I,J)
    noahmp%water%state%IrrigationFracFlood(I,J)                = NoahmpIO%FIFRACT    (I,J)
    noahmp%water%state%IrrigationAmtFlood(I,J)                 = NoahmpIO%IRWATFI    (I,J)
    noahmp%water%state%IrrigationFracMicro(I,J)                = NoahmpIO%MIFRACT    (I,J)
    noahmp%water%state%IrrigationAmtMicro(I,J)                 = NoahmpIO%IRWATMI    (I,J) 
    noahmp%water%state%IrrigationFracSprinkler(I,J)            = NoahmpIO%SIFRACT    (I,J)
    noahmp%water%state%IrrigationAmtSprinkler(I,J)             = NoahmpIO%IRWATSI    (I,J)  
    noahmp%water%state%WaterTableDepth(I,J)                    = NoahmpIO%ZWTXY      (I,J) 
    noahmp%water%state%SoilMoistureToWT(I,J)                   = NoahmpIO%SMCWTDXY   (I,J)
    noahmp%water%state%TileDrainFrac(I,J)                      = NoahmpIO%TD_FRACTION(I,J)
    noahmp%water%state%WaterStorageAquifer(I,J)                = NoahmpIO%WAXY       (I,J)
    noahmp%water%state%WaterStorageSoilAqf(I,J)                = NoahmpIO%WTXY       (I,J)
    noahmp%water%state%WaterStorageLake(I,J)                   = NoahmpIO%WSLAKEXY   (I,J)
    noahmp%water%state%IrrigationFracGrid(I,J)                 = NoahmpIO%IRFRACT    (I,J)
    noahmp%water%state%IrrigationCntSprinkler(I,J)             = NoahmpIO%IRNUMSI    (I,J)     
    noahmp%water%state%IrrigationCntMicro(I,J)                 = NoahmpIO%IRNUMMI    (I,J)
    noahmp%water%state%IrrigationCntFlood(I,J)                 = NoahmpIO%IRNUMFI    (I,J)
    !$acc loop seq
    do LoopInd = -NumSnowLayerMax+1, 0
      noahmp%water%state%SnowIce     (I,LoopInd,J) = NoahmpIO%SNICEXY    (I,LoopInd,J)
      noahmp%water%state%SnowLiqWater(I,LoopInd,J) = NoahmpIO%SNLIQXY    (I,LoopInd,J)
    end do
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
      noahmp%water%state%SoilLiqWater      (I,LoopInd,J) = NoahmpIO%SH2O       (I,LoopInd,J)
      noahmp%water%state%SoilMoisture      (I,LoopInd,J) = NoahmpIO%SMOIS      (I,LoopInd,J)    
      noahmp%water%state%SoilMoistureEqui  (I,LoopInd,J) = NoahmpIO%SMOISEQ    (I,LoopInd,J)

      ! Flux requiring loop over soil layers done here
      noahmp%water%flux%TranspWatLossSoilAcc(I,LoopInd,J)= NoahmpIO%ACC_ETRANIXY(I,LoopInd,J)
    enddo
    noahmp%water%state%RechargeGwDeepWT(I,J)                   = 0.0
    noahmp%water%state%RechargeGwShallowWT(I,J)                = 0.0
    if ( noahmp%config%nmlist%OptWetlandModel > 0 ) then
       noahmp%water%state%SoilSaturateFrac(I,J)                = NoahmpIO%FSATXY     (I,J)
       noahmp%water%state%WaterStorageWetland(I,J)             = NoahmpIO%WSURFXY    (I,J)
    endif
#ifdef WRF_HYDRO
    noahmp%water%state%WaterTableHydro(I,J)                    = NoahmpIO%ZWATBLE2D  (I,J)
    noahmp%water%state%WaterHeadSfc(I,J)                       = NoahmpIO%sfcheadrt  (I,J)
#endif
    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
      !$acc loop seq
      do LoopInd = -NumSnowLayerMax+1, 0
         noahmp%water%state%SnowRadius  (I,LoopInd,J)       = NoahmpIO%SNRDSXY (I,LoopInd,J)
         noahmp%water%state%MassBChydrophi(I,LoopInd,J)     = NoahmpIO%BCPHIXY (I,LoopInd,J)
         noahmp%water%state%MassBChydropho(I,LoopInd,J)     = NoahmpIO%BCPHOXY (I,LoopInd,J)
         noahmp%water%state%MassOChydrophi(I,LoopInd,J)     = NoahmpIO%OCPHIXY (I,LoopInd,J)
         noahmp%water%state%MassOChydropho(I,LoopInd,J)     = NoahmpIO%OCPHOXY (I,LoopInd,J)
         noahmp%water%state%MassDust1(I,LoopInd,J)          = NoahmpIO%DUST1XY (I,LoopInd,J)
         noahmp%water%state%MassDust2(I,LoopInd,J)          = NoahmpIO%DUST2XY (I,LoopInd,J)
         noahmp%water%state%MassDust3(I,LoopInd,J)          = NoahmpIO%DUST3XY (I,LoopInd,J)
         noahmp%water%state%MassDust4(I,LoopInd,J)          = NoahmpIO%DUST4XY (I,LoopInd,J)
         noahmp%water%state%MassDust5(I,LoopInd,J)          = NoahmpIO%DUST5XY (I,LoopInd,J)
         noahmp%water%state%MassConcBChydrophi(I,LoopInd,J) = NoahmpIO%MassConcBCPHIXY (I,LoopInd,J)
         noahmp%water%state%MassConcBChydropho(I,LoopInd,J) = NoahmpIO%MassConcBCPHOXY (I,LoopInd,J)
         noahmp%water%state%MassConcOChydrophi(I,LoopInd,J) = NoahmpIO%MassConcOCPHIXY (I,LoopInd,J)
         noahmp%water%state%MassConcOChydropho(I,LoopInd,J) = NoahmpIO%MassConcOCPHOXY (I,LoopInd,J)
         noahmp%water%state%MassConcDust1(I,LoopInd,J)      = NoahmpIO%MassConcDUST1XY (I,LoopInd,J)
         noahmp%water%state%MassConcDust2(I,LoopInd,J)      = NoahmpIO%MassConcDUST2XY (I,LoopInd,J)
         noahmp%water%state%MassConcDust3(I,LoopInd,J)      = NoahmpIO%MassConcDUST3XY (I,LoopInd,J)
         noahmp%water%state%MassConcDust4(I,LoopInd,J)      = NoahmpIO%MassConcDUST4XY (I,LoopInd,J)
         noahmp%water%state%MassConcDust5(I,LoopInd,J)      = NoahmpIO%MassConcDUST5XY (I,LoopInd,J)
         
         ! water flux variable for SNICAR done here
         noahmp%water%flux%SnowFreezeRate(I,LoopInd,J) = NoahmpIO%SNFRXY(I,LoopInd,J)
      end do
    endif


    ! water flux variables
    noahmp%water%flux%EvapSoilSfcLiqAcc(I,J)                   = NoahmpIO%ACC_QSEVAXY (I,J)
    noahmp%water%flux%SoilSfcInflowAcc(I,J)                    = NoahmpIO%ACC_QINSURXY(I,J)
    noahmp%water%flux%SfcWaterTotChgAcc(I,J)                   = NoahmpIO%ACC_DWATERXY(I,J)
    noahmp%water%flux%PrecipTotAcc(I,J)                        = NoahmpIO%ACC_PRCPXY  (I,J)
    noahmp%water%flux%EvapCanopyNetAcc(I,J)                    = NoahmpIO%ACC_ECANXY  (I,J)
    noahmp%water%flux%TranspirationAcc(I,J)                    = NoahmpIO%ACC_ETRANXY (I,J)
    noahmp%water%flux%EvapGroundNetAcc(I,J)                    = NoahmpIO%ACC_EDIRXY  (I,J)
    noahmp%water%flux%GlacierExcessFlowAcc(I,J)                = NoahmpIO%ACC_GLAFLWXY(I,J)


    ! water parameter variables
    noahmp%water%param%DrainSoilLayerInd(I,J)                  = NoahmpIO%DRAIN_LAYER_OPT_TABLE
    noahmp%water%param%CanopyLiqHoldCap(I,J)                   = NoahmpIO%CH2OP_TABLE(VegType)
    noahmp%water%param%SnowCompactBurdenFac(I,J)               = NoahmpIO%C2_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowCompactAgingFac1(I,J)               = NoahmpIO%C3_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowCompactAgingFac2(I,J)               = NoahmpIO%C4_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowCompactAgingFac3(I,J)               = NoahmpIO%C5_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowCompactAgingMax(I,J)                = NoahmpIO%DM_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowViscosityCoeff(I,J)                 = NoahmpIO%ETA0_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowCompactmAR24(I,J)                   = NoahmpIO%SNOWCOMPACTm_AR24_TABLE
    noahmp%water%param%SnowCompactbAR24(I,J)                   = NoahmpIO%SNOWCOMPACTb_AR24_TABLE
    noahmp%water%param%SnowCompactP1AR24(I,J)                  = NoahmpIO%SNOWCOMPACT_P1_AR24_TABLE
    noahmp%water%param%SnowCompactP2AR24(I,J)                  = NoahmpIO%SNOWCOMPACT_P2_AR24_TABLE
    noahmp%water%param%SnowCompactP3AR24(I,J)                  = NoahmpIO%SNOWCOMPACT_P3_AR24_TABLE
    noahmp%water%param%SnowCoverM1AR25(I,J)                    = NoahmpIO%SCFm1_AR25_TABLE
    noahmp%water%param%SnowCoverM2AR25(I,J)                    = NoahmpIO%SCFm2_AR25_TABLE
    noahmp%water%param%SnowCoverFac1AR25(I,J)                  = NoahmpIO%SCfac1_AR25_TABLE
    noahmp%water%param%SnowCoverFac2AR25(I,J)                  = NoahmpIO%SCfac2_AR25_TABLE
    noahmp%water%param%BurdenFacUpAR24(I,J)                    = NoahmpIO%SNOWCOMPACT_Up_AR24_TABLE
    noahmp%water%param%SnowLiqFracMax(I,J)                     = NoahmpIO%SNLIQMAXFRAC_TABLE
    noahmp%water%param%SnowLiqHoldCap(I,J)                     = NoahmpIO%SSI_TABLE
    noahmp%water%param%SnowLiqReleaseFac(I,J)                  = NoahmpIO%SNOW_RET_FAC_TABLE
    noahmp%water%param%IrriFloodRateFac(I,J)                   = NoahmpIO%FIRTFAC_TABLE
    noahmp%water%param%IrriMicroRate(I,J)                      = NoahmpIO%MICIR_RATE_TABLE
    noahmp%water%param%SoilConductivityRef(I,J)                = NoahmpIO%REFDK_TABLE
    noahmp%water%param%SoilInfilFacRef(I,J)                    = NoahmpIO%REFKDT_TABLE
    noahmp%water%param%GroundFrzCoeff(I,J)                     = NoahmpIO%FRZK_TABLE
    noahmp%water%param%GridTopoIndex(I,J)                      = NoahmpIO%TIMEAN_TABLE
    noahmp%water%param%SoilSfcSatFracMax(I,J)                  = NoahmpIO%FSATMX_TABLE
    noahmp%water%param%SpecYieldGw(I,J)                        = NoahmpIO%ROUS_TABLE
    noahmp%water%param%MicroPoreContent(I,J)                   = NoahmpIO%CMIC_TABLE
    noahmp%water%param%WaterStorageLakeMax(I,J)                = NoahmpIO%WSLMAX_TABLE
    noahmp%water%param%SnoWatEqvMaxGlacier(I,J)                = NoahmpIO%SWEMAXGLA_TABLE
    noahmp%water%param%IrriStopDayBfHarvest(I,J)               = NoahmpIO%IRR_HAR_TABLE
    noahmp%water%param%IrriTriggerLaiMin(I,J)                  = NoahmpIO%IRR_LAI_TABLE
    noahmp%water%param%SoilWatDeficitAllow(I,J)                = NoahmpIO%IRR_MAD_TABLE
    noahmp%water%param%IrriFloodLossFrac(I,J)                  = NoahmpIO%FILOSS_TABLE
    noahmp%water%param%IrriSprinklerRate(I,J)                  = NoahmpIO%SPRIR_RATE_TABLE
    noahmp%water%param%IrriFracThreshold(I,J)                  = NoahmpIO%IRR_FRAC_TABLE
    noahmp%water%param%IrriStopPrecipThr(I,J)                  = NoahmpIO%IR_RAIN_TABLE
    noahmp%water%param%SnowfallDensityMax(I,J)                 = NoahmpIO%SNOWDEN_MAX_TABLE
    noahmp%water%param%SnowMassFullCoverOld(I,J)               = NoahmpIO%SWEMX_TABLE
    noahmp%water%param%SoilMatPotentialWilt(I,J)               = NoahmpIO%PSIWLT_TABLE
    noahmp%water%param%SnowMeltFac(I,J)                        = NoahmpIO%MFSNO_TABLE(VegType)
    noahmp%water%param%SnowCoverFac(I,J)                       = NoahmpIO%SCFFAC_TABLE(VegType)
    noahmp%water%param%InfilFacVic(I,J)                        = NoahmpIO%BVIC_TABLE(SoilType(I,1,J))
    noahmp%water%param%TensionWatDistrInfl(I,J)                = NoahmpIO%AXAJ_TABLE(SoilType(I,1,J))
    noahmp%water%param%TensionWatDistrShp(I,J)                 = NoahmpIO%BXAJ_TABLE(SoilType(I,1,J))
    noahmp%water%param%FreeWatDistrShp(I,J)                    = NoahmpIO%XXAJ_TABLE(SoilType(I,1,J))
    noahmp%water%param%InfilHeteroDynVic(I,J)                  = NoahmpIO%BBVIC_TABLE(SoilType(I,1,J))
    noahmp%water%param%InfilCapillaryDynVic(I,J)               = NoahmpIO%GDVIC_TABLE(SoilType(I,1,J))
    noahmp%water%param%InfilFacDynVic(I,J)                     = NoahmpIO%BDVIC_TABLE(SoilType(I,1,J))
    noahmp%water%param%TileDrainCoeffSp(I,J)                   = NoahmpIO%TD_DC_TABLE(SoilType(I,1,J))
    noahmp%water%param%TileDrainTubeDepth(I,J)                 = NoahmpIO%TD_DEPTH_TABLE(SoilType(I,1,J))
    noahmp%water%param%DrainFacSoilWat(I,J)                    = NoahmpIO%TDSMC_FAC_TABLE(SoilType(I,1,J))
    noahmp%water%param%TileDrainCoeff(I,J)                     = NoahmpIO%TD_DCOEF_TABLE(SoilType(I,1,J))
    noahmp%water%param%DrainDepthToImperv(I,J)                 = NoahmpIO%TD_ADEPTH_TABLE(SoilType(I,1,J))
    noahmp%water%param%LateralWatCondFac(I,J)                  = NoahmpIO%KLAT_FAC_TABLE(SoilType(I,1,J))
    noahmp%water%param%TileDrainDepth(I,J)                     = NoahmpIO%TD_DDRAIN_TABLE(SoilType(I,1,J))
    noahmp%water%param%DrainTubeDist(I,J)                      = NoahmpIO%TD_SPAC_TABLE(SoilType(I,1,J))
    noahmp%water%param%DrainTubeRadius(I,J)                    = NoahmpIO%TD_RADI_TABLE(SoilType(I,1,J))
    noahmp%water%param%DrainWatDepToImperv(I,J)                = NoahmpIO%TD_D_TABLE(SoilType(I,1,J))
    noahmp%water%param%NumSoilLayerRoot(I,J)                   = NoahmpIO%NROOT_TABLE(VegType)
    noahmp%water%param%SoilDrainSlope(I,J)                     = NoahmpIO%SLOPE_TABLE(RunoffSlopeType)
    noahmp%water%param%WetlandCapMax(I,J)                      = NoahmpIO%WCAP_TABLE

    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 )then
       noahmp%water%param%SnowRadiusMin(I,J)                   = NoahmpIO%SnowRadiusMin_TABLE
       noahmp%water%param%FreshSnowRadiusMax(I,J)              = NoahmpIO%FreshSnowRadiusMax_TABLE
       noahmp%water%param%SnowRadiusRefrz(I,J)                 = NoahmpIO%SnowRadiusRefrz_TABLE
       noahmp%water%param%ScavEffMeltScale(I,J)                = NoahmpIO%ScavEffMeltScale_TABLE
       noahmp%water%param%ScavEffMeltBCphi(I,J)                = NoahmpIO%ScavEffMeltBCphi_TABLE
       noahmp%water%param%ScavEffMeltBCpho(I,J)                = NoahmpIO%ScavEffMeltBCpho_TABLE
       noahmp%water%param%ScavEffMeltOCphi(I,J)                = NoahmpIO%ScavEffMeltOCphi_TABLE
       noahmp%water%param%ScavEffMeltOCpho(I,J)                = NoahmpIO%ScavEffMeltOCpho_TABLE
       noahmp%water%param%ScavEffMeltDust1(I,J)                = NoahmpIO%ScavEffMeltDust1_TABLE
       noahmp%water%param%ScavEffMeltDust2(I,J)                = NoahmpIO%ScavEffMeltDust2_TABLE
       noahmp%water%param%ScavEffMeltDust3(I,J)                = NoahmpIO%ScavEffMeltDust3_TABLE
       noahmp%water%param%ScavEffMeltDust4(I,J)                = NoahmpIO%ScavEffMeltDust4_TABLE
       noahmp%water%param%ScavEffMeltDust5(I,J)                = NoahmpIO%ScavEffMeltDust5_TABLE
       noahmp%water%param%SnowRadiusMax(I,J)                   = NoahmpIO%SnowRadiusMax_TABLE
       noahmp%water%param%SnowWetAgeC1Brun89(I,J)              = NoahmpIO%SnowWetAgeC1Brun89_TABLE
       noahmp%water%param%SnowWetAgeC2Brun89(I,J)              = NoahmpIO%SnowWetAgeC2Brun89_TABLE
       noahmp%water%param%SnowAgeScaleFac(I,J)                 = NoahmpIO%SnowAgeScaleFac_TABLE
    endif

    ! soil properties
    do IndexSoilLayer = 1, size(SoilType,2)
       noahmp%water%param%SoilMoistureSat       (I,IndexSoilLayer,J) = NoahmpIO%SMCMAX_TABLE(SoilType(I,IndexSoilLayer,J))
       noahmp%water%param%SoilMoistureWilt      (I,IndexSoilLayer,J) = NoahmpIO%SMCWLT_TABLE(SoilType(I,IndexSoilLayer,J))
       noahmp%water%param%SoilMoistureFieldCap  (I,IndexSoilLayer,J) = NoahmpIO%SMCREF_TABLE(SoilType(I,IndexSoilLayer,J))
       noahmp%water%param%SoilMoistureDry       (I,IndexSoilLayer,J) = NoahmpIO%SMCDRY_TABLE(SoilType(I,IndexSoilLayer,J))
       noahmp%water%param%SoilWatDiffusivitySat (I,IndexSoilLayer,J) = NoahmpIO%DWSAT_TABLE (SoilType(I,IndexSoilLayer,J))
       noahmp%water%param%SoilWatConductivitySat(I,IndexSoilLayer,J) = NoahmpIO%DKSAT_TABLE (SoilType(I,IndexSoilLayer,J))
       noahmp%water%param%SoilExpCoeffB         (I,IndexSoilLayer,J) = NoahmpIO%BEXP_TABLE  (SoilType(I,IndexSoilLayer,J))
       noahmp%water%param%SoilMatPotentialSat   (I,IndexSoilLayer,J) = NoahmpIO%PSISAT_TABLE(SoilType(I,IndexSoilLayer,J))
    enddo
   
    ! spatial varying soil texture and properties directly from input
    if ( noahmp%config%nmlist%OptSoilProperty == 4 ) then
       ! 3D soil properties
       !$acc loop seq
       do IndexSoilLayer = 1, NumSoilLayer
         noahmp%water%param%SoilExpCoeffB(I,IndexSoilLayer,J)          = NoahmpIO%BEXP_3D  (I,IndexSoilLayer,J) ! C-H B exponent
         noahmp%water%param%SoilMoistureDry(I,IndexSoilLayer,J)        = NoahmpIO%SMCDRY_3D(I,IndexSoilLayer,J) ! Soil Moisture Limit: Dry
         noahmp%water%param%SoilMoistureWilt(I,IndexSoilLayer,J)       = NoahmpIO%SMCWLT_3D(I,IndexSoilLayer,J) ! Soil Moisture Limit: Wilt
         noahmp%water%param%SoilMoistureFieldCap(I,IndexSoilLayer,J)   = NoahmpIO%SMCREF_3D(I,IndexSoilLayer,J) ! Soil Moisture Limit: Reference
         noahmp%water%param%SoilMoistureSat(I,IndexSoilLayer,J)        = NoahmpIO%SMCMAX_3D(I,IndexSoilLayer,J) ! Soil Moisture Limit: Max
         noahmp%water%param%SoilWatConductivitySat(I,IndexSoilLayer,J) = NoahmpIO%DKSAT_3D (I,IndexSoilLayer,J) ! Saturated Soil Conductivity
         noahmp%water%param%SoilWatDiffusivitySat(I,IndexSoilLayer,J)  = NoahmpIO%DWSAT_3D (I,IndexSoilLayer,J) ! Saturated Soil Diffusivity
         noahmp%water%param%SoilMatPotentialSat(I,IndexSoilLayer,J)    = NoahmpIO%PSISAT_3D(I,IndexSoilLayer,J) ! Saturated Matric Potential
      enddo
       noahmp%water%param%SoilConductivityRef(I,J)    = NoahmpIO%REFDK_2D (I,J)                ! Reference Soil Conductivity
       noahmp%water%param%SoilInfilFacRef(I,J)        = NoahmpIO%REFKDT_2D(I,J)                ! Soil Infiltration Parameter
       ! 2D additional runoff6~8 parameters
       noahmp%water%param%InfilFacVic(I,J)            = NoahmpIO%BVIC_2D (I,J)                 ! VIC model infiltration parameter
       noahmp%water%param%TensionWatDistrInfl(I,J)    = NoahmpIO%AXAJ_2D (I,J)                 ! Xinanjiang: Tension water distribution inflection parameter
       noahmp%water%param%TensionWatDistrShp(I,J)     = NoahmpIO%BXAJ_2D (I,J)                 ! Xinanjiang: Tension water distribution shape parameter
       noahmp%water%param%FreeWatDistrShp(I,J)        = NoahmpIO%XXAJ_2D (I,J)                 ! Xinanjiang: Free water distribution shape parameter
       noahmp%water%param%InfilFacDynVic(I,J)         = NoahmpIO%BDVIC_2D(I,J)                 ! VIC model infiltration parameter
       noahmp%water%param%InfilCapillaryDynVic(I,J)   = NoahmpIO%GDVIC_2D(I,J)                 ! Mean Capillary Drive for infiltration models
       noahmp%water%param%InfilHeteroDynVic(I,J)      = NoahmpIO%BBVIC_2D(I,J)                 ! DVIC heterogeniety parameter for infiltraton
       ! 2D irrigation params
       noahmp%water%param%IrriFracThreshold(I,J)      = NoahmpIO%IRR_FRAC_2D  (I,J)            ! irrigation Fraction
       noahmp%water%param%IrriStopDayBfHarvest(I,J)   = NoahmpIO%IRR_HAR_2D   (I,J)            ! number of days before harvest date to stop irrigation 
       noahmp%water%param%IrriTriggerLaiMin(I,J)      = NoahmpIO%IRR_LAI_2D   (I,J)            ! Minimum lai to trigger irrigation
       noahmp%water%param%SoilWatDeficitAllow(I,J)    = NoahmpIO%IRR_MAD_2D   (I,J)            ! management allowable deficit (0-1)
       noahmp%water%param%IrriFloodLossFrac(I,J)      = NoahmpIO%FILOSS_2D    (I,J)            ! fraction of flood irrigation loss (0-1) 
       noahmp%water%param%IrriSprinklerRate(I,J)      = NoahmpIO%SPRIR_RATE_2D(I,J)            ! mm/h, sprinkler irrigation rate
       noahmp%water%param%IrriMicroRate(I,J)          = NoahmpIO%MICIR_RATE_2D(I,J)            ! mm/h, micro irrigation rate
       noahmp%water%param%IrriFloodRateFac(I,J)       = NoahmpIO%FIRTFAC_2D   (I,J)            ! flood application rate factor
       noahmp%water%param%IrriStopPrecipThr(I,J)      = NoahmpIO%IR_RAIN_2D   (I,J)            ! maximum precipitation to stop irrigation trigger
       ! 2D tile drainage parameters
       noahmp%water%param%LateralWatCondFac(I,J)      = NoahmpIO%KLAT_FAC (I,J)                ! factor multiplier to hydraulic conductivity
       noahmp%water%param%DrainFacSoilWat(I,J)        = NoahmpIO%TDSMC_FAC(I,J)                ! factor multiplier to field capacity
       noahmp%water%param%TileDrainCoeffSp(I,J)       = NoahmpIO%TD_DC    (I,J)                ! drainage coefficient for simple
       noahmp%water%param%TileDrainCoeff(I,J)         = NoahmpIO%TD_DCOEF (I,J)                ! drainge coefficient for Hooghoudt 
       noahmp%water%param%TileDrainDepth(I,J)         = NoahmpIO%TD_DDRAIN(I,J)                ! depth of drain
       noahmp%water%param%DrainTubeRadius(I,J)        = NoahmpIO%TD_RADI  (I,J)                ! tile tube radius
       noahmp%water%param%DrainTubeDist(I,J)          = NoahmpIO%TD_SPAC  (I,J)                ! tile spacing
    endif

    ! spatial varying wetland parameters from input
    if ( noahmp%config%nmlist%OptWetlandModel == 2 ) then
       noahmp%water%param%SoilSfcSatFracMax(I,J)      = NoahmpIO%FSATMX(I,J)
       noahmp%water%param%WetlandCapMax(I,J)          = NoahmpIO%WCAP(I,J)
    endif 

    ! derived water parameters
    noahmp%water%param%SoilInfilMaxCoeff(I,J)  = noahmp%water%param%SoilInfilFacRef(I,J) *           &
                                            noahmp%water%param%SoilWatConductivitySat(I,1,J) / &
                                            noahmp%water%param%SoilConductivityRef(I,J)
    if ( FlagUrban .eqv. .true. ) then
       !$acc loop seq
       do IndexSoilLayer = 1, NumSoilLayer
         noahmp%water%param%SoilMoistureSat(I,IndexSoilLayer,J)      = 0.45
         noahmp%water%param%SoilMoistureFieldCap(I,IndexSoilLayer,J) = 0.42
         noahmp%water%param%SoilMoistureWilt(I,IndexSoilLayer,J)     = 0.40
         noahmp%water%param%SoilMoistureDry(I,IndexSoilLayer,J)      = 0.40
       enddo
    endif

    if ( SoilType(I,1,J) /= 14 ) then
       noahmp%water%param%SoilImpervFracCoeff(I,J) = noahmp%water%param%GroundFrzCoeff(I,J) *       &
                                                ((noahmp%water%param%SoilMoistureSat(I,1,J) / &
                                                 noahmp%water%param%SoilMoistureFieldCap(I,1,J)) * (0.412/0.468))
    endif

   !  noahmp%water%state%SnowIceFracPrev(I,J) = 0.0
    noahmp%water%state%SnowIceFracPrev(I,NumSnowLayerNeg+1:0,J) = NoahmpIO%SNICEXY(I,NumSnowLayerNeg+1:0,J) /  & 
                                                              (NoahmpIO%SNICEXY(I,NumSnowLayerNeg+1:0,J) + &
                                                               NoahmpIO%SNLIQXY(I,NumSnowLayerNeg+1:0,J))

    if ( (noahmp%config%nmlist%OptSoilProperty == 3) .and. (.not. noahmp%config%domain%FlagUrban(I,J)) ) then
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
         SoilSand(LoopInd) = 0.01 * NoahmpIO%soilcomp(I,LoopInd,J)
         SoilClay(LoopInd) = 0.01 * NoahmpIO%soilcomp(I,NumSoilLayer+LoopInd,J)
         SoilOrg(LoopInd)  = 0.0
       end do
       if (noahmp%config%nmlist%OptPedotransfer == 1) &
          call PedoTransferSR2006(NoahmpIO,noahmp,SoilSand,SoilClay,SoilOrg,I,J)
    endif

    end associate

      enddo
   enddo

    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 )then
       !$acc loop gang vector collapse(3)
       do I = 1, noahmp%config%domain%NumTempSnwAgeSnicar
       do LoopInd = 1, noahmp%config%domain%NumTempGradSnwAgeSnicar
       do J = 1, noahmp%config%domain%NumDensitySnwAgeSnicar
          noahmp%water%param%snowage_tau(I,LoopInd,J)     = NoahmpIO%snowage_tau(I,LoopInd,J)
          noahmp%water%param%snowage_kappa(I,LoopInd,J)   = NoahmpIO%snowage_kappa(I,LoopInd,J)
          noahmp%water%param%snowage_drdt0(I,LoopInd,J)   = NoahmpIO%snowage_drdt0(I,LoopInd,J)
       end do
       end do
       end do
    endif

   end subroutine WaterVarInTransfer

end module WaterVarInTransferMod
