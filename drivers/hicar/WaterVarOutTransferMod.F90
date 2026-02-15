module WaterVarOutTransferMod

!!! Transfer column (1-D) Noah-MP water variables to 2D NoahmpIO for output

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type
  use NoahmpVarType

  implicit none

contains

!=== Transfer model states to output =====

  subroutine WaterVarOutTransfer(noahmp, NoahmpIO)

    implicit none

    type(noahmp_type),   intent(inout) :: noahmp
    type(NoahmpIO_type), intent(inout) :: NoahmpIO

    integer :: I, J, LoopInd
! -------------------------------------------------------------------------
    associate(                                                         &
              NumSnowLayerMax => noahmp%config%domain%NumSnowLayerMax ,&
              NumSoilLayer    => noahmp%config%domain%NumSoilLayer     &
             )
! -------------------------------------------------------------------------

    !$acc parallel loop collapse(2) present(noahmp, NoahmpIO) private(NumSnowLayerMax, NumSoilLayer)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE 

        if (NoahmpIO%XLAND(I,J) - 1.5 >= 0.0) cycle ! Do out write output for open water points

    ! special treatment for glacier point output
    if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) then ! land ice point
       noahmp%water%state%SnowCoverFrac(I,J)      = 1.0
       noahmp%water%flux%EvapCanopyNet(I,J)       = 0.0
       noahmp%water%flux%Transpiration(I,J)       = 0.0
       noahmp%water%flux%InterceptCanopySnow(I,J) = 0.0
       noahmp%water%flux%InterceptCanopyRain(I,J) = 0.0
       noahmp%water%flux%DripCanopySnow(I,J)      = 0.0
       noahmp%water%flux%DripCanopyRain(I,J)      = 0.0
       noahmp%water%flux%ThroughfallSnow(I,J)     = noahmp%water%flux%SnowfallRefHeight(I,J)
       noahmp%water%flux%ThroughfallRain(I,J)     = noahmp%water%flux%RainfallRefHeight(I,J)
       noahmp%water%flux%SublimCanopyIce(I,J)     = 0.0
       noahmp%water%flux%FrostCanopyIce(I,J)      = 0.0
       noahmp%water%flux%FreezeCanopyLiq(I,J)     = 0.0
       noahmp%water%flux%MeltCanopyIce(I,J)       = 0.0
       noahmp%water%flux%EvapCanopyLiq(I,J)       = 0.0
       noahmp%water%flux%DewCanopyLiq(I,J)        = 0.0
       noahmp%water%state%CanopyIce(I,J)          = 0.0
       noahmp%water%state%CanopyLiqWater(I,J)     = 0.0
       noahmp%water%flux%TileDrain(I,J)           = 0.0
       noahmp%water%flux%RunoffSurface(I,J)       = noahmp%water%flux%RunoffSurface(I,J) * noahmp%config%domain%MainTimeStep
       noahmp%water%flux%RunoffSubsurface(I,J)    = noahmp%water%flux%RunoffSubsurface(I,J) * noahmp%config%domain%MainTimeStep
       NoahmpIO%QFX(I,J)                     = noahmp%water%flux%EvapGroundNet(I,J)
    endif

    if ( noahmp%config%domain%IndicatorIceSfc(I,J) == 0 ) then ! land soil point
       NoahmpIO%QFX(I,J) = noahmp%water%flux%EvapCanopyNet(I,J) + noahmp%water%flux%EvapGroundNet(I,J) + &
                           noahmp%water%flux%Transpiration(I,J) + noahmp%water%flux%EvapIrriSprinkler(I,J)
    endif

    NoahmpIO%SMSTAV      (I,J) = 0.0  ! [maintained as Noah consistency] water
    NoahmpIO%SMSTOT      (I,J) = 0.0  ! [maintained as Noah consistency] water
    NoahmpIO%SFCRUNOFF   (I,J) = NoahmpIO%SFCRUNOFF(I,J) + noahmp%water%flux%RunoffSurface(I,J)
    NoahmpIO%UDRUNOFF    (I,J) = NoahmpIO%UDRUNOFF (I,J) + noahmp%water%flux%RunoffSubsurface(I,J)
    NoahmpIO%QTDRAIN     (I,J) = NoahmpIO%QTDRAIN  (I,J) + noahmp%water%flux%TileDrain(I,J)
    NoahmpIO%SNOWC       (I,J) = noahmp%water%state%SnowCoverFrac(I,J)
    NoahmpIO%SNOW        (I,J) = noahmp%water%state%SnowWaterEquiv(I,J)
    NoahmpIO%SNOWH       (I,J) = noahmp%water%state%SnowDepth(I,J)
    NoahmpIO%CANWAT      (I,J) = noahmp%water%state%CanopyLiqWater(I,J) + noahmp%water%state%CanopyIce(I,J)
    NoahmpIO%ACSNOW      (I,J) = NoahmpIO%ACSNOW(I,J) + NoahmpIO%RAINBL (I,J) * noahmp%water%state%FrozenPrecipFrac(I,J)
    NoahmpIO%ACSNOM      (I,J) = NoahmpIO%ACSNOM(I,J) + noahmp%water%flux%MeltGroundSnow(I,J) * NoahmpIO%DTBL
    NoahmpIO%CANLIQXY    (I,J) = noahmp%water%state%CanopyLiqWater(I,J)
    NoahmpIO%CANICEXY    (I,J) = noahmp%water%state%CanopyIce(I,J)
    NoahmpIO%FWETXY      (I,J) = noahmp%water%state%CanopyWetFrac(I,J)
    NoahmpIO%SNEQVOXY    (I,J) = noahmp%water%state%SnowWaterEquivPrev(I,J)
    NoahmpIO%QSNOWXY     (I,J) = noahmp%water%flux%SnowfallGround(I,J)
    NoahmpIO%QRAINXY     (I,J) = noahmp%water%flux%RainfallGround(I,J)
    NoahmpIO%WSLAKEXY    (I,J) = noahmp%water%state%WaterStorageLake(I,J)
    NoahmpIO%ZWTXY       (I,J) = noahmp%water%state%WaterTableDepth(I,J)
    NoahmpIO%WAXY        (I,J) = noahmp%water%state%WaterStorageAquifer(I,J)
    NoahmpIO%WTXY        (I,J) = noahmp%water%state%WaterStorageSoilAqf(I,J)
    NoahmpIO%RUNSFXY     (I,J) = noahmp%water%flux%RunoffSurface(I,J)
    NoahmpIO%RUNSBXY     (I,J) = noahmp%water%flux%RunoffSubsurface(I,J)
    NoahmpIO%ECANXY      (I,J) = noahmp%water%flux%EvapCanopyNet(I,J)
    NoahmpIO%EDIRXY      (I,J) = noahmp%water%flux%EvapGroundNet(I,J)
    NoahmpIO%ETRANXY     (I,J) = noahmp%water%flux%Transpiration(I,J)
    NoahmpIO%QINTSXY     (I,J) = noahmp%water%flux%InterceptCanopySnow(I,J)
    NoahmpIO%QINTRXY     (I,J) = noahmp%water%flux%InterceptCanopyRain(I,J)
    NoahmpIO%QDRIPSXY    (I,J) = noahmp%water%flux%DripCanopySnow(I,J)
    NoahmpIO%QDRIPRXY    (I,J) = noahmp%water%flux%DripCanopyRain(I,J)
    NoahmpIO%QTHROSXY    (I,J) = noahmp%water%flux%ThroughfallSnow(I,J)
    NoahmpIO%QTHRORXY    (I,J) = noahmp%water%flux%ThroughfallRain(I,J)
    NoahmpIO%QSNSUBXY    (I,J) = noahmp%water%flux%SublimSnowSfcIce(I,J)
    NoahmpIO%QSNFROXY    (I,J) = noahmp%water%flux%FrostSnowSfcIce(I,J)
    NoahmpIO%QSUBCXY     (I,J) = noahmp%water%flux%SublimCanopyIce(I,J)
    NoahmpIO%QFROCXY     (I,J) = noahmp%water%flux%FrostCanopyIce(I,J)
    NoahmpIO%QEVACXY     (I,J) = noahmp%water%flux%EvapCanopyLiq(I,J)
    NoahmpIO%QDEWCXY     (I,J) = noahmp%water%flux%DewCanopyLiq(I,J)
    NoahmpIO%QFRZCXY     (I,J) = noahmp%water%flux%FreezeCanopyLiq(I,J)
    NoahmpIO%QMELTCXY    (I,J) = noahmp%water%flux%MeltCanopyIce(I,J)
    NoahmpIO%QSNBOTXY    (I,J) = noahmp%water%flux%SnowBotOutflow(I,J)
    NoahmpIO%QMELTXY     (I,J) = noahmp%water%flux%MeltGroundSnow(I,J)
    NoahmpIO%PONDINGXY   (I,J) = noahmp%water%state%PondSfcThinSnwTrans(I,J) + &
                                 noahmp%water%state%PondSfcThinSnwComb(I,J) + noahmp%water%state%PondSfcThinSnwMelt(I,J)
    NoahmpIO%FPICEXY     (I,J) = noahmp%water%state%FrozenPrecipFrac(I,J)
    NoahmpIO%RAINLSM     (I,J) = noahmp%water%flux%RainfallRefHeight(I,J)
    NoahmpIO%SNOWLSM     (I,J) = noahmp%water%flux%SnowfallRefHeight(I,J)
    NoahmpIO%ACC_QINSURXY(I,J) = noahmp%water%flux%SoilSfcInflowAcc(I,J)
    NoahmpIO%ACC_QSEVAXY (I,J) = noahmp%water%flux%EvapSoilSfcLiqAcc(I,J)
    NoahmpIO%ACC_DWATERXY(I,J) = noahmp%water%flux%SfcWaterTotChgAcc(I,J)
    NoahmpIO%ACC_PRCPXY  (I,J) = noahmp%water%flux%PrecipTotAcc(I,J)
    NoahmpIO%ACC_ECANXY  (I,J) = noahmp%water%flux%EvapCanopyNetAcc(I,J)
    NoahmpIO%ACC_ETRANXY (I,J) = noahmp%water%flux%TranspirationAcc(I,J)
    NoahmpIO%ACC_EDIRXY  (I,J) = noahmp%water%flux%EvapGroundNetAcc(I,J)
    NoahmpIO%ACC_GLAFLWXY(I,J) = noahmp%water%flux%GlacierExcessFlowAcc(I,J)
    NoahmpIO%RECHXY      (I,J) = NoahmpIO%RECHXY(I,J) + (noahmp%water%state%RechargeGwShallowWT(I,J)*1.0e3)
    NoahmpIO%DEEPRECHXY  (I,J) = NoahmpIO%DEEPRECHXY(I,J) + noahmp%water%state%RechargeGwDeepWT(I,J)
    NoahmpIO%SMCWTDXY    (I,J) = noahmp%water%state%SoilMoistureToWT(I,J)
      !$acc loop seq
      do LoopInd = 1, NumSoilLayer
      NoahmpIO%SMOIS       (I,LoopInd,J)       = noahmp%water%state%SoilMoisture(I,LoopInd,J)
      NoahmpIO%SH2O        (I,LoopInd,J)       = noahmp%water%state%SoilLiqWater(I,LoopInd,J)
      NoahmpIO%ACC_ETRANIXY(I,LoopInd,J)       = noahmp%water%flux%TranspWatLossSoilAcc(I,LoopInd,J)
     end do
     !$acc loop seq
     do LoopInd = -NumSnowLayerMax+1, 0
      NoahmpIO%SNICEXY     (I,LoopInd,J) = noahmp%water%state%SnowIce(I,LoopInd,J)
      NoahmpIO%SNLIQXY     (I,LoopInd,J) = noahmp%water%state%SnowLiqWater(I,LoopInd,J)
      end do
    !SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       !$acc loop seq
       do LoopInd = -NumSnowLayerMax+1, 0
        NoahmpIO%SNRDSXY(I,LoopInd,J) = noahmp%water%state%SnowRadius(I,LoopInd,J)
        NoahmpIO%SNFRXY (I,LoopInd,J) = noahmp%water%flux%SnowFreezeRate(I,LoopInd,J)
        NoahmpIO%BCPHIXY(I,LoopInd,J) = noahmp%water%state%MassBChydrophi(I,LoopInd,J)
        NoahmpIO%BCPHOXY(I,LoopInd,J) = noahmp%water%state%MassBChydropho(I,LoopInd,J)
        NoahmpIO%OCPHIXY(I,LoopInd,J) = noahmp%water%state%MassOChydrophi(I,LoopInd,J)
        NoahmpIO%OCPHOXY(I,LoopInd,J) = noahmp%water%state%MassOChydropho(I,LoopInd,J)
        NoahmpIO%DUST1XY(I,LoopInd,J) = noahmp%water%state%MassDust1(I,LoopInd,J)
        NoahmpIO%DUST2XY(I,LoopInd,J) = noahmp%water%state%MassDust2(I,LoopInd,J)
        NoahmpIO%DUST3XY(I,LoopInd,J) = noahmp%water%state%MassDust3(I,LoopInd,J)
        NoahmpIO%DUST4XY(I,LoopInd,J) = noahmp%water%state%MassDust4(I,LoopInd,J)
        NoahmpIO%DUST5XY(I,LoopInd,J) = noahmp%water%state%MassDust5(I,LoopInd,J)
        NoahmpIO%MassConcBCPHIXY(I,LoopInd,J) = noahmp%water%state%MassConcBChydrophi(I,LoopInd,J)
        NoahmpIO%MassConcBCPHOXY(I,LoopInd,J) = noahmp%water%state%MassConcBChydropho(I,LoopInd,J)
        NoahmpIO%MassConcOCPHIXY(I,LoopInd,J) = noahmp%water%state%MassConcOChydrophi(I,LoopInd,J)
        NoahmpIO%MassConcOCPHOXY(I,LoopInd,J) = noahmp%water%state%MassConcOChydropho(I,LoopInd,J)
        NoahmpIO%MassConcDUST1XY(I,LoopInd,J) = noahmp%water%state%MassConcDust1(I,LoopInd,J)
        NoahmpIO%MassConcDUST2XY(I,LoopInd,J) = noahmp%water%state%MassConcDust2(I,LoopInd,J)
        NoahmpIO%MassConcDUST3XY(I,LoopInd,J) = noahmp%water%state%MassConcDust3(I,LoopInd,J)
        NoahmpIO%MassConcDUST4XY(I,LoopInd,J) = noahmp%water%state%MassConcDust4(I,LoopInd,J)
        NoahmpIO%MassConcDUST5XY(I,LoopInd,J) = noahmp%water%state%MassConcDust5(I,LoopInd,J)
      end do
    endif

    ! irrigation
    NoahmpIO%IRNUMSI(I,J) = noahmp%water%state%IrrigationCntSprinkler(I,J)
    NoahmpIO%IRNUMMI(I,J) = noahmp%water%state%IrrigationCntMicro(I,J)
    NoahmpIO%IRNUMFI(I,J) = noahmp%water%state%IrrigationCntFlood(I,J)
    NoahmpIO%IRWATSI(I,J) = noahmp%water%state%IrrigationAmtSprinkler(I,J)
    NoahmpIO%IRWATMI(I,J) = noahmp%water%state%IrrigationAmtMicro(I,J)
    NoahmpIO%IRWATFI(I,J) = noahmp%water%state%IrrigationAmtFlood(I,J)
    NoahmpIO%IRSIVOL(I,J) = NoahmpIO%IRSIVOL(I,J) + (noahmp%water%flux%IrrigationRateSprinkler(I,J)*1000.0)
    NoahmpIO%IRMIVOL(I,J) = NoahmpIO%IRMIVOL(I,J) + (noahmp%water%flux%IrrigationRateMicro(I,J)*1000.0)
    NoahmpIO%IRFIVOL(I,J) = NoahmpIO%IRFIVOL(I,J) + (noahmp%water%flux%IrrigationRateFlood(I,J)*1000.0)
    NoahmpIO%IRELOSS(I,J) = NoahmpIO%IRELOSS(I,J) + (noahmp%water%flux%EvapIrriSprinkler(I,J)*NoahmpIO%DTBL)

    ! wetland (Zhang2022)
    if ( noahmp%config%nmlist%OptWetlandModel > 0 ) then
       NoahmpIO%WSURFXY(I,J) = noahmp%water%state%WaterStorageWetland(I,J)
       NoahmpIO%FSATXY (I,J) = noahmp%water%state%SoilSaturateFrac(I,J)
    endif

#ifdef WRF_HYDRO
    NoahmpIO%infxsrt   (I,J) = max(noahmp%water%flux%RunoffSurface(I,J), 0.0)               ! mm, surface runoff
    NoahmpIO%soldrain  (I,J) = max(noahmp%water%flux%RunoffSubsurface(I,J), 0.0)            ! mm, underground runoff
    NoahmpIO%qtiledrain(I,J) = max(noahmp%water%flux%TileDrain(I,J), 0.0)                   ! mm, tile drainage
#endif
      end do
    end do
    !$acc end parallel loop
    end associate

  end subroutine WaterVarOutTransfer

end module WaterVarOutTransferMod
