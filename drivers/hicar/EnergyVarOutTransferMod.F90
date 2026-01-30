module EnergyVarOutTransferMod

!!! Transfer column (1-D) Noah-MP Energy variables to 2D NoahmpIO for output

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

  subroutine EnergyVarOutTransfer(noahmp, NoahmpIO)

    implicit none

    type(NoahmpIO_type), intent(inout) :: NoahmpIO
    type(noahmp_type),   intent(inout) :: noahmp

    ! local variables
    integer                          :: I, J
    integer                          :: LoopInd                   ! snow/soil layer loop index
    real(kind=kind_noahmp)           :: LeafAreaIndSunlit         ! sunlit leaf area index [m2/m2]
    real(kind=kind_noahmp)           :: LeafAreaIndShade          ! shaded leaf area index [m2/m2]
    real(kind=kind_noahmp)           :: ResistanceLeafBoundary    ! leaf boundary layer resistance [s/m]
    real(kind=kind_noahmp)           :: ThicknessSnowSoilLayer    ! temporary snow/soil layer thickness [m]

    !$acc parallel loop collapse(2) present(noahmp, NoahmpIO)
      do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
!-----------------------------------------------------------------------
    associate(                                                         &
              NumSoilLayer    => noahmp%config%domain%NumSoilLayer    ,&
              NumSnowLayerMax => noahmp%config%domain%NumSnowLayerMax ,&
              NumSnowLayerNeg => noahmp%config%domain%NumSnowLayerNeg(I,J) ,&
              NumSwRadBand    => noahmp%config%domain%NumSwRadBand    ,&
              IndicatorIceSfc => noahmp%config%domain%IndicatorIceSfc(I,J)  &
             )
!-----------------------------------------------------------------------

    ! special treatment for glacier point output
    if ( IndicatorIceSfc == -1 ) then ! land ice point
       noahmp%energy%state%VegFrac(I,J)             = 0.0
       noahmp%energy%state%RoughLenMomSfcToAtm(I,J) = 0.002
       noahmp%energy%flux%RadSwAbsVeg(I,J)          = 0.0
       noahmp%energy%flux%RadLwNetCanopy(I,J)       = 0.0
       noahmp%energy%flux%RadLwNetVegGrd(I,J)       = 0.0
       noahmp%energy%flux%HeatSensibleCanopy(I,J)   = 0.0
       noahmp%energy%flux%HeatSensibleVegGrd(I,J)   = 0.0
       noahmp%energy%flux%HeatLatentVegGrd(I,J)     = 0.0
       noahmp%energy%flux%HeatGroundVegGrd(I,J)     = 0.0
       noahmp%energy%flux%HeatCanStorageChg(I,J)    = 0.0
       noahmp%energy%flux%HeatLatentCanTransp(I,J)  = 0.0
       noahmp%energy%flux%HeatLatentCanEvap(I,J)    = 0.0
       noahmp%energy%flux%HeatPrecipAdvCanopy(I,J)  = 0.0
       noahmp%energy%flux%HeatPrecipAdvVegGrd(I,J)  = 0.0
       noahmp%energy%flux%HeatLatentCanopy(I,J)     = 0.0
       noahmp%energy%flux%HeatLatentTransp(I,J)     = 0.0
       noahmp%energy%flux%RadLwNetBareGrd(I,J)      = noahmp%energy%flux%RadLwNetSfc(I,J)
       noahmp%energy%flux%HeatSensibleBareGrd(I,J)  = noahmp%energy%flux%HeatSensibleSfc(I,J)
       noahmp%energy%flux%HeatLatentBareGrd(I,J)    = noahmp%energy%flux%HeatLatentGrd(I,J)
       noahmp%energy%flux%HeatGroundBareGrd(I,J)    = noahmp%energy%flux%HeatGroundTot(I,J)
       noahmp%energy%state%TemperatureGrdBare(I,J)  = noahmp%energy%state%TemperatureGrd(I,J)
       noahmp%energy%state%ExchCoeffShBare(I,J)     = noahmp%energy%state%ExchCoeffShSfc(I,J)
       NoahmpIO%LH(I,J)                        = noahmp%energy%flux%HeatLatentGrd(I,J)
    endif

    if ( IndicatorIceSfc == 0 ) then ! land soil point
       NoahmpIO%LH(I,J) = noahmp%energy%flux%HeatLatentGrd(I,J) + noahmp%energy%flux%HeatLatentCanopy(I,J) + &
                          noahmp%energy%flux%HeatLatentTransp(I,J) + noahmp%energy%flux%HeatLatentIrriEvap(I,J) 
    endif

    ! energy flux variables
    NoahmpIO%HFX        (I,J) = noahmp%energy%flux%HeatSensibleSfc(I,J)
    NoahmpIO%GRDFLX     (I,J) = noahmp%energy%flux%HeatGroundTot(I,J)
    NoahmpIO%FSAXY      (I,J) = noahmp%energy%flux%RadSwAbsSfc(I,J)
    NoahmpIO%FIRAXY     (I,J) = noahmp%energy%flux%RadLwNetSfc(I,J)
    NoahmpIO%APARXY     (I,J) = noahmp%energy%flux%RadPhotoActAbsCan(I,J)
    NoahmpIO%SAVXY      (I,J) = noahmp%energy%flux%RadSwAbsVeg(I,J)
    NoahmpIO%SAGXY      (I,J) = noahmp%energy%flux%RadSwAbsGrd(I,J)
    NoahmpIO%IRCXY      (I,J) = noahmp%energy%flux%RadLwNetCanopy(I,J)
    NoahmpIO%IRGXY      (I,J) = noahmp%energy%flux%RadLwNetVegGrd(I,J)
    NoahmpIO%SHCXY      (I,J) = noahmp%energy%flux%HeatSensibleCanopy(I,J)
    NoahmpIO%SHGXY      (I,J) = noahmp%energy%flux%HeatSensibleVegGrd(I,J)
    NoahmpIO%EVGXY      (I,J) = noahmp%energy%flux%HeatLatentVegGrd(I,J)
    NoahmpIO%GHVXY      (I,J) = noahmp%energy%flux%HeatGroundVegGrd(I,J)
    NoahmpIO%IRBXY      (I,J) = noahmp%energy%flux%RadLwNetBareGrd(I,J)
    NoahmpIO%SHBXY      (I,J) = noahmp%energy%flux%HeatSensibleBareGrd(I,J)
    NoahmpIO%EVBXY      (I,J) = noahmp%energy%flux%HeatLatentBareGrd(I,J)
    NoahmpIO%GHBXY      (I,J) = noahmp%energy%flux%HeatGroundBareGrd(I,J)
    NoahmpIO%TRXY       (I,J) = noahmp%energy%flux%HeatLatentCanTransp(I,J)
    NoahmpIO%EVCXY      (I,J) = noahmp%energy%flux%HeatLatentCanEvap(I,J)
    NoahmpIO%CANHSXY    (I,J) = noahmp%energy%flux%HeatCanStorageChg(I,J)
    NoahmpIO%PAHXY      (I,J) = noahmp%energy%flux%HeatPrecipAdvSfc(I,J)
    NoahmpIO%PAHGXY     (I,J) = noahmp%energy%flux%HeatPrecipAdvVegGrd(I,J)
    NoahmpIO%PAHVXY     (I,J) = noahmp%energy%flux%HeatPrecipAdvCanopy(I,J)
    NoahmpIO%PAHBXY     (I,J) = noahmp%energy%flux%HeatPrecipAdvBareGrd(I,J)
    NoahmpIO%ACC_SSOILXY(I,J) = noahmp%energy%flux%HeatGroundTotAcc(I,J)
    NoahmpIO%EFLXBXY    (I,J) = noahmp%energy%flux%HeatFromSoilBot(I,J)

    ! energy state variables
    NoahmpIO%TSK     (I,J) = noahmp%energy%state%TemperatureRadSfc(I,J)
    NoahmpIO%EMISS   (I,J) = noahmp%energy%state%EmissivitySfc(I,J)
    NoahmpIO%QSFC    (I,J) = noahmp%energy%state%SpecHumiditySfcMean(I,J)
    NoahmpIO%TVXY    (I,J) = noahmp%energy%state%TemperatureCanopy(I,J)
    NoahmpIO%TGXY    (I,J) = noahmp%energy%state%TemperatureGrd(I,J)
    NoahmpIO%EAHXY   (I,J) = noahmp%energy%state%PressureVaporCanAir(I,J)
    NoahmpIO%TAHXY   (I,J) = noahmp%energy%state%TemperatureCanopyAir(I,J)
    NoahmpIO%CMXY    (I,J) = noahmp%energy%state%ExchCoeffMomSfc(I,J)
    NoahmpIO%CHXY    (I,J) = noahmp%energy%state%ExchCoeffShSfc(I,J)
    NoahmpIO%ALBOLDXY(I,J) = noahmp%energy%state%AlbedoSnowPrev(I,J)
    NoahmpIO%LAI     (I,J) = noahmp%energy%state%LeafAreaIndex(I,J)
    NoahmpIO%XSAIXY  (I,J) = noahmp%energy%state%StemAreaIndex(I,J)
    NoahmpIO%TAUSSXY (I,J) = noahmp%energy%state%SnowAgeNondim(I,J)
    NoahmpIO%Z0      (I,J) = noahmp%energy%state%RoughLenMomSfcToAtm(I,J)
    NoahmpIO%T2MVXY  (I,J) = noahmp%energy%state%TemperatureAir2mVeg(I,J)
    NoahmpIO%T2MBXY  (I,J) = noahmp%energy%state%TemperatureAir2mBare(I,J)
    NoahmpIO%TRADXY  (I,J) = noahmp%energy%state%TemperatureRadSfc(I,J)
    NoahmpIO%FVEGXY  (I,J) = noahmp%energy%state%VegFrac(I,J)
    NoahmpIO%RSSUNXY (I,J) = noahmp%energy%state%ResistanceStomataSunlit(I,J)
    NoahmpIO%RSSHAXY (I,J) = noahmp%energy%state%ResistanceStomataShade(I,J)
    NoahmpIO%BGAPXY  (I,J) = noahmp%energy%state%GapBtwCanopy(I,J)
    NoahmpIO%WGAPXY  (I,J) = noahmp%energy%state%GapInCanopy(I,J)
    NoahmpIO%TGVXY   (I,J) = noahmp%energy%state%TemperatureGrdVeg(I,J)
    NoahmpIO%TGBXY   (I,J) = noahmp%energy%state%TemperatureGrdBare(I,J)
    NoahmpIO%CHVXY   (I,J) = noahmp%energy%state%ExchCoeffShAbvCan(I,J)
    NoahmpIO%CHBXY   (I,J) = noahmp%energy%state%ExchCoeffShBare(I,J)
    NoahmpIO%CHLEAFXY(I,J) = noahmp%energy%state%ExchCoeffShLeaf(I,J)
    NoahmpIO%CHUCXY  (I,J) = noahmp%energy%state%ExchCoeffShUndCan(I,J)
    NoahmpIO%CHV2XY  (I,J) = noahmp%energy%state%ExchCoeffSh2mVeg(I,J)
    NoahmpIO%CHB2XY  (I,J) = noahmp%energy%state%ExchCoeffSh2mBare(I,J)
    NoahmpIO%Q2MVXY  (I,J) = noahmp%energy%state%SpecHumidity2mVeg(I,J) /(1.0-noahmp%energy%state%SpecHumidity2mVeg(I,J))  ! spec humidity to mixing ratio
    NoahmpIO%Q2MBXY  (I,J) = noahmp%energy%state%SpecHumidity2mBare(I,J)/(1.0-noahmp%energy%state%SpecHumidity2mBare(I,J))
    NoahmpIO%ALBEDO  (I,J) = noahmp%energy%state%AlbedoSfc(I,J)
    NoahmpIO%IRRSPLH (I,J) = NoahmpIO%IRRSPLH(I,J) + &
                             (noahmp%energy%flux%HeatLatentIrriEvap(I,J) * noahmp%config%domain%MainTimeStep)
    NoahmpIO%TSLB    (I,1:NumSoilLayer,J)       = noahmp%energy%state%TemperatureSoilSnow(I,1:NumSoilLayer,J)
    NoahmpIO%TSNOXY  (I,-NumSnowLayerMax+1:0,J) = noahmp%energy%state%TemperatureSoilSnow(I,-NumSnowLayerMax+1:0,J)

    NoahmpIO%ALBSOILDIRXY(I,1:NumSwRadBand,J) = noahmp%energy%state%AlbedoSoilDir(I,1:NumSwRadBand,J)
    NoahmpIO%ALBSOILDIFXY(I,1:NumSwRadBand,J) = noahmp%energy%state%AlbedoSoilDif(I,1:NumSwRadBand,J)
    NoahmpIO%ALBSFCDIRXY (I,1:NumSwRadBand,J) = noahmp%energy%state%AlbedoSfcDir (I,1:NumSwRadBand,J)
    NoahmpIO%ALBSFCDIFXY (I,1:NumSwRadBand,J) = noahmp%energy%state%AlbedoSfcDif (I,1:NumSwRadBand,J)
    NoahmpIO%ALBSNOWDIRXY(I,1:NumSwRadBand,J) = noahmp%energy%state%AlbedoSnowDir(I,1:NumSwRadBand,J)
    NoahmpIO%ALBSNOWDIFXY(I,1:NumSwRadBand,J) = noahmp%energy%state%AlbedoSnowDif(I,1:NumSwRadBand,J)

    ! New Calculation of total Canopy/Stomatal Conductance Based on Bonan et al. (2011), Inverse of Canopy Resistance (below)
    LeafAreaIndSunlit      = max(noahmp%energy%state%LeafAreaIndSunlit(I,J), 0.0)
    LeafAreaIndShade       = max(noahmp%energy%state%LeafAreaIndShade(I,J), 0.0)
    ResistanceLeafBoundary = max(noahmp%energy%state%ResistanceLeafBoundary(I,J), 0.0)
    if ( (noahmp%energy%state%ResistanceStomataSunlit(I,J) <= 0.0) .or. (noahmp%energy%state%ResistanceStomataShade(I,J) <= 0.0) .or. &
         (LeafAreaIndSunlit == 0.0) .or. (LeafAreaIndShade == 0.0)       .or. &
         (noahmp%energy%state%ResistanceStomataSunlit(I,J) == undefined_real) .or. &
         (noahmp%energy%state%ResistanceStomataShade(I,J) == undefined_real) ) then
       NoahmpIO%RS   (I,J) = 0.0
    else
       NoahmpIO%RS   (I,J) = ((1.0 / (noahmp%energy%state%ResistanceStomataSunlit(I,J) + ResistanceLeafBoundary) * &
                              noahmp%energy%state%LeafAreaIndSunlit(I,J)) + &
                             ((1.0 / (noahmp%energy%state%ResistanceStomataShade(I,J) + ResistanceLeafBoundary)) * &
                              noahmp%energy%state%LeafAreaIndShade(I,J)))
       NoahmpIO%RS   (I,J) = 1.0 / NoahmpIO%RS (I,J) ! Resistance
    endif

    ! calculation of snow and soil energy storage
    NoahmpIO%SNOWENERGY(I,J) = 0.0
    NoahmpIO%SOILENERGY(I,J) = 0.0
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       if ( LoopInd == NumSnowLayerNeg+1 ) then
          ThicknessSnowSoilLayer = -noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd,J)
       else
          ThicknessSnowSoilLayer = noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd-1,J) - &
                                   noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd,J)
       endif
       if ( LoopInd >= 1 ) then
          NoahmpIO%SOILENERGY(I,J) = NoahmpIO%SOILENERGY(I,J) + ThicknessSnowSoilLayer * &
                                     noahmp%energy%state%HeatCapacSoilSnow(I,LoopInd,J) * &
                                     (noahmp%energy%state%TemperatureSoilSnow(I,LoopInd,J) - 273.16) * 0.001
       else
          NoahmpIO%SNOWENERGY(I,J) = NoahmpIO%SNOWENERGY(I,J) + ThicknessSnowSoilLayer * &
                                     noahmp%energy%state%HeatCapacSoilSnow(I,LoopInd,J) * &
                                     (noahmp%energy%state%TemperatureSoilSnow(I,LoopInd,J) - 273.16) * 0.001
       endif
    enddo

    end associate

      end do
      end do
      !$acc end parallel loop
  end subroutine EnergyVarOutTransfer

end module EnergyVarOutTransferMod
