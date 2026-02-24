module EnergyVarInitMod

!!! Initialize column (1-D) Noah-MP energy variables
!!! Energy variables should be first defined in EnergyVarType.F90

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
  subroutine EnergyVarInitDefault(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer :: I, J      ! grid indices
    integer :: LoopInd, k   ! loop index

    ! Domain and layer bounds for allocations
    associate(                                                                    &
              ITS                   => noahmp%config%domain%ITS                  ,&
              ITE                   => noahmp%config%domain%ITE                  ,&
              JTS                   => noahmp%config%domain%JTS                  ,&
              JTE                   => noahmp%config%domain%JTE                  ,&
              NumSnowLayerMax       => noahmp%config%domain%NumSnowLayerMax      ,&
              NumSoilLayer          => noahmp%config%domain%NumSoilLayer         ,&
              NumSwRadBand          => noahmp%config%domain%NumSwRadBand         ,&
              NumSnicarRadBand      => noahmp%config%domain%NumSnicarRadBand     ,&
              NumRadiusSnwMieSnicar => noahmp%config%domain%NumRadiusSnwMieSnicar &
             )

    ! Allocate 2D energy flux arrays and transfer to GPU
    if ( .not. allocated(noahmp%energy%flux%HeatLatentCanopy) ) then
       allocate( noahmp%energy%flux%HeatLatentCanopy(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentTransp) ) then
       allocate( noahmp%energy%flux%HeatLatentTransp(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentGrd) ) then
       allocate( noahmp%energy%flux%HeatLatentGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentIrriEvap) ) then
       allocate( noahmp%energy%flux%HeatLatentIrriEvap(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatPrecipAdvCanopy) ) then
       allocate( noahmp%energy%flux%HeatPrecipAdvCanopy(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatPrecipAdvVegGrd) ) then
       allocate( noahmp%energy%flux%HeatPrecipAdvVegGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatPrecipAdvBareGrd) ) then
       allocate( noahmp%energy%flux%HeatPrecipAdvBareGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatPrecipAdvSfc) ) then
       allocate( noahmp%energy%flux%HeatPrecipAdvSfc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatSensibleCanopy) ) then
       allocate( noahmp%energy%flux%HeatSensibleCanopy(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentCanEvap) ) then
       allocate( noahmp%energy%flux%HeatLatentCanEvap(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatSensibleVegGrd) ) then
       allocate( noahmp%energy%flux%HeatSensibleVegGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatSensibleSfc) ) then
       allocate( noahmp%energy%flux%HeatSensibleSfc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentVegGrd) ) then
       allocate( noahmp%energy%flux%HeatLatentVegGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentCanTransp) ) then
       allocate( noahmp%energy%flux%HeatLatentCanTransp(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatGroundVegGrd) ) then
       allocate( noahmp%energy%flux%HeatGroundVegGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatSensibleBareGrd) ) then
       allocate( noahmp%energy%flux%HeatSensibleBareGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentBareGrd) ) then
       allocate( noahmp%energy%flux%HeatLatentBareGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatGroundBareGrd) ) then
       allocate( noahmp%energy%flux%HeatGroundBareGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatGroundTot) ) then
       allocate( noahmp%energy%flux%HeatGroundTot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatGroundTotMean) ) then
       allocate( noahmp%energy%flux%HeatGroundTotMean(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatFromSoilBot) ) then
       allocate( noahmp%energy%flux%HeatFromSoilBot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatCanStorageChg) ) then
       allocate( noahmp%energy%flux%HeatCanStorageChg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatGroundTotAcc) ) then
       allocate( noahmp%energy%flux%HeatGroundTotAcc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadPhotoActAbsSunlit) ) then
       allocate( noahmp%energy%flux%RadPhotoActAbsSunlit(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadPhotoActAbsShade) ) then
       allocate( noahmp%energy%flux%RadPhotoActAbsShade(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwAbsVeg) ) then
       allocate( noahmp%energy%flux%RadSwAbsVeg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwAbsGrd) ) then
       allocate( noahmp%energy%flux%RadSwAbsGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwAbsSfc) ) then
       allocate( noahmp%energy%flux%RadSwAbsSfc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflSfc) ) then
       allocate( noahmp%energy%flux%RadSwReflSfc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflVeg) ) then
       allocate( noahmp%energy%flux%RadSwReflVeg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflGrd) ) then
       allocate( noahmp%energy%flux%RadSwReflGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadLwNetCanopy) ) then
       allocate( noahmp%energy%flux%RadLwNetCanopy(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadLwNetSfc) ) then
       allocate( noahmp%energy%flux%RadLwNetSfc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadPhotoActAbsCan) ) then
       allocate( noahmp%energy%flux%RadPhotoActAbsCan(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadLwEmitSfc) ) then
       allocate( noahmp%energy%flux%RadLwEmitSfc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadLwNetVegGrd) ) then
       allocate( noahmp%energy%flux%RadLwNetVegGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadLwNetBareGrd) ) then
       allocate( noahmp%energy%flux%RadLwNetBareGrd(ITS:ITE,JTS:JTE) )
    endif

    ! Allocate 3D energy state arrays and transfer to GPU

    if (.not. allocated(noahmp%energy%state%FlagFrozenCanopy) ) then
       allocate( noahmp%energy%state%FlagFrozenCanopy(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%FlagFrozenGround) ) then
       allocate( noahmp%energy%state%FlagFrozenGround(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%LeafAreaIndEff) ) then
       allocate( noahmp%energy%state%LeafAreaIndEff(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%StemAreaIndEff) ) then
       allocate( noahmp%energy%state%StemAreaIndEff(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%LeafAreaIndex) ) then
       allocate( noahmp%energy%state%LeafAreaIndex(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%StemAreaIndex) ) then
       allocate( noahmp%energy%state%StemAreaIndex(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%VegAreaIndEff) ) then
       allocate( noahmp%energy%state%VegAreaIndEff(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%VegFrac) ) then
       allocate( noahmp%energy%state%VegFrac(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureGrd) ) then
       allocate( noahmp%energy%state%TemperatureGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureCanopy) ) then
       allocate( noahmp%energy%state%TemperatureCanopy(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureSfc) ) then
       allocate( noahmp%energy%state%TemperatureSfc(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureRootZone) ) then
       allocate( noahmp%energy%state%TemperatureRootZone(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%PressureVaporRefHeight) ) then
       allocate( noahmp%energy%state%PressureVaporRefHeight(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%SnowAgeFac) ) then
       allocate( noahmp%energy%state%SnowAgeFac(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%SnowAgeNondim) ) then
       allocate( noahmp%energy%state%SnowAgeNondim(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%AlbedoSnowPrev) ) then
       allocate( noahmp%energy%state%AlbedoSnowPrev(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%VegAreaProjDir) ) then
       allocate( noahmp%energy%state%VegAreaProjDir(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%GapBtwCanopy) ) then
       allocate( noahmp%energy%state%GapBtwCanopy(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%GapInCanopy) ) then
       allocate( noahmp%energy%state%GapInCanopy(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%GapCanopyDif) ) then
       allocate( noahmp%energy%state%GapCanopyDif(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%GapCanopyDir) ) then
       allocate( noahmp%energy%state%GapCanopyDir(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%CanopySunlitFrac) ) then
       allocate( noahmp%energy%state%CanopySunlitFrac(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%CanopyShadeFrac) ) then
       allocate( noahmp%energy%state%CanopyShadeFrac(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%LeafAreaIndSunlit) ) then
       allocate( noahmp%energy%state%LeafAreaIndSunlit(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%LeafAreaIndShade) ) then
       allocate( noahmp%energy%state%LeafAreaIndShade(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatCanopy) ) then
       allocate( noahmp%energy%state%VapPresSatCanopy(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatGrdVeg) ) then
       allocate( noahmp%energy%state%VapPresSatGrdVeg(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatGrdBare) ) then
       allocate( noahmp%energy%state%VapPresSatGrdBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatCanTempD) ) then
       allocate( noahmp%energy%state%VapPresSatCanTempD(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatGrdVegTempD) ) then
       allocate( noahmp%energy%state%VapPresSatGrdVegTempD(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatGrdBareTempD) ) then
       allocate( noahmp%energy%state%VapPresSatGrdBareTempD(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%PressureVaporCanAir) ) then
       allocate( noahmp%energy%state%PressureVaporCanAir(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%PressureAtmosCO2) ) then
       allocate( noahmp%energy%state%PressureAtmosCO2(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%PressureAtmosO2) ) then
       allocate( noahmp%energy%state%PressureAtmosO2(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceStomataSunlit) ) then
       allocate( noahmp%energy%state%ResistanceStomataSunlit(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceStomataShade) ) then
       allocate( noahmp%energy%state%ResistanceStomataShade(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%DensityAirRefHeight) ) then
       allocate( noahmp%energy%state%DensityAirRefHeight(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureCanopyAir) ) then
       allocate( noahmp%energy%state%TemperatureCanopyAir(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ZeroPlaneDispSfc) ) then
       allocate( noahmp%energy%state%ZeroPlaneDispSfc(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ZeroPlaneDispGrd) ) then
       allocate( noahmp%energy%state%ZeroPlaneDispGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenMomGrd) ) then
       allocate( noahmp%energy%state%RoughLenMomGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenMomSfc) ) then
       allocate( noahmp%energy%state%RoughLenMomSfc(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenShCanopy) ) then
       allocate( noahmp%energy%state%RoughLenShCanopy(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenShVegGrd) ) then
       allocate( noahmp%energy%state%RoughLenShVegGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenShBareGrd) ) then
       allocate( noahmp%energy%state%RoughLenShBareGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%CanopyHeight) ) then
       allocate( noahmp%energy%state%CanopyHeight(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%WindSpdCanopyTop) ) then
       allocate( noahmp%energy%state%WindSpdCanopyTop(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%FrictionVelVeg) ) then
       allocate( noahmp%energy%state%FrictionVelVeg(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%FrictionVelBare) ) then
       allocate( noahmp%energy%state%FrictionVelBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%WindExtCoeffCanopy) ) then
       allocate( noahmp%energy%state%WindExtCoeffCanopy(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabParaUndCan) ) then
       allocate( noahmp%energy%state%MoStabParaUndCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabParaAbvCan) ) then
       allocate( noahmp%energy%state%MoStabParaAbvCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabParaBare) ) then
       allocate( noahmp%energy%state%MoStabParaBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabParaVeg2m) ) then
       allocate( noahmp%energy%state%MoStabParaVeg2m(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabParaBare2m) ) then
       allocate( noahmp%energy%state%MoStabParaBare2m(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoLengthUndCan) ) then
       allocate( noahmp%energy%state%MoLengthUndCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoLengthAbvCan) ) then
       allocate( noahmp%energy%state%MoLengthAbvCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoLengthBare) ) then
       allocate( noahmp%energy%state%MoLengthBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrShUndCan) ) then
       allocate( noahmp%energy%state%MoStabCorrShUndCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrMomAbvCan) ) then
       allocate( noahmp%energy%state%MoStabCorrMomAbvCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrShAbvCan) ) then
       allocate( noahmp%energy%state%MoStabCorrShAbvCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrMomVeg2m) ) then
       allocate( noahmp%energy%state%MoStabCorrMomVeg2m(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrShVeg2m) ) then
       allocate( noahmp%energy%state%MoStabCorrShVeg2m(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrShBare) ) then
       allocate( noahmp%energy%state%MoStabCorrShBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrMomBare) ) then
       allocate( noahmp%energy%state%MoStabCorrMomBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrMomBare2m) ) then
       allocate( noahmp%energy%state%MoStabCorrMomBare2m(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrShBare2m) ) then
       allocate( noahmp%energy%state%MoStabCorrShBare2m(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffMomSfc) ) then
       allocate( noahmp%energy%state%ExchCoeffMomSfc(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffMomAbvCan) ) then
       allocate( noahmp%energy%state%ExchCoeffMomAbvCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffMomBare) ) then
       allocate( noahmp%energy%state%ExchCoeffMomBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffShSfc) ) then
       allocate( noahmp%energy%state%ExchCoeffShSfc(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffShAbvCan) ) then
       allocate( noahmp%energy%state%ExchCoeffShAbvCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffShBare) ) then
       allocate( noahmp%energy%state%ExchCoeffShBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffSh2mVegMo) ) then
       allocate( noahmp%energy%state%ExchCoeffSh2mVegMo(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffSh2mBareMo) ) then
       allocate( noahmp%energy%state%ExchCoeffSh2mBareMo(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffSh2mVeg) ) then
       allocate( noahmp%energy%state%ExchCoeffSh2mVeg(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffLhAbvCan) ) then
       allocate( noahmp%energy%state%ExchCoeffLhAbvCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffLhTransp) ) then
       allocate( noahmp%energy%state%ExchCoeffLhTransp(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffLhEvap) ) then
       allocate( noahmp%energy%state%ExchCoeffLhEvap(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffLhUndCan) ) then
       allocate( noahmp%energy%state%ExchCoeffLhUndCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceMomUndCan) ) then
       allocate( noahmp%energy%state%ResistanceMomUndCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceShUndCan) ) then
       allocate( noahmp%energy%state%ResistanceShUndCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceLhUndCan) ) then
       allocate( noahmp%energy%state%ResistanceLhUndCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceMomAbvCan) ) then
       allocate( noahmp%energy%state%ResistanceMomAbvCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceShAbvCan) ) then
       allocate( noahmp%energy%state%ResistanceShAbvCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceLhAbvCan) ) then
       allocate( noahmp%energy%state%ResistanceLhAbvCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceMomBareGrd) ) then
       allocate( noahmp%energy%state%ResistanceMomBareGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceShBareGrd) ) then
       allocate( noahmp%energy%state%ResistanceShBareGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceLhBareGrd) ) then
       allocate( noahmp%energy%state%ResistanceLhBareGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceLeafBoundary) ) then
       allocate( noahmp%energy%state%ResistanceLeafBoundary(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperaturePotRefHeight) ) then
       allocate( noahmp%energy%state%TemperaturePotRefHeight(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%WindSpdRefHeight) ) then
       allocate( noahmp%energy%state%WindSpdRefHeight(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%FrictionVelVertVeg) ) then
       allocate( noahmp%energy%state%FrictionVelVertVeg(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%FrictionVelVertBare) ) then
       allocate( noahmp%energy%state%FrictionVelVertBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%EmissivityVeg) ) then
       allocate( noahmp%energy%state%EmissivityVeg(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%EmissivityGrd) ) then
       allocate( noahmp%energy%state%EmissivityGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceGrdEvap) ) then
       allocate( noahmp%energy%state%ResistanceGrdEvap(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%PsychConstCanopy) ) then
       allocate( noahmp%energy%state%PsychConstCanopy(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%LatHeatVapCanopy) ) then
       allocate( noahmp%energy%state%LatHeatVapCanopy(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%PsychConstGrd) ) then
       allocate( noahmp%energy%state%PsychConstGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%LatHeatVapGrd) ) then
       allocate( noahmp%energy%state%LatHeatVapGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%RelHumidityGrd) ) then
       allocate( noahmp%energy%state%RelHumidityGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%SpecHumiditySfc) ) then
       allocate( noahmp%energy%state%SpecHumiditySfc(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%SpecHumiditySfcMean) ) then
       allocate( noahmp%energy%state%SpecHumiditySfcMean(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%SpecHumidity2mVeg) ) then
       allocate( noahmp%energy%state%SpecHumidity2mVeg(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%SpecHumidity2mBare) ) then
       allocate( noahmp%energy%state%SpecHumidity2mBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%SpecHumidity2m) ) then
       allocate( noahmp%energy%state%SpecHumidity2m(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureGrdVeg) ) then
       allocate( noahmp%energy%state%TemperatureGrdVeg(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureGrdBare) ) then
       allocate( noahmp%energy%state%TemperatureGrdBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressEwVeg) ) then
       allocate( noahmp%energy%state%WindStressEwVeg(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressNsVeg) ) then
       allocate( noahmp%energy%state%WindStressNsVeg(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressEwBare) ) then
       allocate( noahmp%energy%state%WindStressEwBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressNsBare) ) then
       allocate( noahmp%energy%state%WindStressNsBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressEwSfc) ) then
       allocate( noahmp%energy%state%WindStressEwSfc(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressNsSfc) ) then
       allocate( noahmp%energy%state%WindStressNsSfc(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureAir2mVeg) ) then
       allocate( noahmp%energy%state%TemperatureAir2mVeg(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureAir2mBare) ) then
       allocate( noahmp%energy%state%TemperatureAir2mBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureAir2m) ) then
       allocate( noahmp%energy%state%TemperatureAir2m(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffShLeaf) ) then
       allocate( noahmp%energy%state%ExchCoeffShLeaf(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffShUndCan) ) then
       allocate( noahmp%energy%state%ExchCoeffShUndCan(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffSh2mBare) ) then
       allocate( noahmp%energy%state%ExchCoeffSh2mBare(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%RefHeightAboveGrd) ) then
       allocate( noahmp%energy%state%RefHeightAboveGrd(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%CanopyFracSnowBury) ) then
       allocate( noahmp%energy%state%CanopyFracSnowBury(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%DepthSoilTempBotToSno) ) then
       allocate( noahmp%energy%state%DepthSoilTempBotToSno(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenMomSfcToAtm) ) then
       allocate( noahmp%energy%state%RoughLenMomSfcToAtm(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureRadSfc) ) then
       allocate( noahmp%energy%state%TemperatureRadSfc(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%EmissivitySfc) ) then
       allocate( noahmp%energy%state%EmissivitySfc(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%AlbedoSfc) ) then
       allocate( noahmp%energy%state%AlbedoSfc(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%EnergyBalanceError) ) then
       allocate( noahmp%energy%state%EnergyBalanceError(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%RadSwBalanceError) ) then
       allocate( noahmp%energy%state%RadSwBalanceError(ITS:ITE,JTS:JTE) )
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureSoilSnow) ) then
       allocate( noahmp%energy%state%TemperatureSoilSnow(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%ThermConductSoilSnow) ) then
       allocate( noahmp%energy%state%ThermConductSoilSnow(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%HeatCapacSoilSnow) ) then
       allocate( noahmp%energy%state%HeatCapacSoilSnow(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%PhaseChgFacSoilSnow) ) then
       allocate( noahmp%energy%state%PhaseChgFacSoilSnow(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%HeatCapacVolSnow) ) then
       allocate( noahmp%energy%state%HeatCapacVolSnow(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%ThermConductSnow) ) then
       allocate( noahmp%energy%state%ThermConductSnow(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%HeatCapacVolSoil) ) then
       allocate( noahmp%energy%state%HeatCapacVolSoil(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%ThermConductSoil) ) then
       allocate( noahmp%energy%state%ThermConductSoil(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%HeatCapacGlaIce) ) then
       allocate( noahmp%energy%state%HeatCapacGlaIce(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%ThermConductGlaIce) ) then
       allocate( noahmp%energy%state%ThermConductGlaIce(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSnowDir) ) then
       allocate( noahmp%energy%state%AlbedoSnowDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSnowDif) ) then
       allocate( noahmp%energy%state%AlbedoSnowDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSoilDir) ) then
       allocate( noahmp%energy%state%AlbedoSoilDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSoilDif) ) then
       allocate( noahmp%energy%state%AlbedoSoilDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoGrdDir) ) then
       allocate( noahmp%energy%state%AlbedoGrdDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoGrdDif) ) then
       allocate( noahmp%energy%state%AlbedoGrdDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%ReflectanceVeg) ) then
       allocate( noahmp%energy%state%ReflectanceVeg(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%TransmittanceVeg) ) then
       allocate( noahmp%energy%state%TransmittanceVeg(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSfcDir) ) then
       allocate( noahmp%energy%state%AlbedoSfcDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSfcDif) ) then
       allocate( noahmp%energy%state%AlbedoSfcDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif

    ! Allocate 3D energy flux arrays and transfer to GPU
    if ( .not. allocated(noahmp%energy%flux%RadSwAbsVegDir) ) then
       allocate( noahmp%energy%flux%RadSwAbsVegDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwAbsVegDif) ) then
       allocate( noahmp%energy%flux%RadSwAbsVegDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDirTranGrdDir) ) then
       allocate( noahmp%energy%flux%RadSwDirTranGrdDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDirTranGrdDif) ) then
       allocate( noahmp%energy%flux%RadSwDirTranGrdDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDifTranGrdDir) ) then
       allocate( noahmp%energy%flux%RadSwDifTranGrdDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDifTranGrdDif) ) then
       allocate( noahmp%energy%flux%RadSwDifTranGrdDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflVegDir) ) then
       allocate( noahmp%energy%flux%RadSwReflVegDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflVegDif) ) then
       allocate( noahmp%energy%flux%RadSwReflVegDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflGrdDir) ) then
       allocate( noahmp%energy%flux%RadSwReflGrdDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflGrdDif) ) then
       allocate( noahmp%energy%flux%RadSwReflGrdDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDownDir) ) then
       allocate( noahmp%energy%flux%RadSwDownDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDownDif) ) then
       allocate( noahmp%energy%flux%RadSwDownDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwPenetrateGrd) ) then
       allocate( noahmp%energy%flux%RadSwPenetrateGrd(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
    endif

    ! SNICAR flux arrays
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( .not. allocated(noahmp%energy%flux%FracRadSwAbsSnowDir) ) then
          allocate( noahmp%energy%flux%FracRadSwAbsSnowDir(ITS:ITE,-NumSnowLayerMax+1:1,1:NumSwRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%flux%FracRadSwAbsSnowDif) ) then
          allocate( noahmp%energy%flux%FracRadSwAbsSnowDif(ITS:ITE,-NumSnowLayerMax+1:1,1:NumSwRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%flux%RadSwAbsSnowSoilLayer) ) then
          allocate( noahmp%energy%flux%RadSwAbsSnowSoilLayer(ITS:ITE,-NumSnowLayerMax+1:1,JTS:JTE) )
       endif
    endif

    ! Allocate 3D energy parameter arrays and transfer to GPU
    if ( .not. allocated(noahmp%energy%param%LeafAreaIndexMon) ) then
       allocate( noahmp%energy%param%LeafAreaIndexMon(ITS:ITE,1:12,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%StemAreaIndexMon) ) then
       allocate( noahmp%energy%param%StemAreaIndexMon(ITS:ITE,1:12,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SoilQuartzFrac) ) then
       allocate( noahmp%energy%param%SoilQuartzFrac(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%AlbedoSoilSat) ) then
       allocate( noahmp%energy%param%AlbedoSoilSat(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%AlbedoSoilDry) ) then
       allocate( noahmp%energy%param%AlbedoSoilDry(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%AlbedoLakeFrz) ) then
       allocate( noahmp%energy%param%AlbedoLakeFrz(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%ScatterCoeffSnow) ) then
       allocate( noahmp%energy%param%ScatterCoeffSnow(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%ReflectanceLeaf) ) then
       allocate( noahmp%energy%param%ReflectanceLeaf(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%ReflectanceStem) ) then
       allocate( noahmp%energy%param%ReflectanceStem(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%TransmittanceLeaf) ) then
       allocate( noahmp%energy%param%TransmittanceLeaf(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%TransmittanceStem) ) then
       allocate( noahmp%energy%param%TransmittanceStem(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%EmissivitySoilLake) ) then
       allocate( noahmp%energy%param%EmissivitySoilLake(ITS:ITE,1:2,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%AlbedoLandIce) ) then
       allocate( noahmp%energy%param%AlbedoLandIce(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
    endif

    ! SNICAR parameter arrays - these are lookup tables, not spatially varying
    ! Keep as 1D/2D since they are the same for all grid points
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( .not. allocated(noahmp%energy%param%RadSwWgtDir) ) then
          allocate( noahmp%energy%param%RadSwWgtDir(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%RadSwWgtDif) ) then
          allocate( noahmp%energy%param%RadSwWgtDif(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbSnwRadDir) ) then
          allocate( noahmp%energy%param%SsAlbSnwRadDir(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmSnwRadDir) ) then
          allocate( noahmp%energy%param%AsyPrmSnwRadDir(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassSnwRadDir) ) then
          allocate( noahmp%energy%param%ExtCffMassSnwRadDir(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbSnwRadDif) ) then
          allocate( noahmp%energy%param%SsAlbSnwRadDif(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmSnwRadDif) ) then
          allocate( noahmp%energy%param%AsyPrmSnwRadDif(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassSnwRadDif) ) then
          allocate( noahmp%energy%param%ExtCffMassSnwRadDif(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbBCphi) ) then
          allocate( noahmp%energy%param%SsAlbBCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmBCphi) ) then
          allocate( noahmp%energy%param%AsyPrmBCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassBCphi) ) then
          allocate( noahmp%energy%param%ExtCffMassBCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbBCpho) ) then
          allocate( noahmp%energy%param%SsAlbBCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmBCpho) ) then
          allocate( noahmp%energy%param%AsyPrmBCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassBCpho) ) then
          allocate( noahmp%energy%param%ExtCffMassBCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbOCphi) ) then
          allocate( noahmp%energy%param%SsAlbOCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmOCphi) ) then
          allocate( noahmp%energy%param%AsyPrmOCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassOCphi) ) then
          allocate( noahmp%energy%param%ExtCffMassOCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbOCpho) ) then
          allocate( noahmp%energy%param%SsAlbOCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmOCpho) ) then
          allocate( noahmp%energy%param%AsyPrmOCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassOCpho) ) then
          allocate( noahmp%energy%param%ExtCffMassOCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbDustB1) ) then
          allocate( noahmp%energy%param%SsAlbDustB1(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmDustB1) ) then
          allocate( noahmp%energy%param%AsyPrmDustB1(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassDustB1) ) then
          allocate( noahmp%energy%param%ExtCffMassDustB1(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbDustB2) ) then
          allocate( noahmp%energy%param%SsAlbDustB2(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmDustB2) ) then
          allocate( noahmp%energy%param%AsyPrmDustB2(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassDustB2) ) then
          allocate( noahmp%energy%param%ExtCffMassDustB2(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbDustB3) ) then
          allocate( noahmp%energy%param%SsAlbDustB3(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmDustB3) ) then
          allocate( noahmp%energy%param%AsyPrmDustB3(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassDustB3) ) then
          allocate( noahmp%energy%param%ExtCffMassDustB3(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbDustB4) ) then
          allocate( noahmp%energy%param%SsAlbDustB4(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmDustB4) ) then
          allocate( noahmp%energy%param%AsyPrmDustB4(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassDustB4) ) then
          allocate( noahmp%energy%param%ExtCffMassDustB4(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbDustB5) ) then
          allocate( noahmp%energy%param%SsAlbDustB5(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmDustB5) ) then
          allocate( noahmp%energy%param%AsyPrmDustB5(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassDustB5) ) then
          allocate( noahmp%energy%param%ExtCffMassDustB5(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
       endif
    endif

    ! Allocate 2D energy param arrays and transfer to GPU
    if ( .not. allocated(noahmp%energy%param%TreeCrownRadius) ) then
       allocate( noahmp%energy%param%TreeCrownRadius(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%HeightCanopyTop) ) then
       allocate( noahmp%energy%param%HeightCanopyTop(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%HeightCanopyBot) ) then
       allocate( noahmp%energy%param%HeightCanopyBot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%RoughLenMomVeg) ) then
       allocate( noahmp%energy%param%RoughLenMomVeg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%TreeDensity) ) then
       allocate( noahmp%energy%param%TreeDensity(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%CanopyOrientIndex) ) then
       allocate( noahmp%energy%param%CanopyOrientIndex(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%UpscatterCoeffSnowDir) ) then
       allocate( noahmp%energy%param%UpscatterCoeffSnowDir(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%UpscatterCoeffSnowDif) ) then
       allocate( noahmp%energy%param%UpscatterCoeffSnowDif(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SoilHeatCapacity) ) then
       allocate( noahmp%energy%param%SoilHeatCapacity(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SnowAgeFacBats) ) then
       allocate( noahmp%energy%param%SnowAgeFacBats(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SnowGrowVapFacBats) ) then
       allocate( noahmp%energy%param%SnowGrowVapFacBats(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SnowSootFacBats) ) then
       allocate( noahmp%energy%param%SnowSootFacBats(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SnowGrowFrzFacBats) ) then
       allocate( noahmp%energy%param%SnowGrowFrzFacBats(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SolarZenithAdjBats) ) then
       allocate( noahmp%energy%param%SolarZenithAdjBats(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%FreshSnoAlbVisBats) ) then
       allocate( noahmp%energy%param%FreshSnoAlbVisBats(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%FreshSnoAlbNirBats) ) then
       allocate( noahmp%energy%param%FreshSnoAlbNirBats(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SnoAgeFacDifVisBats) ) then
       allocate( noahmp%energy%param%SnoAgeFacDifVisBats(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SnoAgeFacDifNirBats) ) then
       allocate( noahmp%energy%param%SnoAgeFacDifNirBats(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SzaFacDirVisBats) ) then
       allocate( noahmp%energy%param%SzaFacDirVisBats(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SzaFacDirNirBats) ) then
       allocate( noahmp%energy%param%SzaFacDirNirBats(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SnowAlbRefClass) ) then
       allocate( noahmp%energy%param%SnowAlbRefClass(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SnowAgeFacClass) ) then
       allocate( noahmp%energy%param%SnowAgeFacClass(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%SnowAlbFreshClass) ) then
       allocate( noahmp%energy%param%SnowAlbFreshClass(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%ConductanceLeafMin) ) then
       allocate( noahmp%energy%param%ConductanceLeafMin(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%Co2MmConst25C) ) then
       allocate( noahmp%energy%param%Co2MmConst25C(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%O2MmConst25C) ) then
       allocate( noahmp%energy%param%O2MmConst25C(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%Co2MmConstQ10) ) then
       allocate( noahmp%energy%param%Co2MmConstQ10(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%O2MmConstQ10) ) then
       allocate( noahmp%energy%param%O2MmConstQ10(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%RadiationStressFac) ) then
       allocate( noahmp%energy%param%RadiationStressFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%ResistanceStomataMin) ) then
       allocate( noahmp%energy%param%ResistanceStomataMin(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%ResistanceStomataMax) ) then
       allocate( noahmp%energy%param%ResistanceStomataMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%AirTempOptimTransp) ) then
       allocate( noahmp%energy%param%AirTempOptimTransp(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%VaporPresDeficitFac) ) then
       allocate( noahmp%energy%param%VaporPresDeficitFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%LeafDimLength) ) then
       allocate( noahmp%energy%param%LeafDimLength(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%ZilitinkevichCoeff) ) then
       allocate( noahmp%energy%param%ZilitinkevichCoeff(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%EmissivitySnow) ) then
       allocate( noahmp%energy%param%EmissivitySnow(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%CanopyWindExtFac) ) then
       allocate( noahmp%energy%param%CanopyWindExtFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%RoughLenMomSnow) ) then
       allocate( noahmp%energy%param%RoughLenMomSnow(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%RoughLenMomSoil) ) then
       allocate( noahmp%energy%param%RoughLenMomSoil(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%RoughLenMomLake) ) then
       allocate( noahmp%energy%param%RoughLenMomLake(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%EmissivityIceSfc) ) then
       allocate( noahmp%energy%param%EmissivityIceSfc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%ResistanceSoilExp) ) then
       allocate( noahmp%energy%param%ResistanceSoilExp(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%ResistanceSnowSfc) ) then
       allocate( noahmp%energy%param%ResistanceSnowSfc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%VegFracGreen) ) then
       allocate( noahmp%energy%param%VegFracGreen(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%VegFracAnnMax) ) then
       allocate( noahmp%energy%param%VegFracAnnMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%energy%param%HeatCapacCanFac) ) then
       allocate( noahmp%energy%param%HeatCapacCanFac(ITS:ITE,JTS:JTE) )
    endif

    end associate

    ! Now initialize all 2D and 3D arrays in parallel loop

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%energy%flux%HeatLatentCanopy,   &
    !$acc   noahmp%energy%flux%HeatLatentTransp,   &
    !$acc   noahmp%energy%flux%HeatLatentGrd,   &
    !$acc   noahmp%energy%flux%HeatLatentIrriEvap,   &
    !$acc   noahmp%energy%flux%HeatPrecipAdvCanopy,   &
    !$acc   noahmp%energy%flux%HeatPrecipAdvVegGrd,   &
    !$acc   noahmp%energy%flux%HeatPrecipAdvBareGrd,   &
    !$acc   noahmp%energy%flux%HeatPrecipAdvSfc,   &
    !$acc   noahmp%energy%flux%HeatSensibleCanopy,   &
    !$acc   noahmp%energy%flux%HeatLatentCanEvap,   &
    !$acc   noahmp%energy%flux%HeatSensibleVegGrd,   &
    !$acc   noahmp%energy%flux%HeatSensibleSfc,   &
    !$acc   noahmp%energy%flux%HeatLatentVegGrd,   &
    !$acc   noahmp%energy%flux%HeatLatentCanTransp,   &
    !$acc   noahmp%energy%flux%HeatGroundVegGrd,   &
    !$acc   noahmp%energy%flux%HeatSensibleBareGrd,   &
    !$acc   noahmp%energy%flux%HeatLatentBareGrd,   &
    !$acc   noahmp%energy%flux%HeatGroundBareGrd,   &
    !$acc   noahmp%energy%flux%HeatGroundTot,   &
    !$acc   noahmp%energy%flux%HeatGroundTotMean,   &
    !$acc   noahmp%energy%flux%HeatFromSoilBot,   &
    !$acc   noahmp%energy%flux%HeatCanStorageChg,   &
    !$acc   noahmp%energy%flux%HeatGroundTotAcc,   &
    !$acc   noahmp%energy%flux%RadPhotoActAbsSunlit,   &
    !$acc   noahmp%energy%flux%RadPhotoActAbsShade,   &
    !$acc   noahmp%energy%flux%RadSwAbsVeg,   &
    !$acc   noahmp%energy%flux%RadSwAbsGrd,   &
    !$acc   noahmp%energy%flux%RadSwAbsSfc,   &
    !$acc   noahmp%energy%flux%RadSwReflSfc,   &
    !$acc   noahmp%energy%flux%RadSwReflVeg    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%energy%flux%RadSwReflGrd,   &
    !$acc   noahmp%energy%flux%RadLwNetCanopy,   &
    !$acc   noahmp%energy%flux%RadLwNetSfc,   &
    !$acc   noahmp%energy%flux%RadPhotoActAbsCan,   &
    !$acc   noahmp%energy%flux%RadLwEmitSfc,   &
    !$acc   noahmp%energy%flux%RadLwNetVegGrd,   &
    !$acc   noahmp%energy%flux%RadLwNetBareGrd,   &
    !$acc   noahmp%energy%state%FlagFrozenCanopy,   &
    !$acc   noahmp%energy%state%FlagFrozenGround,   &
    !$acc   noahmp%energy%state%LeafAreaIndEff,   &
    !$acc   noahmp%energy%state%StemAreaIndEff,   &
    !$acc   noahmp%energy%state%LeafAreaIndex,   &
    !$acc   noahmp%energy%state%StemAreaIndex,   &
    !$acc   noahmp%energy%state%VegAreaIndEff,   &
    !$acc   noahmp%energy%state%VegFrac,   &
    !$acc   noahmp%energy%state%TemperatureGrd,   &
    !$acc   noahmp%energy%state%TemperatureCanopy,   &
    !$acc   noahmp%energy%state%TemperatureSfc,   &
    !$acc   noahmp%energy%state%TemperatureRootZone,   &
    !$acc   noahmp%energy%state%PressureVaporRefHeight,   &
    !$acc   noahmp%energy%state%SnowAgeFac,   &
    !$acc   noahmp%energy%state%SnowAgeNondim,   &
    !$acc   noahmp%energy%state%AlbedoSnowPrev,   &
    !$acc   noahmp%energy%state%VegAreaProjDir,   &
    !$acc   noahmp%energy%state%GapBtwCanopy,   &
    !$acc   noahmp%energy%state%GapInCanopy,   &
    !$acc   noahmp%energy%state%GapCanopyDif,   &
    !$acc   noahmp%energy%state%GapCanopyDir,   &
    !$acc   noahmp%energy%state%CanopySunlitFrac,   &
    !$acc   noahmp%energy%state%CanopyShadeFrac    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%energy%state%LeafAreaIndSunlit,   &
    !$acc   noahmp%energy%state%LeafAreaIndShade,   &
    !$acc   noahmp%energy%state%VapPresSatCanopy,   &
    !$acc   noahmp%energy%state%VapPresSatGrdVeg,   &
    !$acc   noahmp%energy%state%VapPresSatGrdBare,   &
    !$acc   noahmp%energy%state%VapPresSatCanTempD,   &
    !$acc   noahmp%energy%state%VapPresSatGrdVegTempD,   &
    !$acc   noahmp%energy%state%VapPresSatGrdBareTempD,   &
    !$acc   noahmp%energy%state%PressureVaporCanAir,   &
    !$acc   noahmp%energy%state%PressureAtmosCO2,   &
    !$acc   noahmp%energy%state%PressureAtmosO2,   &
    !$acc   noahmp%energy%state%ResistanceStomataSunlit,   &
    !$acc   noahmp%energy%state%ResistanceStomataShade,   &
    !$acc   noahmp%energy%state%DensityAirRefHeight,   &
    !$acc   noahmp%energy%state%TemperatureCanopyAir,   &
    !$acc   noahmp%energy%state%ZeroPlaneDispSfc,   &
    !$acc   noahmp%energy%state%ZeroPlaneDispGrd,   &
    !$acc   noahmp%energy%state%RoughLenMomGrd,   &
    !$acc   noahmp%energy%state%RoughLenMomSfc,   &
    !$acc   noahmp%energy%state%RoughLenShCanopy,   &
    !$acc   noahmp%energy%state%RoughLenShVegGrd,   &
    !$acc   noahmp%energy%state%RoughLenShBareGrd,   &
    !$acc   noahmp%energy%state%CanopyHeight,   &
    !$acc   noahmp%energy%state%WindSpdCanopyTop,   &
    !$acc   noahmp%energy%state%FrictionVelVeg,   &
    !$acc   noahmp%energy%state%FrictionVelBare,   &
    !$acc   noahmp%energy%state%WindExtCoeffCanopy,   &
    !$acc   noahmp%energy%state%MoStabParaUndCan,   &
    !$acc   noahmp%energy%state%MoStabParaAbvCan,   &
    !$acc   noahmp%energy%state%MoStabParaBare    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%energy%state%MoStabParaVeg2m,   &
    !$acc   noahmp%energy%state%MoStabParaBare2m,   &
    !$acc   noahmp%energy%state%MoLengthUndCan,   &
    !$acc   noahmp%energy%state%MoLengthAbvCan,   &
    !$acc   noahmp%energy%state%MoLengthBare,   &
    !$acc   noahmp%energy%state%MoStabCorrShUndCan,   &
    !$acc   noahmp%energy%state%MoStabCorrMomAbvCan,   &
    !$acc   noahmp%energy%state%MoStabCorrShAbvCan,   &
    !$acc   noahmp%energy%state%MoStabCorrMomVeg2m,   &
    !$acc   noahmp%energy%state%MoStabCorrShVeg2m,   &
    !$acc   noahmp%energy%state%MoStabCorrShBare,   &
    !$acc   noahmp%energy%state%MoStabCorrMomBare,   &
    !$acc   noahmp%energy%state%MoStabCorrMomBare2m,   &
    !$acc   noahmp%energy%state%MoStabCorrShBare2m,   &
    !$acc   noahmp%energy%state%ExchCoeffMomSfc,   &
    !$acc   noahmp%energy%state%ExchCoeffMomAbvCan,   &
    !$acc   noahmp%energy%state%ExchCoeffMomBare,   &
    !$acc   noahmp%energy%state%ExchCoeffShSfc,   &
    !$acc   noahmp%energy%state%ExchCoeffShAbvCan,   &
    !$acc   noahmp%energy%state%ExchCoeffShBare,   &
    !$acc   noahmp%energy%state%ExchCoeffSh2mVegMo,   &
    !$acc   noahmp%energy%state%ExchCoeffSh2mBareMo,   &
    !$acc   noahmp%energy%state%ExchCoeffSh2mVeg,   &
    !$acc   noahmp%energy%state%ExchCoeffLhAbvCan,   &
    !$acc   noahmp%energy%state%ExchCoeffLhTransp,   &
    !$acc   noahmp%energy%state%ExchCoeffLhEvap,   &
    !$acc   noahmp%energy%state%ExchCoeffLhUndCan,   &
    !$acc   noahmp%energy%state%ResistanceMomUndCan,   &
    !$acc   noahmp%energy%state%ResistanceShUndCan,   &
    !$acc   noahmp%energy%state%ResistanceLhUndCan    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%energy%state%ResistanceMomAbvCan,   &
    !$acc   noahmp%energy%state%ResistanceShAbvCan,   &
    !$acc   noahmp%energy%state%ResistanceLhAbvCan,   &
    !$acc   noahmp%energy%state%ResistanceMomBareGrd,   &
    !$acc   noahmp%energy%state%ResistanceShBareGrd,   &
    !$acc   noahmp%energy%state%ResistanceLhBareGrd,   &
    !$acc   noahmp%energy%state%ResistanceLeafBoundary,   &
    !$acc   noahmp%energy%state%TemperaturePotRefHeight,   &
    !$acc   noahmp%energy%state%WindSpdRefHeight,   &
    !$acc   noahmp%energy%state%FrictionVelVertVeg,   &
    !$acc   noahmp%energy%state%FrictionVelVertBare,   &
    !$acc   noahmp%energy%state%EmissivityVeg,   &
    !$acc   noahmp%energy%state%EmissivityGrd,   &
    !$acc   noahmp%energy%state%ResistanceGrdEvap,   &
    !$acc   noahmp%energy%state%PsychConstCanopy,   &
    !$acc   noahmp%energy%state%LatHeatVapCanopy,   &
    !$acc   noahmp%energy%state%PsychConstGrd,   &
    !$acc   noahmp%energy%state%LatHeatVapGrd,   &
    !$acc   noahmp%energy%state%RelHumidityGrd,   &
    !$acc   noahmp%energy%state%SpecHumiditySfc,   &
    !$acc   noahmp%energy%state%SpecHumiditySfcMean,   &
    !$acc   noahmp%energy%state%SpecHumidity2mVeg,   &
    !$acc   noahmp%energy%state%SpecHumidity2mBare,   &
    !$acc   noahmp%energy%state%SpecHumidity2m,   &
    !$acc   noahmp%energy%state%TemperatureGrdVeg,   &
    !$acc   noahmp%energy%state%TemperatureGrdBare,   &
    !$acc   noahmp%energy%state%WindStressEwVeg,   &
    !$acc   noahmp%energy%state%WindStressNsVeg,   &
    !$acc   noahmp%energy%state%WindStressEwBare,   &
    !$acc   noahmp%energy%state%WindStressNsBare    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%energy%state%WindStressEwSfc,   &
    !$acc   noahmp%energy%state%WindStressNsSfc,   &
    !$acc   noahmp%energy%state%TemperatureAir2mVeg,   &
    !$acc   noahmp%energy%state%TemperatureAir2mBare,   &
    !$acc   noahmp%energy%state%TemperatureAir2m,   &
    !$acc   noahmp%energy%state%ExchCoeffShLeaf,   &
    !$acc   noahmp%energy%state%ExchCoeffShUndCan,   &
    !$acc   noahmp%energy%state%ExchCoeffSh2mBare,   &
    !$acc   noahmp%energy%state%RefHeightAboveGrd,   &
    !$acc   noahmp%energy%state%CanopyFracSnowBury,   &
    !$acc   noahmp%energy%state%DepthSoilTempBotToSno,   &
    !$acc   noahmp%energy%state%RoughLenMomSfcToAtm,   &
    !$acc   noahmp%energy%state%TemperatureRadSfc,   &
    !$acc   noahmp%energy%state%EmissivitySfc,   &
    !$acc   noahmp%energy%state%AlbedoSfc,   &
    !$acc   noahmp%energy%state%EnergyBalanceError,   &
    !$acc   noahmp%energy%state%RadSwBalanceError,   &
    !$acc   noahmp%energy%state%TemperatureSoilSnow,   &
    !$acc   noahmp%energy%state%ThermConductSoilSnow,   &
    !$acc   noahmp%energy%state%HeatCapacSoilSnow,   &
    !$acc   noahmp%energy%state%PhaseChgFacSoilSnow,   &
    !$acc   noahmp%energy%state%HeatCapacVolSnow,   &
    !$acc   noahmp%energy%state%ThermConductSnow,   &
    !$acc   noahmp%energy%state%HeatCapacVolSoil,   &
    !$acc   noahmp%energy%state%ThermConductSoil,   &
    !$acc   noahmp%energy%state%HeatCapacGlaIce,   &
    !$acc   noahmp%energy%state%ThermConductGlaIce,   &
    !$acc   noahmp%energy%state%AlbedoSnowDir,   &
    !$acc   noahmp%energy%state%AlbedoSnowDif,   &
    !$acc   noahmp%energy%state%AlbedoSoilDir    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%energy%state%AlbedoSoilDif,   &
    !$acc   noahmp%energy%state%AlbedoGrdDir,   &
    !$acc   noahmp%energy%state%AlbedoGrdDif,   &
    !$acc   noahmp%energy%state%ReflectanceVeg,   &
    !$acc   noahmp%energy%state%TransmittanceVeg,   &
    !$acc   noahmp%energy%state%AlbedoSfcDir,   &
    !$acc   noahmp%energy%state%AlbedoSfcDif,   &
    !$acc   noahmp%energy%flux%RadSwAbsVegDir,   &
    !$acc   noahmp%energy%flux%RadSwAbsVegDif,   &
    !$acc   noahmp%energy%flux%RadSwDirTranGrdDir,   &
    !$acc   noahmp%energy%flux%RadSwDirTranGrdDif,   &
    !$acc   noahmp%energy%flux%RadSwDifTranGrdDir,   &
    !$acc   noahmp%energy%flux%RadSwDifTranGrdDif,   &
    !$acc   noahmp%energy%flux%RadSwReflVegDir,   &
    !$acc   noahmp%energy%flux%RadSwReflVegDif,   &
    !$acc   noahmp%energy%flux%RadSwReflGrdDir,   &
    !$acc   noahmp%energy%flux%RadSwReflGrdDif,   &
    !$acc   noahmp%energy%flux%RadSwDownDir,   &
    !$acc   noahmp%energy%flux%RadSwDownDif,   &
    !$acc   noahmp%energy%flux%RadSwPenetrateGrd,   &
    !$acc   noahmp%energy%param%LeafAreaIndexMon,   &
    !$acc   noahmp%energy%param%StemAreaIndexMon,   &
    !$acc   noahmp%energy%param%SoilQuartzFrac,   &
    !$acc   noahmp%energy%param%AlbedoSoilSat,   &
    !$acc   noahmp%energy%param%AlbedoSoilDry,   &
    !$acc   noahmp%energy%param%AlbedoLakeFrz,   &
    !$acc   noahmp%energy%param%ScatterCoeffSnow,   &
    !$acc   noahmp%energy%param%ReflectanceLeaf,   &
    !$acc   noahmp%energy%param%ReflectanceStem,   &
    !$acc   noahmp%energy%param%TransmittanceLeaf    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%energy%param%TransmittanceStem,   &
    !$acc   noahmp%energy%param%EmissivitySoilLake,   &
    !$acc   noahmp%energy%param%AlbedoLandIce,   &
    !$acc   noahmp%energy%param%TreeCrownRadius,   &
    !$acc   noahmp%energy%param%HeightCanopyTop,   &
    !$acc   noahmp%energy%param%HeightCanopyBot,   &
    !$acc   noahmp%energy%param%RoughLenMomVeg,   &
    !$acc   noahmp%energy%param%TreeDensity,   &
    !$acc   noahmp%energy%param%CanopyOrientIndex,   &
    !$acc   noahmp%energy%param%UpscatterCoeffSnowDir,   &
    !$acc   noahmp%energy%param%UpscatterCoeffSnowDif,   &
    !$acc   noahmp%energy%param%SoilHeatCapacity,   &
    !$acc   noahmp%energy%param%SnowAgeFacBats,   &
    !$acc   noahmp%energy%param%SnowGrowVapFacBats,   &
    !$acc   noahmp%energy%param%SnowSootFacBats,   &
    !$acc   noahmp%energy%param%SnowGrowFrzFacBats,   &
    !$acc   noahmp%energy%param%SolarZenithAdjBats,   &
    !$acc   noahmp%energy%param%FreshSnoAlbVisBats,   &
    !$acc   noahmp%energy%param%FreshSnoAlbNirBats,   &
    !$acc   noahmp%energy%param%SnoAgeFacDifVisBats,   &
    !$acc   noahmp%energy%param%SnoAgeFacDifNirBats,   &
    !$acc   noahmp%energy%param%SzaFacDirVisBats,   &
    !$acc   noahmp%energy%param%SzaFacDirNirBats,   &
    !$acc   noahmp%energy%param%SnowAlbRefClass,   &
    !$acc   noahmp%energy%param%SnowAgeFacClass,   &
    !$acc   noahmp%energy%param%SnowAlbFreshClass,   &
    !$acc   noahmp%energy%param%ConductanceLeafMin,   &
    !$acc   noahmp%energy%param%Co2MmConst25C,   &
    !$acc   noahmp%energy%param%O2MmConst25C,   &
    !$acc   noahmp%energy%param%Co2MmConstQ10    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%energy%param%O2MmConstQ10,   &
    !$acc   noahmp%energy%param%RadiationStressFac,   &
    !$acc   noahmp%energy%param%ResistanceStomataMin,   &
    !$acc   noahmp%energy%param%ResistanceStomataMax,   &
    !$acc   noahmp%energy%param%AirTempOptimTransp,   &
    !$acc   noahmp%energy%param%VaporPresDeficitFac,   &
    !$acc   noahmp%energy%param%LeafDimLength,   &
    !$acc   noahmp%energy%param%ZilitinkevichCoeff,   &
    !$acc   noahmp%energy%param%EmissivitySnow,   &
    !$acc   noahmp%energy%param%CanopyWindExtFac,   &
    !$acc   noahmp%energy%param%RoughLenMomSnow,   &
    !$acc   noahmp%energy%param%RoughLenMomSoil,   &
    !$acc   noahmp%energy%param%RoughLenMomLake,   &
    !$acc   noahmp%energy%param%EmissivityIceSfc,   &
    !$acc   noahmp%energy%param%ResistanceSoilExp,   &
    !$acc   noahmp%energy%param%ResistanceSnowSfc,   &
    !$acc   noahmp%energy%param%VegFracGreen,   &
    !$acc   noahmp%energy%param%VegFracAnnMax,   &
    !$acc   noahmp%energy%param%HeatCapacCanFac    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
    !$acc enter data create(              &
    !$acc   noahmp%energy%flux%FracRadSwAbsSnowDir,   &
    !$acc   noahmp%energy%flux%FracRadSwAbsSnowDif,   &
    !$acc   noahmp%energy%flux%RadSwAbsSnowSoilLayer,   &
    !$acc   noahmp%energy%param%RadSwWgtDir,   &
    !$acc   noahmp%energy%param%RadSwWgtDif,   &
    !$acc   noahmp%energy%param%SsAlbSnwRadDir,   &
    !$acc   noahmp%energy%param%AsyPrmSnwRadDir,   &
    !$acc   noahmp%energy%param%ExtCffMassSnwRadDir,   &
    !$acc   noahmp%energy%param%SsAlbSnwRadDif,   &
    !$acc   noahmp%energy%param%AsyPrmSnwRadDif,   &
    !$acc   noahmp%energy%param%ExtCffMassSnwRadDif,   &
    !$acc   noahmp%energy%param%SsAlbBCphi,   &
    !$acc   noahmp%energy%param%AsyPrmBCphi,   &
    !$acc   noahmp%energy%param%ExtCffMassBCphi,   &
    !$acc   noahmp%energy%param%SsAlbBCpho,   &
    !$acc   noahmp%energy%param%AsyPrmBCpho,   &
    !$acc   noahmp%energy%param%ExtCffMassBCpho,   &
    !$acc   noahmp%energy%param%SsAlbOCphi,   &
    !$acc   noahmp%energy%param%AsyPrmOCphi,   &
    !$acc   noahmp%energy%param%ExtCffMassOCphi,   &
    !$acc   noahmp%energy%param%SsAlbOCpho,   &
    !$acc   noahmp%energy%param%AsyPrmOCpho,   &
    !$acc   noahmp%energy%param%ExtCffMassOCpho,   &
    !$acc   noahmp%energy%param%SsAlbDustB1,   &
    !$acc   noahmp%energy%param%AsyPrmDustB1,   &
    !$acc   noahmp%energy%param%ExtCffMassDustB1,   &
    !$acc   noahmp%energy%param%SsAlbDustB2,   &
    !$acc   noahmp%energy%param%AsyPrmDustB2,   &
    !$acc   noahmp%energy%param%ExtCffMassDustB2,   &
    !$acc   noahmp%energy%param%SsAlbDustB3    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)
    !$acc enter data create(              &
    !$acc   noahmp%energy%param%AsyPrmDustB3,   &
    !$acc   noahmp%energy%param%ExtCffMassDustB3,   &
    !$acc   noahmp%energy%param%SsAlbDustB4,   &
    !$acc   noahmp%energy%param%AsyPrmDustB4,   &
    !$acc   noahmp%energy%param%ExtCffMassDustB4,   &
    !$acc   noahmp%energy%param%SsAlbDustB5,   &
    !$acc   noahmp%energy%param%AsyPrmDustB5,   &
    !$acc   noahmp%energy%param%ExtCffMassDustB5    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)
    endif

    !$acc wait(NOAHMP_ACC_QUEUE)

    ! SNICAR lookup table initialization (moved here after device creates)
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       !$acc parallel loop gang vector collapse(2) default(present) private(k, LoopInd)
       do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
        !$acc loop seq
        do LoopInd = 1, noahmp%config%domain%NumSnicarRadBand
        !$acc loop seq
         do k = 1, noahmp%config%domain%NumRadiusSnwMieSnicar
           noahmp%energy%param%SsAlbSnwRadDir       (I,k,LoopInd,J) = undefined_real
           noahmp%energy%param%AsyPrmSnwRadDir      (I,k,LoopInd,J) = undefined_real
           noahmp%energy%param%ExtCffMassSnwRadDir  (I,k,LoopInd,J) = undefined_real
           noahmp%energy%param%SsAlbSnwRadDif       (I,k,LoopInd,J) = undefined_real
           noahmp%energy%param%AsyPrmSnwRadDif      (I,k,LoopInd,J) = undefined_real
           noahmp%energy%param%ExtCffMassSnwRadDif  (I,k,LoopInd,J) = undefined_real
         end do
         noahmp%energy%param%RadSwWgtDir           (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%RadSwWgtDif           (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%SsAlbBCphi            (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%AsyPrmBCphi           (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%ExtCffMassBCphi       (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%SsAlbBCpho            (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%AsyPrmBCpho           (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%ExtCffMassBCpho       (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%SsAlbOCphi            (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%AsyPrmOCphi           (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%ExtCffMassOCphi       (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%SsAlbOCpho            (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%AsyPrmOCpho           (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%ExtCffMassOCpho       (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%SsAlbDustB1           (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%AsyPrmDustB1          (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%ExtCffMassDustB1      (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%SsAlbDustB2           (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%AsyPrmDustB2          (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%ExtCffMassDustB2      (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%SsAlbDustB3           (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%AsyPrmDustB3          (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%ExtCffMassDustB3      (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%SsAlbDustB4           (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%AsyPrmDustB4          (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%ExtCffMassDustB4      (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%SsAlbDustB5           (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%AsyPrmDustB5          (I,LoopInd,J)       = undefined_real
         noahmp%energy%param%ExtCffMassDustB5      (I,LoopInd,J)       = undefined_real
         end do
      enddo
    enddo
    endif

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        ! Initialize 2D energy state scalars
        noahmp%energy%state%FlagFrozenCanopy(I,J)        = .false.
        noahmp%energy%state%FlagFrozenGround(I,J)        = .false.
        noahmp%energy%state%LeafAreaIndEff(I,J)          = undefined_real
        noahmp%energy%state%StemAreaIndEff(I,J)          = undefined_real
        noahmp%energy%state%LeafAreaIndex(I,J)           = undefined_real
        noahmp%energy%state%StemAreaIndex(I,J)           = undefined_real
        noahmp%energy%state%VegAreaIndEff(I,J)           = undefined_real
        noahmp%energy%state%VegFrac(I,J)                 = undefined_real
        noahmp%energy%state%PressureVaporRefHeight(I,J)  = undefined_real
        noahmp%energy%state%SnowAgeFac(I,J)              = undefined_real
        noahmp%energy%state%SnowAgeNondim(I,J)           = undefined_real
        noahmp%energy%state%AlbedoSnowPrev(I,J)          = undefined_real
        noahmp%energy%state%VegAreaProjDir(I,J)          = undefined_real
        noahmp%energy%state%GapBtwCanopy(I,J)            = undefined_real
        noahmp%energy%state%GapInCanopy(I,J)             = undefined_real
        noahmp%energy%state%GapCanopyDif(I,J)            = undefined_real
        noahmp%energy%state%GapCanopyDir(I,J)            = undefined_real
        noahmp%energy%state%CanopySunlitFrac(I,J)        = undefined_real
        noahmp%energy%state%CanopyShadeFrac(I,J)         = undefined_real
        noahmp%energy%state%LeafAreaIndSunlit(I,J)       = undefined_real
        noahmp%energy%state%LeafAreaIndShade(I,J)        = undefined_real
        noahmp%energy%state%VapPresSatCanopy(I,J)        = undefined_real
        noahmp%energy%state%VapPresSatGrdVeg(I,J)        = undefined_real
        noahmp%energy%state%VapPresSatGrdBare(I,J)       = undefined_real
        noahmp%energy%state%VapPresSatCanTempD(I,J)      = undefined_real
        noahmp%energy%state%VapPresSatGrdVegTempD(I,J)   = undefined_real
        noahmp%energy%state%VapPresSatGrdBareTempD(I,J)  = undefined_real
        noahmp%energy%state%PressureVaporCanAir(I,J)     = undefined_real
        noahmp%energy%state%PressureAtmosCO2(I,J)        = undefined_real
        noahmp%energy%state%PressureAtmosO2(I,J)         = undefined_real
        noahmp%energy%state%ResistanceStomataSunlit(I,J) = undefined_real
        noahmp%energy%state%ResistanceStomataShade(I,J)  = undefined_real
        noahmp%energy%state%DensityAirRefHeight(I,J)     = undefined_real
        noahmp%energy%state%TemperatureCanopyAir(I,J)    = undefined_real
        noahmp%energy%state%ZeroPlaneDispSfc(I,J)        = undefined_real
        noahmp%energy%state%ZeroPlaneDispGrd(I,J)        = undefined_real
        noahmp%energy%state%RoughLenMomGrd(I,J)          = undefined_real
        noahmp%energy%state%RoughLenMomSfc(I,J)          = undefined_real
        noahmp%energy%state%CanopyHeight(I,J)            = undefined_real
        noahmp%energy%state%WindSpdCanopyTop(I,J)        = undefined_real
        noahmp%energy%state%RoughLenShCanopy(I,J)        = undefined_real
        noahmp%energy%state%RoughLenShVegGrd(I,J)        = undefined_real
        noahmp%energy%state%RoughLenShBareGrd(I,J)       = undefined_real
        noahmp%energy%state%FrictionVelVeg(I,J)          = undefined_real
        noahmp%energy%state%FrictionVelBare(I,J)         = undefined_real
        noahmp%energy%state%WindExtCoeffCanopy(I,J)      = undefined_real
        noahmp%energy%state%MoStabParaUndCan(I,J)        = undefined_real
        noahmp%energy%state%MoStabParaAbvCan(I,J)        = undefined_real
        noahmp%energy%state%MoStabParaBare(I,J)          = undefined_real
        noahmp%energy%state%MoStabParaVeg2m(I,J)         = undefined_real
        noahmp%energy%state%MoStabParaBare2m(I,J)        = undefined_real
        noahmp%energy%state%MoLengthUndCan(I,J)          = undefined_real
        noahmp%energy%state%MoLengthAbvCan(I,J)          = undefined_real
        noahmp%energy%state%MoLengthBare(I,J)            = undefined_real
        noahmp%energy%state%MoStabCorrShUndCan(I,J)      = undefined_real
        noahmp%energy%state%MoStabCorrMomAbvCan(I,J)     = undefined_real
        noahmp%energy%state%MoStabCorrShAbvCan(I,J)      = undefined_real
        noahmp%energy%state%MoStabCorrMomVeg2m(I,J)      = undefined_real
        noahmp%energy%state%MoStabCorrShVeg2m(I,J)       = undefined_real
        noahmp%energy%state%MoStabCorrShBare(I,J)        = undefined_real
        noahmp%energy%state%MoStabCorrMomBare(I,J)       = undefined_real
        noahmp%energy%state%MoStabCorrMomBare2m(I,J)     = undefined_real
        noahmp%energy%state%MoStabCorrShBare2m(I,J)      = undefined_real
        noahmp%energy%state%ExchCoeffMomSfc(I,J)         = undefined_real
        noahmp%energy%state%ExchCoeffMomAbvCan(I,J)      = undefined_real
        noahmp%energy%state%ExchCoeffMomBare(I,J)        = undefined_real
        noahmp%energy%state%ExchCoeffShSfc(I,J)          = undefined_real
        noahmp%energy%state%ExchCoeffShBare(I,J)         = undefined_real
        noahmp%energy%state%ExchCoeffShAbvCan(I,J)       = undefined_real
        noahmp%energy%state%ExchCoeffShLeaf(I,J)         = undefined_real
        noahmp%energy%state%ExchCoeffShUndCan(I,J)       = undefined_real
        noahmp%energy%state%ExchCoeffSh2mVegMo(I,J)      = undefined_real
        noahmp%energy%state%ExchCoeffSh2mBareMo(I,J)     = undefined_real
        noahmp%energy%state%ExchCoeffSh2mVeg(I,J)        = undefined_real
        noahmp%energy%state%ExchCoeffSh2mBare(I,J)       = undefined_real
        noahmp%energy%state%ExchCoeffLhAbvCan(I,J)       = undefined_real
        noahmp%energy%state%ExchCoeffLhTransp(I,J)       = undefined_real
        noahmp%energy%state%ExchCoeffLhEvap(I,J)         = undefined_real
        noahmp%energy%state%ExchCoeffLhUndCan(I,J)       = undefined_real
        noahmp%energy%state%ResistanceMomUndCan(I,J)     = undefined_real
        noahmp%energy%state%ResistanceShUndCan(I,J)      = undefined_real
        noahmp%energy%state%ResistanceLhUndCan(I,J)      = undefined_real
        noahmp%energy%state%ResistanceMomAbvCan(I,J)     = undefined_real
        noahmp%energy%state%ResistanceShAbvCan(I,J)      = undefined_real
        noahmp%energy%state%ResistanceLhAbvCan(I,J)      = undefined_real
        noahmp%energy%state%ResistanceMomBareGrd(I,J)    = undefined_real
        noahmp%energy%state%ResistanceShBareGrd(I,J)     = undefined_real
        noahmp%energy%state%ResistanceLhBareGrd(I,J)     = undefined_real
        noahmp%energy%state%ResistanceLeafBoundary(I,J)  = undefined_real
        noahmp%energy%state%TemperaturePotRefHeight(I,J) = undefined_real
        noahmp%energy%state%WindSpdRefHeight(I,J)        = undefined_real
        noahmp%energy%state%FrictionVelVertVeg(I,J)      = undefined_real
        noahmp%energy%state%FrictionVelVertBare(I,J)     = undefined_real
        noahmp%energy%state%EmissivityVeg(I,J)           = undefined_real
        noahmp%energy%state%EmissivityGrd(I,J)           = undefined_real
        noahmp%energy%state%ResistanceGrdEvap(I,J)       = undefined_real
        noahmp%energy%state%PsychConstCanopy(I,J)        = undefined_real
        noahmp%energy%state%LatHeatVapCanopy(I,J)        = undefined_real
        noahmp%energy%state%PsychConstGrd(I,J)           = undefined_real
        noahmp%energy%state%LatHeatVapGrd(I,J)           = undefined_real
        noahmp%energy%state%RelHumidityGrd(I,J)          = undefined_real
        noahmp%energy%state%SpecHumiditySfcMean(I,J)     = undefined_real
        noahmp%energy%state%SpecHumiditySfc(I,J)         = undefined_real
        noahmp%energy%state%SpecHumidity2mVeg(I,J)       = undefined_real
        noahmp%energy%state%SpecHumidity2mBare(I,J)      = undefined_real
        noahmp%energy%state%SpecHumidity2m(I,J)          = undefined_real
        noahmp%energy%state%TemperatureSfc(I,J)          = undefined_real
        noahmp%energy%state%TemperatureGrd(I,J)          = undefined_real
        noahmp%energy%state%TemperatureCanopy(I,J)       = undefined_real
        noahmp%energy%state%TemperatureGrdVeg(I,J)       = undefined_real
        noahmp%energy%state%TemperatureGrdBare(I,J)      = undefined_real
        noahmp%energy%state%TemperatureRootZone(I,J)     = undefined_real
        noahmp%energy%state%WindStressEwVeg(I,J)         = undefined_real
        noahmp%energy%state%WindStressNsVeg(I,J)         = undefined_real
        noahmp%energy%state%WindStressEwBare(I,J)        = undefined_real
        noahmp%energy%state%WindStressNsBare(I,J)        = undefined_real
        noahmp%energy%state%WindStressEwSfc(I,J)         = undefined_real
        noahmp%energy%state%WindStressNsSfc(I,J)         = undefined_real
        noahmp%energy%state%TemperatureAir2mVeg(I,J)     = undefined_real
        noahmp%energy%state%TemperatureAir2mBare(I,J)    = undefined_real
        noahmp%energy%state%TemperatureAir2m(I,J)        = undefined_real
        noahmp%energy%state%CanopyFracSnowBury(I,J)      = undefined_real
        noahmp%energy%state%DepthSoilTempBotToSno(I,J)   = undefined_real
        noahmp%energy%state%RoughLenMomSfcToAtm(I,J)     = undefined_real
        noahmp%energy%state%TemperatureRadSfc(I,J)       = undefined_real
        noahmp%energy%state%EmissivitySfc(I,J)           = undefined_real
        noahmp%energy%state%AlbedoSfc(I,J)               = undefined_real
        noahmp%energy%state%EnergyBalanceError(I,J)      = undefined_real
        noahmp%energy%state%RadSwBalanceError(I,J)       = undefined_real
        noahmp%energy%state%RefHeightAboveGrd(I,J)       = undefined_real

        ! Initialize 2D energy flux scalars
        noahmp%energy%flux%HeatLatentCanopy(I,J)         = undefined_real
        noahmp%energy%flux%HeatLatentTransp(I,J)         = undefined_real
        noahmp%energy%flux%HeatLatentGrd(I,J)            = undefined_real
        noahmp%energy%flux%HeatPrecipAdvCanopy(I,J)      = undefined_real
        noahmp%energy%flux%HeatPrecipAdvVegGrd(I,J)      = undefined_real
        noahmp%energy%flux%HeatPrecipAdvBareGrd(I,J)     = undefined_real
        noahmp%energy%flux%HeatPrecipAdvSfc(I,J)         = undefined_real
        noahmp%energy%flux%RadPhotoActAbsSunlit(I,J)     = undefined_real
        noahmp%energy%flux%RadPhotoActAbsShade(I,J)      = undefined_real
        noahmp%energy%flux%RadSwAbsVeg(I,J)              = undefined_real
        noahmp%energy%flux%RadSwAbsGrd(I,J)              = undefined_real
        noahmp%energy%flux%RadSwAbsSfc(I,J)              = undefined_real
        noahmp%energy%flux%RadSwReflSfc(I,J)             = undefined_real
        noahmp%energy%flux%RadSwReflVeg(I,J)             = undefined_real
        noahmp%energy%flux%RadSwReflGrd(I,J)             = undefined_real
        noahmp%energy%flux%RadLwNetCanopy(I,J)           = undefined_real
        noahmp%energy%flux%HeatSensibleCanopy(I,J)       = undefined_real
        noahmp%energy%flux%HeatLatentCanEvap(I,J)        = undefined_real
        noahmp%energy%flux%RadLwNetVegGrd(I,J)           = undefined_real
        noahmp%energy%flux%HeatSensibleVegGrd(I,J)       = undefined_real
        noahmp%energy%flux%HeatLatentVegGrd(I,J)         = undefined_real
        noahmp%energy%flux%HeatLatentCanTransp(I,J)      = undefined_real
        noahmp%energy%flux%HeatGroundVegGrd(I,J)         = undefined_real
        noahmp%energy%flux%RadLwNetBareGrd(I,J)          = undefined_real
        noahmp%energy%flux%HeatSensibleBareGrd(I,J)      = undefined_real
        noahmp%energy%flux%HeatLatentBareGrd(I,J)        = undefined_real
        noahmp%energy%flux%HeatGroundBareGrd(I,J)        = undefined_real
        noahmp%energy%flux%HeatGroundTot(I,J)            = undefined_real
        noahmp%energy%flux%HeatFromSoilBot(I,J)          = undefined_real
        noahmp%energy%flux%RadLwNetSfc(I,J)              = undefined_real
        noahmp%energy%flux%HeatSensibleSfc(I,J)          = undefined_real
        noahmp%energy%flux%RadPhotoActAbsCan(I,J)        = undefined_real
        noahmp%energy%flux%RadLwEmitSfc(I,J)             = undefined_real
        noahmp%energy%flux%HeatCanStorageChg(I,J)        = undefined_real
        noahmp%energy%flux%HeatGroundTotAcc(I,J)         = undefined_real
        noahmp%energy%flux%HeatGroundTotMean(I,J)        = undefined_real
        noahmp%energy%flux%HeatLatentIrriEvap(I,J)       = 0.0

        ! Initialize 2D energy parameter scalars
        noahmp%energy%param%TreeCrownRadius(I,J)         = undefined_real
        noahmp%energy%param%HeightCanopyTop(I,J)         = undefined_real
        noahmp%energy%param%HeightCanopyBot(I,J)         = undefined_real
        noahmp%energy%param%RoughLenMomVeg(I,J)          = undefined_real
        noahmp%energy%param%TreeDensity(I,J)             = undefined_real
        noahmp%energy%param%CanopyOrientIndex(I,J)       = undefined_real
        noahmp%energy%param%UpscatterCoeffSnowDir(I,J)   = undefined_real
        noahmp%energy%param%UpscatterCoeffSnowDif(I,J)   = undefined_real
        noahmp%energy%param%SoilHeatCapacity(I,J)        = undefined_real
        noahmp%energy%param%SnowAgeFacBats(I,J)          = undefined_real
        noahmp%energy%param%SnowGrowVapFacBats(I,J)      = undefined_real
        noahmp%energy%param%SnowSootFacBats(I,J)         = undefined_real
        noahmp%energy%param%SnowGrowFrzFacBats(I,J)      = undefined_real
        noahmp%energy%param%SolarZenithAdjBats(I,J)      = undefined_real
        noahmp%energy%param%FreshSnoAlbVisBats(I,J)      = undefined_real
        noahmp%energy%param%FreshSnoAlbNirBats(I,J)      = undefined_real
        noahmp%energy%param%SnoAgeFacDifVisBats(I,J)     = undefined_real
        noahmp%energy%param%SnoAgeFacDifNirBats(I,J)     = undefined_real
        noahmp%energy%param%SzaFacDirVisBats(I,J)        = undefined_real
        noahmp%energy%param%SzaFacDirNirBats(I,J)        = undefined_real
        noahmp%energy%param%SnowAlbRefClass(I,J)         = undefined_real
        noahmp%energy%param%SnowAgeFacClass(I,J)         = undefined_real
        noahmp%energy%param%SnowAlbFreshClass(I,J)       = undefined_real
        noahmp%energy%param%ConductanceLeafMin(I,J)      = undefined_real
        noahmp%energy%param%Co2MmConst25C(I,J)           = undefined_real
        noahmp%energy%param%O2MmConst25C(I,J)            = undefined_real
        noahmp%energy%param%Co2MmConstQ10(I,J)           = undefined_real
        noahmp%energy%param%O2MmConstQ10(I,J)            = undefined_real
        noahmp%energy%param%RadiationStressFac(I,J)      = undefined_real
        noahmp%energy%param%ResistanceStomataMin(I,J)    = undefined_real
        noahmp%energy%param%ResistanceStomataMax(I,J)    = undefined_real
        noahmp%energy%param%AirTempOptimTransp(I,J)      = undefined_real
        noahmp%energy%param%VaporPresDeficitFac(I,J)     = undefined_real
        noahmp%energy%param%LeafDimLength(I,J)           = undefined_real
        noahmp%energy%param%ZilitinkevichCoeff(I,J)      = undefined_real
        noahmp%energy%param%EmissivitySnow(I,J)          = undefined_real
        noahmp%energy%param%CanopyWindExtFac(I,J)        = undefined_real
        noahmp%energy%param%RoughLenMomSnow(I,J)         = undefined_real
        noahmp%energy%param%RoughLenMomSoil(I,J)         = undefined_real
        noahmp%energy%param%RoughLenMomLake(I,J)         = undefined_real
        noahmp%energy%param%EmissivityIceSfc(I,J)        = undefined_real
        noahmp%energy%param%ResistanceSoilExp(I,J)       = undefined_real
        noahmp%energy%param%ResistanceSnowSfc(I,J)       = undefined_real
        noahmp%energy%param%VegFracAnnMax(I,J)           = undefined_real
        noahmp%energy%param%VegFracGreen(I,J)            = undefined_real
        noahmp%energy%param%HeatCapacCanFac(I,J)         = undefined_real

        ! Initialize 3D energy state arrays
        !$acc loop seq
        do LoopInd = -noahmp%config%domain%NumSnowLayerMax+1, noahmp%config%domain%NumSoilLayer
           noahmp%energy%state%TemperatureSoilSnow(I,LoopInd,J)  = undefined_real
           noahmp%energy%state%ThermConductSoilSnow(I,LoopInd,J) = undefined_real
           noahmp%energy%state%HeatCapacSoilSnow(I,LoopInd,J)    = undefined_real
           noahmp%energy%state%PhaseChgFacSoilSnow(I,LoopInd,J)  = undefined_real
           noahmp%energy%flux%RadSwPenetrateGrd(I,LoopInd,J)     = undefined_real
        enddo

        !$acc loop seq
        do LoopInd = -noahmp%config%domain%NumSnowLayerMax+1, 0
           noahmp%energy%state%HeatCapacVolSnow(I,LoopInd,J)  = undefined_real
           noahmp%energy%state%ThermConductSnow(I,LoopInd,J)  = undefined_real
        enddo

        !$acc loop seq
        do LoopInd = 1, noahmp%config%domain%NumSoilLayer
           noahmp%energy%state%HeatCapacVolSoil(I,LoopInd,J)   = undefined_real
           noahmp%energy%state%ThermConductSoil(I,LoopInd,J)   = undefined_real
           noahmp%energy%state%HeatCapacGlaIce(I,LoopInd,J)    = undefined_real
           noahmp%energy%state%ThermConductGlaIce(I,LoopInd,J) = undefined_real
           noahmp%energy%param%SoilQuartzFrac(I,LoopInd,J)     = undefined_real
        enddo

        !$acc loop seq
        do LoopInd = 1, noahmp%config%domain%NumSwRadBand
           noahmp%energy%state%AlbedoSnowDir(I,LoopInd,J)    = undefined_real
           noahmp%energy%state%AlbedoSnowDif(I,LoopInd,J)    = undefined_real
           noahmp%energy%state%AlbedoSoilDir(I,LoopInd,J)    = 0.0
           noahmp%energy%state%AlbedoSoilDif(I,LoopInd,J)    = 0.0
           noahmp%energy%state%AlbedoGrdDir(I,LoopInd,J)     = undefined_real
           noahmp%energy%state%AlbedoGrdDif(I,LoopInd,J)     = undefined_real
           noahmp%energy%state%ReflectanceVeg(I,LoopInd,J)   = undefined_real
           noahmp%energy%state%TransmittanceVeg(I,LoopInd,J) = undefined_real
           noahmp%energy%state%AlbedoSfcDir(I,LoopInd,J)     = undefined_real
           noahmp%energy%state%AlbedoSfcDif(I,LoopInd,J)     = undefined_real
           noahmp%energy%flux%RadSwAbsVegDir(I,LoopInd,J)    = undefined_real
           noahmp%energy%flux%RadSwAbsVegDif(I,LoopInd,J)    = undefined_real
           noahmp%energy%flux%RadSwDirTranGrdDir(I,LoopInd,J)= undefined_real
           noahmp%energy%flux%RadSwDirTranGrdDif(I,LoopInd,J)= undefined_real
           noahmp%energy%flux%RadSwDifTranGrdDir(I,LoopInd,J)= undefined_real
           noahmp%energy%flux%RadSwDifTranGrdDif(I,LoopInd,J)= undefined_real
           noahmp%energy%flux%RadSwReflVegDir(I,LoopInd,J)   = undefined_real
           noahmp%energy%flux%RadSwReflVegDif(I,LoopInd,J)   = undefined_real
           noahmp%energy%flux%RadSwReflGrdDir(I,LoopInd,J)   = undefined_real
           noahmp%energy%flux%RadSwReflGrdDif(I,LoopInd,J)   = undefined_real
           noahmp%energy%flux%RadSwDownDir(I,LoopInd,J)      = undefined_real
           noahmp%energy%flux%RadSwDownDif(I,LoopInd,J)      = undefined_real
           noahmp%energy%param%AlbedoSoilSat(I,LoopInd,J)    = undefined_real
           noahmp%energy%param%AlbedoSoilDry(I,LoopInd,J)    = undefined_real
           noahmp%energy%param%AlbedoLakeFrz(I,LoopInd,J)    = undefined_real
           noahmp%energy%param%ScatterCoeffSnow(I,LoopInd,J) = undefined_real
           noahmp%energy%param%ReflectanceLeaf(I,LoopInd,J)  = undefined_real
           noahmp%energy%param%ReflectanceStem(I,LoopInd,J)  = undefined_real
           noahmp%energy%param%TransmittanceLeaf(I,LoopInd,J)= undefined_real
           noahmp%energy%param%TransmittanceStem(I,LoopInd,J)= undefined_real
           noahmp%energy%param%AlbedoLandIce(I,LoopInd,J)    = undefined_real
        enddo

        !$acc loop seq
        do LoopInd = 1, 12
           noahmp%energy%param%LeafAreaIndexMon(I,LoopInd,J) = undefined_real
           noahmp%energy%param%StemAreaIndexMon(I,LoopInd,J) = undefined_real
        enddo

        !$acc loop seq
        do LoopInd = 1, 2
           noahmp%energy%param%EmissivitySoilLake(I,LoopInd,J) = undefined_real
        enddo

      end do
    end do
    !$acc end parallel loop

  end subroutine EnergyVarInitDefault


  subroutine EnergyVarExitDevice(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    !$acc exit data delete(               &
    !$acc   noahmp%energy%flux%HeatLatentCanopy,   &
    !$acc   noahmp%energy%flux%HeatLatentTransp,   &
    !$acc   noahmp%energy%flux%HeatLatentGrd,   &
    !$acc   noahmp%energy%flux%HeatLatentIrriEvap,   &
    !$acc   noahmp%energy%flux%HeatPrecipAdvCanopy,   &
    !$acc   noahmp%energy%flux%HeatPrecipAdvVegGrd,   &
    !$acc   noahmp%energy%flux%HeatPrecipAdvBareGrd,   &
    !$acc   noahmp%energy%flux%HeatPrecipAdvSfc,   &
    !$acc   noahmp%energy%flux%HeatSensibleCanopy,   &
    !$acc   noahmp%energy%flux%HeatLatentCanEvap,   &
    !$acc   noahmp%energy%flux%HeatSensibleVegGrd,   &
    !$acc   noahmp%energy%flux%HeatSensibleSfc,   &
    !$acc   noahmp%energy%flux%HeatLatentVegGrd,   &
    !$acc   noahmp%energy%flux%HeatLatentCanTransp,   &
    !$acc   noahmp%energy%flux%HeatGroundVegGrd,   &
    !$acc   noahmp%energy%flux%HeatSensibleBareGrd,   &
    !$acc   noahmp%energy%flux%HeatLatentBareGrd,   &
    !$acc   noahmp%energy%flux%HeatGroundBareGrd,   &
    !$acc   noahmp%energy%flux%HeatGroundTot,   &
    !$acc   noahmp%energy%flux%HeatGroundTotMean,   &
    !$acc   noahmp%energy%flux%HeatFromSoilBot,   &
    !$acc   noahmp%energy%flux%HeatCanStorageChg,   &
    !$acc   noahmp%energy%flux%HeatGroundTotAcc,   &
    !$acc   noahmp%energy%flux%RadPhotoActAbsSunlit,   &
    !$acc   noahmp%energy%flux%RadPhotoActAbsShade,   &
    !$acc   noahmp%energy%flux%RadSwAbsVeg,   &
    !$acc   noahmp%energy%flux%RadSwAbsGrd,   &
    !$acc   noahmp%energy%flux%RadSwAbsSfc,   &
    !$acc   noahmp%energy%flux%RadSwReflSfc,   &
    !$acc   noahmp%energy%flux%RadSwReflVeg    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%energy%flux%RadSwReflGrd,   &
    !$acc   noahmp%energy%flux%RadLwNetCanopy,   &
    !$acc   noahmp%energy%flux%RadLwNetSfc,   &
    !$acc   noahmp%energy%flux%RadPhotoActAbsCan,   &
    !$acc   noahmp%energy%flux%RadLwEmitSfc,   &
    !$acc   noahmp%energy%flux%RadLwNetVegGrd,   &
    !$acc   noahmp%energy%flux%RadLwNetBareGrd,   &
    !$acc   noahmp%energy%state%FlagFrozenCanopy,   &
    !$acc   noahmp%energy%state%FlagFrozenGround,   &
    !$acc   noahmp%energy%state%LeafAreaIndEff,   &
    !$acc   noahmp%energy%state%StemAreaIndEff,   &
    !$acc   noahmp%energy%state%LeafAreaIndex,   &
    !$acc   noahmp%energy%state%StemAreaIndex,   &
    !$acc   noahmp%energy%state%VegAreaIndEff,   &
    !$acc   noahmp%energy%state%VegFrac,   &
    !$acc   noahmp%energy%state%TemperatureGrd,   &
    !$acc   noahmp%energy%state%TemperatureCanopy,   &
    !$acc   noahmp%energy%state%TemperatureSfc,   &
    !$acc   noahmp%energy%state%TemperatureRootZone,   &
    !$acc   noahmp%energy%state%PressureVaporRefHeight,   &
    !$acc   noahmp%energy%state%SnowAgeFac,   &
    !$acc   noahmp%energy%state%SnowAgeNondim,   &
    !$acc   noahmp%energy%state%AlbedoSnowPrev,   &
    !$acc   noahmp%energy%state%VegAreaProjDir,   &
    !$acc   noahmp%energy%state%GapBtwCanopy,   &
    !$acc   noahmp%energy%state%GapInCanopy,   &
    !$acc   noahmp%energy%state%GapCanopyDif,   &
    !$acc   noahmp%energy%state%GapCanopyDir,   &
    !$acc   noahmp%energy%state%CanopySunlitFrac,   &
    !$acc   noahmp%energy%state%CanopyShadeFrac    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%energy%state%LeafAreaIndSunlit,   &
    !$acc   noahmp%energy%state%LeafAreaIndShade,   &
    !$acc   noahmp%energy%state%VapPresSatCanopy,   &
    !$acc   noahmp%energy%state%VapPresSatGrdVeg,   &
    !$acc   noahmp%energy%state%VapPresSatGrdBare,   &
    !$acc   noahmp%energy%state%VapPresSatCanTempD,   &
    !$acc   noahmp%energy%state%VapPresSatGrdVegTempD,   &
    !$acc   noahmp%energy%state%VapPresSatGrdBareTempD,   &
    !$acc   noahmp%energy%state%PressureVaporCanAir,   &
    !$acc   noahmp%energy%state%PressureAtmosCO2,   &
    !$acc   noahmp%energy%state%PressureAtmosO2,   &
    !$acc   noahmp%energy%state%ResistanceStomataSunlit,   &
    !$acc   noahmp%energy%state%ResistanceStomataShade,   &
    !$acc   noahmp%energy%state%DensityAirRefHeight,   &
    !$acc   noahmp%energy%state%TemperatureCanopyAir,   &
    !$acc   noahmp%energy%state%ZeroPlaneDispSfc,   &
    !$acc   noahmp%energy%state%ZeroPlaneDispGrd,   &
    !$acc   noahmp%energy%state%RoughLenMomGrd,   &
    !$acc   noahmp%energy%state%RoughLenMomSfc,   &
    !$acc   noahmp%energy%state%RoughLenShCanopy,   &
    !$acc   noahmp%energy%state%RoughLenShVegGrd,   &
    !$acc   noahmp%energy%state%RoughLenShBareGrd,   &
    !$acc   noahmp%energy%state%CanopyHeight,   &
    !$acc   noahmp%energy%state%WindSpdCanopyTop,   &
    !$acc   noahmp%energy%state%FrictionVelVeg,   &
    !$acc   noahmp%energy%state%FrictionVelBare,   &
    !$acc   noahmp%energy%state%WindExtCoeffCanopy,   &
    !$acc   noahmp%energy%state%MoStabParaUndCan,   &
    !$acc   noahmp%energy%state%MoStabParaAbvCan,   &
    !$acc   noahmp%energy%state%MoStabParaBare    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%energy%state%MoStabParaVeg2m,   &
    !$acc   noahmp%energy%state%MoStabParaBare2m,   &
    !$acc   noahmp%energy%state%MoLengthUndCan,   &
    !$acc   noahmp%energy%state%MoLengthAbvCan,   &
    !$acc   noahmp%energy%state%MoLengthBare,   &
    !$acc   noahmp%energy%state%MoStabCorrShUndCan,   &
    !$acc   noahmp%energy%state%MoStabCorrMomAbvCan,   &
    !$acc   noahmp%energy%state%MoStabCorrShAbvCan,   &
    !$acc   noahmp%energy%state%MoStabCorrMomVeg2m,   &
    !$acc   noahmp%energy%state%MoStabCorrShVeg2m,   &
    !$acc   noahmp%energy%state%MoStabCorrShBare,   &
    !$acc   noahmp%energy%state%MoStabCorrMomBare,   &
    !$acc   noahmp%energy%state%MoStabCorrMomBare2m,   &
    !$acc   noahmp%energy%state%MoStabCorrShBare2m,   &
    !$acc   noahmp%energy%state%ExchCoeffMomSfc,   &
    !$acc   noahmp%energy%state%ExchCoeffMomAbvCan,   &
    !$acc   noahmp%energy%state%ExchCoeffMomBare,   &
    !$acc   noahmp%energy%state%ExchCoeffShSfc,   &
    !$acc   noahmp%energy%state%ExchCoeffShAbvCan,   &
    !$acc   noahmp%energy%state%ExchCoeffShBare,   &
    !$acc   noahmp%energy%state%ExchCoeffSh2mVegMo,   &
    !$acc   noahmp%energy%state%ExchCoeffSh2mBareMo,   &
    !$acc   noahmp%energy%state%ExchCoeffSh2mVeg,   &
    !$acc   noahmp%energy%state%ExchCoeffLhAbvCan,   &
    !$acc   noahmp%energy%state%ExchCoeffLhTransp,   &
    !$acc   noahmp%energy%state%ExchCoeffLhEvap,   &
    !$acc   noahmp%energy%state%ExchCoeffLhUndCan,   &
    !$acc   noahmp%energy%state%ResistanceMomUndCan,   &
    !$acc   noahmp%energy%state%ResistanceShUndCan,   &
    !$acc   noahmp%energy%state%ResistanceLhUndCan    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%energy%state%ResistanceMomAbvCan,   &
    !$acc   noahmp%energy%state%ResistanceShAbvCan,   &
    !$acc   noahmp%energy%state%ResistanceLhAbvCan,   &
    !$acc   noahmp%energy%state%ResistanceMomBareGrd,   &
    !$acc   noahmp%energy%state%ResistanceShBareGrd,   &
    !$acc   noahmp%energy%state%ResistanceLhBareGrd,   &
    !$acc   noahmp%energy%state%ResistanceLeafBoundary,   &
    !$acc   noahmp%energy%state%TemperaturePotRefHeight,   &
    !$acc   noahmp%energy%state%WindSpdRefHeight,   &
    !$acc   noahmp%energy%state%FrictionVelVertVeg,   &
    !$acc   noahmp%energy%state%FrictionVelVertBare,   &
    !$acc   noahmp%energy%state%EmissivityVeg,   &
    !$acc   noahmp%energy%state%EmissivityGrd,   &
    !$acc   noahmp%energy%state%ResistanceGrdEvap,   &
    !$acc   noahmp%energy%state%PsychConstCanopy,   &
    !$acc   noahmp%energy%state%LatHeatVapCanopy,   &
    !$acc   noahmp%energy%state%PsychConstGrd,   &
    !$acc   noahmp%energy%state%LatHeatVapGrd,   &
    !$acc   noahmp%energy%state%RelHumidityGrd,   &
    !$acc   noahmp%energy%state%SpecHumiditySfc,   &
    !$acc   noahmp%energy%state%SpecHumiditySfcMean,   &
    !$acc   noahmp%energy%state%SpecHumidity2mVeg,   &
    !$acc   noahmp%energy%state%SpecHumidity2mBare,   &
    !$acc   noahmp%energy%state%SpecHumidity2m,   &
    !$acc   noahmp%energy%state%TemperatureGrdVeg,   &
    !$acc   noahmp%energy%state%TemperatureGrdBare,   &
    !$acc   noahmp%energy%state%WindStressEwVeg,   &
    !$acc   noahmp%energy%state%WindStressNsVeg,   &
    !$acc   noahmp%energy%state%WindStressEwBare,   &
    !$acc   noahmp%energy%state%WindStressNsBare    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%energy%state%WindStressEwSfc,   &
    !$acc   noahmp%energy%state%WindStressNsSfc,   &
    !$acc   noahmp%energy%state%TemperatureAir2mVeg,   &
    !$acc   noahmp%energy%state%TemperatureAir2mBare,   &
    !$acc   noahmp%energy%state%TemperatureAir2m,   &
    !$acc   noahmp%energy%state%ExchCoeffShLeaf,   &
    !$acc   noahmp%energy%state%ExchCoeffShUndCan,   &
    !$acc   noahmp%energy%state%ExchCoeffSh2mBare,   &
    !$acc   noahmp%energy%state%RefHeightAboveGrd,   &
    !$acc   noahmp%energy%state%CanopyFracSnowBury,   &
    !$acc   noahmp%energy%state%DepthSoilTempBotToSno,   &
    !$acc   noahmp%energy%state%RoughLenMomSfcToAtm,   &
    !$acc   noahmp%energy%state%TemperatureRadSfc,   &
    !$acc   noahmp%energy%state%EmissivitySfc,   &
    !$acc   noahmp%energy%state%AlbedoSfc,   &
    !$acc   noahmp%energy%state%EnergyBalanceError,   &
    !$acc   noahmp%energy%state%RadSwBalanceError,   &
    !$acc   noahmp%energy%state%TemperatureSoilSnow,   &
    !$acc   noahmp%energy%state%ThermConductSoilSnow,   &
    !$acc   noahmp%energy%state%HeatCapacSoilSnow,   &
    !$acc   noahmp%energy%state%PhaseChgFacSoilSnow,   &
    !$acc   noahmp%energy%state%HeatCapacVolSnow,   &
    !$acc   noahmp%energy%state%ThermConductSnow,   &
    !$acc   noahmp%energy%state%HeatCapacVolSoil,   &
    !$acc   noahmp%energy%state%ThermConductSoil,   &
    !$acc   noahmp%energy%state%HeatCapacGlaIce,   &
    !$acc   noahmp%energy%state%ThermConductGlaIce,   &
    !$acc   noahmp%energy%state%AlbedoSnowDir,   &
    !$acc   noahmp%energy%state%AlbedoSnowDif,   &
    !$acc   noahmp%energy%state%AlbedoSoilDir    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%energy%state%AlbedoSoilDif,   &
    !$acc   noahmp%energy%state%AlbedoGrdDir,   &
    !$acc   noahmp%energy%state%AlbedoGrdDif,   &
    !$acc   noahmp%energy%state%ReflectanceVeg,   &
    !$acc   noahmp%energy%state%TransmittanceVeg,   &
    !$acc   noahmp%energy%state%AlbedoSfcDir,   &
    !$acc   noahmp%energy%state%AlbedoSfcDif,   &
    !$acc   noahmp%energy%flux%RadSwAbsVegDir,   &
    !$acc   noahmp%energy%flux%RadSwAbsVegDif,   &
    !$acc   noahmp%energy%flux%RadSwDirTranGrdDir,   &
    !$acc   noahmp%energy%flux%RadSwDirTranGrdDif,   &
    !$acc   noahmp%energy%flux%RadSwDifTranGrdDir,   &
    !$acc   noahmp%energy%flux%RadSwDifTranGrdDif,   &
    !$acc   noahmp%energy%flux%RadSwReflVegDir,   &
    !$acc   noahmp%energy%flux%RadSwReflVegDif,   &
    !$acc   noahmp%energy%flux%RadSwReflGrdDir,   &
    !$acc   noahmp%energy%flux%RadSwReflGrdDif,   &
    !$acc   noahmp%energy%flux%RadSwDownDir,   &
    !$acc   noahmp%energy%flux%RadSwDownDif,   &
    !$acc   noahmp%energy%flux%RadSwPenetrateGrd,   &
    !$acc   noahmp%energy%param%LeafAreaIndexMon,   &
    !$acc   noahmp%energy%param%StemAreaIndexMon,   &
    !$acc   noahmp%energy%param%SoilQuartzFrac,   &
    !$acc   noahmp%energy%param%AlbedoSoilSat,   &
    !$acc   noahmp%energy%param%AlbedoSoilDry,   &
    !$acc   noahmp%energy%param%AlbedoLakeFrz,   &
    !$acc   noahmp%energy%param%ScatterCoeffSnow,   &
    !$acc   noahmp%energy%param%ReflectanceLeaf,   &
    !$acc   noahmp%energy%param%ReflectanceStem,   &
    !$acc   noahmp%energy%param%TransmittanceLeaf    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%energy%param%TransmittanceStem,   &
    !$acc   noahmp%energy%param%EmissivitySoilLake,   &
    !$acc   noahmp%energy%param%AlbedoLandIce,   &
    !$acc   noahmp%energy%param%TreeCrownRadius,   &
    !$acc   noahmp%energy%param%HeightCanopyTop,   &
    !$acc   noahmp%energy%param%HeightCanopyBot,   &
    !$acc   noahmp%energy%param%RoughLenMomVeg,   &
    !$acc   noahmp%energy%param%TreeDensity,   &
    !$acc   noahmp%energy%param%CanopyOrientIndex,   &
    !$acc   noahmp%energy%param%UpscatterCoeffSnowDir,   &
    !$acc   noahmp%energy%param%UpscatterCoeffSnowDif,   &
    !$acc   noahmp%energy%param%SoilHeatCapacity,   &
    !$acc   noahmp%energy%param%SnowAgeFacBats,   &
    !$acc   noahmp%energy%param%SnowGrowVapFacBats,   &
    !$acc   noahmp%energy%param%SnowSootFacBats,   &
    !$acc   noahmp%energy%param%SnowGrowFrzFacBats,   &
    !$acc   noahmp%energy%param%SolarZenithAdjBats,   &
    !$acc   noahmp%energy%param%FreshSnoAlbVisBats,   &
    !$acc   noahmp%energy%param%FreshSnoAlbNirBats,   &
    !$acc   noahmp%energy%param%SnoAgeFacDifVisBats,   &
    !$acc   noahmp%energy%param%SnoAgeFacDifNirBats,   &
    !$acc   noahmp%energy%param%SzaFacDirVisBats,   &
    !$acc   noahmp%energy%param%SzaFacDirNirBats,   &
    !$acc   noahmp%energy%param%SnowAlbRefClass,   &
    !$acc   noahmp%energy%param%SnowAgeFacClass,   &
    !$acc   noahmp%energy%param%SnowAlbFreshClass,   &
    !$acc   noahmp%energy%param%ConductanceLeafMin,   &
    !$acc   noahmp%energy%param%Co2MmConst25C,   &
    !$acc   noahmp%energy%param%O2MmConst25C,   &
    !$acc   noahmp%energy%param%Co2MmConstQ10    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%energy%param%O2MmConstQ10,   &
    !$acc   noahmp%energy%param%RadiationStressFac,   &
    !$acc   noahmp%energy%param%ResistanceStomataMin,   &
    !$acc   noahmp%energy%param%ResistanceStomataMax,   &
    !$acc   noahmp%energy%param%AirTempOptimTransp,   &
    !$acc   noahmp%energy%param%VaporPresDeficitFac,   &
    !$acc   noahmp%energy%param%LeafDimLength,   &
    !$acc   noahmp%energy%param%ZilitinkevichCoeff,   &
    !$acc   noahmp%energy%param%EmissivitySnow,   &
    !$acc   noahmp%energy%param%CanopyWindExtFac,   &
    !$acc   noahmp%energy%param%RoughLenMomSnow,   &
    !$acc   noahmp%energy%param%RoughLenMomSoil,   &
    !$acc   noahmp%energy%param%RoughLenMomLake,   &
    !$acc   noahmp%energy%param%EmissivityIceSfc,   &
    !$acc   noahmp%energy%param%ResistanceSoilExp,   &
    !$acc   noahmp%energy%param%ResistanceSnowSfc,   &
    !$acc   noahmp%energy%param%VegFracGreen,   &
    !$acc   noahmp%energy%param%VegFracAnnMax,   &
    !$acc   noahmp%energy%param%HeatCapacCanFac    &
    !$acc   )

    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
    !$acc exit data delete(               &
    !$acc   noahmp%energy%flux%FracRadSwAbsSnowDir,   &
    !$acc   noahmp%energy%flux%FracRadSwAbsSnowDif,   &
    !$acc   noahmp%energy%flux%RadSwAbsSnowSoilLayer,   &
    !$acc   noahmp%energy%param%RadSwWgtDir,   &
    !$acc   noahmp%energy%param%RadSwWgtDif,   &
    !$acc   noahmp%energy%param%SsAlbSnwRadDir,   &
    !$acc   noahmp%energy%param%AsyPrmSnwRadDir,   &
    !$acc   noahmp%energy%param%ExtCffMassSnwRadDir,   &
    !$acc   noahmp%energy%param%SsAlbSnwRadDif,   &
    !$acc   noahmp%energy%param%AsyPrmSnwRadDif,   &
    !$acc   noahmp%energy%param%ExtCffMassSnwRadDif,   &
    !$acc   noahmp%energy%param%SsAlbBCphi,   &
    !$acc   noahmp%energy%param%AsyPrmBCphi,   &
    !$acc   noahmp%energy%param%ExtCffMassBCphi,   &
    !$acc   noahmp%energy%param%SsAlbBCpho,   &
    !$acc   noahmp%energy%param%AsyPrmBCpho,   &
    !$acc   noahmp%energy%param%ExtCffMassBCpho,   &
    !$acc   noahmp%energy%param%SsAlbOCphi,   &
    !$acc   noahmp%energy%param%AsyPrmOCphi,   &
    !$acc   noahmp%energy%param%ExtCffMassOCphi,   &
    !$acc   noahmp%energy%param%SsAlbOCpho,   &
    !$acc   noahmp%energy%param%AsyPrmOCpho,   &
    !$acc   noahmp%energy%param%ExtCffMassOCpho,   &
    !$acc   noahmp%energy%param%SsAlbDustB1,   &
    !$acc   noahmp%energy%param%AsyPrmDustB1,   &
    !$acc   noahmp%energy%param%ExtCffMassDustB1,   &
    !$acc   noahmp%energy%param%SsAlbDustB2,   &
    !$acc   noahmp%energy%param%AsyPrmDustB2,   &
    !$acc   noahmp%energy%param%ExtCffMassDustB2,   &
    !$acc   noahmp%energy%param%SsAlbDustB3    &
    !$acc   )
    !$acc exit data delete(               &
    !$acc   noahmp%energy%param%AsyPrmDustB3,   &
    !$acc   noahmp%energy%param%ExtCffMassDustB3,   &
    !$acc   noahmp%energy%param%SsAlbDustB4,   &
    !$acc   noahmp%energy%param%AsyPrmDustB4,   &
    !$acc   noahmp%energy%param%ExtCffMassDustB4,   &
    !$acc   noahmp%energy%param%SsAlbDustB5,   &
    !$acc   noahmp%energy%param%AsyPrmDustB5,   &
    !$acc   noahmp%energy%param%ExtCffMassDustB5    &
    !$acc   )
    endif

  end subroutine EnergyVarExitDevice
end module EnergyVarInitMod
