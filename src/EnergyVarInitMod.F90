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
       !$acc enter data create(noahmp%energy%flux%HeatLatentCanopy)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentTransp) ) then
       allocate( noahmp%energy%flux%HeatLatentTransp(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatLatentTransp)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentGrd) ) then
       allocate( noahmp%energy%flux%HeatLatentGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatLatentGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentIrriEvap) ) then
       allocate( noahmp%energy%flux%HeatLatentIrriEvap(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatLatentIrriEvap)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatPrecipAdvCanopy) ) then
       allocate( noahmp%energy%flux%HeatPrecipAdvCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatPrecipAdvCanopy)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatPrecipAdvVegGrd) ) then
       allocate( noahmp%energy%flux%HeatPrecipAdvVegGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatPrecipAdvVegGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatPrecipAdvBareGrd) ) then
       allocate( noahmp%energy%flux%HeatPrecipAdvBareGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatPrecipAdvBareGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatPrecipAdvSfc) ) then
       allocate( noahmp%energy%flux%HeatPrecipAdvSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatPrecipAdvSfc)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatSensibleCanopy) ) then
       allocate( noahmp%energy%flux%HeatSensibleCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatSensibleCanopy)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentCanEvap) ) then
       allocate( noahmp%energy%flux%HeatLatentCanEvap(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatLatentCanEvap)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatSensibleVegGrd) ) then
       allocate( noahmp%energy%flux%HeatSensibleVegGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatSensibleVegGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatSensibleSfc) ) then
       allocate( noahmp%energy%flux%HeatSensibleSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatSensibleSfc)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentVegGrd) ) then
       allocate( noahmp%energy%flux%HeatLatentVegGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatLatentVegGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentCanTransp) ) then
       allocate( noahmp%energy%flux%HeatLatentCanTransp(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatLatentCanTransp)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatGroundVegGrd) ) then
       allocate( noahmp%energy%flux%HeatGroundVegGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatGroundVegGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatSensibleBareGrd) ) then
       allocate( noahmp%energy%flux%HeatSensibleBareGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatSensibleBareGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatLatentBareGrd) ) then
       allocate( noahmp%energy%flux%HeatLatentBareGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatLatentBareGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatGroundBareGrd) ) then
       allocate( noahmp%energy%flux%HeatGroundBareGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatGroundBareGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatGroundTot) ) then
       allocate( noahmp%energy%flux%HeatGroundTot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatGroundTot)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatGroundTotMean) ) then
       allocate( noahmp%energy%flux%HeatGroundTotMean(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatGroundTotMean)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatFromSoilBot) ) then
       allocate( noahmp%energy%flux%HeatFromSoilBot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatFromSoilBot)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatCanStorageChg) ) then
       allocate( noahmp%energy%flux%HeatCanStorageChg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatCanStorageChg)
    endif
    if ( .not. allocated(noahmp%energy%flux%HeatGroundTotAcc) ) then
       allocate( noahmp%energy%flux%HeatGroundTotAcc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%HeatGroundTotAcc)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadPhotoActAbsSunlit) ) then
       allocate( noahmp%energy%flux%RadPhotoActAbsSunlit(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadPhotoActAbsSunlit)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadPhotoActAbsShade) ) then
       allocate( noahmp%energy%flux%RadPhotoActAbsShade(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadPhotoActAbsShade)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwAbsVeg) ) then
       allocate( noahmp%energy%flux%RadSwAbsVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwAbsVeg)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwAbsGrd) ) then
       allocate( noahmp%energy%flux%RadSwAbsGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwAbsGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwAbsSfc) ) then
       allocate( noahmp%energy%flux%RadSwAbsSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwAbsSfc)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflSfc) ) then
       allocate( noahmp%energy%flux%RadSwReflSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwReflSfc)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflVeg) ) then
       allocate( noahmp%energy%flux%RadSwReflVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwReflVeg)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflGrd) ) then
       allocate( noahmp%energy%flux%RadSwReflGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwReflGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadLwNetCanopy) ) then
       allocate( noahmp%energy%flux%RadLwNetCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadLwNetCanopy)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadLwNetSfc) ) then
       allocate( noahmp%energy%flux%RadLwNetSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadLwNetSfc)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadPhotoActAbsCan) ) then
       allocate( noahmp%energy%flux%RadPhotoActAbsCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadPhotoActAbsCan)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadLwEmitSfc) ) then
       allocate( noahmp%energy%flux%RadLwEmitSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadLwEmitSfc)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadLwNetVegGrd) ) then
       allocate( noahmp%energy%flux%RadLwNetVegGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadLwNetVegGrd)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadLwNetBareGrd) ) then
       allocate( noahmp%energy%flux%RadLwNetBareGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadLwNetBareGrd)
    endif

    ! Allocate 3D energy state arrays and transfer to GPU

    if (.not. allocated(noahmp%energy%state%FlagFrozenCanopy) ) then
       allocate( noahmp%energy%state%FlagFrozenCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%FlagFrozenCanopy)
    end if

    if ( .not. allocated(noahmp%energy%state%FlagFrozenGround) ) then
       allocate( noahmp%energy%state%FlagFrozenGround(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%FlagFrozenGround)
    end if

    if ( .not. allocated(noahmp%energy%state%LeafAreaIndEff) ) then
       allocate( noahmp%energy%state%LeafAreaIndEff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%LeafAreaIndEff)
    end if

    if ( .not. allocated(noahmp%energy%state%StemAreaIndEff) ) then
       allocate( noahmp%energy%state%StemAreaIndEff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%StemAreaIndEff)
    end if

    if ( .not. allocated(noahmp%energy%state%LeafAreaIndex) ) then
       allocate( noahmp%energy%state%LeafAreaIndex(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%LeafAreaIndex)
    end if

    if ( .not. allocated(noahmp%energy%state%StemAreaIndex) ) then
       allocate( noahmp%energy%state%StemAreaIndex(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%StemAreaIndex)
    end if

    if ( .not. allocated(noahmp%energy%state%VegAreaIndEff) ) then
       allocate( noahmp%energy%state%VegAreaIndEff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%VegAreaIndEff)
    end if

    if ( .not. allocated(noahmp%energy%state%VegFrac) ) then
       allocate( noahmp%energy%state%VegFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%VegFrac)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureGrd) ) then
       allocate( noahmp%energy%state%TemperatureGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureCanopy) ) then
       allocate( noahmp%energy%state%TemperatureCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureCanopy)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureSfc) ) then
       allocate( noahmp%energy%state%TemperatureSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureSfc)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureRootZone) ) then
       allocate( noahmp%energy%state%TemperatureRootZone(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureRootZone)
    end if

    if ( .not. allocated(noahmp%energy%state%PressureVaporRefHeight) ) then
       allocate( noahmp%energy%state%PressureVaporRefHeight(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%PressureVaporRefHeight)
    end if

    if ( .not. allocated(noahmp%energy%state%SnowAgeFac) ) then
       allocate( noahmp%energy%state%SnowAgeFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%SnowAgeFac)
    end if

    if ( .not. allocated(noahmp%energy%state%SnowAgeNondim) ) then
       allocate( noahmp%energy%state%SnowAgeNondim(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%SnowAgeNondim)
    end if

    if ( .not. allocated(noahmp%energy%state%AlbedoSnowPrev) ) then
       allocate( noahmp%energy%state%AlbedoSnowPrev(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%AlbedoSnowPrev)
    end if

    if ( .not. allocated(noahmp%energy%state%VegAreaProjDir) ) then
       allocate( noahmp%energy%state%VegAreaProjDir(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%VegAreaProjDir)
    end if

    if ( .not. allocated(noahmp%energy%state%GapBtwCanopy) ) then
       allocate( noahmp%energy%state%GapBtwCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%GapBtwCanopy)
    end if

    if ( .not. allocated(noahmp%energy%state%GapInCanopy) ) then
       allocate( noahmp%energy%state%GapInCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%GapInCanopy)
    end if

    if ( .not. allocated(noahmp%energy%state%GapCanopyDif) ) then
       allocate( noahmp%energy%state%GapCanopyDif(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%GapCanopyDif)
    end if

    if ( .not. allocated(noahmp%energy%state%GapCanopyDir) ) then
       allocate( noahmp%energy%state%GapCanopyDir(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%GapCanopyDir)
    end if

    if ( .not. allocated(noahmp%energy%state%CanopySunlitFrac) ) then
       allocate( noahmp%energy%state%CanopySunlitFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%CanopySunlitFrac)
    end if

    if ( .not. allocated(noahmp%energy%state%CanopyShadeFrac) ) then
       allocate( noahmp%energy%state%CanopyShadeFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%CanopyShadeFrac)
    end if

    if ( .not. allocated(noahmp%energy%state%LeafAreaIndSunlit) ) then
       allocate( noahmp%energy%state%LeafAreaIndSunlit(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%LeafAreaIndSunlit)
    end if

    if ( .not. allocated(noahmp%energy%state%LeafAreaIndShade) ) then
       allocate( noahmp%energy%state%LeafAreaIndShade(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%LeafAreaIndShade)
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatCanopy) ) then
       allocate( noahmp%energy%state%VapPresSatCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%VapPresSatCanopy)
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatGrdVeg) ) then
       allocate( noahmp%energy%state%VapPresSatGrdVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%VapPresSatGrdVeg)
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatGrdBare) ) then
       allocate( noahmp%energy%state%VapPresSatGrdBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%VapPresSatGrdBare)
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatCanTempD) ) then
       allocate( noahmp%energy%state%VapPresSatCanTempD(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%VapPresSatCanTempD)
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatGrdVegTempD) ) then
       allocate( noahmp%energy%state%VapPresSatGrdVegTempD(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%VapPresSatGrdVegTempD)
    end if

    if ( .not. allocated(noahmp%energy%state%VapPresSatGrdBareTempD) ) then
       allocate( noahmp%energy%state%VapPresSatGrdBareTempD(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%VapPresSatGrdBareTempD)
    end if

    if ( .not. allocated(noahmp%energy%state%PressureVaporCanAir) ) then
       allocate( noahmp%energy%state%PressureVaporCanAir(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%PressureVaporCanAir)
    end if

    if ( .not. allocated(noahmp%energy%state%PressureAtmosCO2) ) then
       allocate( noahmp%energy%state%PressureAtmosCO2(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%PressureAtmosCO2)
    end if

    if ( .not. allocated(noahmp%energy%state%PressureAtmosO2) ) then
       allocate( noahmp%energy%state%PressureAtmosO2(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%PressureAtmosO2)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceStomataSunlit) ) then
       allocate( noahmp%energy%state%ResistanceStomataSunlit(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceStomataSunlit)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceStomataShade) ) then
       allocate( noahmp%energy%state%ResistanceStomataShade(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceStomataShade)
    end if

    if ( .not. allocated(noahmp%energy%state%DensityAirRefHeight) ) then
       allocate( noahmp%energy%state%DensityAirRefHeight(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%DensityAirRefHeight)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureCanopyAir) ) then
       allocate( noahmp%energy%state%TemperatureCanopyAir(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureCanopyAir)
    end if

    if ( .not. allocated(noahmp%energy%state%ZeroPlaneDispSfc) ) then
       allocate( noahmp%energy%state%ZeroPlaneDispSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ZeroPlaneDispSfc)
    end if

    if ( .not. allocated(noahmp%energy%state%ZeroPlaneDispGrd) ) then
       allocate( noahmp%energy%state%ZeroPlaneDispGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ZeroPlaneDispGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenMomGrd) ) then
       allocate( noahmp%energy%state%RoughLenMomGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%RoughLenMomGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenMomSfc) ) then
       allocate( noahmp%energy%state%RoughLenMomSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%RoughLenMomSfc)
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenShCanopy) ) then
       allocate( noahmp%energy%state%RoughLenShCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%RoughLenShCanopy)
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenShVegGrd) ) then
       allocate( noahmp%energy%state%RoughLenShVegGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%RoughLenShVegGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenShBareGrd) ) then
       allocate( noahmp%energy%state%RoughLenShBareGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%RoughLenShBareGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%CanopyHeight) ) then
       allocate( noahmp%energy%state%CanopyHeight(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%CanopyHeight)
    end if

    if ( .not. allocated(noahmp%energy%state%WindSpdCanopyTop) ) then
       allocate( noahmp%energy%state%WindSpdCanopyTop(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%WindSpdCanopyTop)
    end if

    if ( .not. allocated(noahmp%energy%state%FrictionVelVeg) ) then
       allocate( noahmp%energy%state%FrictionVelVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%FrictionVelVeg)
    end if

    if ( .not. allocated(noahmp%energy%state%FrictionVelBare) ) then
       allocate( noahmp%energy%state%FrictionVelBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%FrictionVelBare)
    end if

    if ( .not. allocated(noahmp%energy%state%WindExtCoeffCanopy) ) then
       allocate( noahmp%energy%state%WindExtCoeffCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%WindExtCoeffCanopy)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabParaUndCan) ) then
       allocate( noahmp%energy%state%MoStabParaUndCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabParaUndCan)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabParaAbvCan) ) then
       allocate( noahmp%energy%state%MoStabParaAbvCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabParaAbvCan)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabParaBare) ) then
       allocate( noahmp%energy%state%MoStabParaBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabParaBare)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabParaVeg2m) ) then
       allocate( noahmp%energy%state%MoStabParaVeg2m(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabParaVeg2m)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabParaBare2m) ) then
       allocate( noahmp%energy%state%MoStabParaBare2m(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabParaBare2m)
    end if

    if ( .not. allocated(noahmp%energy%state%MoLengthUndCan) ) then
       allocate( noahmp%energy%state%MoLengthUndCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoLengthUndCan)
    end if

    if ( .not. allocated(noahmp%energy%state%MoLengthAbvCan) ) then
       allocate( noahmp%energy%state%MoLengthAbvCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoLengthAbvCan)
    end if

    if ( .not. allocated(noahmp%energy%state%MoLengthBare) ) then
       allocate( noahmp%energy%state%MoLengthBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoLengthBare)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrShUndCan) ) then
       allocate( noahmp%energy%state%MoStabCorrShUndCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabCorrShUndCan)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrMomAbvCan) ) then
       allocate( noahmp%energy%state%MoStabCorrMomAbvCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabCorrMomAbvCan)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrShAbvCan) ) then
       allocate( noahmp%energy%state%MoStabCorrShAbvCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabCorrShAbvCan)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrMomVeg2m) ) then
       allocate( noahmp%energy%state%MoStabCorrMomVeg2m(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabCorrMomVeg2m)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrShVeg2m) ) then
       allocate( noahmp%energy%state%MoStabCorrShVeg2m(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabCorrShVeg2m)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrShBare) ) then
       allocate( noahmp%energy%state%MoStabCorrShBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabCorrShBare)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrMomBare) ) then
       allocate( noahmp%energy%state%MoStabCorrMomBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabCorrMomBare)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrMomBare2m) ) then
       allocate( noahmp%energy%state%MoStabCorrMomBare2m(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabCorrMomBare2m)
    end if

    if ( .not. allocated(noahmp%energy%state%MoStabCorrShBare2m) ) then
       allocate( noahmp%energy%state%MoStabCorrShBare2m(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%MoStabCorrShBare2m)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffMomSfc) ) then
       allocate( noahmp%energy%state%ExchCoeffMomSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffMomSfc)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffMomAbvCan) ) then
       allocate( noahmp%energy%state%ExchCoeffMomAbvCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffMomAbvCan)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffMomBare) ) then
       allocate( noahmp%energy%state%ExchCoeffMomBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffMomBare)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffShSfc) ) then
       allocate( noahmp%energy%state%ExchCoeffShSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffShSfc)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffShAbvCan) ) then
       allocate( noahmp%energy%state%ExchCoeffShAbvCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffShAbvCan)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffShBare) ) then
       allocate( noahmp%energy%state%ExchCoeffShBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffShBare)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffSh2mVegMo) ) then
       allocate( noahmp%energy%state%ExchCoeffSh2mVegMo(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffSh2mVegMo)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffSh2mBareMo) ) then
       allocate( noahmp%energy%state%ExchCoeffSh2mBareMo(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffSh2mBareMo)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffSh2mVeg) ) then
       allocate( noahmp%energy%state%ExchCoeffSh2mVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffSh2mVeg)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffLhAbvCan) ) then
       allocate( noahmp%energy%state%ExchCoeffLhAbvCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffLhAbvCan)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffLhTransp) ) then
       allocate( noahmp%energy%state%ExchCoeffLhTransp(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffLhTransp)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffLhEvap) ) then
       allocate( noahmp%energy%state%ExchCoeffLhEvap(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffLhEvap)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffLhUndCan) ) then
       allocate( noahmp%energy%state%ExchCoeffLhUndCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffLhUndCan)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceMomUndCan) ) then
       allocate( noahmp%energy%state%ResistanceMomUndCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceMomUndCan)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceShUndCan) ) then
       allocate( noahmp%energy%state%ResistanceShUndCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceShUndCan)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceLhUndCan) ) then
       allocate( noahmp%energy%state%ResistanceLhUndCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceLhUndCan)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceMomAbvCan) ) then
       allocate( noahmp%energy%state%ResistanceMomAbvCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceMomAbvCan)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceShAbvCan) ) then
       allocate( noahmp%energy%state%ResistanceShAbvCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceShAbvCan)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceLhAbvCan) ) then
       allocate( noahmp%energy%state%ResistanceLhAbvCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceLhAbvCan)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceMomBareGrd) ) then
       allocate( noahmp%energy%state%ResistanceMomBareGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceMomBareGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceShBareGrd) ) then
       allocate( noahmp%energy%state%ResistanceShBareGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceShBareGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceLhBareGrd) ) then
       allocate( noahmp%energy%state%ResistanceLhBareGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceLhBareGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceLeafBoundary) ) then
       allocate( noahmp%energy%state%ResistanceLeafBoundary(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceLeafBoundary)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperaturePotRefHeight) ) then
       allocate( noahmp%energy%state%TemperaturePotRefHeight(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperaturePotRefHeight)
    end if

    if ( .not. allocated(noahmp%energy%state%WindSpdRefHeight) ) then
       allocate( noahmp%energy%state%WindSpdRefHeight(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%WindSpdRefHeight)
    end if

    if ( .not. allocated(noahmp%energy%state%FrictionVelVertVeg) ) then
       allocate( noahmp%energy%state%FrictionVelVertVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%FrictionVelVertVeg)
    end if

    if ( .not. allocated(noahmp%energy%state%FrictionVelVertBare) ) then
       allocate( noahmp%energy%state%FrictionVelVertBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%FrictionVelVertBare)
    end if

    if ( .not. allocated(noahmp%energy%state%EmissivityVeg) ) then
       allocate( noahmp%energy%state%EmissivityVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%EmissivityVeg)
    end if

    if ( .not. allocated(noahmp%energy%state%EmissivityGrd) ) then
       allocate( noahmp%energy%state%EmissivityGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%EmissivityGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%ResistanceGrdEvap) ) then
       allocate( noahmp%energy%state%ResistanceGrdEvap(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ResistanceGrdEvap)
    end if

    if ( .not. allocated(noahmp%energy%state%PsychConstCanopy) ) then
       allocate( noahmp%energy%state%PsychConstCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%PsychConstCanopy)
    end if

    if ( .not. allocated(noahmp%energy%state%LatHeatVapCanopy) ) then
       allocate( noahmp%energy%state%LatHeatVapCanopy(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%LatHeatVapCanopy)
    end if

    if ( .not. allocated(noahmp%energy%state%PsychConstGrd) ) then
       allocate( noahmp%energy%state%PsychConstGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%PsychConstGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%LatHeatVapGrd) ) then
       allocate( noahmp%energy%state%LatHeatVapGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%LatHeatVapGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%RelHumidityGrd) ) then
       allocate( noahmp%energy%state%RelHumidityGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%RelHumidityGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%SpecHumiditySfc) ) then
       allocate( noahmp%energy%state%SpecHumiditySfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%SpecHumiditySfc)
    end if

    if ( .not. allocated(noahmp%energy%state%SpecHumiditySfcMean) ) then
       allocate( noahmp%energy%state%SpecHumiditySfcMean(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%SpecHumiditySfcMean)
    end if

    if ( .not. allocated(noahmp%energy%state%SpecHumidity2mVeg) ) then
       allocate( noahmp%energy%state%SpecHumidity2mVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%SpecHumidity2mVeg)
    end if

    if ( .not. allocated(noahmp%energy%state%SpecHumidity2mBare) ) then
       allocate( noahmp%energy%state%SpecHumidity2mBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%SpecHumidity2mBare)
    end if

    if ( .not. allocated(noahmp%energy%state%SpecHumidity2m) ) then
       allocate( noahmp%energy%state%SpecHumidity2m(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%SpecHumidity2m)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureGrdVeg) ) then
       allocate( noahmp%energy%state%TemperatureGrdVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureGrdVeg)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureGrdBare) ) then
       allocate( noahmp%energy%state%TemperatureGrdBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureGrdBare)
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressEwVeg) ) then
       allocate( noahmp%energy%state%WindStressEwVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%WindStressEwVeg)
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressNsVeg) ) then
       allocate( noahmp%energy%state%WindStressNsVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%WindStressNsVeg)
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressEwBare) ) then
       allocate( noahmp%energy%state%WindStressEwBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%WindStressEwBare)
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressNsBare) ) then
       allocate( noahmp%energy%state%WindStressNsBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%WindStressNsBare)
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressEwSfc) ) then
       allocate( noahmp%energy%state%WindStressEwSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%WindStressEwSfc)
    end if

    if ( .not. allocated(noahmp%energy%state%WindStressNsSfc) ) then
       allocate( noahmp%energy%state%WindStressNsSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%WindStressNsSfc)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureAir2mVeg) ) then
       allocate( noahmp%energy%state%TemperatureAir2mVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureAir2mVeg)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureAir2mBare) ) then
       allocate( noahmp%energy%state%TemperatureAir2mBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureAir2mBare)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureAir2m) ) then
       allocate( noahmp%energy%state%TemperatureAir2m(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureAir2m)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffShLeaf) ) then
       allocate( noahmp%energy%state%ExchCoeffShLeaf(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffShLeaf)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffShUndCan) ) then
       allocate( noahmp%energy%state%ExchCoeffShUndCan(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffShUndCan)
    end if

    if ( .not. allocated(noahmp%energy%state%ExchCoeffSh2mBare) ) then
       allocate( noahmp%energy%state%ExchCoeffSh2mBare(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ExchCoeffSh2mBare)
    end if

    if ( .not. allocated(noahmp%energy%state%RefHeightAboveGrd) ) then
       allocate( noahmp%energy%state%RefHeightAboveGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%RefHeightAboveGrd)
    end if

    if ( .not. allocated(noahmp%energy%state%CanopyFracSnowBury) ) then
       allocate( noahmp%energy%state%CanopyFracSnowBury(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%CanopyFracSnowBury)
    end if

    if ( .not. allocated(noahmp%energy%state%DepthSoilTempBotToSno) ) then
       allocate( noahmp%energy%state%DepthSoilTempBotToSno(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%DepthSoilTempBotToSno)
    end if

    if ( .not. allocated(noahmp%energy%state%RoughLenMomSfcToAtm) ) then
       allocate( noahmp%energy%state%RoughLenMomSfcToAtm(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%RoughLenMomSfcToAtm)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureRadSfc) ) then
       allocate( noahmp%energy%state%TemperatureRadSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureRadSfc)
    end if

    if ( .not. allocated(noahmp%energy%state%EmissivitySfc) ) then
       allocate( noahmp%energy%state%EmissivitySfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%EmissivitySfc)
    end if

    if ( .not. allocated(noahmp%energy%state%AlbedoSfc) ) then
       allocate( noahmp%energy%state%AlbedoSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%AlbedoSfc)
    end if

    if ( .not. allocated(noahmp%energy%state%EnergyBalanceError) ) then
       allocate( noahmp%energy%state%EnergyBalanceError(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%EnergyBalanceError)
    end if

    if ( .not. allocated(noahmp%energy%state%RadSwBalanceError) ) then
       allocate( noahmp%energy%state%RadSwBalanceError(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%RadSwBalanceError)
    end if

    if ( .not. allocated(noahmp%energy%state%TemperatureSoilSnow) ) then
       allocate( noahmp%energy%state%TemperatureSoilSnow(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TemperatureSoilSnow)
    endif
    if ( .not. allocated(noahmp%energy%state%ThermConductSoilSnow) ) then
       allocate( noahmp%energy%state%ThermConductSoilSnow(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ThermConductSoilSnow)
    endif
    if ( .not. allocated(noahmp%energy%state%HeatCapacSoilSnow) ) then
       allocate( noahmp%energy%state%HeatCapacSoilSnow(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%HeatCapacSoilSnow)
    endif
    if ( .not. allocated(noahmp%energy%state%PhaseChgFacSoilSnow) ) then
       allocate( noahmp%energy%state%PhaseChgFacSoilSnow(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%PhaseChgFacSoilSnow)
    endif
    if ( .not. allocated(noahmp%energy%state%HeatCapacVolSnow) ) then
       allocate( noahmp%energy%state%HeatCapacVolSnow(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%HeatCapacVolSnow)
    endif
    if ( .not. allocated(noahmp%energy%state%ThermConductSnow) ) then
       allocate( noahmp%energy%state%ThermConductSnow(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ThermConductSnow)
    endif
    if ( .not. allocated(noahmp%energy%state%HeatCapacVolSoil) ) then
       allocate( noahmp%energy%state%HeatCapacVolSoil(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%HeatCapacVolSoil)
    endif
    if ( .not. allocated(noahmp%energy%state%ThermConductSoil) ) then
       allocate( noahmp%energy%state%ThermConductSoil(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ThermConductSoil)
    endif
    if ( .not. allocated(noahmp%energy%state%HeatCapacGlaIce) ) then
       allocate( noahmp%energy%state%HeatCapacGlaIce(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%HeatCapacGlaIce)
    endif
    if ( .not. allocated(noahmp%energy%state%ThermConductGlaIce) ) then
       allocate( noahmp%energy%state%ThermConductGlaIce(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ThermConductGlaIce)
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSnowDir) ) then
       allocate( noahmp%energy%state%AlbedoSnowDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%AlbedoSnowDir)
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSnowDif) ) then
       allocate( noahmp%energy%state%AlbedoSnowDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%AlbedoSnowDif)
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSoilDir) ) then
       allocate( noahmp%energy%state%AlbedoSoilDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%AlbedoSoilDir)
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSoilDif) ) then
       allocate( noahmp%energy%state%AlbedoSoilDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%AlbedoSoilDif)
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoGrdDir) ) then
       allocate( noahmp%energy%state%AlbedoGrdDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%AlbedoGrdDir)
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoGrdDif) ) then
       allocate( noahmp%energy%state%AlbedoGrdDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%AlbedoGrdDif)
    endif
    if ( .not. allocated(noahmp%energy%state%ReflectanceVeg) ) then
       allocate( noahmp%energy%state%ReflectanceVeg(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%ReflectanceVeg)
    endif
    if ( .not. allocated(noahmp%energy%state%TransmittanceVeg) ) then
       allocate( noahmp%energy%state%TransmittanceVeg(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%TransmittanceVeg)
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSfcDir) ) then
       allocate( noahmp%energy%state%AlbedoSfcDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%AlbedoSfcDir)
    endif
    if ( .not. allocated(noahmp%energy%state%AlbedoSfcDif) ) then
       allocate( noahmp%energy%state%AlbedoSfcDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%state%AlbedoSfcDif)
    endif

    ! Allocate 3D energy flux arrays and transfer to GPU
    if ( .not. allocated(noahmp%energy%flux%RadSwAbsVegDir) ) then
       allocate( noahmp%energy%flux%RadSwAbsVegDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwAbsVegDir)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwAbsVegDif) ) then
       allocate( noahmp%energy%flux%RadSwAbsVegDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwAbsVegDif)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDirTranGrdDir) ) then
       allocate( noahmp%energy%flux%RadSwDirTranGrdDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwDirTranGrdDir)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDirTranGrdDif) ) then
       allocate( noahmp%energy%flux%RadSwDirTranGrdDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwDirTranGrdDif)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDifTranGrdDir) ) then
       allocate( noahmp%energy%flux%RadSwDifTranGrdDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwDifTranGrdDir)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDifTranGrdDif) ) then
       allocate( noahmp%energy%flux%RadSwDifTranGrdDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwDifTranGrdDif)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflVegDir) ) then
       allocate( noahmp%energy%flux%RadSwReflVegDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwReflVegDir)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflVegDif) ) then
       allocate( noahmp%energy%flux%RadSwReflVegDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwReflVegDif)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflGrdDir) ) then
       allocate( noahmp%energy%flux%RadSwReflGrdDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwReflGrdDir)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwReflGrdDif) ) then
       allocate( noahmp%energy%flux%RadSwReflGrdDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwReflGrdDif)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDownDir) ) then
       allocate( noahmp%energy%flux%RadSwDownDir(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwDownDir)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwDownDif) ) then
       allocate( noahmp%energy%flux%RadSwDownDif(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwDownDif)
    endif
    if ( .not. allocated(noahmp%energy%flux%RadSwPenetrateGrd) ) then
       allocate( noahmp%energy%flux%RadSwPenetrateGrd(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%energy%flux%RadSwPenetrateGrd)
    endif

    ! SNICAR flux arrays
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( .not. allocated(noahmp%energy%flux%FracRadSwAbsSnowDir) ) then
          allocate( noahmp%energy%flux%FracRadSwAbsSnowDir(ITS:ITE,-NumSnowLayerMax+1:1,1:NumSwRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%flux%FracRadSwAbsSnowDir)
       endif
       if ( .not. allocated(noahmp%energy%flux%FracRadSwAbsSnowDif) ) then
          allocate( noahmp%energy%flux%FracRadSwAbsSnowDif(ITS:ITE,-NumSnowLayerMax+1:1,1:NumSwRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%flux%FracRadSwAbsSnowDif)
       endif
       if ( .not. allocated(noahmp%energy%flux%RadSwAbsSnowSoilLayer) ) then
          allocate( noahmp%energy%flux%RadSwAbsSnowSoilLayer(ITS:ITE,-NumSnowLayerMax+1:1,JTS:JTE) )
          !$acc enter data create(noahmp%energy%flux%RadSwAbsSnowSoilLayer)
       endif
    endif

    ! Allocate 3D energy parameter arrays and transfer to GPU
    if ( .not. allocated(noahmp%energy%param%LeafAreaIndexMon) ) then
       allocate( noahmp%energy%param%LeafAreaIndexMon(ITS:ITE,1:12,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%LeafAreaIndexMon)
    endif
    if ( .not. allocated(noahmp%energy%param%StemAreaIndexMon) ) then
       allocate( noahmp%energy%param%StemAreaIndexMon(ITS:ITE,1:12,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%StemAreaIndexMon)
    endif
    if ( .not. allocated(noahmp%energy%param%SoilQuartzFrac) ) then
       allocate( noahmp%energy%param%SoilQuartzFrac(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SoilQuartzFrac)
    endif
    if ( .not. allocated(noahmp%energy%param%AlbedoSoilSat) ) then
       allocate( noahmp%energy%param%AlbedoSoilSat(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%AlbedoSoilSat)
    endif
    if ( .not. allocated(noahmp%energy%param%AlbedoSoilDry) ) then
       allocate( noahmp%energy%param%AlbedoSoilDry(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%AlbedoSoilDry)
    endif
    if ( .not. allocated(noahmp%energy%param%AlbedoLakeFrz) ) then
       allocate( noahmp%energy%param%AlbedoLakeFrz(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%AlbedoLakeFrz)
    endif
    if ( .not. allocated(noahmp%energy%param%ScatterCoeffSnow) ) then
       allocate( noahmp%energy%param%ScatterCoeffSnow(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%ScatterCoeffSnow)
    endif
    if ( .not. allocated(noahmp%energy%param%ReflectanceLeaf) ) then
       allocate( noahmp%energy%param%ReflectanceLeaf(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%ReflectanceLeaf)
    endif
    if ( .not. allocated(noahmp%energy%param%ReflectanceStem) ) then
       allocate( noahmp%energy%param%ReflectanceStem(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%ReflectanceStem)
    endif
    if ( .not. allocated(noahmp%energy%param%TransmittanceLeaf) ) then
       allocate( noahmp%energy%param%TransmittanceLeaf(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%TransmittanceLeaf)
    endif
    if ( .not. allocated(noahmp%energy%param%TransmittanceStem) ) then
       allocate( noahmp%energy%param%TransmittanceStem(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%TransmittanceStem)
    endif
    if ( .not. allocated(noahmp%energy%param%EmissivitySoilLake) ) then
       allocate( noahmp%energy%param%EmissivitySoilLake(ITS:ITE,1:2,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%EmissivitySoilLake)
    endif
    if ( .not. allocated(noahmp%energy%param%AlbedoLandIce) ) then
       allocate( noahmp%energy%param%AlbedoLandIce(ITS:ITE,1:NumSwRadBand,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%AlbedoLandIce)
    endif

    ! SNICAR parameter arrays - these are lookup tables, not spatially varying
    ! Keep as 1D/2D since they are the same for all grid points
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( .not. allocated(noahmp%energy%param%RadSwWgtDir) ) then
          allocate( noahmp%energy%param%RadSwWgtDir(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%RadSwWgtDir)
       endif
       if ( .not. allocated(noahmp%energy%param%RadSwWgtDif) ) then
          allocate( noahmp%energy%param%RadSwWgtDif(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%RadSwWgtDif)
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbSnwRadDir) ) then
          allocate( noahmp%energy%param%SsAlbSnwRadDir(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%SsAlbSnwRadDir)
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmSnwRadDir) ) then
          allocate( noahmp%energy%param%AsyPrmSnwRadDir(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%AsyPrmSnwRadDir)
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassSnwRadDir) ) then
          allocate( noahmp%energy%param%ExtCffMassSnwRadDir(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%ExtCffMassSnwRadDir)
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbSnwRadDif) ) then
          allocate( noahmp%energy%param%SsAlbSnwRadDif(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%SsAlbSnwRadDif)
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmSnwRadDif) ) then
          allocate( noahmp%energy%param%AsyPrmSnwRadDif(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%AsyPrmSnwRadDif)
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassSnwRadDif) ) then
          allocate( noahmp%energy%param%ExtCffMassSnwRadDif(ITS:ITE,1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%ExtCffMassSnwRadDif)
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbBCphi) ) then
          allocate( noahmp%energy%param%SsAlbBCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%SsAlbBCphi)
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmBCphi) ) then
          allocate( noahmp%energy%param%AsyPrmBCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%AsyPrmBCphi)
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassBCphi) ) then
          allocate( noahmp%energy%param%ExtCffMassBCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%ExtCffMassBCphi)
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbBCpho) ) then
          allocate( noahmp%energy%param%SsAlbBCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%SsAlbBCpho)
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmBCpho) ) then
          allocate( noahmp%energy%param%AsyPrmBCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%AsyPrmBCpho)
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassBCpho) ) then
          allocate( noahmp%energy%param%ExtCffMassBCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%ExtCffMassBCpho)
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbOCphi) ) then
          allocate( noahmp%energy%param%SsAlbOCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%SsAlbOCphi)
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmOCphi) ) then
          allocate( noahmp%energy%param%AsyPrmOCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%AsyPrmOCphi)
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassOCphi) ) then
          allocate( noahmp%energy%param%ExtCffMassOCphi(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%ExtCffMassOCphi)
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbOCpho) ) then
          allocate( noahmp%energy%param%SsAlbOCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%SsAlbOCpho)
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmOCpho) ) then
          allocate( noahmp%energy%param%AsyPrmOCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%AsyPrmOCpho)
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassOCpho) ) then
          allocate( noahmp%energy%param%ExtCffMassOCpho(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%ExtCffMassOCpho)
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbDustB1) ) then
          allocate( noahmp%energy%param%SsAlbDustB1(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%SsAlbDustB1)
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmDustB1) ) then
          allocate( noahmp%energy%param%AsyPrmDustB1(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%AsyPrmDustB1)
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassDustB1) ) then
          allocate( noahmp%energy%param%ExtCffMassDustB1(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%ExtCffMassDustB1)
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbDustB2) ) then
          allocate( noahmp%energy%param%SsAlbDustB2(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%SsAlbDustB2)
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmDustB2) ) then
          allocate( noahmp%energy%param%AsyPrmDustB2(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%AsyPrmDustB2)
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassDustB2) ) then
          allocate( noahmp%energy%param%ExtCffMassDustB2(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%ExtCffMassDustB2)
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbDustB3) ) then
          allocate( noahmp%energy%param%SsAlbDustB3(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%SsAlbDustB3)
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmDustB3) ) then
          allocate( noahmp%energy%param%AsyPrmDustB3(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%AsyPrmDustB3)
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassDustB3) ) then
          allocate( noahmp%energy%param%ExtCffMassDustB3(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%ExtCffMassDustB3)
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbDustB4) ) then
          allocate( noahmp%energy%param%SsAlbDustB4(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%SsAlbDustB4)
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmDustB4) ) then
          allocate( noahmp%energy%param%AsyPrmDustB4(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%AsyPrmDustB4)
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassDustB4) ) then
          allocate( noahmp%energy%param%ExtCffMassDustB4(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%ExtCffMassDustB4)
       endif
       if ( .not. allocated(noahmp%energy%param%SsAlbDustB5) ) then
          allocate( noahmp%energy%param%SsAlbDustB5(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%SsAlbDustB5)
       endif
       if ( .not. allocated(noahmp%energy%param%AsyPrmDustB5) ) then
          allocate( noahmp%energy%param%AsyPrmDustB5(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%AsyPrmDustB5)
       endif
       if ( .not. allocated(noahmp%energy%param%ExtCffMassDustB5) ) then
          allocate( noahmp%energy%param%ExtCffMassDustB5(ITS:ITE,1:NumSnicarRadBand,JTS:JTE) )
          !$acc enter data create(noahmp%energy%param%ExtCffMassDustB5)
       endif

       ! Initialize SNICAR lookup tables (not spatially varying)
       !$acc parallel loop gang vector collapse(2) present(noahmp)
       do I = ITS, ITE
         do J = JTS, JTE
        !$acc loop seq
        do LoopInd = 1, NumSnicarRadBand
        !$acc loop seq
         do k = 1, NumRadiusSnwMieSnicar
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

    ! Allocate 2D energy param arrays and transfer to GPU
    if ( .not. allocated(noahmp%energy%param%TreeCrownRadius) ) then
       allocate( noahmp%energy%param%TreeCrownRadius(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%TreeCrownRadius)
    endif
    if ( .not. allocated(noahmp%energy%param%HeightCanopyTop) ) then
       allocate( noahmp%energy%param%HeightCanopyTop(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%HeightCanopyTop)
    endif
    if ( .not. allocated(noahmp%energy%param%HeightCanopyBot) ) then
       allocate( noahmp%energy%param%HeightCanopyBot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%HeightCanopyBot)
    endif
    if ( .not. allocated(noahmp%energy%param%RoughLenMomVeg) ) then
       allocate( noahmp%energy%param%RoughLenMomVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%RoughLenMomVeg)
    endif
    if ( .not. allocated(noahmp%energy%param%TreeDensity) ) then
       allocate( noahmp%energy%param%TreeDensity(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%TreeDensity)
    endif
    if ( .not. allocated(noahmp%energy%param%CanopyOrientIndex) ) then
       allocate( noahmp%energy%param%CanopyOrientIndex(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%CanopyOrientIndex)
    endif
    if ( .not. allocated(noahmp%energy%param%UpscatterCoeffSnowDir) ) then
       allocate( noahmp%energy%param%UpscatterCoeffSnowDir(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%UpscatterCoeffSnowDir)
    endif
    if ( .not. allocated(noahmp%energy%param%UpscatterCoeffSnowDif) ) then
       allocate( noahmp%energy%param%UpscatterCoeffSnowDif(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%UpscatterCoeffSnowDif)
    endif
    if ( .not. allocated(noahmp%energy%param%SoilHeatCapacity) ) then
       allocate( noahmp%energy%param%SoilHeatCapacity(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SoilHeatCapacity)
    endif
    if ( .not. allocated(noahmp%energy%param%SnowAgeFacBats) ) then
       allocate( noahmp%energy%param%SnowAgeFacBats(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SnowAgeFacBats)
    endif
    if ( .not. allocated(noahmp%energy%param%SnowGrowVapFacBats) ) then
       allocate( noahmp%energy%param%SnowGrowVapFacBats(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SnowGrowVapFacBats)
    endif
    if ( .not. allocated(noahmp%energy%param%SnowSootFacBats) ) then
       allocate( noahmp%energy%param%SnowSootFacBats(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SnowSootFacBats)
    endif
    if ( .not. allocated(noahmp%energy%param%SnowGrowFrzFacBats) ) then
       allocate( noahmp%energy%param%SnowGrowFrzFacBats(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SnowGrowFrzFacBats)
    endif
    if ( .not. allocated(noahmp%energy%param%SolarZenithAdjBats) ) then
       allocate( noahmp%energy%param%SolarZenithAdjBats(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SolarZenithAdjBats)
    endif
    if ( .not. allocated(noahmp%energy%param%FreshSnoAlbVisBats) ) then
       allocate( noahmp%energy%param%FreshSnoAlbVisBats(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%FreshSnoAlbVisBats)
    endif
    if ( .not. allocated(noahmp%energy%param%FreshSnoAlbNirBats) ) then
       allocate( noahmp%energy%param%FreshSnoAlbNirBats(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%FreshSnoAlbNirBats)
    endif
    if ( .not. allocated(noahmp%energy%param%SnoAgeFacDifVisBats) ) then
       allocate( noahmp%energy%param%SnoAgeFacDifVisBats(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SnoAgeFacDifVisBats)
    endif
    if ( .not. allocated(noahmp%energy%param%SnoAgeFacDifNirBats) ) then
       allocate( noahmp%energy%param%SnoAgeFacDifNirBats(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SnoAgeFacDifNirBats)
    endif
    if ( .not. allocated(noahmp%energy%param%SzaFacDirVisBats) ) then
       allocate( noahmp%energy%param%SzaFacDirVisBats(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SzaFacDirVisBats)
    endif
    if ( .not. allocated(noahmp%energy%param%SzaFacDirNirBats) ) then
       allocate( noahmp%energy%param%SzaFacDirNirBats(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SzaFacDirNirBats)
    endif
    if ( .not. allocated(noahmp%energy%param%SnowAlbRefClass) ) then
       allocate( noahmp%energy%param%SnowAlbRefClass(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SnowAlbRefClass)
    endif
    if ( .not. allocated(noahmp%energy%param%SnowAgeFacClass) ) then
       allocate( noahmp%energy%param%SnowAgeFacClass(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SnowAgeFacClass)
    endif
    if ( .not. allocated(noahmp%energy%param%SnowAlbFreshClass) ) then
       allocate( noahmp%energy%param%SnowAlbFreshClass(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%SnowAlbFreshClass)
    endif
    if ( .not. allocated(noahmp%energy%param%ConductanceLeafMin) ) then
       allocate( noahmp%energy%param%ConductanceLeafMin(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%ConductanceLeafMin)
    endif
    if ( .not. allocated(noahmp%energy%param%Co2MmConst25C) ) then
       allocate( noahmp%energy%param%Co2MmConst25C(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%Co2MmConst25C)
    endif
    if ( .not. allocated(noahmp%energy%param%O2MmConst25C) ) then
       allocate( noahmp%energy%param%O2MmConst25C(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%O2MmConst25C)
    endif
    if ( .not. allocated(noahmp%energy%param%Co2MmConstQ10) ) then
       allocate( noahmp%energy%param%Co2MmConstQ10(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%Co2MmConstQ10)
    endif
    if ( .not. allocated(noahmp%energy%param%O2MmConstQ10) ) then
       allocate( noahmp%energy%param%O2MmConstQ10(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%O2MmConstQ10)
    endif
    if ( .not. allocated(noahmp%energy%param%RadiationStressFac) ) then
       allocate( noahmp%energy%param%RadiationStressFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%RadiationStressFac)
    endif
    if ( .not. allocated(noahmp%energy%param%ResistanceStomataMin) ) then
       allocate( noahmp%energy%param%ResistanceStomataMin(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%ResistanceStomataMin)
    endif
    if ( .not. allocated(noahmp%energy%param%ResistanceStomataMax) ) then
       allocate( noahmp%energy%param%ResistanceStomataMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%ResistanceStomataMax)
    endif
    if ( .not. allocated(noahmp%energy%param%AirTempOptimTransp) ) then
       allocate( noahmp%energy%param%AirTempOptimTransp(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%AirTempOptimTransp)
    endif
    if ( .not. allocated(noahmp%energy%param%VaporPresDeficitFac) ) then
       allocate( noahmp%energy%param%VaporPresDeficitFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%VaporPresDeficitFac)
    endif
    if ( .not. allocated(noahmp%energy%param%LeafDimLength) ) then
       allocate( noahmp%energy%param%LeafDimLength(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%LeafDimLength)
    endif
    if ( .not. allocated(noahmp%energy%param%ZilitinkevichCoeff) ) then
       allocate( noahmp%energy%param%ZilitinkevichCoeff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%ZilitinkevichCoeff)
    endif
    if ( .not. allocated(noahmp%energy%param%EmissivitySnow) ) then
       allocate( noahmp%energy%param%EmissivitySnow(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%EmissivitySnow)
    endif
    if ( .not. allocated(noahmp%energy%param%CanopyWindExtFac) ) then
       allocate( noahmp%energy%param%CanopyWindExtFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%CanopyWindExtFac)
    endif
    if ( .not. allocated(noahmp%energy%param%RoughLenMomSnow) ) then
       allocate( noahmp%energy%param%RoughLenMomSnow(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%RoughLenMomSnow)
    endif
    if ( .not. allocated(noahmp%energy%param%RoughLenMomSoil) ) then
       allocate( noahmp%energy%param%RoughLenMomSoil(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%RoughLenMomSoil)
    endif
    if ( .not. allocated(noahmp%energy%param%RoughLenMomLake) ) then
       allocate( noahmp%energy%param%RoughLenMomLake(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%RoughLenMomLake)
    endif
    if ( .not. allocated(noahmp%energy%param%EmissivityIceSfc) ) then
       allocate( noahmp%energy%param%EmissivityIceSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%EmissivityIceSfc)
    endif
    if ( .not. allocated(noahmp%energy%param%ResistanceSoilExp) ) then
       allocate( noahmp%energy%param%ResistanceSoilExp(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%ResistanceSoilExp)
    endif
    if ( .not. allocated(noahmp%energy%param%ResistanceSnowSfc) ) then
       allocate( noahmp%energy%param%ResistanceSnowSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%ResistanceSnowSfc)
    endif
    if ( .not. allocated(noahmp%energy%param%VegFracGreen) ) then
       allocate( noahmp%energy%param%VegFracGreen(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%VegFracGreen)
    endif
    if ( .not. allocated(noahmp%energy%param%VegFracAnnMax) ) then
       allocate( noahmp%energy%param%VegFracAnnMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%VegFracAnnMax)
    endif
    if ( .not. allocated(noahmp%energy%param%HeatCapacCanFac) ) then
       allocate( noahmp%energy%param%HeatCapacCanFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%energy%param%HeatCapacCanFac)
    endif

    end associate

    ! Now initialize all 2D and 3D arrays in parallel loop
    !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
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

end module EnergyVarInitMod
