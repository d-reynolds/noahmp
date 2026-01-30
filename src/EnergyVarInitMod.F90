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

    ! Allocate 3D energy state arrays and transfer to GPU
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
        !$acc loop seq
        do LoopInd = 1, NumSnicarRadBand
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
