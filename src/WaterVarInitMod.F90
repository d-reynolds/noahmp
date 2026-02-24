module WaterVarInitMod

!!! Initialize column (1-D) Noah-MP water variables
!!! Water variables should be first defined in WaterVarType.F90

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
  subroutine WaterVarInitDefault(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer :: I, J      ! grid indices
    integer :: LoopInd   ! loop index

    ! Domain and layer bounds for allocations
    associate(                                                                         &
              ITS                     => noahmp%config%domain%ITS                     ,&
              ITE                     => noahmp%config%domain%ITE                     ,&
              JTS                     => noahmp%config%domain%JTS                     ,&
              JTE                     => noahmp%config%domain%JTE                     ,&
              NumSnowLayerMax         => noahmp%config%domain%NumSnowLayerMax         ,&
              NumSoilLayer            => noahmp%config%domain%NumSoilLayer            ,&
              NumDensitySnwAgeSnicar  => noahmp%config%domain%NumDensitySnwAgeSnicar  ,&
              NumTempGradSnwAgeSnicar => noahmp%config%domain%NumTempGradSnwAgeSnicar ,&
              NumTempSnwAgeSnicar     => noahmp%config%domain%NumTempSnwAgeSnicar      &
             )

    ! Allocate 3D water state arrays and transfer to GPU
    if ( .not. allocated(noahmp%water%state%IndexPhaseChange) ) then
       allocate( noahmp%water%state%IndexPhaseChange(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilSupercoolWater) ) then
       allocate( noahmp%water%state%SoilSupercoolWater(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowIce) ) then
       allocate( noahmp%water%state%SnowIce(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowLiqWater) ) then
       allocate( noahmp%water%state%SnowLiqWater(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowIceVol) ) then
       allocate( noahmp%water%state%SnowIceVol(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowLiqWaterVol) ) then
       allocate( noahmp%water%state%SnowLiqWaterVol(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowIceFracPrev) ) then
       allocate( noahmp%water%state%SnowIceFracPrev(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowIceFrac) ) then
       allocate( noahmp%water%state%SnowIceFrac(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowEffPorosity) ) then
       allocate( noahmp%water%state%SnowEffPorosity(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilLiqWater) ) then
       allocate( noahmp%water%state%SoilLiqWater(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilIce) ) then
       allocate( noahmp%water%state%SoilIce(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilMoisture) ) then
       allocate( noahmp%water%state%SoilMoisture(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilImpervFrac) ) then
       allocate( noahmp%water%state%SoilImpervFrac(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilWatConductivity) ) then
       allocate( noahmp%water%state%SoilWatConductivity(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilWatDiffusivity) ) then
       allocate( noahmp%water%state%SoilWatDiffusivity(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilEffPorosity) ) then
       allocate( noahmp%water%state%SoilEffPorosity(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilIceFrac) ) then
       allocate( noahmp%water%state%SoilIceFrac(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilMoistureEqui) ) then
       allocate( noahmp%water%state%SoilMoistureEqui(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilTranspFac) ) then
       allocate( noahmp%water%state%SoilTranspFac(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilMatPotential) ) then
       allocate( noahmp%water%state%SoilMatPotential(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif

    ! SNICAR state arrays
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( .not. allocated(noahmp%water%state%SnowRadius) ) then
          allocate( noahmp%water%state%SnowRadius(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassBChydropho) ) then
          allocate( noahmp%water%state%MassBChydropho(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassBChydrophi) ) then
          allocate( noahmp%water%state%MassBChydrophi(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassOChydropho) ) then
          allocate( noahmp%water%state%MassOChydropho(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassOChydrophi) ) then
          allocate( noahmp%water%state%MassOChydrophi(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassDust1) ) then
          allocate( noahmp%water%state%MassDust1(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassDust2) ) then
          allocate( noahmp%water%state%MassDust2(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassDust3) ) then
          allocate( noahmp%water%state%MassDust3(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassDust4) ) then
          allocate( noahmp%water%state%MassDust4(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassDust5) ) then
          allocate( noahmp%water%state%MassDust5(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassConcBChydropho) ) then
          allocate( noahmp%water%state%MassConcBChydropho(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassConcBChydrophi) ) then
          allocate( noahmp%water%state%MassConcBChydrophi(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassConcOChydropho) ) then
          allocate( noahmp%water%state%MassConcOChydropho(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassConcOChydrophi) ) then
          allocate( noahmp%water%state%MassConcOChydrophi(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassConcDust1) ) then
          allocate( noahmp%water%state%MassConcDust1(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassConcDust2) ) then
          allocate( noahmp%water%state%MassConcDust2(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassConcDust3) ) then
          allocate( noahmp%water%state%MassConcDust3(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassConcDust4) ) then
          allocate( noahmp%water%state%MassConcDust4(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
       if ( .not. allocated(noahmp%water%state%MassConcDust5) ) then
          allocate( noahmp%water%state%MassConcDust5(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
    endif

    ! Allocate 3D water flux arrays and transfer to GPU
    if ( .not. allocated(noahmp%water%flux%CompactionSnowAging) ) then
       allocate( noahmp%water%flux%CompactionSnowAging(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%CompactionSnowBurden) ) then
       allocate( noahmp%water%flux%CompactionSnowBurden(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%CompactionSnowMelt) ) then
       allocate( noahmp%water%flux%CompactionSnowMelt(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%CompactionSnowTot) ) then
       allocate( noahmp%water%flux%CompactionSnowTot(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%TranspWatLossSoil) ) then
       allocate( noahmp%water%flux%TranspWatLossSoil(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%TranspWatLossSoilAcc) ) then
       allocate( noahmp%water%flux%TranspWatLossSoilAcc(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%TranspWatLossSoilMean) ) then
       allocate( noahmp%water%flux%TranspWatLossSoilMean(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%OutflowSnowLayer) ) then
       allocate( noahmp%water%flux%OutflowSnowLayer(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
    endif

    ! SNICAR flux arrays
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( .not. allocated(noahmp%water%flux%SnowFreezeRate) ) then
          allocate( noahmp%water%flux%SnowFreezeRate(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       endif
    endif

    ! Allocate 3D water parameter arrays and transfer to GPU
    if ( .not. allocated(noahmp%water%param%SoilMoistureSat) ) then
       allocate( noahmp%water%param%SoilMoistureSat(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilMoistureWilt) ) then
       allocate( noahmp%water%param%SoilMoistureWilt(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilMoistureFieldCap) ) then
       allocate( noahmp%water%param%SoilMoistureFieldCap(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilMoistureDry) ) then
       allocate( noahmp%water%param%SoilMoistureDry(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilWatDiffusivitySat) ) then
       allocate( noahmp%water%param%SoilWatDiffusivitySat(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilWatConductivitySat) ) then
       allocate( noahmp%water%param%SoilWatConductivitySat(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilExpCoeffB) ) then
       allocate( noahmp%water%param%SoilExpCoeffB(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilMatPotentialSat) ) then
       allocate( noahmp%water%param%SoilMatPotentialSat(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
    endif

    ! SNICAR parameter lookup tables - keep as 3D (not spatially varying)
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( .not. allocated(noahmp%water%param%snowage_tau) ) then
          allocate( noahmp%water%param%snowage_tau(NumDensitySnwAgeSnicar,NumTempGradSnwAgeSnicar,NumTempSnwAgeSnicar) )
       endif
       if ( .not. allocated(noahmp%water%param%snowage_kappa) ) then
          allocate( noahmp%water%param%snowage_kappa(NumDensitySnwAgeSnicar,NumTempGradSnwAgeSnicar,NumTempSnwAgeSnicar) )
       endif
       if ( .not. allocated(noahmp%water%param%snowage_drdt0) ) then
          allocate( noahmp%water%param%snowage_drdt0(NumDensitySnwAgeSnicar,NumTempGradSnwAgeSnicar,NumTempSnwAgeSnicar) )
       endif
    endif

    ! Allocate 2D water state arrays and transfer to GPU
    if ( .not. allocated(noahmp%water%state%IrrigationCntSprinkler) ) then
       allocate( noahmp%water%state%IrrigationCntSprinkler(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationCntMicro) ) then
       allocate( noahmp%water%state%IrrigationCntMicro(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationCntFlood) ) then
       allocate( noahmp%water%state%IrrigationCntFlood(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%CanopyTotalWater) ) then
       allocate( noahmp%water%state%CanopyTotalWater(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%CanopyWetFrac) ) then
       allocate( noahmp%water%state%CanopyWetFrac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowfallDensity) ) then
       allocate( noahmp%water%state%SnowfallDensity(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%CanopyLiqWater) ) then
       allocate( noahmp%water%state%CanopyLiqWater(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%CanopyIce) ) then
       allocate( noahmp%water%state%CanopyIce(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%CanopyIceMax) ) then
       allocate( noahmp%water%state%CanopyIceMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%CanopyLiqWaterMax) ) then
       allocate( noahmp%water%state%CanopyLiqWaterMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowDepth) ) then
       allocate( noahmp%water%state%SnowDepth(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowWaterEquiv) ) then
       allocate( noahmp%water%state%SnowWaterEquiv(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowWaterEquivPrev) ) then
       allocate( noahmp%water%state%SnowWaterEquivPrev(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%PondSfcThinSnwMelt) ) then
       allocate( noahmp%water%state%PondSfcThinSnwMelt(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%PondSfcThinSnwComb) ) then
       allocate( noahmp%water%state%PondSfcThinSnwComb(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%PondSfcThinSnwTrans) ) then
       allocate( noahmp%water%state%PondSfcThinSnwTrans(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationFracFlood) ) then
       allocate( noahmp%water%state%IrrigationFracFlood(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationAmtFlood) ) then
       allocate( noahmp%water%state%IrrigationAmtFlood(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationFracMicro) ) then
       allocate( noahmp%water%state%IrrigationFracMicro(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationAmtMicro) ) then
       allocate( noahmp%water%state%IrrigationAmtMicro(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationFracSprinkler) ) then
       allocate( noahmp%water%state%IrrigationFracSprinkler(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationAmtSprinkler) ) then
       allocate( noahmp%water%state%IrrigationAmtSprinkler(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%WaterTableDepth) ) then
       allocate( noahmp%water%state%WaterTableDepth(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilIceMax) ) then
       allocate( noahmp%water%state%SoilIceMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilLiqWaterMin) ) then
       allocate( noahmp%water%state%SoilLiqWaterMin(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilSaturateFrac) ) then
       allocate( noahmp%water%state%SoilSaturateFrac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilImpervFracMax) ) then
       allocate( noahmp%water%state%SoilImpervFracMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilMoistureToWT) ) then
       allocate( noahmp%water%state%SoilMoistureToWT(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%RechargeGwDeepWT) ) then
       allocate( noahmp%water%state%RechargeGwDeepWT(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%RechargeGwShallowWT) ) then
       allocate( noahmp%water%state%RechargeGwShallowWT(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilSaturationExcess) ) then
       allocate( noahmp%water%state%SoilSaturationExcess(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%WaterTableHydro) ) then
       allocate( noahmp%water%state%WaterTableHydro(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%TileDrainFrac) ) then
       allocate( noahmp%water%state%TileDrainFrac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageAquifer) ) then
       allocate( noahmp%water%state%WaterStorageAquifer(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageSoilAqf) ) then
       allocate( noahmp%water%state%WaterStorageSoilAqf(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageLake) ) then
       allocate( noahmp%water%state%WaterStorageLake(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageWetland) ) then
       allocate( noahmp%water%state%WaterStorageWetland(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%WaterHeadSfc) ) then
       allocate( noahmp%water%state%WaterHeadSfc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationFracGrid) ) then
       allocate( noahmp%water%state%IrrigationFracGrid(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%PrecipAreaFrac) ) then
       allocate( noahmp%water%state%PrecipAreaFrac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowCoverFrac) ) then
       allocate( noahmp%water%state%SnowCoverFrac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilTranspFacAcc) ) then
       allocate( noahmp%water%state%SoilTranspFacAcc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%FrozenPrecipFrac) ) then
       allocate( noahmp%water%state%FrozenPrecipFrac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilWaterRootZone) ) then
       allocate( noahmp%water%state%SoilWaterRootZone(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SoilWaterStress) ) then
       allocate( noahmp%water%state%SoilWaterStress(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageTotBeg) ) then
       allocate( noahmp%water%state%WaterStorageTotBeg(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%WaterBalanceError) ) then
       allocate( noahmp%water%state%WaterBalanceError(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageTotEnd) ) then
       allocate( noahmp%water%state%WaterStorageTotEnd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%state%SnowRadiusFresh) ) then
       allocate( noahmp%water%state%SnowRadiusFresh(ITS:ITE,JTS:JTE) )
    endif

    ! Allocate 2D water flux arrays and transfer to GPU
    if ( .not. allocated(noahmp%water%flux%RainfallRefHeight) ) then
       allocate( noahmp%water%flux%RainfallRefHeight(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%SnowfallRefHeight) ) then
       allocate( noahmp%water%flux%SnowfallRefHeight(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%PrecipTotRefHeight) ) then
       allocate( noahmp%water%flux%PrecipTotRefHeight(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%PrecipConvTotRefHeight) ) then
       allocate( noahmp%water%flux%PrecipConvTotRefHeight(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%PrecipLargeSclRefHeight) ) then
       allocate( noahmp%water%flux%PrecipLargeSclRefHeight(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%EvapCanopyNet) ) then
       allocate( noahmp%water%flux%EvapCanopyNet(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%Transpiration) ) then
       allocate( noahmp%water%flux%Transpiration(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%EvapCanopyLiq) ) then
       allocate( noahmp%water%flux%EvapCanopyLiq(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%DewCanopyLiq) ) then
       allocate( noahmp%water%flux%DewCanopyLiq(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%FrostCanopyIce) ) then
       allocate( noahmp%water%flux%FrostCanopyIce(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%SublimCanopyIce) ) then
       allocate( noahmp%water%flux%SublimCanopyIce(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%MeltCanopyIce) ) then
       allocate( noahmp%water%flux%MeltCanopyIce(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%FreezeCanopyLiq) ) then
       allocate( noahmp%water%flux%FreezeCanopyLiq(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%SnowfallGround) ) then
       allocate( noahmp%water%flux%SnowfallGround(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%SnowDepthIncr) ) then
       allocate( noahmp%water%flux%SnowDepthIncr(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%FrostSnowSfcIce) ) then
       allocate( noahmp%water%flux%FrostSnowSfcIce(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%SublimSnowSfcIce) ) then
       allocate( noahmp%water%flux%SublimSnowSfcIce(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%RainfallGround) ) then
       allocate( noahmp%water%flux%RainfallGround(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%SnowBotOutflow) ) then
       allocate( noahmp%water%flux%SnowBotOutflow(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%GlacierExcessFlow) ) then
       allocate( noahmp%water%flux%GlacierExcessFlow(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%IrrigationRateFlood) ) then
       allocate( noahmp%water%flux%IrrigationRateFlood(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%IrrigationRateMicro) ) then
       allocate( noahmp%water%flux%IrrigationRateMicro(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%IrrigationRateSprinkler) ) then
       allocate( noahmp%water%flux%IrrigationRateSprinkler(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%IrriEvapLossSprinkler) ) then
       allocate( noahmp%water%flux%IrriEvapLossSprinkler(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%SoilSfcInflow) ) then
       allocate( noahmp%water%flux%SoilSfcInflow(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%RunoffSurface) ) then
       allocate( noahmp%water%flux%RunoffSurface(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%RunoffSubsurface) ) then
       allocate( noahmp%water%flux%RunoffSubsurface(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%InfilRateSfc) ) then
       allocate( noahmp%water%flux%InfilRateSfc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%EvapSoilSfcLiq) ) then
       allocate( noahmp%water%flux%EvapSoilSfcLiq(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%DrainSoilBot) ) then
       allocate( noahmp%water%flux%DrainSoilBot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%TileDrain) ) then
       allocate( noahmp%water%flux%TileDrain(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%RechargeGw) ) then
       allocate( noahmp%water%flux%RechargeGw(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%DischargeGw) ) then
       allocate( noahmp%water%flux%DischargeGw(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%VaporizeGrd) ) then
       allocate( noahmp%water%flux%VaporizeGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%CondenseVapGrd) ) then
       allocate( noahmp%water%flux%CondenseVapGrd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%DewSoilSfcLiq) ) then
       allocate( noahmp%water%flux%DewSoilSfcLiq(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%EvapIrriSprinkler) ) then
       allocate( noahmp%water%flux%EvapIrriSprinkler(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%InterceptCanopyRain) ) then
       allocate( noahmp%water%flux%InterceptCanopyRain(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%DripCanopyRain) ) then
       allocate( noahmp%water%flux%DripCanopyRain(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%ThroughfallRain) ) then
       allocate( noahmp%water%flux%ThroughfallRain(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%InterceptCanopySnow) ) then
       allocate( noahmp%water%flux%InterceptCanopySnow(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%DripCanopySnow) ) then
       allocate( noahmp%water%flux%DripCanopySnow(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%ThroughfallSnow) ) then
       allocate( noahmp%water%flux%ThroughfallSnow(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%EvapGroundNet) ) then
       allocate( noahmp%water%flux%EvapGroundNet(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%MeltGroundSnow) ) then
       allocate( noahmp%water%flux%MeltGroundSnow(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%WaterToAtmosTotal) ) then
       allocate( noahmp%water%flux%WaterToAtmosTotal(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%EvapSoilSfcLiqAcc) ) then
       allocate( noahmp%water%flux%EvapSoilSfcLiqAcc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%SoilSfcInflowAcc) ) then
       allocate( noahmp%water%flux%SoilSfcInflowAcc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%SfcWaterTotChgAcc) ) then
       allocate( noahmp%water%flux%SfcWaterTotChgAcc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%PrecipTotAcc) ) then
       allocate( noahmp%water%flux%PrecipTotAcc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%EvapCanopyNetAcc) ) then
       allocate( noahmp%water%flux%EvapCanopyNetAcc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%TranspirationAcc) ) then
       allocate( noahmp%water%flux%TranspirationAcc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%EvapGroundNetAcc) ) then
       allocate( noahmp%water%flux%EvapGroundNetAcc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%GlacierExcessFlowAcc) ) then
       allocate( noahmp%water%flux%GlacierExcessFlowAcc(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%EvapSoilSfcLiqMean) ) then
       allocate( noahmp%water%flux%EvapSoilSfcLiqMean(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%flux%SoilSfcInflowMean) ) then
       allocate( noahmp%water%flux%SoilSfcInflowMean(ITS:ITE,JTS:JTE) )
    endif

    ! Allocate 2D water param arrays and transfer to GPU
    if ( .not. allocated(noahmp%water%param%SnowCoverFac) ) then
       allocate( noahmp%water%param%SnowCoverFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%DrainSoilLayerInd) ) then
       allocate( noahmp%water%param%DrainSoilLayerInd(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%TileDrainTubeDepth) ) then
       allocate( noahmp%water%param%TileDrainTubeDepth(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%NumSoilLayerRoot) ) then
       allocate( noahmp%water%param%NumSoilLayerRoot(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%IrriStopDayBfHarvest) ) then
       allocate( noahmp%water%param%IrriStopDayBfHarvest(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%CanopyLiqHoldCap) ) then
       allocate( noahmp%water%param%CanopyLiqHoldCap(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactBurdenFac) ) then
       allocate( noahmp%water%param%SnowCompactBurdenFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactAgingFac1) ) then
       allocate( noahmp%water%param%SnowCompactAgingFac1(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactAgingFac2) ) then
       allocate( noahmp%water%param%SnowCompactAgingFac2(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactAgingFac3) ) then
       allocate( noahmp%water%param%SnowCompactAgingFac3(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactAgingMax) ) then
       allocate( noahmp%water%param%SnowCompactAgingMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowViscosityCoeff) ) then
       allocate( noahmp%water%param%SnowViscosityCoeff(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactmAR24) ) then
       allocate( noahmp%water%param%SnowCompactmAR24(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactbAR24) ) then
       allocate( noahmp%water%param%SnowCompactbAR24(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactP1AR24) ) then
       allocate( noahmp%water%param%SnowCompactP1AR24(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactP2AR24) ) then
       allocate( noahmp%water%param%SnowCompactP2AR24(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactP3AR24) ) then
       allocate( noahmp%water%param%SnowCompactP3AR24(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%BurdenFacUpAR24) ) then
       allocate( noahmp%water%param%BurdenFacUpAR24(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCoverM1AR25) ) then
       allocate( noahmp%water%param%SnowCoverM1AR25(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCoverM2AR25) ) then
       allocate( noahmp%water%param%SnowCoverM2AR25(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCoverFac1AR25) ) then
       allocate( noahmp%water%param%SnowCoverFac1AR25(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowCoverFac2AR25) ) then
       allocate( noahmp%water%param%SnowCoverFac2AR25(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowLiqFracMax) ) then
       allocate( noahmp%water%param%SnowLiqFracMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowLiqHoldCap) ) then
       allocate( noahmp%water%param%SnowLiqHoldCap(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowLiqReleaseFac) ) then
       allocate( noahmp%water%param%SnowLiqReleaseFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%IrriFloodRateFac) ) then
       allocate( noahmp%water%param%IrriFloodRateFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%IrriMicroRate) ) then
       allocate( noahmp%water%param%IrriMicroRate(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilInfilMaxCoeff) ) then
       allocate( noahmp%water%param%SoilInfilMaxCoeff(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilImpervFracCoeff) ) then
       allocate( noahmp%water%param%SoilImpervFracCoeff(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%InfilFacVic) ) then
       allocate( noahmp%water%param%InfilFacVic(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%TensionWatDistrInfl) ) then
       allocate( noahmp%water%param%TensionWatDistrInfl(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%TensionWatDistrShp) ) then
       allocate( noahmp%water%param%TensionWatDistrShp(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%FreeWatDistrShp) ) then
       allocate( noahmp%water%param%FreeWatDistrShp(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%InfilHeteroDynVic) ) then
       allocate( noahmp%water%param%InfilHeteroDynVic(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%InfilCapillaryDynVic) ) then
       allocate( noahmp%water%param%InfilCapillaryDynVic(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%InfilFacDynVic) ) then
       allocate( noahmp%water%param%InfilFacDynVic(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilDrainSlope) ) then
       allocate( noahmp%water%param%SoilDrainSlope(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%TileDrainCoeffSp) ) then
       allocate( noahmp%water%param%TileDrainCoeffSp(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%DrainFacSoilWat) ) then
       allocate( noahmp%water%param%DrainFacSoilWat(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%TileDrainCoeff) ) then
       allocate( noahmp%water%param%TileDrainCoeff(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%DrainDepthToImperv) ) then
       allocate( noahmp%water%param%DrainDepthToImperv(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%LateralWatCondFac) ) then
       allocate( noahmp%water%param%LateralWatCondFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%TileDrainDepth) ) then
       allocate( noahmp%water%param%TileDrainDepth(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%DrainTubeDist) ) then
       allocate( noahmp%water%param%DrainTubeDist(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%DrainTubeRadius) ) then
       allocate( noahmp%water%param%DrainTubeRadius(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%DrainWatDepToImperv) ) then
       allocate( noahmp%water%param%DrainWatDepToImperv(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%RunoffDecayFac) ) then
       allocate( noahmp%water%param%RunoffDecayFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%BaseflowCoeff) ) then
       allocate( noahmp%water%param%BaseflowCoeff(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%GridTopoIndex) ) then
       allocate( noahmp%water%param%GridTopoIndex(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilSfcSatFracMax) ) then
       allocate( noahmp%water%param%SoilSfcSatFracMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SpecYieldGw) ) then
       allocate( noahmp%water%param%SpecYieldGw(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%MicroPoreContent) ) then
       allocate( noahmp%water%param%MicroPoreContent(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%WaterStorageLakeMax) ) then
       allocate( noahmp%water%param%WaterStorageLakeMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnoWatEqvMaxGlacier) ) then
       allocate( noahmp%water%param%SnoWatEqvMaxGlacier(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilConductivityRef) ) then
       allocate( noahmp%water%param%SoilConductivityRef(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilInfilFacRef) ) then
       allocate( noahmp%water%param%SoilInfilFacRef(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%GroundFrzCoeff) ) then
       allocate( noahmp%water%param%GroundFrzCoeff(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%IrriTriggerLaiMin) ) then
       allocate( noahmp%water%param%IrriTriggerLaiMin(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilWatDeficitAllow) ) then
       allocate( noahmp%water%param%SoilWatDeficitAllow(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%IrriFloodLossFrac) ) then
       allocate( noahmp%water%param%IrriFloodLossFrac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%IrriSprinklerRate) ) then
       allocate( noahmp%water%param%IrriSprinklerRate(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%IrriFracThreshold) ) then
       allocate( noahmp%water%param%IrriFracThreshold(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%IrriStopPrecipThr) ) then
       allocate( noahmp%water%param%IrriStopPrecipThr(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowfallDensityMax) ) then
       allocate( noahmp%water%param%SnowfallDensityMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowMassFullCoverOld) ) then
       allocate( noahmp%water%param%SnowMassFullCoverOld(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SoilMatPotentialWilt) ) then
       allocate( noahmp%water%param%SoilMatPotentialWilt(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowMeltFac) ) then
       allocate( noahmp%water%param%SnowMeltFac(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%WetlandCapMax) ) then
       allocate( noahmp%water%param%WetlandCapMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowRadiusMin) ) then
       allocate( noahmp%water%param%SnowRadiusMin(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%FreshSnowRadiusMax) ) then
       allocate( noahmp%water%param%FreshSnowRadiusMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowRadiusRefrz) ) then
       allocate( noahmp%water%param%SnowRadiusRefrz(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltScale) ) then
       allocate( noahmp%water%param%ScavEffMeltScale(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltBCphi) ) then
       allocate( noahmp%water%param%ScavEffMeltBCphi(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltBCpho) ) then
       allocate( noahmp%water%param%ScavEffMeltBCpho(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltOCphi) ) then
       allocate( noahmp%water%param%ScavEffMeltOCphi(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltOCpho) ) then
       allocate( noahmp%water%param%ScavEffMeltOCpho(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltDust1) ) then
       allocate( noahmp%water%param%ScavEffMeltDust1(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltDust2) ) then
       allocate( noahmp%water%param%ScavEffMeltDust2(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltDust3) ) then
       allocate( noahmp%water%param%ScavEffMeltDust3(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltDust4) ) then
       allocate( noahmp%water%param%ScavEffMeltDust4(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltDust5) ) then
       allocate( noahmp%water%param%ScavEffMeltDust5(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowRadiusMax) ) then
       allocate( noahmp%water%param%SnowRadiusMax(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowWetAgeC1Brun89) ) then
       allocate( noahmp%water%param%SnowWetAgeC1Brun89(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowWetAgeC2Brun89) ) then
       allocate( noahmp%water%param%SnowWetAgeC2Brun89(ITS:ITE,JTS:JTE) )
    endif
    if ( .not. allocated(noahmp%water%param%SnowAgeScaleFac) ) then
       allocate( noahmp%water%param%SnowAgeScaleFac(ITS:ITE,JTS:JTE) )
    endif

    end associate

    ! Initialize all 2D and 3D arrays in parallel loop

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%water%state%IndexPhaseChange,   &
    !$acc   noahmp%water%state%SoilSupercoolWater,   &
    !$acc   noahmp%water%state%SnowIce,   &
    !$acc   noahmp%water%state%SnowLiqWater,   &
    !$acc   noahmp%water%state%SnowIceVol,   &
    !$acc   noahmp%water%state%SnowLiqWaterVol,   &
    !$acc   noahmp%water%state%SnowIceFracPrev,   &
    !$acc   noahmp%water%state%SnowIceFrac,   &
    !$acc   noahmp%water%state%SnowEffPorosity,   &
    !$acc   noahmp%water%state%SoilLiqWater,   &
    !$acc   noahmp%water%state%SoilIce,   &
    !$acc   noahmp%water%state%SoilMoisture,   &
    !$acc   noahmp%water%state%SoilImpervFrac,   &
    !$acc   noahmp%water%state%SoilWatConductivity,   &
    !$acc   noahmp%water%state%SoilWatDiffusivity,   &
    !$acc   noahmp%water%state%SoilEffPorosity,   &
    !$acc   noahmp%water%state%SoilIceFrac,   &
    !$acc   noahmp%water%state%SoilMoistureEqui,   &
    !$acc   noahmp%water%state%SoilTranspFac,   &
    !$acc   noahmp%water%state%SoilMatPotential,   &
    !$acc   noahmp%water%flux%CompactionSnowAging,   &
    !$acc   noahmp%water%flux%CompactionSnowBurden,   &
    !$acc   noahmp%water%flux%CompactionSnowMelt,   &
    !$acc   noahmp%water%flux%CompactionSnowTot,   &
    !$acc   noahmp%water%flux%TranspWatLossSoil,   &
    !$acc   noahmp%water%flux%TranspWatLossSoilAcc,   &
    !$acc   noahmp%water%flux%TranspWatLossSoilMean,   &
    !$acc   noahmp%water%flux%OutflowSnowLayer,   &
    !$acc   noahmp%water%param%SoilMoistureSat,   &
    !$acc   noahmp%water%param%SoilMoistureWilt    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%water%param%SoilMoistureFieldCap,   &
    !$acc   noahmp%water%param%SoilMoistureDry,   &
    !$acc   noahmp%water%param%SoilWatDiffusivitySat,   &
    !$acc   noahmp%water%param%SoilWatConductivitySat,   &
    !$acc   noahmp%water%param%SoilExpCoeffB,   &
    !$acc   noahmp%water%param%SoilMatPotentialSat,   &
    !$acc   noahmp%water%state%IrrigationCntSprinkler,   &
    !$acc   noahmp%water%state%IrrigationCntMicro,   &
    !$acc   noahmp%water%state%IrrigationCntFlood,   &
    !$acc   noahmp%water%state%CanopyTotalWater,   &
    !$acc   noahmp%water%state%CanopyWetFrac,   &
    !$acc   noahmp%water%state%SnowfallDensity,   &
    !$acc   noahmp%water%state%CanopyLiqWater,   &
    !$acc   noahmp%water%state%CanopyIce,   &
    !$acc   noahmp%water%state%CanopyIceMax,   &
    !$acc   noahmp%water%state%CanopyLiqWaterMax,   &
    !$acc   noahmp%water%state%SnowDepth,   &
    !$acc   noahmp%water%state%SnowWaterEquiv,   &
    !$acc   noahmp%water%state%SnowWaterEquivPrev,   &
    !$acc   noahmp%water%state%PondSfcThinSnwMelt,   &
    !$acc   noahmp%water%state%PondSfcThinSnwComb,   &
    !$acc   noahmp%water%state%PondSfcThinSnwTrans,   &
    !$acc   noahmp%water%state%IrrigationFracFlood,   &
    !$acc   noahmp%water%state%IrrigationAmtFlood,   &
    !$acc   noahmp%water%state%IrrigationFracMicro,   &
    !$acc   noahmp%water%state%IrrigationAmtMicro,   &
    !$acc   noahmp%water%state%IrrigationFracSprinkler,   &
    !$acc   noahmp%water%state%IrrigationAmtSprinkler,   &
    !$acc   noahmp%water%state%WaterTableDepth,   &
    !$acc   noahmp%water%state%SoilIceMax    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%water%state%SoilLiqWaterMin,   &
    !$acc   noahmp%water%state%SoilSaturateFrac,   &
    !$acc   noahmp%water%state%SoilImpervFracMax,   &
    !$acc   noahmp%water%state%SoilMoistureToWT,   &
    !$acc   noahmp%water%state%RechargeGwDeepWT,   &
    !$acc   noahmp%water%state%RechargeGwShallowWT,   &
    !$acc   noahmp%water%state%SoilSaturationExcess,   &
    !$acc   noahmp%water%state%WaterTableHydro,   &
    !$acc   noahmp%water%state%TileDrainFrac,   &
    !$acc   noahmp%water%state%WaterStorageAquifer,   &
    !$acc   noahmp%water%state%WaterStorageSoilAqf,   &
    !$acc   noahmp%water%state%WaterStorageLake,   &
    !$acc   noahmp%water%state%WaterStorageWetland,   &
    !$acc   noahmp%water%state%WaterHeadSfc,   &
    !$acc   noahmp%water%state%IrrigationFracGrid,   &
    !$acc   noahmp%water%state%PrecipAreaFrac,   &
    !$acc   noahmp%water%state%SnowCoverFrac,   &
    !$acc   noahmp%water%state%SoilTranspFacAcc,   &
    !$acc   noahmp%water%state%FrozenPrecipFrac,   &
    !$acc   noahmp%water%state%SoilWaterRootZone,   &
    !$acc   noahmp%water%state%SoilWaterStress,   &
    !$acc   noahmp%water%state%WaterStorageTotBeg,   &
    !$acc   noahmp%water%state%WaterBalanceError,   &
    !$acc   noahmp%water%state%WaterStorageTotEnd,   &
    !$acc   noahmp%water%state%SnowRadiusFresh,   &
    !$acc   noahmp%water%flux%RainfallRefHeight,   &
    !$acc   noahmp%water%flux%SnowfallRefHeight,   &
    !$acc   noahmp%water%flux%PrecipTotRefHeight,   &
    !$acc   noahmp%water%flux%PrecipConvTotRefHeight,   &
    !$acc   noahmp%water%flux%PrecipLargeSclRefHeight    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%water%flux%EvapCanopyNet,   &
    !$acc   noahmp%water%flux%Transpiration,   &
    !$acc   noahmp%water%flux%EvapCanopyLiq,   &
    !$acc   noahmp%water%flux%DewCanopyLiq,   &
    !$acc   noahmp%water%flux%FrostCanopyIce,   &
    !$acc   noahmp%water%flux%SublimCanopyIce,   &
    !$acc   noahmp%water%flux%MeltCanopyIce,   &
    !$acc   noahmp%water%flux%FreezeCanopyLiq,   &
    !$acc   noahmp%water%flux%SnowfallGround,   &
    !$acc   noahmp%water%flux%SnowDepthIncr,   &
    !$acc   noahmp%water%flux%FrostSnowSfcIce,   &
    !$acc   noahmp%water%flux%SublimSnowSfcIce,   &
    !$acc   noahmp%water%flux%RainfallGround,   &
    !$acc   noahmp%water%flux%SnowBotOutflow,   &
    !$acc   noahmp%water%flux%GlacierExcessFlow,   &
    !$acc   noahmp%water%flux%IrrigationRateFlood,   &
    !$acc   noahmp%water%flux%IrrigationRateMicro,   &
    !$acc   noahmp%water%flux%IrrigationRateSprinkler,   &
    !$acc   noahmp%water%flux%IrriEvapLossSprinkler,   &
    !$acc   noahmp%water%flux%SoilSfcInflow,   &
    !$acc   noahmp%water%flux%RunoffSurface,   &
    !$acc   noahmp%water%flux%RunoffSubsurface,   &
    !$acc   noahmp%water%flux%InfilRateSfc,   &
    !$acc   noahmp%water%flux%EvapSoilSfcLiq,   &
    !$acc   noahmp%water%flux%DrainSoilBot,   &
    !$acc   noahmp%water%flux%TileDrain,   &
    !$acc   noahmp%water%flux%RechargeGw,   &
    !$acc   noahmp%water%flux%DischargeGw,   &
    !$acc   noahmp%water%flux%VaporizeGrd,   &
    !$acc   noahmp%water%flux%CondenseVapGrd    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%water%flux%DewSoilSfcLiq,   &
    !$acc   noahmp%water%flux%EvapIrriSprinkler,   &
    !$acc   noahmp%water%flux%InterceptCanopyRain,   &
    !$acc   noahmp%water%flux%DripCanopyRain,   &
    !$acc   noahmp%water%flux%ThroughfallRain,   &
    !$acc   noahmp%water%flux%InterceptCanopySnow,   &
    !$acc   noahmp%water%flux%DripCanopySnow,   &
    !$acc   noahmp%water%flux%ThroughfallSnow,   &
    !$acc   noahmp%water%flux%EvapGroundNet,   &
    !$acc   noahmp%water%flux%MeltGroundSnow,   &
    !$acc   noahmp%water%flux%WaterToAtmosTotal,   &
    !$acc   noahmp%water%flux%EvapSoilSfcLiqAcc,   &
    !$acc   noahmp%water%flux%SoilSfcInflowAcc,   &
    !$acc   noahmp%water%flux%SfcWaterTotChgAcc,   &
    !$acc   noahmp%water%flux%PrecipTotAcc,   &
    !$acc   noahmp%water%flux%EvapCanopyNetAcc,   &
    !$acc   noahmp%water%flux%TranspirationAcc,   &
    !$acc   noahmp%water%flux%EvapGroundNetAcc,   &
    !$acc   noahmp%water%flux%GlacierExcessFlowAcc,   &
    !$acc   noahmp%water%flux%EvapSoilSfcLiqMean,   &
    !$acc   noahmp%water%flux%SoilSfcInflowMean,   &
    !$acc   noahmp%water%param%SnowCoverFac,   &
    !$acc   noahmp%water%param%DrainSoilLayerInd,   &
    !$acc   noahmp%water%param%TileDrainTubeDepth,   &
    !$acc   noahmp%water%param%NumSoilLayerRoot,   &
    !$acc   noahmp%water%param%IrriStopDayBfHarvest,   &
    !$acc   noahmp%water%param%CanopyLiqHoldCap,   &
    !$acc   noahmp%water%param%SnowCompactBurdenFac,   &
    !$acc   noahmp%water%param%SnowCompactAgingFac1,   &
    !$acc   noahmp%water%param%SnowCompactAgingFac2    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%water%param%SnowCompactAgingFac3,   &
    !$acc   noahmp%water%param%SnowCompactAgingMax,   &
    !$acc   noahmp%water%param%SnowViscosityCoeff,   &
    !$acc   noahmp%water%param%SnowCompactmAR24,   &
    !$acc   noahmp%water%param%SnowCompactbAR24,   &
    !$acc   noahmp%water%param%SnowCompactP1AR24,   &
    !$acc   noahmp%water%param%SnowCompactP2AR24,   &
    !$acc   noahmp%water%param%SnowCompactP3AR24,   &
    !$acc   noahmp%water%param%BurdenFacUpAR24,   &
    !$acc   noahmp%water%param%SnowCoverM1AR25,   &
    !$acc   noahmp%water%param%SnowCoverM2AR25,   &
    !$acc   noahmp%water%param%SnowCoverFac1AR25,   &
    !$acc   noahmp%water%param%SnowCoverFac2AR25,   &
    !$acc   noahmp%water%param%SnowLiqFracMax,   &
    !$acc   noahmp%water%param%SnowLiqHoldCap,   &
    !$acc   noahmp%water%param%SnowLiqReleaseFac,   &
    !$acc   noahmp%water%param%IrriFloodRateFac,   &
    !$acc   noahmp%water%param%IrriMicroRate,   &
    !$acc   noahmp%water%param%SoilInfilMaxCoeff,   &
    !$acc   noahmp%water%param%SoilImpervFracCoeff,   &
    !$acc   noahmp%water%param%InfilFacVic,   &
    !$acc   noahmp%water%param%TensionWatDistrInfl,   &
    !$acc   noahmp%water%param%TensionWatDistrShp,   &
    !$acc   noahmp%water%param%FreeWatDistrShp,   &
    !$acc   noahmp%water%param%InfilHeteroDynVic,   &
    !$acc   noahmp%water%param%InfilCapillaryDynVic,   &
    !$acc   noahmp%water%param%InfilFacDynVic,   &
    !$acc   noahmp%water%param%SoilDrainSlope,   &
    !$acc   noahmp%water%param%TileDrainCoeffSp,   &
    !$acc   noahmp%water%param%DrainFacSoilWat    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%water%param%TileDrainCoeff,   &
    !$acc   noahmp%water%param%DrainDepthToImperv,   &
    !$acc   noahmp%water%param%LateralWatCondFac,   &
    !$acc   noahmp%water%param%TileDrainDepth,   &
    !$acc   noahmp%water%param%DrainTubeDist,   &
    !$acc   noahmp%water%param%DrainTubeRadius,   &
    !$acc   noahmp%water%param%DrainWatDepToImperv,   &
    !$acc   noahmp%water%param%RunoffDecayFac,   &
    !$acc   noahmp%water%param%BaseflowCoeff,   &
    !$acc   noahmp%water%param%GridTopoIndex,   &
    !$acc   noahmp%water%param%SoilSfcSatFracMax,   &
    !$acc   noahmp%water%param%SpecYieldGw,   &
    !$acc   noahmp%water%param%MicroPoreContent,   &
    !$acc   noahmp%water%param%WaterStorageLakeMax,   &
    !$acc   noahmp%water%param%SnoWatEqvMaxGlacier,   &
    !$acc   noahmp%water%param%SoilConductivityRef,   &
    !$acc   noahmp%water%param%SoilInfilFacRef,   &
    !$acc   noahmp%water%param%GroundFrzCoeff,   &
    !$acc   noahmp%water%param%IrriTriggerLaiMin,   &
    !$acc   noahmp%water%param%SoilWatDeficitAllow,   &
    !$acc   noahmp%water%param%IrriFloodLossFrac,   &
    !$acc   noahmp%water%param%IrriSprinklerRate,   &
    !$acc   noahmp%water%param%IrriFracThreshold,   &
    !$acc   noahmp%water%param%IrriStopPrecipThr,   &
    !$acc   noahmp%water%param%SnowfallDensityMax,   &
    !$acc   noahmp%water%param%SnowMassFullCoverOld,   &
    !$acc   noahmp%water%param%SoilMatPotentialWilt,   &
    !$acc   noahmp%water%param%SnowMeltFac,   &
    !$acc   noahmp%water%param%WetlandCapMax,   &
    !$acc   noahmp%water%param%SnowRadiusMin    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    ! Batched async device create
    !$acc enter data create(              &
    !$acc   noahmp%water%param%FreshSnowRadiusMax,   &
    !$acc   noahmp%water%param%SnowRadiusRefrz,   &
    !$acc   noahmp%water%param%ScavEffMeltScale,   &
    !$acc   noahmp%water%param%ScavEffMeltBCphi,   &
    !$acc   noahmp%water%param%ScavEffMeltBCpho,   &
    !$acc   noahmp%water%param%ScavEffMeltOCphi,   &
    !$acc   noahmp%water%param%ScavEffMeltOCpho,   &
    !$acc   noahmp%water%param%ScavEffMeltDust1,   &
    !$acc   noahmp%water%param%ScavEffMeltDust2,   &
    !$acc   noahmp%water%param%ScavEffMeltDust3,   &
    !$acc   noahmp%water%param%ScavEffMeltDust4,   &
    !$acc   noahmp%water%param%ScavEffMeltDust5,   &
    !$acc   noahmp%water%param%SnowRadiusMax,   &
    !$acc   noahmp%water%param%SnowWetAgeC1Brun89,   &
    !$acc   noahmp%water%param%SnowWetAgeC2Brun89,   &
    !$acc   noahmp%water%param%SnowAgeScaleFac    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
    !$acc enter data create(              &
    !$acc   noahmp%water%state%SnowRadius,   &
    !$acc   noahmp%water%state%MassBChydropho,   &
    !$acc   noahmp%water%state%MassBChydrophi,   &
    !$acc   noahmp%water%state%MassOChydropho,   &
    !$acc   noahmp%water%state%MassOChydrophi,   &
    !$acc   noahmp%water%state%MassDust1,   &
    !$acc   noahmp%water%state%MassDust2,   &
    !$acc   noahmp%water%state%MassDust3,   &
    !$acc   noahmp%water%state%MassDust4,   &
    !$acc   noahmp%water%state%MassDust5,   &
    !$acc   noahmp%water%state%MassConcBChydropho,   &
    !$acc   noahmp%water%state%MassConcBChydrophi,   &
    !$acc   noahmp%water%state%MassConcOChydropho,   &
    !$acc   noahmp%water%state%MassConcOChydrophi,   &
    !$acc   noahmp%water%state%MassConcDust1,   &
    !$acc   noahmp%water%state%MassConcDust2,   &
    !$acc   noahmp%water%state%MassConcDust3,   &
    !$acc   noahmp%water%state%MassConcDust4,   &
    !$acc   noahmp%water%state%MassConcDust5,   &
    !$acc   noahmp%water%flux%SnowFreezeRate,   &
    !$acc   noahmp%water%param%snowage_tau,   &
    !$acc   noahmp%water%param%snowage_kappa,   &
    !$acc   noahmp%water%param%snowage_drdt0    &
    !$acc   ) async(NOAHMP_ACC_QUEUE)
    endif

    !$acc wait(NOAHMP_ACC_QUEUE)

    ! Initialize SNICAR snowage lookup tables (device memory now ready)
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       !$acc parallel loop vector collapse(3)
       do J = 1, noahmp%config%domain%NumTempSnwAgeSnicar
         do LoopInd = 1, noahmp%config%domain%NumTempGradSnwAgeSnicar
           do I = 1, noahmp%config%domain%NumDensitySnwAgeSnicar
             noahmp%water%param%snowage_tau(I,LoopInd,J)  = undefined_real
             noahmp%water%param%snowage_kappa(I,LoopInd,J)= undefined_real
             noahmp%water%param%snowage_drdt0(I,LoopInd,J)= undefined_real
           end do
         end do
       end do
    endif

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        ! Initialize 2D water state scalars
        noahmp%water%state%IrrigationCntSprinkler(I,J)  = undefined_int
        noahmp%water%state%IrrigationCntMicro(I,J)      = undefined_int
        noahmp%water%state%IrrigationCntFlood(I,J)      = undefined_int
        noahmp%water%state%IrrigationFracFlood(I,J)     = undefined_real
        noahmp%water%state%IrrigationAmtFlood(I,J)      = undefined_real
        noahmp%water%state%IrrigationFracMicro(I,J)     = undefined_real
        noahmp%water%state%IrrigationAmtMicro(I,J)      = undefined_real
        noahmp%water%state%IrrigationFracSprinkler(I,J) = undefined_real
        noahmp%water%state%IrrigationAmtSprinkler(I,J)  = undefined_real
        noahmp%water%state%IrrigationFracGrid(I,J)      = undefined_real
        noahmp%water%state%CanopyLiqWater(I,J)          = undefined_real
        noahmp%water%state%CanopyIce(I,J)               = undefined_real
        noahmp%water%state%CanopyTotalWater(I,J)        = undefined_real
        noahmp%water%state%CanopyWetFrac(I,J)           = undefined_real
        noahmp%water%state%CanopyIceMax(I,J)            = undefined_real
        noahmp%water%state%CanopyLiqWaterMax(I,J)       = undefined_real
        noahmp%water%state%SnowfallDensity(I,J)         = undefined_real
        noahmp%water%state%SnowDepth(I,J)               = undefined_real
        noahmp%water%state%SnowWaterEquiv(I,J)          = undefined_real
        noahmp%water%state%SnowWaterEquivPrev(I,J)      = undefined_real
        noahmp%water%state%SnowCoverFrac(I,J)           = undefined_real
        noahmp%water%state%PondSfcThinSnwMelt(I,J)      = undefined_real
        noahmp%water%state%PondSfcThinSnwComb(I,J)      = undefined_real
        noahmp%water%state%PondSfcThinSnwTrans(I,J)     = undefined_real
        noahmp%water%state%SoilIceMax(I,J)              = undefined_real
        noahmp%water%state%SoilLiqWaterMin(I,J)         = undefined_real
        noahmp%water%state%SoilSaturateFrac(I,J)        = undefined_real
        noahmp%water%state%SoilImpervFracMax(I,J)       = undefined_real
        noahmp%water%state%SoilMoistureToWT(I,J)        = undefined_real
        noahmp%water%state%SoilTranspFacAcc(I,J)        = undefined_real
        noahmp%water%state%SoilWaterRootZone(I,J)       = undefined_real
        noahmp%water%state%SoilWaterStress(I,J)         = undefined_real
        noahmp%water%state%SoilSaturationExcess(I,J)    = undefined_real
        noahmp%water%state%RechargeGwDeepWT(I,J)        = undefined_real
        noahmp%water%state%RechargeGwShallowWT(I,J)     = undefined_real
        noahmp%water%state%WaterTableHydro(I,J)         = undefined_real
        noahmp%water%state%WaterTableDepth(I,J)         = undefined_real
        noahmp%water%state%WaterStorageAquifer(I,J)     = undefined_real
        noahmp%water%state%WaterStorageSoilAqf(I,J)     = undefined_real
        noahmp%water%state%WaterStorageLake(I,J)        = undefined_real
        noahmp%water%state%WaterStorageTotBeg(I,J)      = undefined_real
        noahmp%water%state%WaterBalanceError(I,J)       = undefined_real
        noahmp%water%state%WaterStorageTotEnd(I,J)      = undefined_real
        noahmp%water%state%WaterHeadSfc(I,J)            = undefined_real
        noahmp%water%state%PrecipAreaFrac(I,J)          = undefined_real
        noahmp%water%state%TileDrainFrac(I,J)           = undefined_real
        noahmp%water%state%FrozenPrecipFrac(I,J)        = undefined_real
        noahmp%water%state%WaterStorageWetland(I,J)     = undefined_real

        ! Initialize 2D water flux scalars
        noahmp%water%flux%PrecipTotRefHeight(I,J)       = undefined_real
        noahmp%water%flux%RainfallRefHeight(I,J)        = undefined_real
        noahmp%water%flux%SnowfallRefHeight(I,J)        = undefined_real
        noahmp%water%flux%PrecipConvTotRefHeight(I,J)   = undefined_real
        noahmp%water%flux%PrecipLargeSclRefHeight(I,J)  = undefined_real
        noahmp%water%flux%EvapCanopyNet(I,J)            = undefined_real
        noahmp%water%flux%Transpiration(I,J)            = undefined_real
        noahmp%water%flux%EvapCanopyLiq(I,J)            = undefined_real
        noahmp%water%flux%DewCanopyLiq(I,J)             = undefined_real
        noahmp%water%flux%FrostCanopyIce(I,J)           = undefined_real
        noahmp%water%flux%SublimCanopyIce(I,J)          = undefined_real
        noahmp%water%flux%MeltCanopyIce(I,J)            = undefined_real
        noahmp%water%flux%FreezeCanopyLiq(I,J)          = undefined_real
        noahmp%water%flux%SnowfallGround(I,J)           = undefined_real
        noahmp%water%flux%SnowDepthIncr(I,J)            = undefined_real
        noahmp%water%flux%FrostSnowSfcIce(I,J)          = undefined_real
        noahmp%water%flux%SublimSnowSfcIce(I,J)         = undefined_real
        noahmp%water%flux%RainfallGround(I,J)           = undefined_real
        noahmp%water%flux%SnowBotOutflow(I,J)           = undefined_real
        noahmp%water%flux%GlacierExcessFlow(I,J)        = undefined_real
        noahmp%water%flux%SoilSfcInflow(I,J)            = undefined_real
        noahmp%water%flux%RunoffSurface(I,J)            = undefined_real
        noahmp%water%flux%RunoffSubsurface(I,J)         = undefined_real
        noahmp%water%flux%InfilRateSfc(I,J)             = undefined_real
        noahmp%water%flux%EvapSoilSfcLiq(I,J)           = undefined_real
        noahmp%water%flux%DrainSoilBot(I,J)             = undefined_real
        noahmp%water%flux%RechargeGw(I,J)               = undefined_real
        noahmp%water%flux%DischargeGw(I,J)              = undefined_real
        noahmp%water%flux%VaporizeGrd(I,J)              = undefined_real
        noahmp%water%flux%CondenseVapGrd(I,J)           = undefined_real
        noahmp%water%flux%DewSoilSfcLiq(I,J)            = undefined_real
        noahmp%water%flux%InterceptCanopyRain(I,J)      = undefined_real
        noahmp%water%flux%DripCanopyRain(I,J)           = undefined_real
        noahmp%water%flux%ThroughfallRain(I,J)          = undefined_real
        noahmp%water%flux%InterceptCanopySnow(I,J)      = undefined_real
        noahmp%water%flux%DripCanopySnow(I,J)           = undefined_real
        noahmp%water%flux%ThroughfallSnow(I,J)          = undefined_real
        noahmp%water%flux%EvapGroundNet(I,J)            = undefined_real
        noahmp%water%flux%MeltGroundSnow(I,J)           = undefined_real
        noahmp%water%flux%WaterToAtmosTotal(I,J)        = undefined_real
        noahmp%water%flux%EvapSoilSfcLiqAcc(I,J)        = undefined_real
        noahmp%water%flux%SoilSfcInflowAcc(I,J)         = undefined_real
        noahmp%water%flux%GlacierExcessFlowAcc(I,J)     = undefined_real
        noahmp%water%flux%SfcWaterTotChgAcc(I,J)        = undefined_real
        noahmp%water%flux%PrecipTotAcc(I,J)             = undefined_real
        noahmp%water%flux%EvapCanopyNetAcc(I,J)         = undefined_real
        noahmp%water%flux%TranspirationAcc(I,J)         = undefined_real
        noahmp%water%flux%EvapGroundNetAcc(I,J)         = undefined_real
        noahmp%water%flux%EvapSoilSfcLiqMean(I,J)       = undefined_real
        noahmp%water%flux%SoilSfcInflowMean(I,J)        = undefined_real
        noahmp%water%flux%IrrigationRateFlood(I,J)      = 0.0
        noahmp%water%flux%IrrigationRateMicro(I,J)      = 0.0
        noahmp%water%flux%IrrigationRateSprinkler(I,J)  = 0.0
        noahmp%water%flux%IrriEvapLossSprinkler(I,J)    = 0.0
        noahmp%water%flux%EvapIrriSprinkler(I,J)        = 0.0
        noahmp%water%flux%TileDrain(I,J)                = 0.0

        ! Initialize 2D water parameter scalars
        noahmp%water%param%DrainSoilLayerInd(I,J)       = undefined_int
        noahmp%water%param%TileDrainTubeDepth(I,J)      = undefined_int
        noahmp%water%param%NumSoilLayerRoot(I,J)        = undefined_int
        noahmp%water%param%IrriStopDayBfHarvest(I,J)    = undefined_int
        noahmp%water%param%CanopyLiqHoldCap(I,J)        = undefined_real
        noahmp%water%param%SnowCompactBurdenFac(I,J)    = undefined_real
        noahmp%water%param%SnowCompactAgingFac1(I,J)    = undefined_real
        noahmp%water%param%SnowCompactAgingFac2(I,J)    = undefined_real
        noahmp%water%param%SnowCompactAgingFac3(I,J)    = undefined_real
        noahmp%water%param%SnowCompactAgingMax(I,J)     = undefined_real
        noahmp%water%param%SnowViscosityCoeff(I,J)      = undefined_real
        noahmp%water%param%SnowCompactmAR24(I,J)        = undefined_real
        noahmp%water%param%SnowCompactbAR24(I,J)        = undefined_real
        noahmp%water%param%SnowCompactP1AR24(I,J)       = undefined_real
        noahmp%water%param%SnowCompactP2AR24(I,J)       = undefined_real
        noahmp%water%param%SnowCompactP3AR24(I,J)       = undefined_real
        noahmp%water%param%SnowCoverM1AR25(I,J)         = undefined_real
        noahmp%water%param%SnowCoverM2AR25(I,J)         = undefined_real
        noahmp%water%param%SnowCoverFac1AR25(I,J)       = undefined_real
        noahmp%water%param%SnowCoverFac2AR25(I,J)       = undefined_real
        noahmp%water%param%BurdenFacUpAR24(I,J)         = undefined_real
        noahmp%water%param%SnowLiqFracMax(I,J)          = undefined_real
        noahmp%water%param%SnowLiqHoldCap(I,J)          = undefined_real
        noahmp%water%param%SnowLiqReleaseFac(I,J)       = undefined_real
        noahmp%water%param%IrriFloodRateFac(I,J)        = undefined_real
        noahmp%water%param%IrriMicroRate(I,J)           = undefined_real
        noahmp%water%param%SoilInfilMaxCoeff(I,J)       = undefined_real
        noahmp%water%param%SoilImpervFracCoeff(I,J)     = undefined_real
        noahmp%water%param%InfilFacVic(I,J)             = undefined_real
        noahmp%water%param%TensionWatDistrInfl(I,J)     = undefined_real
        noahmp%water%param%TensionWatDistrShp(I,J)      = undefined_real
        noahmp%water%param%FreeWatDistrShp(I,J)         = undefined_real
        noahmp%water%param%InfilHeteroDynVic(I,J)       = undefined_real
        noahmp%water%param%InfilCapillaryDynVic(I,J)    = undefined_real
        noahmp%water%param%InfilFacDynVic(I,J)          = undefined_real
        noahmp%water%param%SoilDrainSlope(I,J)          = undefined_real
        noahmp%water%param%TileDrainCoeffSp(I,J)        = undefined_real
        noahmp%water%param%DrainFacSoilWat(I,J)         = undefined_real
        noahmp%water%param%TileDrainCoeff(I,J)          = undefined_real
        noahmp%water%param%DrainDepthToImperv(I,J)      = undefined_real
        noahmp%water%param%LateralWatCondFac(I,J)       = undefined_real
        noahmp%water%param%TileDrainDepth(I,J)          = undefined_real
        noahmp%water%param%DrainTubeDist(I,J)           = undefined_real
        noahmp%water%param%DrainTubeRadius(I,J)         = undefined_real
        noahmp%water%param%DrainWatDepToImperv(I,J)     = undefined_real
        noahmp%water%param%RunoffDecayFac(I,J)          = undefined_real
        noahmp%water%param%BaseflowCoeff(I,J)           = undefined_real
        noahmp%water%param%GridTopoIndex(I,J)           = undefined_real
        noahmp%water%param%SoilSfcSatFracMax(I,J)       = undefined_real
        noahmp%water%param%SpecYieldGw(I,J)             = undefined_real
        noahmp%water%param%MicroPoreContent(I,J)        = undefined_real
        noahmp%water%param%WaterStorageLakeMax(I,J)     = undefined_real
        noahmp%water%param%SnoWatEqvMaxGlacier(I,J)     = undefined_real
        noahmp%water%param%SoilConductivityRef(I,J)     = undefined_real
        noahmp%water%param%SoilInfilFacRef(I,J)         = undefined_real
        noahmp%water%param%GroundFrzCoeff(I,J)          = undefined_real
        noahmp%water%param%IrriTriggerLaiMin(I,J)       = undefined_real
        noahmp%water%param%SoilWatDeficitAllow(I,J)     = undefined_real
        noahmp%water%param%IrriFloodLossFrac(I,J)       = undefined_real
        noahmp%water%param%IrriSprinklerRate(I,J)       = undefined_real
        noahmp%water%param%IrriFracThreshold(I,J)       = undefined_real
        noahmp%water%param%IrriStopPrecipThr(I,J)       = undefined_real
        noahmp%water%param%SnowfallDensityMax(I,J)      = undefined_real
        noahmp%water%param%SnowMassFullCoverOld(I,J)    = undefined_real
        noahmp%water%param%SoilMatPotentialWilt(I,J)    = undefined_real
        noahmp%water%param%SnowMeltFac(I,J)             = undefined_real
        noahmp%water%param%SnowCoverFac(I,J)            = undefined_real
        noahmp%water%param%WetlandCapMax(I,J)           = undefined_real

        ! Initialize SNICAR 2D scalars
        if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
           noahmp%water%state%SnowRadiusFresh(I,J)      = undefined_real
           noahmp%water%param%SnowRadiusMin(I,J)        = undefined_real
           noahmp%water%param%FreshSnowRadiusMax(I,J)   = undefined_real
           noahmp%water%param%SnowRadiusRefrz(I,J)      = undefined_real
           noahmp%water%param%ScavEffMeltScale(I,J)     = undefined_real
           noahmp%water%param%ScavEffMeltBCphi(I,J)     = undefined_real
           noahmp%water%param%ScavEffMeltBCpho(I,J)     = undefined_real
           noahmp%water%param%ScavEffMeltOCphi(I,J)     = undefined_real
           noahmp%water%param%ScavEffMeltOCpho(I,J)     = undefined_real
           noahmp%water%param%ScavEffMeltDust1(I,J)     = undefined_real
           noahmp%water%param%ScavEffMeltDust2(I,J)     = undefined_real
           noahmp%water%param%ScavEffMeltDust3(I,J)     = undefined_real
           noahmp%water%param%ScavEffMeltDust4(I,J)     = undefined_real
           noahmp%water%param%ScavEffMeltDust5(I,J)     = undefined_real
           noahmp%water%param%SnowRadiusMax(I,J)        = undefined_real
           noahmp%water%param%SnowWetAgeC1Brun89(I,J)   = undefined_real
           noahmp%water%param%SnowWetAgeC2Brun89(I,J)   = undefined_real
           noahmp%water%param%SnowAgeScaleFac(I,J)      = undefined_real
        endif

        ! Initialize 3D water state arrays
        !$acc loop seq
        do LoopInd = -noahmp%config%domain%NumSnowLayerMax+1, noahmp%config%domain%NumSoilLayer
           noahmp%water%state%IndexPhaseChange(I,LoopInd,J)   = undefined_int
           noahmp%water%state%SoilSupercoolWater(I,LoopInd,J) = undefined_real
        enddo

        !$acc loop seq
        do LoopInd = -noahmp%config%domain%NumSnowLayerMax+1, 0
           noahmp%water%state%SnowIce(I,LoopInd,J)         = undefined_real
           noahmp%water%state%SnowLiqWater(I,LoopInd,J)    = undefined_real
           noahmp%water%state%SnowIceVol(I,LoopInd,J)      = undefined_real
           noahmp%water%state%SnowLiqWaterVol(I,LoopInd,J) = undefined_real
           noahmp%water%state%SnowIceFracPrev(I,LoopInd,J) = undefined_real
           noahmp%water%state%SnowIceFrac(I,LoopInd,J)     = undefined_real
           noahmp%water%state%SnowEffPorosity(I,LoopInd,J) = undefined_real
           noahmp%water%flux%CompactionSnowAging(I,LoopInd,J)  = undefined_real
           noahmp%water%flux%CompactionSnowBurden(I,LoopInd,J) = undefined_real
           noahmp%water%flux%CompactionSnowMelt(I,LoopInd,J)   = undefined_real
           noahmp%water%flux%CompactionSnowTot(I,LoopInd,J)    = undefined_real
           noahmp%water%flux%OutflowSnowLayer(I,LoopInd,J)     = undefined_real
        enddo

        !$acc loop seq
        do LoopInd = 1, noahmp%config%domain%NumSoilLayer
           noahmp%water%state%SoilLiqWater(I,LoopInd,J)       = undefined_real
           noahmp%water%state%SoilIce(I,LoopInd,J)            = undefined_real
           noahmp%water%state%SoilMoisture(I,LoopInd,J)       = undefined_real
           noahmp%water%state%SoilImpervFrac(I,LoopInd,J)     = undefined_real
           noahmp%water%state%SoilWatConductivity(I,LoopInd,J)= undefined_real
           noahmp%water%state%SoilWatDiffusivity(I,LoopInd,J) = undefined_real
           noahmp%water%state%SoilEffPorosity(I,LoopInd,J)    = undefined_real
           noahmp%water%state%SoilIceFrac(I,LoopInd,J)        = undefined_real
           noahmp%water%state%SoilMoistureEqui(I,LoopInd,J)   = undefined_real
           noahmp%water%state%SoilTranspFac(I,LoopInd,J)      = undefined_real
           noahmp%water%state%SoilMatPotential(I,LoopInd,J)   = undefined_real
           noahmp%water%flux%TranspWatLossSoil(I,LoopInd,J)     = undefined_real
           noahmp%water%flux%TranspWatLossSoilAcc(I,LoopInd,J)  = undefined_real
           noahmp%water%flux%TranspWatLossSoilMean(I,LoopInd,J) = undefined_real
           noahmp%water%param%SoilMoistureSat(I,LoopInd,J)       = undefined_real
           noahmp%water%param%SoilMoistureWilt(I,LoopInd,J)      = undefined_real
           noahmp%water%param%SoilMoistureFieldCap(I,LoopInd,J)  = undefined_real
           noahmp%water%param%SoilMoistureDry(I,LoopInd,J)       = undefined_real
           noahmp%water%param%SoilWatDiffusivitySat(I,LoopInd,J) = undefined_real
           noahmp%water%param%SoilWatConductivitySat(I,LoopInd,J)= undefined_real
           noahmp%water%param%SoilExpCoeffB(I,LoopInd,J)         = undefined_real
           noahmp%water%param%SoilMatPotentialSat(I,LoopInd,J)   = undefined_real
        enddo

        ! Initialize SNICAR 3D state arrays
        if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
           !$acc loop seq
           do LoopInd = -noahmp%config%domain%NumSnowLayerMax+1, 0
              noahmp%water%state%SnowRadius(I,LoopInd,J)         = undefined_real
              noahmp%water%state%MassBChydropho(I,LoopInd,J)     = undefined_real
              noahmp%water%state%MassBChydrophi(I,LoopInd,J)     = undefined_real
              noahmp%water%state%MassOChydropho(I,LoopInd,J)     = undefined_real
              noahmp%water%state%MassOChydrophi(I,LoopInd,J)     = undefined_real
              noahmp%water%state%MassDust1(I,LoopInd,J)          = undefined_real
              noahmp%water%state%MassDust2(I,LoopInd,J)          = undefined_real
              noahmp%water%state%MassDust3(I,LoopInd,J)          = undefined_real
              noahmp%water%state%MassDust4(I,LoopInd,J)          = undefined_real
              noahmp%water%state%MassDust5(I,LoopInd,J)          = undefined_real
              noahmp%water%state%MassConcBChydropho(I,LoopInd,J) = undefined_real
              noahmp%water%state%MassConcBChydrophi(I,LoopInd,J) = undefined_real
              noahmp%water%state%MassConcOChydropho(I,LoopInd,J) = undefined_real
              noahmp%water%state%MassConcOChydrophi(I,LoopInd,J) = undefined_real
              noahmp%water%state%MassConcDust1(I,LoopInd,J)      = undefined_real
              noahmp%water%state%MassConcDust2(I,LoopInd,J)      = undefined_real
              noahmp%water%state%MassConcDust3(I,LoopInd,J)      = undefined_real
              noahmp%water%state%MassConcDust4(I,LoopInd,J)      = undefined_real
              noahmp%water%state%MassConcDust5(I,LoopInd,J)      = undefined_real
              noahmp%water%flux%SnowFreezeRate(I,LoopInd,J)      = undefined_real
           enddo
        endif

      end do
    end do
    !$acc end parallel loop

  end subroutine WaterVarInitDefault


  subroutine WaterVarExitDevice(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    !$acc exit data delete(               &
    !$acc   noahmp%water%state%IndexPhaseChange,   &
    !$acc   noahmp%water%state%SoilSupercoolWater,   &
    !$acc   noahmp%water%state%SnowIce,   &
    !$acc   noahmp%water%state%SnowLiqWater,   &
    !$acc   noahmp%water%state%SnowIceVol,   &
    !$acc   noahmp%water%state%SnowLiqWaterVol,   &
    !$acc   noahmp%water%state%SnowIceFracPrev,   &
    !$acc   noahmp%water%state%SnowIceFrac,   &
    !$acc   noahmp%water%state%SnowEffPorosity,   &
    !$acc   noahmp%water%state%SoilLiqWater,   &
    !$acc   noahmp%water%state%SoilIce,   &
    !$acc   noahmp%water%state%SoilMoisture,   &
    !$acc   noahmp%water%state%SoilImpervFrac,   &
    !$acc   noahmp%water%state%SoilWatConductivity,   &
    !$acc   noahmp%water%state%SoilWatDiffusivity,   &
    !$acc   noahmp%water%state%SoilEffPorosity,   &
    !$acc   noahmp%water%state%SoilIceFrac,   &
    !$acc   noahmp%water%state%SoilMoistureEqui,   &
    !$acc   noahmp%water%state%SoilTranspFac,   &
    !$acc   noahmp%water%state%SoilMatPotential,   &
    !$acc   noahmp%water%flux%CompactionSnowAging,   &
    !$acc   noahmp%water%flux%CompactionSnowBurden,   &
    !$acc   noahmp%water%flux%CompactionSnowMelt,   &
    !$acc   noahmp%water%flux%CompactionSnowTot,   &
    !$acc   noahmp%water%flux%TranspWatLossSoil,   &
    !$acc   noahmp%water%flux%TranspWatLossSoilAcc,   &
    !$acc   noahmp%water%flux%TranspWatLossSoilMean,   &
    !$acc   noahmp%water%flux%OutflowSnowLayer,   &
    !$acc   noahmp%water%param%SoilMoistureSat,   &
    !$acc   noahmp%water%param%SoilMoistureWilt    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%water%param%SoilMoistureFieldCap,   &
    !$acc   noahmp%water%param%SoilMoistureDry,   &
    !$acc   noahmp%water%param%SoilWatDiffusivitySat,   &
    !$acc   noahmp%water%param%SoilWatConductivitySat,   &
    !$acc   noahmp%water%param%SoilExpCoeffB,   &
    !$acc   noahmp%water%param%SoilMatPotentialSat,   &
    !$acc   noahmp%water%state%IrrigationCntSprinkler,   &
    !$acc   noahmp%water%state%IrrigationCntMicro,   &
    !$acc   noahmp%water%state%IrrigationCntFlood,   &
    !$acc   noahmp%water%state%CanopyTotalWater,   &
    !$acc   noahmp%water%state%CanopyWetFrac,   &
    !$acc   noahmp%water%state%SnowfallDensity,   &
    !$acc   noahmp%water%state%CanopyLiqWater,   &
    !$acc   noahmp%water%state%CanopyIce,   &
    !$acc   noahmp%water%state%CanopyIceMax,   &
    !$acc   noahmp%water%state%CanopyLiqWaterMax,   &
    !$acc   noahmp%water%state%SnowDepth,   &
    !$acc   noahmp%water%state%SnowWaterEquiv,   &
    !$acc   noahmp%water%state%SnowWaterEquivPrev,   &
    !$acc   noahmp%water%state%PondSfcThinSnwMelt,   &
    !$acc   noahmp%water%state%PondSfcThinSnwComb,   &
    !$acc   noahmp%water%state%PondSfcThinSnwTrans,   &
    !$acc   noahmp%water%state%IrrigationFracFlood,   &
    !$acc   noahmp%water%state%IrrigationAmtFlood,   &
    !$acc   noahmp%water%state%IrrigationFracMicro,   &
    !$acc   noahmp%water%state%IrrigationAmtMicro,   &
    !$acc   noahmp%water%state%IrrigationFracSprinkler,   &
    !$acc   noahmp%water%state%IrrigationAmtSprinkler,   &
    !$acc   noahmp%water%state%WaterTableDepth,   &
    !$acc   noahmp%water%state%SoilIceMax    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%water%state%SoilLiqWaterMin,   &
    !$acc   noahmp%water%state%SoilSaturateFrac,   &
    !$acc   noahmp%water%state%SoilImpervFracMax,   &
    !$acc   noahmp%water%state%SoilMoistureToWT,   &
    !$acc   noahmp%water%state%RechargeGwDeepWT,   &
    !$acc   noahmp%water%state%RechargeGwShallowWT,   &
    !$acc   noahmp%water%state%SoilSaturationExcess,   &
    !$acc   noahmp%water%state%WaterTableHydro,   &
    !$acc   noahmp%water%state%TileDrainFrac,   &
    !$acc   noahmp%water%state%WaterStorageAquifer,   &
    !$acc   noahmp%water%state%WaterStorageSoilAqf,   &
    !$acc   noahmp%water%state%WaterStorageLake,   &
    !$acc   noahmp%water%state%WaterStorageWetland,   &
    !$acc   noahmp%water%state%WaterHeadSfc,   &
    !$acc   noahmp%water%state%IrrigationFracGrid,   &
    !$acc   noahmp%water%state%PrecipAreaFrac,   &
    !$acc   noahmp%water%state%SnowCoverFrac,   &
    !$acc   noahmp%water%state%SoilTranspFacAcc,   &
    !$acc   noahmp%water%state%FrozenPrecipFrac,   &
    !$acc   noahmp%water%state%SoilWaterRootZone,   &
    !$acc   noahmp%water%state%SoilWaterStress,   &
    !$acc   noahmp%water%state%WaterStorageTotBeg,   &
    !$acc   noahmp%water%state%WaterBalanceError,   &
    !$acc   noahmp%water%state%WaterStorageTotEnd,   &
    !$acc   noahmp%water%state%SnowRadiusFresh,   &
    !$acc   noahmp%water%flux%RainfallRefHeight,   &
    !$acc   noahmp%water%flux%SnowfallRefHeight,   &
    !$acc   noahmp%water%flux%PrecipTotRefHeight,   &
    !$acc   noahmp%water%flux%PrecipConvTotRefHeight,   &
    !$acc   noahmp%water%flux%PrecipLargeSclRefHeight    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%water%flux%EvapCanopyNet,   &
    !$acc   noahmp%water%flux%Transpiration,   &
    !$acc   noahmp%water%flux%EvapCanopyLiq,   &
    !$acc   noahmp%water%flux%DewCanopyLiq,   &
    !$acc   noahmp%water%flux%FrostCanopyIce,   &
    !$acc   noahmp%water%flux%SublimCanopyIce,   &
    !$acc   noahmp%water%flux%MeltCanopyIce,   &
    !$acc   noahmp%water%flux%FreezeCanopyLiq,   &
    !$acc   noahmp%water%flux%SnowfallGround,   &
    !$acc   noahmp%water%flux%SnowDepthIncr,   &
    !$acc   noahmp%water%flux%FrostSnowSfcIce,   &
    !$acc   noahmp%water%flux%SublimSnowSfcIce,   &
    !$acc   noahmp%water%flux%RainfallGround,   &
    !$acc   noahmp%water%flux%SnowBotOutflow,   &
    !$acc   noahmp%water%flux%GlacierExcessFlow,   &
    !$acc   noahmp%water%flux%IrrigationRateFlood,   &
    !$acc   noahmp%water%flux%IrrigationRateMicro,   &
    !$acc   noahmp%water%flux%IrrigationRateSprinkler,   &
    !$acc   noahmp%water%flux%IrriEvapLossSprinkler,   &
    !$acc   noahmp%water%flux%SoilSfcInflow,   &
    !$acc   noahmp%water%flux%RunoffSurface,   &
    !$acc   noahmp%water%flux%RunoffSubsurface,   &
    !$acc   noahmp%water%flux%InfilRateSfc,   &
    !$acc   noahmp%water%flux%EvapSoilSfcLiq,   &
    !$acc   noahmp%water%flux%DrainSoilBot,   &
    !$acc   noahmp%water%flux%TileDrain,   &
    !$acc   noahmp%water%flux%RechargeGw,   &
    !$acc   noahmp%water%flux%DischargeGw,   &
    !$acc   noahmp%water%flux%VaporizeGrd,   &
    !$acc   noahmp%water%flux%CondenseVapGrd    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%water%flux%DewSoilSfcLiq,   &
    !$acc   noahmp%water%flux%EvapIrriSprinkler,   &
    !$acc   noahmp%water%flux%InterceptCanopyRain,   &
    !$acc   noahmp%water%flux%DripCanopyRain,   &
    !$acc   noahmp%water%flux%ThroughfallRain,   &
    !$acc   noahmp%water%flux%InterceptCanopySnow,   &
    !$acc   noahmp%water%flux%DripCanopySnow,   &
    !$acc   noahmp%water%flux%ThroughfallSnow,   &
    !$acc   noahmp%water%flux%EvapGroundNet,   &
    !$acc   noahmp%water%flux%MeltGroundSnow,   &
    !$acc   noahmp%water%flux%WaterToAtmosTotal,   &
    !$acc   noahmp%water%flux%EvapSoilSfcLiqAcc,   &
    !$acc   noahmp%water%flux%SoilSfcInflowAcc,   &
    !$acc   noahmp%water%flux%SfcWaterTotChgAcc,   &
    !$acc   noahmp%water%flux%PrecipTotAcc,   &
    !$acc   noahmp%water%flux%EvapCanopyNetAcc,   &
    !$acc   noahmp%water%flux%TranspirationAcc,   &
    !$acc   noahmp%water%flux%EvapGroundNetAcc,   &
    !$acc   noahmp%water%flux%GlacierExcessFlowAcc,   &
    !$acc   noahmp%water%flux%EvapSoilSfcLiqMean,   &
    !$acc   noahmp%water%flux%SoilSfcInflowMean,   &
    !$acc   noahmp%water%param%SnowCoverFac,   &
    !$acc   noahmp%water%param%DrainSoilLayerInd,   &
    !$acc   noahmp%water%param%TileDrainTubeDepth,   &
    !$acc   noahmp%water%param%NumSoilLayerRoot,   &
    !$acc   noahmp%water%param%IrriStopDayBfHarvest,   &
    !$acc   noahmp%water%param%CanopyLiqHoldCap,   &
    !$acc   noahmp%water%param%SnowCompactBurdenFac,   &
    !$acc   noahmp%water%param%SnowCompactAgingFac1,   &
    !$acc   noahmp%water%param%SnowCompactAgingFac2    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%water%param%SnowCompactAgingFac3,   &
    !$acc   noahmp%water%param%SnowCompactAgingMax,   &
    !$acc   noahmp%water%param%SnowViscosityCoeff,   &
    !$acc   noahmp%water%param%SnowCompactmAR24,   &
    !$acc   noahmp%water%param%SnowCompactbAR24,   &
    !$acc   noahmp%water%param%SnowCompactP1AR24,   &
    !$acc   noahmp%water%param%SnowCompactP2AR24,   &
    !$acc   noahmp%water%param%SnowCompactP3AR24,   &
    !$acc   noahmp%water%param%BurdenFacUpAR24,   &
    !$acc   noahmp%water%param%SnowCoverM1AR25,   &
    !$acc   noahmp%water%param%SnowCoverM2AR25,   &
    !$acc   noahmp%water%param%SnowCoverFac1AR25,   &
    !$acc   noahmp%water%param%SnowCoverFac2AR25,   &
    !$acc   noahmp%water%param%SnowLiqFracMax,   &
    !$acc   noahmp%water%param%SnowLiqHoldCap,   &
    !$acc   noahmp%water%param%SnowLiqReleaseFac,   &
    !$acc   noahmp%water%param%IrriFloodRateFac,   &
    !$acc   noahmp%water%param%IrriMicroRate,   &
    !$acc   noahmp%water%param%SoilInfilMaxCoeff,   &
    !$acc   noahmp%water%param%SoilImpervFracCoeff,   &
    !$acc   noahmp%water%param%InfilFacVic,   &
    !$acc   noahmp%water%param%TensionWatDistrInfl,   &
    !$acc   noahmp%water%param%TensionWatDistrShp,   &
    !$acc   noahmp%water%param%FreeWatDistrShp,   &
    !$acc   noahmp%water%param%InfilHeteroDynVic,   &
    !$acc   noahmp%water%param%InfilCapillaryDynVic,   &
    !$acc   noahmp%water%param%InfilFacDynVic,   &
    !$acc   noahmp%water%param%SoilDrainSlope,   &
    !$acc   noahmp%water%param%TileDrainCoeffSp,   &
    !$acc   noahmp%water%param%DrainFacSoilWat    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%water%param%TileDrainCoeff,   &
    !$acc   noahmp%water%param%DrainDepthToImperv,   &
    !$acc   noahmp%water%param%LateralWatCondFac,   &
    !$acc   noahmp%water%param%TileDrainDepth,   &
    !$acc   noahmp%water%param%DrainTubeDist,   &
    !$acc   noahmp%water%param%DrainTubeRadius,   &
    !$acc   noahmp%water%param%DrainWatDepToImperv,   &
    !$acc   noahmp%water%param%RunoffDecayFac,   &
    !$acc   noahmp%water%param%BaseflowCoeff,   &
    !$acc   noahmp%water%param%GridTopoIndex,   &
    !$acc   noahmp%water%param%SoilSfcSatFracMax,   &
    !$acc   noahmp%water%param%SpecYieldGw,   &
    !$acc   noahmp%water%param%MicroPoreContent,   &
    !$acc   noahmp%water%param%WaterStorageLakeMax,   &
    !$acc   noahmp%water%param%SnoWatEqvMaxGlacier,   &
    !$acc   noahmp%water%param%SoilConductivityRef,   &
    !$acc   noahmp%water%param%SoilInfilFacRef,   &
    !$acc   noahmp%water%param%GroundFrzCoeff,   &
    !$acc   noahmp%water%param%IrriTriggerLaiMin,   &
    !$acc   noahmp%water%param%SoilWatDeficitAllow,   &
    !$acc   noahmp%water%param%IrriFloodLossFrac,   &
    !$acc   noahmp%water%param%IrriSprinklerRate,   &
    !$acc   noahmp%water%param%IrriFracThreshold,   &
    !$acc   noahmp%water%param%IrriStopPrecipThr,   &
    !$acc   noahmp%water%param%SnowfallDensityMax,   &
    !$acc   noahmp%water%param%SnowMassFullCoverOld,   &
    !$acc   noahmp%water%param%SoilMatPotentialWilt,   &
    !$acc   noahmp%water%param%SnowMeltFac,   &
    !$acc   noahmp%water%param%WetlandCapMax,   &
    !$acc   noahmp%water%param%SnowRadiusMin    &
    !$acc   )

    !$acc exit data delete(               &
    !$acc   noahmp%water%param%FreshSnowRadiusMax,   &
    !$acc   noahmp%water%param%SnowRadiusRefrz,   &
    !$acc   noahmp%water%param%ScavEffMeltScale,   &
    !$acc   noahmp%water%param%ScavEffMeltBCphi,   &
    !$acc   noahmp%water%param%ScavEffMeltBCpho,   &
    !$acc   noahmp%water%param%ScavEffMeltOCphi,   &
    !$acc   noahmp%water%param%ScavEffMeltOCpho,   &
    !$acc   noahmp%water%param%ScavEffMeltDust1,   &
    !$acc   noahmp%water%param%ScavEffMeltDust2,   &
    !$acc   noahmp%water%param%ScavEffMeltDust3,   &
    !$acc   noahmp%water%param%ScavEffMeltDust4,   &
    !$acc   noahmp%water%param%ScavEffMeltDust5,   &
    !$acc   noahmp%water%param%SnowRadiusMax,   &
    !$acc   noahmp%water%param%SnowWetAgeC1Brun89,   &
    !$acc   noahmp%water%param%SnowWetAgeC2Brun89,   &
    !$acc   noahmp%water%param%SnowAgeScaleFac    &
    !$acc   )

    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
    !$acc exit data delete(               &
    !$acc   noahmp%water%state%SnowRadius,   &
    !$acc   noahmp%water%state%MassBChydropho,   &
    !$acc   noahmp%water%state%MassBChydrophi,   &
    !$acc   noahmp%water%state%MassOChydropho,   &
    !$acc   noahmp%water%state%MassOChydrophi,   &
    !$acc   noahmp%water%state%MassDust1,   &
    !$acc   noahmp%water%state%MassDust2,   &
    !$acc   noahmp%water%state%MassDust3,   &
    !$acc   noahmp%water%state%MassDust4,   &
    !$acc   noahmp%water%state%MassDust5,   &
    !$acc   noahmp%water%state%MassConcBChydropho,   &
    !$acc   noahmp%water%state%MassConcBChydrophi,   &
    !$acc   noahmp%water%state%MassConcOChydropho,   &
    !$acc   noahmp%water%state%MassConcOChydrophi,   &
    !$acc   noahmp%water%state%MassConcDust1,   &
    !$acc   noahmp%water%state%MassConcDust2,   &
    !$acc   noahmp%water%state%MassConcDust3,   &
    !$acc   noahmp%water%state%MassConcDust4,   &
    !$acc   noahmp%water%state%MassConcDust5,   &
    !$acc   noahmp%water%flux%SnowFreezeRate,   &
    !$acc   noahmp%water%param%snowage_tau,   &
    !$acc   noahmp%water%param%snowage_kappa,   &
    !$acc   noahmp%water%param%snowage_drdt0    &
    !$acc   )
    endif

  end subroutine WaterVarExitDevice
end module WaterVarInitMod
