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
       !$acc enter data create(noahmp%water%state%IndexPhaseChange)
    endif
    if ( .not. allocated(noahmp%water%state%SoilSupercoolWater) ) then
       allocate( noahmp%water%state%SoilSupercoolWater(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilSupercoolWater)
    endif
    if ( .not. allocated(noahmp%water%state%SnowIce) ) then
       allocate( noahmp%water%state%SnowIce(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowIce)
    endif
    if ( .not. allocated(noahmp%water%state%SnowLiqWater) ) then
       allocate( noahmp%water%state%SnowLiqWater(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowLiqWater)
    endif
    if ( .not. allocated(noahmp%water%state%SnowIceVol) ) then
       allocate( noahmp%water%state%SnowIceVol(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowIceVol)
    endif
    if ( .not. allocated(noahmp%water%state%SnowLiqWaterVol) ) then
       allocate( noahmp%water%state%SnowLiqWaterVol(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowLiqWaterVol)
    endif
    if ( .not. allocated(noahmp%water%state%SnowIceFracPrev) ) then
       allocate( noahmp%water%state%SnowIceFracPrev(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowIceFracPrev)
    endif
    if ( .not. allocated(noahmp%water%state%SnowIceFrac) ) then
       allocate( noahmp%water%state%SnowIceFrac(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowIceFrac)
    endif
    if ( .not. allocated(noahmp%water%state%SnowEffPorosity) ) then
       allocate( noahmp%water%state%SnowEffPorosity(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowEffPorosity)
    endif
    if ( .not. allocated(noahmp%water%state%SoilLiqWater) ) then
       allocate( noahmp%water%state%SoilLiqWater(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilLiqWater)
    endif
    if ( .not. allocated(noahmp%water%state%SoilIce) ) then
       allocate( noahmp%water%state%SoilIce(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilIce)
    endif
    if ( .not. allocated(noahmp%water%state%SoilMoisture) ) then
       allocate( noahmp%water%state%SoilMoisture(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilMoisture)
    endif
    if ( .not. allocated(noahmp%water%state%SoilImpervFrac) ) then
       allocate( noahmp%water%state%SoilImpervFrac(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilImpervFrac)
    endif
    if ( .not. allocated(noahmp%water%state%SoilWatConductivity) ) then
       allocate( noahmp%water%state%SoilWatConductivity(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilWatConductivity)
    endif
    if ( .not. allocated(noahmp%water%state%SoilWatDiffusivity) ) then
       allocate( noahmp%water%state%SoilWatDiffusivity(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilWatDiffusivity)
    endif
    if ( .not. allocated(noahmp%water%state%SoilEffPorosity) ) then
       allocate( noahmp%water%state%SoilEffPorosity(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilEffPorosity)
    endif
    if ( .not. allocated(noahmp%water%state%SoilIceFrac) ) then
       allocate( noahmp%water%state%SoilIceFrac(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilIceFrac)
    endif
    if ( .not. allocated(noahmp%water%state%SoilMoistureEqui) ) then
       allocate( noahmp%water%state%SoilMoistureEqui(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilMoistureEqui)
    endif
    if ( .not. allocated(noahmp%water%state%SoilTranspFac) ) then
       allocate( noahmp%water%state%SoilTranspFac(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilTranspFac)
    endif
    if ( .not. allocated(noahmp%water%state%SoilMatPotential) ) then
       allocate( noahmp%water%state%SoilMatPotential(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilMatPotential)
    endif

    ! SNICAR state arrays
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( .not. allocated(noahmp%water%state%SnowRadius) ) then
          allocate( noahmp%water%state%SnowRadius(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%SnowRadius)
       endif
       if ( .not. allocated(noahmp%water%state%MassBChydropho) ) then
          allocate( noahmp%water%state%MassBChydropho(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassBChydropho)
       endif
       if ( .not. allocated(noahmp%water%state%MassBChydrophi) ) then
          allocate( noahmp%water%state%MassBChydrophi(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassBChydrophi)
       endif
       if ( .not. allocated(noahmp%water%state%MassOChydropho) ) then
          allocate( noahmp%water%state%MassOChydropho(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassOChydropho)
       endif
       if ( .not. allocated(noahmp%water%state%MassOChydrophi) ) then
          allocate( noahmp%water%state%MassOChydrophi(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassOChydrophi)
       endif
       if ( .not. allocated(noahmp%water%state%MassDust1) ) then
          allocate( noahmp%water%state%MassDust1(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassDust1)
       endif
       if ( .not. allocated(noahmp%water%state%MassDust2) ) then
          allocate( noahmp%water%state%MassDust2(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassDust2)
       endif
       if ( .not. allocated(noahmp%water%state%MassDust3) ) then
          allocate( noahmp%water%state%MassDust3(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassDust3)
       endif
       if ( .not. allocated(noahmp%water%state%MassDust4) ) then
          allocate( noahmp%water%state%MassDust4(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassDust4)
       endif
       if ( .not. allocated(noahmp%water%state%MassDust5) ) then
          allocate( noahmp%water%state%MassDust5(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassDust5)
       endif
       if ( .not. allocated(noahmp%water%state%MassConcBChydropho) ) then
          allocate( noahmp%water%state%MassConcBChydropho(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassConcBChydropho)
       endif
       if ( .not. allocated(noahmp%water%state%MassConcBChydrophi) ) then
          allocate( noahmp%water%state%MassConcBChydrophi(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassConcBChydrophi)
       endif
       if ( .not. allocated(noahmp%water%state%MassConcOChydropho) ) then
          allocate( noahmp%water%state%MassConcOChydropho(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassConcOChydropho)
       endif
       if ( .not. allocated(noahmp%water%state%MassConcOChydrophi) ) then
          allocate( noahmp%water%state%MassConcOChydrophi(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassConcOChydrophi)
       endif
       if ( .not. allocated(noahmp%water%state%MassConcDust1) ) then
          allocate( noahmp%water%state%MassConcDust1(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassConcDust1)
       endif
       if ( .not. allocated(noahmp%water%state%MassConcDust2) ) then
          allocate( noahmp%water%state%MassConcDust2(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassConcDust2)
       endif
       if ( .not. allocated(noahmp%water%state%MassConcDust3) ) then
          allocate( noahmp%water%state%MassConcDust3(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassConcDust3)
       endif
       if ( .not. allocated(noahmp%water%state%MassConcDust4) ) then
          allocate( noahmp%water%state%MassConcDust4(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassConcDust4)
       endif
       if ( .not. allocated(noahmp%water%state%MassConcDust5) ) then
          allocate( noahmp%water%state%MassConcDust5(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%state%MassConcDust5)
       endif
    endif

    ! Allocate 3D water flux arrays and transfer to GPU
    if ( .not. allocated(noahmp%water%flux%CompactionSnowAging) ) then
       allocate( noahmp%water%flux%CompactionSnowAging(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%CompactionSnowAging)
    endif
    if ( .not. allocated(noahmp%water%flux%CompactionSnowBurden) ) then
       allocate( noahmp%water%flux%CompactionSnowBurden(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%CompactionSnowBurden)
    endif
    if ( .not. allocated(noahmp%water%flux%CompactionSnowMelt) ) then
       allocate( noahmp%water%flux%CompactionSnowMelt(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%CompactionSnowMelt)
    endif
    if ( .not. allocated(noahmp%water%flux%CompactionSnowTot) ) then
       allocate( noahmp%water%flux%CompactionSnowTot(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%CompactionSnowTot)
    endif
    if ( .not. allocated(noahmp%water%flux%TranspWatLossSoil) ) then
       allocate( noahmp%water%flux%TranspWatLossSoil(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%TranspWatLossSoil)
    endif
    if ( .not. allocated(noahmp%water%flux%TranspWatLossSoilAcc) ) then
       allocate( noahmp%water%flux%TranspWatLossSoilAcc(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%TranspWatLossSoilAcc)
    endif
    if ( .not. allocated(noahmp%water%flux%TranspWatLossSoilMean) ) then
       allocate( noahmp%water%flux%TranspWatLossSoilMean(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%TranspWatLossSoilMean)
    endif
    if ( .not. allocated(noahmp%water%flux%OutflowSnowLayer) ) then
       allocate( noahmp%water%flux%OutflowSnowLayer(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%OutflowSnowLayer)
    endif

    ! SNICAR flux arrays
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( .not. allocated(noahmp%water%flux%SnowFreezeRate) ) then
          allocate( noahmp%water%flux%SnowFreezeRate(ITS:ITE,-NumSnowLayerMax+1:0,JTS:JTE) )
          !$acc enter data create(noahmp%water%flux%SnowFreezeRate)
       endif
    endif

    ! Allocate 3D water parameter arrays and transfer to GPU
    if ( .not. allocated(noahmp%water%param%SoilMoistureSat) ) then
       allocate( noahmp%water%param%SoilMoistureSat(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilMoistureSat)
    endif
    if ( .not. allocated(noahmp%water%param%SoilMoistureWilt) ) then
       allocate( noahmp%water%param%SoilMoistureWilt(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilMoistureWilt)
    endif
    if ( .not. allocated(noahmp%water%param%SoilMoistureFieldCap) ) then
       allocate( noahmp%water%param%SoilMoistureFieldCap(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilMoistureFieldCap)
    endif
    if ( .not. allocated(noahmp%water%param%SoilMoistureDry) ) then
       allocate( noahmp%water%param%SoilMoistureDry(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilMoistureDry)
    endif
    if ( .not. allocated(noahmp%water%param%SoilWatDiffusivitySat) ) then
       allocate( noahmp%water%param%SoilWatDiffusivitySat(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilWatDiffusivitySat)
    endif
    if ( .not. allocated(noahmp%water%param%SoilWatConductivitySat) ) then
       allocate( noahmp%water%param%SoilWatConductivitySat(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilWatConductivitySat)
    endif
    if ( .not. allocated(noahmp%water%param%SoilExpCoeffB) ) then
       allocate( noahmp%water%param%SoilExpCoeffB(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilExpCoeffB)
    endif
    if ( .not. allocated(noahmp%water%param%SoilMatPotentialSat) ) then
       allocate( noahmp%water%param%SoilMatPotentialSat(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilMatPotentialSat)
    endif

    ! SNICAR parameter lookup tables - keep as 3D (not spatially varying)
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( .not. allocated(noahmp%water%param%snowage_tau) ) then
          allocate( noahmp%water%param%snowage_tau(NumDensitySnwAgeSnicar,NumTempGradSnwAgeSnicar,NumTempSnwAgeSnicar) )
          !$acc enter data create(noahmp%water%param%snowage_tau)
       endif
       if ( .not. allocated(noahmp%water%param%snowage_kappa) ) then
          allocate( noahmp%water%param%snowage_kappa(NumDensitySnwAgeSnicar,NumTempGradSnwAgeSnicar,NumTempSnwAgeSnicar) )
          !$acc enter data create(noahmp%water%param%snowage_kappa)
       endif
       if ( .not. allocated(noahmp%water%param%snowage_drdt0) ) then
          allocate( noahmp%water%param%snowage_drdt0(NumDensitySnwAgeSnicar,NumTempGradSnwAgeSnicar,NumTempSnwAgeSnicar) )
          !$acc enter data create(noahmp%water%param%snowage_drdt0)
       endif
       !$acc parallel loop vector collapse(3)
       do J = 1,size(noahmp%water%param%snowage_tau,3)
         do LoopInd = 1, size(noahmp%water%param%snowage_tau,2)
           do I = 1, size(noahmp%water%param%snowage_tau,1)
             noahmp%water%param%snowage_tau(I,LoopInd,J)  = undefined_real
             noahmp%water%param%snowage_kappa(I,LoopInd,J)= undefined_real
             noahmp%water%param%snowage_drdt0(I,LoopInd,J)= undefined_real
           end do
         end do
       end do
    endif

    end associate

    ! Initialize all 2D and 3D arrays in parallel loop
    !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
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

end module WaterVarInitMod
