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

    ! Allocate 2D water state arrays and transfer to GPU
    if ( .not. allocated(noahmp%water%state%IrrigationCntSprinkler) ) then
       allocate( noahmp%water%state%IrrigationCntSprinkler(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%IrrigationCntSprinkler)
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationCntMicro) ) then
       allocate( noahmp%water%state%IrrigationCntMicro(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%IrrigationCntMicro)
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationCntFlood) ) then
       allocate( noahmp%water%state%IrrigationCntFlood(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%IrrigationCntFlood)
    endif
    if ( .not. allocated(noahmp%water%state%CanopyTotalWater) ) then
       allocate( noahmp%water%state%CanopyTotalWater(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%CanopyTotalWater)
    endif
    if ( .not. allocated(noahmp%water%state%CanopyWetFrac) ) then
       allocate( noahmp%water%state%CanopyWetFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%CanopyWetFrac)
    endif
    if ( .not. allocated(noahmp%water%state%SnowfallDensity) ) then
       allocate( noahmp%water%state%SnowfallDensity(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowfallDensity)
    endif
    if ( .not. allocated(noahmp%water%state%CanopyLiqWater) ) then
       allocate( noahmp%water%state%CanopyLiqWater(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%CanopyLiqWater)
    endif
    if ( .not. allocated(noahmp%water%state%CanopyIce) ) then
       allocate( noahmp%water%state%CanopyIce(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%CanopyIce)
    endif
    if ( .not. allocated(noahmp%water%state%CanopyIceMax) ) then
       allocate( noahmp%water%state%CanopyIceMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%CanopyIceMax)
    endif
    if ( .not. allocated(noahmp%water%state%CanopyLiqWaterMax) ) then
       allocate( noahmp%water%state%CanopyLiqWaterMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%CanopyLiqWaterMax)
    endif
    if ( .not. allocated(noahmp%water%state%SnowDepth) ) then
       allocate( noahmp%water%state%SnowDepth(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowDepth)
    endif
    if ( .not. allocated(noahmp%water%state%SnowWaterEquiv) ) then
       allocate( noahmp%water%state%SnowWaterEquiv(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowWaterEquiv)
    endif
    if ( .not. allocated(noahmp%water%state%SnowWaterEquivPrev) ) then
       allocate( noahmp%water%state%SnowWaterEquivPrev(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowWaterEquivPrev)
    endif
    if ( .not. allocated(noahmp%water%state%PondSfcThinSnwMelt) ) then
       allocate( noahmp%water%state%PondSfcThinSnwMelt(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%PondSfcThinSnwMelt)
    endif
    if ( .not. allocated(noahmp%water%state%PondSfcThinSnwComb) ) then
       allocate( noahmp%water%state%PondSfcThinSnwComb(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%PondSfcThinSnwComb)
    endif
    if ( .not. allocated(noahmp%water%state%PondSfcThinSnwTrans) ) then
       allocate( noahmp%water%state%PondSfcThinSnwTrans(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%PondSfcThinSnwTrans)
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationFracFlood) ) then
       allocate( noahmp%water%state%IrrigationFracFlood(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%IrrigationFracFlood)
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationAmtFlood) ) then
       allocate( noahmp%water%state%IrrigationAmtFlood(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%IrrigationAmtFlood)
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationFracMicro) ) then
       allocate( noahmp%water%state%IrrigationFracMicro(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%IrrigationFracMicro)
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationAmtMicro) ) then
       allocate( noahmp%water%state%IrrigationAmtMicro(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%IrrigationAmtMicro)
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationFracSprinkler) ) then
       allocate( noahmp%water%state%IrrigationFracSprinkler(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%IrrigationFracSprinkler)
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationAmtSprinkler) ) then
       allocate( noahmp%water%state%IrrigationAmtSprinkler(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%IrrigationAmtSprinkler)
    endif
    if ( .not. allocated(noahmp%water%state%WaterTableDepth) ) then
       allocate( noahmp%water%state%WaterTableDepth(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%WaterTableDepth)
    endif
    if ( .not. allocated(noahmp%water%state%SoilIceMax) ) then
       allocate( noahmp%water%state%SoilIceMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilIceMax)
    endif
    if ( .not. allocated(noahmp%water%state%SoilLiqWaterMin) ) then
       allocate( noahmp%water%state%SoilLiqWaterMin(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilLiqWaterMin)
    endif
    if ( .not. allocated(noahmp%water%state%SoilSaturateFrac) ) then
       allocate( noahmp%water%state%SoilSaturateFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilSaturateFrac)
    endif
    if ( .not. allocated(noahmp%water%state%SoilImpervFracMax) ) then
       allocate( noahmp%water%state%SoilImpervFracMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilImpervFracMax)
    endif
    if ( .not. allocated(noahmp%water%state%SoilMoistureToWT) ) then
       allocate( noahmp%water%state%SoilMoistureToWT(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilMoistureToWT)
    endif
    if ( .not. allocated(noahmp%water%state%RechargeGwDeepWT) ) then
       allocate( noahmp%water%state%RechargeGwDeepWT(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%RechargeGwDeepWT)
    endif
    if ( .not. allocated(noahmp%water%state%RechargeGwShallowWT) ) then
       allocate( noahmp%water%state%RechargeGwShallowWT(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%RechargeGwShallowWT)
    endif
    if ( .not. allocated(noahmp%water%state%SoilSaturationExcess) ) then
       allocate( noahmp%water%state%SoilSaturationExcess(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilSaturationExcess)
    endif
    if ( .not. allocated(noahmp%water%state%WaterTableHydro) ) then
       allocate( noahmp%water%state%WaterTableHydro(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%WaterTableHydro)
    endif
    if ( .not. allocated(noahmp%water%state%TileDrainFrac) ) then
       allocate( noahmp%water%state%TileDrainFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%TileDrainFrac)
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageAquifer) ) then
       allocate( noahmp%water%state%WaterStorageAquifer(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%WaterStorageAquifer)
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageSoilAqf) ) then
       allocate( noahmp%water%state%WaterStorageSoilAqf(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%WaterStorageSoilAqf)
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageLake) ) then
       allocate( noahmp%water%state%WaterStorageLake(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%WaterStorageLake)
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageWetland) ) then
       allocate( noahmp%water%state%WaterStorageWetland(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%WaterStorageWetland)
    endif
    if ( .not. allocated(noahmp%water%state%WaterHeadSfc) ) then
       allocate( noahmp%water%state%WaterHeadSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%WaterHeadSfc)
    endif
    if ( .not. allocated(noahmp%water%state%IrrigationFracGrid) ) then
       allocate( noahmp%water%state%IrrigationFracGrid(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%IrrigationFracGrid)
    endif
    if ( .not. allocated(noahmp%water%state%PrecipAreaFrac) ) then
       allocate( noahmp%water%state%PrecipAreaFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%PrecipAreaFrac)
    endif
    if ( .not. allocated(noahmp%water%state%SnowCoverFrac) ) then
       allocate( noahmp%water%state%SnowCoverFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowCoverFrac)
    endif
    if ( .not. allocated(noahmp%water%state%SoilTranspFacAcc) ) then
       allocate( noahmp%water%state%SoilTranspFacAcc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilTranspFacAcc)
    endif
    if ( .not. allocated(noahmp%water%state%FrozenPrecipFrac) ) then
       allocate( noahmp%water%state%FrozenPrecipFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%FrozenPrecipFrac)
    endif
    if ( .not. allocated(noahmp%water%state%SoilWaterRootZone) ) then
       allocate( noahmp%water%state%SoilWaterRootZone(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilWaterRootZone)
    endif
    if ( .not. allocated(noahmp%water%state%SoilWaterStress) ) then
       allocate( noahmp%water%state%SoilWaterStress(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SoilWaterStress)
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageTotBeg) ) then
       allocate( noahmp%water%state%WaterStorageTotBeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%WaterStorageTotBeg)
    endif
    if ( .not. allocated(noahmp%water%state%WaterBalanceError) ) then
       allocate( noahmp%water%state%WaterBalanceError(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%WaterBalanceError)
    endif
    if ( .not. allocated(noahmp%water%state%WaterStorageTotEnd) ) then
       allocate( noahmp%water%state%WaterStorageTotEnd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%WaterStorageTotEnd)
    endif
    if ( .not. allocated(noahmp%water%state%SnowRadiusFresh) ) then
       allocate( noahmp%water%state%SnowRadiusFresh(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%state%SnowRadiusFresh)
    endif

    ! Allocate 2D water flux arrays and transfer to GPU
    if ( .not. allocated(noahmp%water%flux%RainfallRefHeight) ) then
       allocate( noahmp%water%flux%RainfallRefHeight(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%RainfallRefHeight)
    endif
    if ( .not. allocated(noahmp%water%flux%SnowfallRefHeight) ) then
       allocate( noahmp%water%flux%SnowfallRefHeight(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%SnowfallRefHeight)
    endif
    if ( .not. allocated(noahmp%water%flux%PrecipTotRefHeight) ) then
       allocate( noahmp%water%flux%PrecipTotRefHeight(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%PrecipTotRefHeight)
    endif
    if ( .not. allocated(noahmp%water%flux%PrecipConvTotRefHeight) ) then
       allocate( noahmp%water%flux%PrecipConvTotRefHeight(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%PrecipConvTotRefHeight)
    endif
    if ( .not. allocated(noahmp%water%flux%PrecipLargeSclRefHeight) ) then
       allocate( noahmp%water%flux%PrecipLargeSclRefHeight(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%PrecipLargeSclRefHeight)
    endif
    if ( .not. allocated(noahmp%water%flux%EvapCanopyNet) ) then
       allocate( noahmp%water%flux%EvapCanopyNet(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%EvapCanopyNet)
    endif
    if ( .not. allocated(noahmp%water%flux%Transpiration) ) then
       allocate( noahmp%water%flux%Transpiration(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%Transpiration)
    endif
    if ( .not. allocated(noahmp%water%flux%EvapCanopyLiq) ) then
       allocate( noahmp%water%flux%EvapCanopyLiq(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%EvapCanopyLiq)
    endif
    if ( .not. allocated(noahmp%water%flux%DewCanopyLiq) ) then
       allocate( noahmp%water%flux%DewCanopyLiq(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%DewCanopyLiq)
    endif
    if ( .not. allocated(noahmp%water%flux%FrostCanopyIce) ) then
       allocate( noahmp%water%flux%FrostCanopyIce(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%FrostCanopyIce)
    endif
    if ( .not. allocated(noahmp%water%flux%SublimCanopyIce) ) then
       allocate( noahmp%water%flux%SublimCanopyIce(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%SublimCanopyIce)
    endif
    if ( .not. allocated(noahmp%water%flux%MeltCanopyIce) ) then
       allocate( noahmp%water%flux%MeltCanopyIce(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%MeltCanopyIce)
    endif
    if ( .not. allocated(noahmp%water%flux%FreezeCanopyLiq) ) then
       allocate( noahmp%water%flux%FreezeCanopyLiq(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%FreezeCanopyLiq)
    endif
    if ( .not. allocated(noahmp%water%flux%SnowfallGround) ) then
       allocate( noahmp%water%flux%SnowfallGround(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%SnowfallGround)
    endif
    if ( .not. allocated(noahmp%water%flux%SnowDepthIncr) ) then
       allocate( noahmp%water%flux%SnowDepthIncr(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%SnowDepthIncr)
    endif
    if ( .not. allocated(noahmp%water%flux%FrostSnowSfcIce) ) then
       allocate( noahmp%water%flux%FrostSnowSfcIce(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%FrostSnowSfcIce)
    endif
    if ( .not. allocated(noahmp%water%flux%SublimSnowSfcIce) ) then
       allocate( noahmp%water%flux%SublimSnowSfcIce(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%SublimSnowSfcIce)
    endif
    if ( .not. allocated(noahmp%water%flux%RainfallGround) ) then
       allocate( noahmp%water%flux%RainfallGround(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%RainfallGround)
    endif
    if ( .not. allocated(noahmp%water%flux%SnowBotOutflow) ) then
       allocate( noahmp%water%flux%SnowBotOutflow(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%SnowBotOutflow)
    endif
    if ( .not. allocated(noahmp%water%flux%GlacierExcessFlow) ) then
       allocate( noahmp%water%flux%GlacierExcessFlow(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%GlacierExcessFlow)
    endif
    if ( .not. allocated(noahmp%water%flux%IrrigationRateFlood) ) then
       allocate( noahmp%water%flux%IrrigationRateFlood(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%IrrigationRateFlood)
    endif
    if ( .not. allocated(noahmp%water%flux%IrrigationRateMicro) ) then
       allocate( noahmp%water%flux%IrrigationRateMicro(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%IrrigationRateMicro)
    endif
    if ( .not. allocated(noahmp%water%flux%IrrigationRateSprinkler) ) then
       allocate( noahmp%water%flux%IrrigationRateSprinkler(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%IrrigationRateSprinkler)
    endif
    if ( .not. allocated(noahmp%water%flux%IrriEvapLossSprinkler) ) then
       allocate( noahmp%water%flux%IrriEvapLossSprinkler(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%IrriEvapLossSprinkler)
    endif
    if ( .not. allocated(noahmp%water%flux%SoilSfcInflow) ) then
       allocate( noahmp%water%flux%SoilSfcInflow(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%SoilSfcInflow)
    endif
    if ( .not. allocated(noahmp%water%flux%RunoffSurface) ) then
       allocate( noahmp%water%flux%RunoffSurface(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%RunoffSurface)
    endif
    if ( .not. allocated(noahmp%water%flux%RunoffSubsurface) ) then
       allocate( noahmp%water%flux%RunoffSubsurface(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%RunoffSubsurface)
    endif
    if ( .not. allocated(noahmp%water%flux%InfilRateSfc) ) then
       allocate( noahmp%water%flux%InfilRateSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%InfilRateSfc)
    endif
    if ( .not. allocated(noahmp%water%flux%EvapSoilSfcLiq) ) then
       allocate( noahmp%water%flux%EvapSoilSfcLiq(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%EvapSoilSfcLiq)
    endif
    if ( .not. allocated(noahmp%water%flux%DrainSoilBot) ) then
       allocate( noahmp%water%flux%DrainSoilBot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%DrainSoilBot)
    endif
    if ( .not. allocated(noahmp%water%flux%TileDrain) ) then
       allocate( noahmp%water%flux%TileDrain(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%TileDrain)
    endif
    if ( .not. allocated(noahmp%water%flux%RechargeGw) ) then
       allocate( noahmp%water%flux%RechargeGw(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%RechargeGw)
    endif
    if ( .not. allocated(noahmp%water%flux%DischargeGw) ) then
       allocate( noahmp%water%flux%DischargeGw(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%DischargeGw)
    endif
    if ( .not. allocated(noahmp%water%flux%VaporizeGrd) ) then
       allocate( noahmp%water%flux%VaporizeGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%VaporizeGrd)
    endif
    if ( .not. allocated(noahmp%water%flux%CondenseVapGrd) ) then
       allocate( noahmp%water%flux%CondenseVapGrd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%CondenseVapGrd)
    endif
    if ( .not. allocated(noahmp%water%flux%DewSoilSfcLiq) ) then
       allocate( noahmp%water%flux%DewSoilSfcLiq(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%DewSoilSfcLiq)
    endif
    if ( .not. allocated(noahmp%water%flux%EvapIrriSprinkler) ) then
       allocate( noahmp%water%flux%EvapIrriSprinkler(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%EvapIrriSprinkler)
    endif
    if ( .not. allocated(noahmp%water%flux%InterceptCanopyRain) ) then
       allocate( noahmp%water%flux%InterceptCanopyRain(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%InterceptCanopyRain)
    endif
    if ( .not. allocated(noahmp%water%flux%DripCanopyRain) ) then
       allocate( noahmp%water%flux%DripCanopyRain(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%DripCanopyRain)
    endif
    if ( .not. allocated(noahmp%water%flux%ThroughfallRain) ) then
       allocate( noahmp%water%flux%ThroughfallRain(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%ThroughfallRain)
    endif
    if ( .not. allocated(noahmp%water%flux%InterceptCanopySnow) ) then
       allocate( noahmp%water%flux%InterceptCanopySnow(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%InterceptCanopySnow)
    endif
    if ( .not. allocated(noahmp%water%flux%DripCanopySnow) ) then
       allocate( noahmp%water%flux%DripCanopySnow(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%DripCanopySnow)
    endif
    if ( .not. allocated(noahmp%water%flux%ThroughfallSnow) ) then
       allocate( noahmp%water%flux%ThroughfallSnow(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%ThroughfallSnow)
    endif
    if ( .not. allocated(noahmp%water%flux%EvapGroundNet) ) then
       allocate( noahmp%water%flux%EvapGroundNet(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%EvapGroundNet)
    endif
    if ( .not. allocated(noahmp%water%flux%MeltGroundSnow) ) then
       allocate( noahmp%water%flux%MeltGroundSnow(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%MeltGroundSnow)
    endif
    if ( .not. allocated(noahmp%water%flux%WaterToAtmosTotal) ) then
       allocate( noahmp%water%flux%WaterToAtmosTotal(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%WaterToAtmosTotal)
    endif
    if ( .not. allocated(noahmp%water%flux%EvapSoilSfcLiqAcc) ) then
       allocate( noahmp%water%flux%EvapSoilSfcLiqAcc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%EvapSoilSfcLiqAcc)
    endif
    if ( .not. allocated(noahmp%water%flux%SoilSfcInflowAcc) ) then
       allocate( noahmp%water%flux%SoilSfcInflowAcc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%SoilSfcInflowAcc)
    endif
    if ( .not. allocated(noahmp%water%flux%SfcWaterTotChgAcc) ) then
       allocate( noahmp%water%flux%SfcWaterTotChgAcc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%SfcWaterTotChgAcc)
    endif
    if ( .not. allocated(noahmp%water%flux%PrecipTotAcc) ) then
       allocate( noahmp%water%flux%PrecipTotAcc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%PrecipTotAcc)
    endif
    if ( .not. allocated(noahmp%water%flux%EvapCanopyNetAcc) ) then
       allocate( noahmp%water%flux%EvapCanopyNetAcc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%EvapCanopyNetAcc)
    endif
    if ( .not. allocated(noahmp%water%flux%TranspirationAcc) ) then
       allocate( noahmp%water%flux%TranspirationAcc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%TranspirationAcc)
    endif
    if ( .not. allocated(noahmp%water%flux%EvapGroundNetAcc) ) then
       allocate( noahmp%water%flux%EvapGroundNetAcc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%EvapGroundNetAcc)
    endif
    if ( .not. allocated(noahmp%water%flux%GlacierExcessFlowAcc) ) then
       allocate( noahmp%water%flux%GlacierExcessFlowAcc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%GlacierExcessFlowAcc)
    endif
    if ( .not. allocated(noahmp%water%flux%EvapSoilSfcLiqMean) ) then
       allocate( noahmp%water%flux%EvapSoilSfcLiqMean(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%EvapSoilSfcLiqMean)
    endif
    if ( .not. allocated(noahmp%water%flux%SoilSfcInflowMean) ) then
       allocate( noahmp%water%flux%SoilSfcInflowMean(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%flux%SoilSfcInflowMean)
    endif

    ! Allocate 2D water param arrays and transfer to GPU
    if ( .not. allocated(noahmp%water%param%SnowCoverFac) ) then
       allocate( noahmp%water%param%SnowCoverFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCoverFac)
    endif
    if ( .not. allocated(noahmp%water%param%DrainSoilLayerInd) ) then
       allocate( noahmp%water%param%DrainSoilLayerInd(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%DrainSoilLayerInd)
    endif
    if ( .not. allocated(noahmp%water%param%TileDrainTubeDepth) ) then
       allocate( noahmp%water%param%TileDrainTubeDepth(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%TileDrainTubeDepth)
    endif
    if ( .not. allocated(noahmp%water%param%NumSoilLayerRoot) ) then
       allocate( noahmp%water%param%NumSoilLayerRoot(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%NumSoilLayerRoot)
    endif
    if ( .not. allocated(noahmp%water%param%IrriStopDayBfHarvest) ) then
       allocate( noahmp%water%param%IrriStopDayBfHarvest(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%IrriStopDayBfHarvest)
    endif
    if ( .not. allocated(noahmp%water%param%CanopyLiqHoldCap) ) then
       allocate( noahmp%water%param%CanopyLiqHoldCap(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%CanopyLiqHoldCap)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactBurdenFac) ) then
       allocate( noahmp%water%param%SnowCompactBurdenFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCompactBurdenFac)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactAgingFac1) ) then
       allocate( noahmp%water%param%SnowCompactAgingFac1(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCompactAgingFac1)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactAgingFac2) ) then
       allocate( noahmp%water%param%SnowCompactAgingFac2(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCompactAgingFac2)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactAgingFac3) ) then
       allocate( noahmp%water%param%SnowCompactAgingFac3(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCompactAgingFac3)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactAgingMax) ) then
       allocate( noahmp%water%param%SnowCompactAgingMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCompactAgingMax)
    endif
    if ( .not. allocated(noahmp%water%param%SnowViscosityCoeff) ) then
       allocate( noahmp%water%param%SnowViscosityCoeff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowViscosityCoeff)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactmAR24) ) then
       allocate( noahmp%water%param%SnowCompactmAR24(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCompactmAR24)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactbAR24) ) then
       allocate( noahmp%water%param%SnowCompactbAR24(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCompactbAR24)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactP1AR24) ) then
       allocate( noahmp%water%param%SnowCompactP1AR24(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCompactP1AR24)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactP2AR24) ) then
       allocate( noahmp%water%param%SnowCompactP2AR24(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCompactP2AR24)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCompactP3AR24) ) then
       allocate( noahmp%water%param%SnowCompactP3AR24(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCompactP3AR24)
    endif
    if ( .not. allocated(noahmp%water%param%BurdenFacUpAR24) ) then
       allocate( noahmp%water%param%BurdenFacUpAR24(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%BurdenFacUpAR24)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCoverM1AR25) ) then
       allocate( noahmp%water%param%SnowCoverM1AR25(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCoverM1AR25)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCoverM2AR25) ) then
       allocate( noahmp%water%param%SnowCoverM2AR25(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCoverM2AR25)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCoverFac1AR25) ) then
       allocate( noahmp%water%param%SnowCoverFac1AR25(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCoverFac1AR25)
    endif
    if ( .not. allocated(noahmp%water%param%SnowCoverFac2AR25) ) then
       allocate( noahmp%water%param%SnowCoverFac2AR25(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowCoverFac2AR25)
    endif
    if ( .not. allocated(noahmp%water%param%SnowLiqFracMax) ) then
       allocate( noahmp%water%param%SnowLiqFracMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowLiqFracMax)
    endif
    if ( .not. allocated(noahmp%water%param%SnowLiqHoldCap) ) then
       allocate( noahmp%water%param%SnowLiqHoldCap(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowLiqHoldCap)
    endif
    if ( .not. allocated(noahmp%water%param%SnowLiqReleaseFac) ) then
       allocate( noahmp%water%param%SnowLiqReleaseFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowLiqReleaseFac)
    endif
    if ( .not. allocated(noahmp%water%param%IrriFloodRateFac) ) then
       allocate( noahmp%water%param%IrriFloodRateFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%IrriFloodRateFac)
    endif
    if ( .not. allocated(noahmp%water%param%IrriMicroRate) ) then
       allocate( noahmp%water%param%IrriMicroRate(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%IrriMicroRate)
    endif
    if ( .not. allocated(noahmp%water%param%SoilInfilMaxCoeff) ) then
       allocate( noahmp%water%param%SoilInfilMaxCoeff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilInfilMaxCoeff)
    endif
    if ( .not. allocated(noahmp%water%param%SoilImpervFracCoeff) ) then
       allocate( noahmp%water%param%SoilImpervFracCoeff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilImpervFracCoeff)
    endif
    if ( .not. allocated(noahmp%water%param%InfilFacVic) ) then
       allocate( noahmp%water%param%InfilFacVic(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%InfilFacVic)
    endif
    if ( .not. allocated(noahmp%water%param%TensionWatDistrInfl) ) then
       allocate( noahmp%water%param%TensionWatDistrInfl(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%TensionWatDistrInfl)
    endif
    if ( .not. allocated(noahmp%water%param%TensionWatDistrShp) ) then
       allocate( noahmp%water%param%TensionWatDistrShp(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%TensionWatDistrShp)
    endif
    if ( .not. allocated(noahmp%water%param%FreeWatDistrShp) ) then
       allocate( noahmp%water%param%FreeWatDistrShp(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%FreeWatDistrShp)
    endif
    if ( .not. allocated(noahmp%water%param%InfilHeteroDynVic) ) then
       allocate( noahmp%water%param%InfilHeteroDynVic(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%InfilHeteroDynVic)
    endif
    if ( .not. allocated(noahmp%water%param%InfilCapillaryDynVic) ) then
       allocate( noahmp%water%param%InfilCapillaryDynVic(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%InfilCapillaryDynVic)
    endif
    if ( .not. allocated(noahmp%water%param%InfilFacDynVic) ) then
       allocate( noahmp%water%param%InfilFacDynVic(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%InfilFacDynVic)
    endif
    if ( .not. allocated(noahmp%water%param%SoilDrainSlope) ) then
       allocate( noahmp%water%param%SoilDrainSlope(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilDrainSlope)
    endif
    if ( .not. allocated(noahmp%water%param%TileDrainCoeffSp) ) then
       allocate( noahmp%water%param%TileDrainCoeffSp(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%TileDrainCoeffSp)
    endif
    if ( .not. allocated(noahmp%water%param%DrainFacSoilWat) ) then
       allocate( noahmp%water%param%DrainFacSoilWat(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%DrainFacSoilWat)
    endif
    if ( .not. allocated(noahmp%water%param%TileDrainCoeff) ) then
       allocate( noahmp%water%param%TileDrainCoeff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%TileDrainCoeff)
    endif
    if ( .not. allocated(noahmp%water%param%DrainDepthToImperv) ) then
       allocate( noahmp%water%param%DrainDepthToImperv(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%DrainDepthToImperv)
    endif
    if ( .not. allocated(noahmp%water%param%LateralWatCondFac) ) then
       allocate( noahmp%water%param%LateralWatCondFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%LateralWatCondFac)
    endif
    if ( .not. allocated(noahmp%water%param%TileDrainDepth) ) then
       allocate( noahmp%water%param%TileDrainDepth(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%TileDrainDepth)
    endif
    if ( .not. allocated(noahmp%water%param%DrainTubeDist) ) then
       allocate( noahmp%water%param%DrainTubeDist(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%DrainTubeDist)
    endif
    if ( .not. allocated(noahmp%water%param%DrainTubeRadius) ) then
       allocate( noahmp%water%param%DrainTubeRadius(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%DrainTubeRadius)
    endif
    if ( .not. allocated(noahmp%water%param%DrainWatDepToImperv) ) then
       allocate( noahmp%water%param%DrainWatDepToImperv(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%DrainWatDepToImperv)
    endif
    if ( .not. allocated(noahmp%water%param%RunoffDecayFac) ) then
       allocate( noahmp%water%param%RunoffDecayFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%RunoffDecayFac)
    endif
    if ( .not. allocated(noahmp%water%param%BaseflowCoeff) ) then
       allocate( noahmp%water%param%BaseflowCoeff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%BaseflowCoeff)
    endif
    if ( .not. allocated(noahmp%water%param%GridTopoIndex) ) then
       allocate( noahmp%water%param%GridTopoIndex(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%GridTopoIndex)
    endif
    if ( .not. allocated(noahmp%water%param%SoilSfcSatFracMax) ) then
       allocate( noahmp%water%param%SoilSfcSatFracMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilSfcSatFracMax)
    endif
    if ( .not. allocated(noahmp%water%param%SpecYieldGw) ) then
       allocate( noahmp%water%param%SpecYieldGw(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SpecYieldGw)
    endif
    if ( .not. allocated(noahmp%water%param%MicroPoreContent) ) then
       allocate( noahmp%water%param%MicroPoreContent(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%MicroPoreContent)
    endif
    if ( .not. allocated(noahmp%water%param%WaterStorageLakeMax) ) then
       allocate( noahmp%water%param%WaterStorageLakeMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%WaterStorageLakeMax)
    endif
    if ( .not. allocated(noahmp%water%param%SnoWatEqvMaxGlacier) ) then
       allocate( noahmp%water%param%SnoWatEqvMaxGlacier(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnoWatEqvMaxGlacier)
    endif
    if ( .not. allocated(noahmp%water%param%SoilConductivityRef) ) then
       allocate( noahmp%water%param%SoilConductivityRef(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilConductivityRef)
    endif
    if ( .not. allocated(noahmp%water%param%SoilInfilFacRef) ) then
       allocate( noahmp%water%param%SoilInfilFacRef(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilInfilFacRef)
    endif
    if ( .not. allocated(noahmp%water%param%GroundFrzCoeff) ) then
       allocate( noahmp%water%param%GroundFrzCoeff(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%GroundFrzCoeff)
    endif
    if ( .not. allocated(noahmp%water%param%IrriTriggerLaiMin) ) then
       allocate( noahmp%water%param%IrriTriggerLaiMin(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%IrriTriggerLaiMin)
    endif
    if ( .not. allocated(noahmp%water%param%SoilWatDeficitAllow) ) then
       allocate( noahmp%water%param%SoilWatDeficitAllow(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilWatDeficitAllow)
    endif
    if ( .not. allocated(noahmp%water%param%IrriFloodLossFrac) ) then
       allocate( noahmp%water%param%IrriFloodLossFrac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%IrriFloodLossFrac)
    endif
    if ( .not. allocated(noahmp%water%param%IrriSprinklerRate) ) then
       allocate( noahmp%water%param%IrriSprinklerRate(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%IrriSprinklerRate)
    endif
    if ( .not. allocated(noahmp%water%param%IrriFracThreshold) ) then
       allocate( noahmp%water%param%IrriFracThreshold(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%IrriFracThreshold)
    endif
    if ( .not. allocated(noahmp%water%param%IrriStopPrecipThr) ) then
       allocate( noahmp%water%param%IrriStopPrecipThr(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%IrriStopPrecipThr)
    endif
    if ( .not. allocated(noahmp%water%param%SnowfallDensityMax) ) then
       allocate( noahmp%water%param%SnowfallDensityMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowfallDensityMax)
    endif
    if ( .not. allocated(noahmp%water%param%SnowMassFullCoverOld) ) then
       allocate( noahmp%water%param%SnowMassFullCoverOld(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowMassFullCoverOld)
    endif
    if ( .not. allocated(noahmp%water%param%SoilMatPotentialWilt) ) then
       allocate( noahmp%water%param%SoilMatPotentialWilt(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SoilMatPotentialWilt)
    endif
    if ( .not. allocated(noahmp%water%param%SnowMeltFac) ) then
       allocate( noahmp%water%param%SnowMeltFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowMeltFac)
    endif
    if ( .not. allocated(noahmp%water%param%WetlandCapMax) ) then
       allocate( noahmp%water%param%WetlandCapMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%WetlandCapMax)
    endif
    if ( .not. allocated(noahmp%water%param%SnowRadiusMin) ) then
       allocate( noahmp%water%param%SnowRadiusMin(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowRadiusMin)
    endif
    if ( .not. allocated(noahmp%water%param%FreshSnowRadiusMax) ) then
       allocate( noahmp%water%param%FreshSnowRadiusMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%FreshSnowRadiusMax)
    endif
    if ( .not. allocated(noahmp%water%param%SnowRadiusRefrz) ) then
       allocate( noahmp%water%param%SnowRadiusRefrz(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowRadiusRefrz)
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltScale) ) then
       allocate( noahmp%water%param%ScavEffMeltScale(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%ScavEffMeltScale)
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltBCphi) ) then
       allocate( noahmp%water%param%ScavEffMeltBCphi(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%ScavEffMeltBCphi)
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltBCpho) ) then
       allocate( noahmp%water%param%ScavEffMeltBCpho(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%ScavEffMeltBCpho)
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltOCphi) ) then
       allocate( noahmp%water%param%ScavEffMeltOCphi(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%ScavEffMeltOCphi)
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltOCpho) ) then
       allocate( noahmp%water%param%ScavEffMeltOCpho(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%ScavEffMeltOCpho)
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltDust1) ) then
       allocate( noahmp%water%param%ScavEffMeltDust1(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%ScavEffMeltDust1)
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltDust2) ) then
       allocate( noahmp%water%param%ScavEffMeltDust2(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%ScavEffMeltDust2)
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltDust3) ) then
       allocate( noahmp%water%param%ScavEffMeltDust3(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%ScavEffMeltDust3)
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltDust4) ) then
       allocate( noahmp%water%param%ScavEffMeltDust4(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%ScavEffMeltDust4)
    endif
    if ( .not. allocated(noahmp%water%param%ScavEffMeltDust5) ) then
       allocate( noahmp%water%param%ScavEffMeltDust5(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%ScavEffMeltDust5)
    endif
    if ( .not. allocated(noahmp%water%param%SnowRadiusMax) ) then
       allocate( noahmp%water%param%SnowRadiusMax(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowRadiusMax)
    endif
    if ( .not. allocated(noahmp%water%param%SnowWetAgeC1Brun89) ) then
       allocate( noahmp%water%param%SnowWetAgeC1Brun89(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowWetAgeC1Brun89)
    endif
    if ( .not. allocated(noahmp%water%param%SnowWetAgeC2Brun89) ) then
       allocate( noahmp%water%param%SnowWetAgeC2Brun89(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowWetAgeC2Brun89)
    endif
    if ( .not. allocated(noahmp%water%param%SnowAgeScaleFac) ) then
       allocate( noahmp%water%param%SnowAgeScaleFac(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%water%param%SnowAgeScaleFac)
    endif

    end associate

    ! Initialize all 2D and 3D arrays in parallel loop
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

end module WaterVarInitMod
