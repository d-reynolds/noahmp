module BiochemVarInTransferMod

!!! Transfer input 2-D NoahmpIO Biochemistry variables to 1-D column variable
!!! 1-D variables should be first defined in /src/BiochemVarType.F90
!!! 2-D variables should be first defined in NoahmpIOVarType.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type
  use NoahmpVarType

  implicit none

contains

!=== initialize with input data or table values

  subroutine BiochemVarInTransfer(noahmp, NoahmpIO)

    implicit none

    type(noahmp_type),   intent(inout) :: noahmp
    type(NoahmpIO_type), intent(inout) :: NoahmpIO
    integer :: I, J
    integer :: LoopInd

! -------------------------------------------------------------------------
      !$acc parallel loop collapse(2) present(noahmp, NoahmpIO)
      do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
    associate(                                                   &
              VegType      => noahmp%config%domain%VegType(I,J)      ,&
              CropType     => noahmp%config%domain%CropType(I,J)     ,&
              OptCropModel => noahmp%config%nmlist%OptCropModel  &
             )
! -------------------------------------------------------------------------

    ! biochem state variables
    noahmp%biochem%state%PlantGrowStage(I,J)             = NoahmpIO%PGSXY   (I,J)   
    noahmp%biochem%state%LeafMass(I,J)                   = NoahmpIO%LFMASSXY(I,J)
    noahmp%biochem%state%RootMass(I,J)                   = NoahmpIO%RTMASSXY(I,J)
    noahmp%biochem%state%StemMass(I,J)                   = NoahmpIO%STMASSXY(I,J) 
    noahmp%biochem%state%WoodMass(I,J)                   = NoahmpIO%WOODXY  (I,J) 
    noahmp%biochem%state%CarbonMassDeepSoil(I,J)         = NoahmpIO%STBLCPXY(I,J) 
    noahmp%biochem%state%CarbonMassShallowSoil(I,J)      = NoahmpIO%FASTCPXY(I,J)
    noahmp%biochem%state%GrainMass(I,J)                  = NoahmpIO%GRAINXY (I,J)  
    noahmp%biochem%state%GrowDegreeDay(I,J)              = NoahmpIO%GDDXY   (I,J)  
    noahmp%biochem%state%NitrogenConcFoliage(I,J)        = 1.0  ! for now, set to nitrogen saturation

    ! biochem parameter variables
    noahmp%biochem%param%NitrogenConcFoliageMax(I,J)     = NoahmpIO%FOLNMX_TABLE (VegType)
    noahmp%biochem%param%QuantumEfficiency25C(I,J)       = NoahmpIO%QE25_TABLE   (VegType)
    noahmp%biochem%param%CarboxylRateMax25C(I,J)         = NoahmpIO%VCMX25_TABLE (VegType)
    noahmp%biochem%param%CarboxylRateMaxQ10(I,J)         = NoahmpIO%AVCMX_TABLE  (VegType)
    noahmp%biochem%param%PhotosynPathC3(I,J)             = NoahmpIO%C3PSN_TABLE  (VegType)
    noahmp%biochem%param%SlopeConductToPhotosyn(I,J)     = NoahmpIO%MP_TABLE     (VegType)
    noahmp%biochem%param%RespMaintQ10(I,J)               = NoahmpIO%ARM_TABLE    (VegType)
    noahmp%biochem%param%RespMaintLeaf25C(I,J)           = NoahmpIO%RMF25_TABLE  (VegType)
    noahmp%biochem%param%RespMaintStem25C(I,J)           = NoahmpIO%RMS25_TABLE  (VegType)
    noahmp%biochem%param%RespMaintRoot25C(I,J)           = NoahmpIO%RMR25_TABLE  (VegType)
    noahmp%biochem%param%WoodToRootRatio(I,J)            = NoahmpIO%WRRAT_TABLE  (VegType)
    noahmp%biochem%param%WoodPoolIndex(I,J)              = NoahmpIO%WDPOOL_TABLE (VegType)
    noahmp%biochem%param%TurnoverCoeffLeafVeg(I,J)       = NoahmpIO%LTOVRC_TABLE (VegType)
    noahmp%biochem%param%TemperaureLeafFreeze(I,J)       = NoahmpIO%TDLEF_TABLE  (VegType)
    noahmp%biochem%param%LeafDeathWaterCoeffVeg(I,J)     = NoahmpIO%DILEFW_TABLE (VegType)
    noahmp%biochem%param%LeafDeathTempCoeffVeg(I,J)      = NoahmpIO%DILEFC_TABLE (VegType)
    noahmp%biochem%param%GrowthRespFrac(I,J)             = NoahmpIO%FRAGR_TABLE  (VegType)
    noahmp%biochem%param%MicroRespCoeff(I,J)             = NoahmpIO%MRP_TABLE    (VegType)
    noahmp%biochem%param%TemperatureMinPhotosyn(I,J)     = NoahmpIO%TMIN_TABLE   (VegType)
    noahmp%biochem%param%LeafAreaPerMass1side(I,J)       = NoahmpIO%SLA_TABLE    (VegType)
    noahmp%biochem%param%StemAreaIndexMin(I,J)           = NoahmpIO%XSAMIN_TABLE (VegType)
    noahmp%biochem%param%WoodAllocFac(I,J)               = NoahmpIO%BF_TABLE     (VegType)
    noahmp%biochem%param%WaterStressCoeff(I,J)           = NoahmpIO%WSTRC_TABLE  (VegType)
    noahmp%biochem%param%LeafAreaIndexMin(I,J)           = NoahmpIO%LAIMIN_TABLE (VegType)
    noahmp%biochem%param%TurnoverCoeffRootVeg(I,J)       = NoahmpIO%RTOVRC_TABLE (VegType)
    noahmp%biochem%param%WoodRespCoeff(I,J)              = NoahmpIO%RSWOODC_TABLE(VegType)
    ! crop model specific parameters
    if ( (OptCropModel > 0) .and. (CropType > 0) ) then
       noahmp%biochem%param%DatePlanting(I,J)            = NoahmpIO%PLTDAY_TABLE   (CropType)
       noahmp%biochem%param%DateHarvest(I,J)             = NoahmpIO%HSDAY_TABLE    (CropType)
       noahmp%biochem%param%NitrogenConcFoliageMax(I,J)  = NoahmpIO%FOLNMXI_TABLE  (CropType)
       noahmp%biochem%param%QuantumEfficiency25C(I,J)    = NoahmpIO%QE25I_TABLE    (CropType)
       noahmp%biochem%param%CarboxylRateMax25C(I,J)      = NoahmpIO%VCMX25I_TABLE  (CropType)
       noahmp%biochem%param%CarboxylRateMaxQ10(I,J)      = NoahmpIO%AVCMXI_TABLE   (CropType)
       noahmp%biochem%param%PhotosynPathC3(I,J)          = NoahmpIO%C3PSNI_TABLE   (CropType)
       noahmp%biochem%param%SlopeConductToPhotosyn(I,J)  = NoahmpIO%MPI_TABLE      (CropType)
       noahmp%biochem%param%RespMaintQ10(I,J)            = NoahmpIO%Q10MR_TABLE    (CropType)
       noahmp%biochem%param%RespMaintLeaf25C(I,J)        = NoahmpIO%LFMR25_TABLE   (CropType)
       noahmp%biochem%param%RespMaintStem25C(I,J)        = NoahmpIO%STMR25_TABLE   (CropType)
       noahmp%biochem%param%RespMaintRoot25C(I,J)        = NoahmpIO%RTMR25_TABLE   (CropType)
       noahmp%biochem%param%GrowthRespFrac(I,J)          = NoahmpIO%FRA_GR_TABLE   (CropType)
       noahmp%biochem%param%TemperaureLeafFreeze(I,J)    = NoahmpIO%LEFREEZ_TABLE  (CropType)
       noahmp%biochem%param%LeafAreaPerBiomass(I,J)      = NoahmpIO%BIO2LAI_TABLE  (CropType)
       noahmp%biochem%param%TempBaseGrowDegDay(I,J)      = NoahmpIO%GDDTBASE_TABLE (CropType)
       noahmp%biochem%param%TempMaxGrowDegDay(I,J)       = NoahmpIO%GDDTCUT_TABLE  (CropType)
       noahmp%biochem%param%GrowDegDayEmerg(I,J)         = NoahmpIO%GDDS1_TABLE    (CropType)
       noahmp%biochem%param%GrowDegDayInitVeg(I,J)       = NoahmpIO%GDDS2_TABLE    (CropType)
       noahmp%biochem%param%GrowDegDayPostVeg(I,J)       = NoahmpIO%GDDS3_TABLE    (CropType)
       noahmp%biochem%param%GrowDegDayInitReprod(I,J)    = NoahmpIO%GDDS4_TABLE    (CropType)
       noahmp%biochem%param%GrowDegDayMature(I,J)        = NoahmpIO%GDDS5_TABLE    (CropType)
       noahmp%biochem%param%PhotosynRadFrac(I,J)         = NoahmpIO%I2PAR_TABLE    (CropType)
       noahmp%biochem%param%TempMinCarbonAssim(I,J)      = NoahmpIO%TASSIM0_TABLE  (CropType)
       noahmp%biochem%param%TempMaxCarbonAssim(I,J)      = NoahmpIO%TASSIM1_TABLE  (CropType)
       noahmp%biochem%param%TempMaxCarbonAssimMax(I,J)   = NoahmpIO%TASSIM2_TABLE  (CropType)
       noahmp%biochem%param%CarbonAssimRefMax(I,J)       = NoahmpIO%AREF_TABLE     (CropType)
       noahmp%biochem%param%LightExtCoeff(I,J)           = NoahmpIO%K_TABLE        (CropType)
       noahmp%biochem%param%LightUseEfficiency(I,J)      = NoahmpIO%EPSI_TABLE     (CropType)
       noahmp%biochem%param%CarbonAssimReducFac(I,J)     = NoahmpIO%PSNRF_TABLE    (CropType)
       noahmp%biochem%param%RespMaintGrain25C(I,J)       = NoahmpIO%GRAINMR25_TABLE(CropType)
       !$acc loop seq
       do LoopInd = 1, noahmp%config%domain%NumCropGrowStage
         noahmp%biochem%param%LeafDeathTempCoeffCrop(I,LoopInd,J)  = NoahmpIO%DILE_FC_TABLE  (CropType,LoopInd)
         noahmp%biochem%param%LeafDeathWaterCoeffCrop(I,LoopInd,J) = NoahmpIO%DILE_FW_TABLE  (CropType,LoopInd)
         noahmp%biochem%param%CarbohydrLeafToGrain(I,LoopInd,J)    = NoahmpIO%LFCT_TABLE     (CropType,LoopInd)
         noahmp%biochem%param%CarbohydrStemToGrain(I,LoopInd,J)    = NoahmpIO%STCT_TABLE     (CropType,LoopInd)
         noahmp%biochem%param%CarbohydrRootToGrain(I,LoopInd,J)    = NoahmpIO%RTCT_TABLE     (CropType,LoopInd)
         noahmp%biochem%param%CarbohydrFracToLeaf(I,LoopInd,J)     = NoahmpIO%LFPT_TABLE     (CropType,LoopInd)
         noahmp%biochem%param%CarbohydrFracToStem(I,LoopInd,J)     = NoahmpIO%STPT_TABLE     (CropType,LoopInd)
         noahmp%biochem%param%CarbohydrFracToRoot(I,LoopInd,J)     = NoahmpIO%RTPT_TABLE     (CropType,LoopInd)
         noahmp%biochem%param%CarbohydrFracToGrain(I,LoopInd,J)    = NoahmpIO%GRAINPT_TABLE  (CropType,LoopInd)
         noahmp%biochem%param%TurnoverCoeffLeafCrop(I,LoopInd,J)   = NoahmpIO%LF_OVRC_TABLE  (CropType,LoopInd)
         noahmp%biochem%param%TurnoverCoeffStemCrop(I,LoopInd,J)   = NoahmpIO%ST_OVRC_TABLE  (CropType,LoopInd)
         noahmp%biochem%param%TurnoverCoeffRootCrop(I,LoopInd,J)   = NoahmpIO%RT_OVRC_TABLE  (CropType,LoopInd)
       end do
       if ( OptCropModel == 1 ) then
          if ( (NoahmpIO%PLANTING(I,J)>0) .and. (NoahmpIO%PLANTING(I,J)<367) ) then
             noahmp%biochem%param%DatePlanting(I,J)      = NoahmpIO%PLANTING(I,J)
          endif ! 2D input map exist
          if ( (NoahmpIO%HARVEST(I,J)>0) .and. (NoahmpIO%HARVEST(I,J)<367) ) then
             noahmp%biochem%param%DateHarvest(I,J)       = NoahmpIO%HARVEST(I,J)
          endif ! 2D input map exist
          if ( (NoahmpIO%SEASON_GDD(I,J)>0.0) .and. (NoahmpIO%SEASON_GDD(I,J)<10000.0) ) then
             noahmp%biochem%param%GrowDegDayEmerg(I,J)   = NoahmpIO%SEASON_GDD(I,J) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayEmerg(I,J)
             noahmp%biochem%param%GrowDegDayInitVeg(I,J) = NoahmpIO%SEASON_GDD(I,J) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayInitVeg(I,J)
             noahmp%biochem%param%GrowDegDayPostVeg(I,J) = NoahmpIO%SEASON_GDD(I,J) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayPostVeg(I,J)
             noahmp%biochem%param%GrowDegDayInitReprod(I,J) = NoahmpIO%SEASON_GDD(I,J) / 1770.0 * &
                                                         noahmp%biochem%param%GrowDegDayInitReprod(I,J)
             noahmp%biochem%param%GrowDegDayMature(I,J)  = NoahmpIO%SEASON_GDD(I,J) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayMature(I,J)
          endif ! 2D input map exist
       endif ! OptCropModel == 1
    endif ! activate crop parameters

    if ( noahmp%config%nmlist%OptIrrigation == 2 ) then
       if ( (NoahmpIO%PLANTING(I,J)>0) .and. (NoahmpIO%PLANTING(I,J)<367) ) then
          noahmp%biochem%param%DatePlanting(I,J) = NoahmpIO%PLANTING(I,J)
       endif ! 2D input map exist
       if ( (NoahmpIO%HARVEST(I,J)>0) .and. (NoahmpIO%HARVEST(I,J)<367) ) then
          noahmp%biochem%param%DateHarvest(I,J)  = NoahmpIO%HARVEST (I,J)
       endif ! 2D input map exist
    endif
    
      end associate
    end do
    end do
    !$acc end parallel loop

  end subroutine BiochemVarInTransfer

end module BiochemVarInTransferMod
