module BiochemCropMainMod

!!! Main Biogeochemistry module for dynamic crop (as opposed to natural vegetation)
!!! currently only include carbon processes (RE Dickinson et al.(1998) and Liu et al., 2014))
 
  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use CarbonFluxCropMod,     only : CarbonFluxCrop
  use CropGrowDegreeDayMod,  only : CropGrowDegreeDay
  use CropPhotosynthesisMod, only : CropPhotosynthesis
        
  implicit none
        
contains
     
  subroutine BiochemCropMain(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: CARBON_CROP
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Modified by Xing Liu, 2014
! Refactered code: C. He, P. Valayamkunnath & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------
        
    implicit none
        
    type(noahmp_type), intent(inout) :: noahmp
    
! local variables
    integer                          :: I, J         ! grid indices
    integer                          :: LoopInd      ! loop index
    
!-------------------------------------------------------------------------
   !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! condition to cycle moved from NoahmoMain to here for better GPU performance
    if ( .not. (noahmp%config%domain%FlagDynamicCrop(I,J) .and. (noahmp%config%nmlist%OptCropModel == 1) ) ) cycle
!------------------------------------------------------------------------
    associate(                                                                       &
              VegType                => noahmp%config%domain%VegType(I,J)                ,& ! in,    vegetation type
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,    depth [m] of layer-bottom from soil surface
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    snow/soil layer thickness [m]
              IndexWaterPoint        => noahmp%config%domain%IndexWaterPoint        ,& ! in,    water point flag
              IndexIcePoint          => noahmp%config%domain%IndexIcePoint          ,& ! in,    land ice flag
              IndexBarrenPoint       => noahmp%config%domain%IndexBarrenPoint       ,& ! in,    bare soil flag
              FlagUrban              => noahmp%config%domain%FlagUrban(I,J)              ,& ! in,    urban point flag
              NumSoilLayerRoot       => noahmp%water%param%NumSoilLayerRoot(I,J)         ,& ! in,    number of soil layers with root present
              SoilMoistureSat        => noahmp%water%param%SoilMoistureSat          ,& ! in,    saturated value of soil moisture [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! in,    soil moisture (ice + liq.) [m3/m3]
              SoilTranspFacAcc       => noahmp%water%state%SoilTranspFacAcc(I,J)         ,& ! in,    accumulated soil water transpiration factor (0 to 1)
              LeafMass               => noahmp%biochem%state%LeafMass(I,J)               ,& ! inout, leaf mass [g/m2]
              RootMass               => noahmp%biochem%state%RootMass(I,J)               ,& ! inout, mass of fine roots [g/m2]
              StemMass               => noahmp%biochem%state%StemMass(I,J)               ,& ! inout, stem mass [g/m2]
              WoodMass               => noahmp%biochem%state%WoodMass(I,J)               ,& ! inout, mass of wood (incl. woody roots) [g/m2]
              CarbonMassDeepSoil     => noahmp%biochem%state%CarbonMassDeepSoil(I,J)     ,& ! inout, stable carbon in deep soil [g/m2]
              CarbonMassShallowSoil  => noahmp%biochem%state%CarbonMassShallowSoil(I,J)  ,& ! inout, short-lived carbon in shallow soil [g/m2]
              LeafAreaIndex          => noahmp%energy%state%LeafAreaIndex(I,J)           ,& ! inout, leaf area index
              StemAreaIndex          => noahmp%energy%state%StemAreaIndex(I,J)           ,& ! inout, stem area index
              GrossPriProduction     => noahmp%biochem%flux%GrossPriProduction(I,J)      ,& ! out,   net instantaneous assimilation [g/m2/s C]
              NetPriProductionTot    => noahmp%biochem%flux%NetPriProductionTot(I,J)     ,& ! out,   net primary productivity [g/m2/s C]
              NetEcoExchange         => noahmp%biochem%flux%NetEcoExchange(I,J)          ,& ! out,   net ecosystem exchange [g/m2/s CO2]
              RespirationPlantTot    => noahmp%biochem%flux%RespirationPlantTot(I,J)     ,& ! out,   total plant respiration [g/m2/s C]
              RespirationSoilOrg     => noahmp%biochem%flux%RespirationSoilOrg(I,J)      ,& ! out,   soil organic respiration [g/m2/s C]
              CarbonMassSoilTot      => noahmp%biochem%state%CarbonMassSoilTot(I,J)      ,& ! out,   total soil carbon [g/m2 C]
              CarbonMassLiveTot      => noahmp%biochem%state%CarbonMassLiveTot(I,J)      ,& ! out,   total living carbon ([g/m2 C]
              GrainMass              => noahmp%biochem%state%GrainMass(I,J)              ,& ! out,   mass of grain [g/m2]
              SoilWaterRootZone      => noahmp%water%state%SoilWaterRootZone(I,J)        ,& ! out,   root zone soil water
              SoilWaterStress        => noahmp%water%state%SoilWaterStress(I,J)           & ! out,   water stress coeficient (1.0 for wilting)
             )
!------------------------------------------------------------------------

    ! initialize
    NetEcoExchange      = 0.0
    NetPriProductionTot = 0.0
    GrossPriProduction  = 0.0

    ! no biogeochemistry in non-vegetated points
    if ( (VegType == IndexWaterPoint) .or. (VegType == IndexBarrenPoint) .or. &
         (VegType == IndexIcePoint  ) .or. (FlagUrban .eqv. .true.) ) then
       LeafAreaIndex         = 0.0
       StemAreaIndex         = 0.0
       GrossPriProduction    = 0.0
       NetPriProductionTot   = 0.0
       NetEcoExchange        = 0.0
       RespirationPlantTot   = 0.0
       RespirationSoilOrg    = 0.0
       CarbonMassSoilTot     = 0.0
       CarbonMassLiveTot     = 0.0
       LeafMass              = 0.0
       RootMass              = 0.0
       StemMass              = 0.0
       WoodMass              = 0.0
       CarbonMassDeepSoil    = 0.0
       CarbonMassShallowSoil = 0.0
       GrainMass             = 0.0
       cycle
    endif
    
    ! start biogeochemistry process
    ! water stress
    SoilWaterStress   = 1.0 - SoilTranspFacAcc
    SoilWaterRootZone = 0.0
    !$acc loop seq
    do LoopInd = 1, NumSoilLayerRoot
       SoilWaterRootZone = SoilWaterRootZone + SoilMoisture(I,LoopInd,J) / SoilMoistureSat(I,LoopInd,J) * &
                                               ThicknessSnowSoilLayer(I,LoopInd,J) / (-DepthSoilLayer(I,LoopInd,J))
    enddo

    end associate

      end do
    end do
   !$acc end parallel loop

    ! start crop carbon process
    ! Note: The following CropPhotosynthesis is not used currently. 
    ! Photosynthesis rate is directly from calculations in the energy part (similar to the treatment in CARBON subroutine)    

    !call CropPhotosynthesis(noahmp)
    call CropGrowDegreeDay(noahmp)
    call CarbonFluxCrop(noahmp)
  end subroutine BiochemCropMain
    
end module BiochemCropMainMod
