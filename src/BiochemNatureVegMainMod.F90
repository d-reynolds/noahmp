module BiochemNatureVegMainMod

!!! Main Biogeochemistry module for dynamic generic vegetation (as opposed to explicit crop scheme)
!!! currently only include carbon processes (RE Dickinson et al.(1998) and Guo-Yue Niu(2004))

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use CarbonFluxNatureVegMod,  only : CarbonFluxNatureVeg

  implicit none

contains

  subroutine BiochemNatureVegMain(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: CARBON
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: I, J         ! grid indices
    integer                          :: LoopInd      ! loop index


    associate(                                                                       &
              VegType                => noahmp%config%domain%VegType                ,& ! in,    vegetation type
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,    depth [m] of layer-bottom from soil surface
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    snow/soil layer thickness [m]
              IndexWaterPoint        => noahmp%config%domain%IndexWaterPoint        ,& ! in,    water point flag
              IndexIcePoint          => noahmp%config%domain%IndexIcePoint          ,& ! in,    land ice flag
              IndexBarrenPoint       => noahmp%config%domain%IndexBarrenPoint       ,& ! in,    bare soil flag
              FlagUrban              => noahmp%config%domain%FlagUrban              ,& ! in,    urban point flag
              NumSoilLayerRoot       => noahmp%water%param%NumSoilLayerRoot         ,& ! in,    number of soil layers with root present
              SoilMoistureSat        => noahmp%water%param%SoilMoistureSat          ,& ! in,    saturated value of soil moisture [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! in,    soil moisture (ice + liq.) [m3/m3]
              SoilTranspFacAcc       => noahmp%water%state%SoilTranspFacAcc         ,& ! in,    accumulated soil water transpiration factor (0 to 1)
              LeafAreaPerMass1side   => noahmp%biochem%param%LeafAreaPerMass1side   ,& ! in,    single-side leaf area per Kg [m2/kg]
              LeafMass               => noahmp%biochem%state%LeafMass               ,& ! inout, leaf mass [g/m2]
              RootMass               => noahmp%biochem%state%RootMass               ,& ! inout, mass of fine roots [g/m2]
              StemMass               => noahmp%biochem%state%StemMass               ,& ! inout, stem mass [g/m2]
              WoodMass               => noahmp%biochem%state%WoodMass               ,& ! inout, mass of wood (incl. woody roots) [g/m2]
              CarbonMassDeepSoil     => noahmp%biochem%state%CarbonMassDeepSoil     ,& ! inout, stable carbon in deep soil [g/m2]
              CarbonMassShallowSoil  => noahmp%biochem%state%CarbonMassShallowSoil  ,& ! inout, short-lived carbon in shallow soil [g/m2]
              LeafAreaIndex          => noahmp%energy%state%LeafAreaIndex           ,& ! inout, leaf area index
              StemAreaIndex          => noahmp%energy%state%StemAreaIndex           ,& ! inout, stem area index
              GrossPriProduction     => noahmp%biochem%flux%GrossPriProduction      ,& ! out,   net instantaneous assimilation [g/m2/s C]
              NetPriProductionTot    => noahmp%biochem%flux%NetPriProductionTot     ,& ! out,   net primary productivity [g/m2/s C]
              NetEcoExchange         => noahmp%biochem%flux%NetEcoExchange          ,& ! out,   net ecosystem exchange [g/m2/s CO2]
              RespirationPlantTot    => noahmp%biochem%flux%RespirationPlantTot     ,& ! out,   total plant respiration [g/m2/s C]
              RespirationSoilOrg     => noahmp%biochem%flux%RespirationSoilOrg      ,& ! out,   soil organic respiration [g/m2/s C]
              CarbonMassSoilTot      => noahmp%biochem%state%CarbonMassSoilTot      ,& ! out,   total soil carbon [g/m2 C]
              CarbonMassLiveTot      => noahmp%biochem%state%CarbonMassLiveTot      ,& ! out,   total living carbon ([g/m2 C]
              SoilWaterRootZone      => noahmp%water%state%SoilWaterRootZone        ,& ! out,   root zone soil water
              SoilWaterStress        => noahmp%water%state%SoilWaterStress          ,& ! out,   water stress coeficient (1. for wilting)
              LeafAreaPerMass        => noahmp%biochem%state%LeafAreaPerMass         & ! out,   leaf area per unit mass [m2/g]
             )

   !$acc parallel loop collapse(2) gang vector default(present) &
   !$acc private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! condition to cycle moved from NoahmpMainMod to here
    if ( noahmp%config%domain%FlagDynamicVeg(I,J) .eqv. .false. ) cycle
!-----------------------------------------------------------------------

    ! initialize
    NetEcoExchange(I,J)      = 0.0
    NetPriProductionTot(I,J) = 0.0
    GrossPriProduction(I,J)  = 0.0

    ! no biogeochemistry in non-vegetated points
    if ( (VegType(I,J) == IndexWaterPoint) .or. (VegType(I,J) == IndexBarrenPoint) .or. &
         (VegType(I,J) == IndexIcePoint  ) .or. (FlagUrban(I,J) .eqv. .true.) ) then
       LeafAreaIndex(I,J)         = 0.0
       StemAreaIndex(I,J)         = 0.0
       GrossPriProduction(I,J)    = 0.0
       NetPriProductionTot(I,J)   = 0.0
       NetEcoExchange(I,J)        = 0.0
       RespirationPlantTot(I,J)   = 0.0
       RespirationSoilOrg(I,J)    = 0.0
       CarbonMassSoilTot(I,J)     = 0.0
       CarbonMassLiveTot(I,J)     = 0.0
       LeafMass(I,J)              = 0.0
       RootMass(I,J)              = 0.0
       StemMass(I,J)              = 0.0
       WoodMass(I,J)              = 0.0
       CarbonMassDeepSoil(I,J)    = 0.0
       CarbonMassShallowSoil(I,J) = 0.0
       cycle
    endif

    ! start biogeochemistry process
    LeafAreaPerMass(I,J) = LeafAreaPerMass1side(I,J) / 1000.0   ! m2/kg -> m2/g

    ! water stress
    SoilWaterStress(I,J)   = 1.0 - SoilTranspFacAcc(I,J)
    SoilWaterRootZone(I,J) = 0.0
    !$acc loop seq
    do LoopInd = 1, NumSoilLayerRoot(I,J)
       SoilWaterRootZone(I,J) = SoilWaterRootZone(I,J) + SoilMoisture(I,LoopInd,J) / SoilMoistureSat(I,LoopInd,J) * &
                                               ThicknessSnowSoilLayer(I,LoopInd,J) / (-DepthSoilLayer(I,NumSoilLayerRoot(I,J),J))
    enddo


      end do
    end do
   !$acc end parallel loop

    ! start carbon process
    call CarbonFluxNatureVeg(noahmp)


    end associate

  end subroutine BiochemNatureVegMain

end module BiochemNatureVegMainMod
