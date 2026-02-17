module ShallowWaterTableMmfMod

!!! Diagnoses water table depth and computes recharge when the water table is 
!!! within the resolved soil layers, according to the Miguez-Macho&Fan scheme

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine ShallowWaterTableMMF(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SHALLOWWATERTABLE
! Original code: Miguez-Macho&Fan (Miguez-Macho et al 2007, Fan et al 2007)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd                            ! do-loop index
    integer                          :: IndAbvWatTbl                       ! layer index above water table layer
    integer                          :: IndWatTbl                          ! layer index where the water table layer is
    real(kind=kind_noahmp)           :: WatTblDepthOld                     ! old water table depth
    real(kind=kind_noahmp)           :: ThicknessUpLy                      ! upper layer thickness
    real(kind=kind_noahmp)           :: SoilMoistDeep                      ! deep layer soil moisture
    real(kind=kind_noahmp)           :: DepthSoilLayer0(0:noahmp%config%domain%NumSoilLayer)   ! temporary soil depth
    integer                          :: I, J                              ! grid indices
    associate(                                                                       &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              SoilTimeStep           => noahmp%config%domain%SoilTimeStep           ,& ! in,    noahmp soil timestep [s]
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,    depth of soil layer-bottom [m]
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
              SoilMoistureEqui       => noahmp%water%state%SoilMoistureEqui         ,& ! in,    equilibrium soil water  content [m3/m3]
              SoilMoistureSat        => noahmp%water%param%SoilMoistureSat          ,& ! in,    saturated value of soil moisture [m3/m3]
              SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat      ,& ! in,    saturated soil matric potential [m]
              SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB            ,& ! in,    soil B parameter
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! inout, total soil water content [m3/m3]
              WaterTableDepth        => noahmp%water%state%WaterTableDepth          ,& ! inout, water table depth [m]
              SoilMoistureToWT       => noahmp%water%state%SoilMoistureToWT         ,& ! inout, soil moisture between bottom of soil & water table
              DrainSoilBot           => noahmp%water%flux%DrainSoilBot              ,& ! inout, soil bottom drainage [m/s]
              RechargeGwShallowWT    => noahmp%water%state%RechargeGwShallowWT       & ! out,   groundwater recharge (net vertical flux across water table), positive up
             )

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd, IndAbvWatTbl, IndWatTbl, WatTblDepthOld, ThicknessUpLy) &
    !$acc                                                       private(SoilMoistDeep, DepthSoilLayer0)
      do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
        do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points

    ! initialization
    DepthSoilLayer0(0)              = 0.0
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       DepthSoilLayer0(LoopInd) = DepthSoilLayer(I,LoopInd,J)
    enddo

    ! find the layer where the water table is
    !$acc loop seq
    do LoopInd = NumSoilLayer, 1, -1
       if ( (WaterTableDepth(I,J)+1.0e-6) < DepthSoilLayer0(LoopInd) ) exit
    enddo
    IndAbvWatTbl = LoopInd

    IndWatTbl    = IndAbvWatTbl + 1          ! layer where the water table is
    if ( IndWatTbl <= NumSoilLayer ) then    ! water table depth in the resolved layers
       WatTblDepthOld = WaterTableDepth(I,J)
       if ( SoilMoisture(I,IndWatTbl,J) > SoilMoistureEqui(I,IndWatTbl,J) ) then
          if ( SoilMoisture(I,IndWatTbl,J) == SoilMoistureSat(I,IndWatTbl,J) ) then ! wtd went to the layer above
             WaterTableDepth(I,J)     = DepthSoilLayer0(IndAbvWatTbl)
             RechargeGwShallowWT = -(WatTblDepthOld - WaterTableDepth(I,J)) * &
                                   (SoilMoistureSat(I,IndWatTbl,J) - SoilMoistureEqui(I,IndWatTbl,J))
             IndAbvWatTbl        = IndAbvWatTbl-1
             IndWatTbl           = IndWatTbl-1
             if ( IndWatTbl >= 1 ) then
                if ( SoilMoisture(I,IndWatTbl,J) > SoilMoistureEqui(I,IndWatTbl,J) ) then
                   WatTblDepthOld      = WaterTableDepth(I,J)
                   WaterTableDepth(I,J)     = min((SoilMoisture(I,IndWatTbl,J)*ThicknessSnowSoilLayer(I,IndWatTbl,J) - &
                                              SoilMoistureEqui(I,IndWatTbl,J)*DepthSoilLayer0(IndAbvWatTbl) + &
                                              SoilMoistureSat(I,IndWatTbl,J)*DepthSoilLayer0(IndWatTbl)) /    &
                                             (SoilMoistureSat(I,IndWatTbl,J)-SoilMoistureEqui(I,IndWatTbl,J)),    &
                                             DepthSoilLayer0(IndAbvWatTbl) )
                   RechargeGwShallowWT = RechargeGwShallowWT - (WatTblDepthOld-WaterTableDepth(I,J)) * &
                                         (SoilMoistureSat(I,IndWatTbl,J)-SoilMoistureEqui(I,IndWatTbl,J))
                endif
             endif
          else  ! water table depth stays in the layer
             WaterTableDepth(I,J) = min((SoilMoisture(I,IndWatTbl,J)*ThicknessSnowSoilLayer(I,IndWatTbl,J) - &
                                    SoilMoistureEqui(I,IndWatTbl,J)*DepthSoilLayer0(IndAbvWatTbl) + &
                                    SoilMoistureSat(I,IndWatTbl,J)*DepthSoilLayer0(IndWatTbl) ) /   &
                                   (SoilMoistureSat(I,IndWatTbl,J)-SoilMoistureEqui(I,IndWatTbl,J)),    &
                                   DepthSoilLayer0(IndAbvWatTbl))
             RechargeGwShallowWT = -(WatTblDepthOld-WaterTableDepth(I,J)) * &
                                    (SoilMoistureSat(I,IndWatTbl,J) - SoilMoistureEqui(I,IndWatTbl,J))
          endif
       else   ! water table depth has gone down to the layer below
          WaterTableDepth(I,J)     = DepthSoilLayer0(IndWatTbl)
          RechargeGwShallowWT = -(WatTblDepthOld-WaterTableDepth(I,J)) * &
                                 (SoilMoistureSat(I,IndWatTbl,J) - SoilMoistureEqui(I,IndWatTbl,J))
          IndWatTbl           = IndWatTbl + 1
          IndAbvWatTbl        = IndAbvWatTbl + 1
          ! water table depth crossed to the layer below. Now adjust it there
          if ( IndWatTbl <= NumSoilLayer ) then
             WatTblDepthOld = WaterTableDepth(I,J)
             if ( SoilMoisture(I,IndWatTbl,J) > SoilMoistureEqui(I,IndWatTbl,J) ) then
                WaterTableDepth(I,J)  = min((SoilMoisture(I,IndWatTbl,J)*ThicknessSnowSoilLayer(I,IndWatTbl,J) - &
                                        SoilMoistureEqui(I,IndWatTbl,J)*DepthSoilLayer0(IndAbvWatTbl) + &
                                        SoilMoistureSat(I,IndWatTbl,J)*DepthSoilLayer0(IndWatTbl) ) /   &
                                       (SoilMoistureSat(I,IndWatTbl,J)-SoilMoistureEqui(I,IndWatTbl,J)),    &
                                       DepthSoilLayer0(IndAbvWatTbl))
             else
                WaterTableDepth(I,J)  = DepthSoilLayer0(IndWatTbl)
             endif
             RechargeGwShallowWT = RechargeGwShallowWT - (WatTblDepthOld-WaterTableDepth(I,J)) *         &
                                   (SoilMoistureSat(I,IndWatTbl,J) - SoilMoistureEqui(I,IndWatTbl,J))
          else
             WatTblDepthOld      = WaterTableDepth(I,J)
             ! restore smoi to equilibrium value with water from the ficticious layer below
             ! SoilMoistureToWT  = SoilMoistureToWT - (SoilMoistureEqui(I,NumSoilLayer,J)-SoilMoisture(I,NumSoilLayer,J))
             ! DrainSoilBot      = DrainSoilBot - 1000 * &
             !                     (SoilMoistureEqui(I,NumSoilLayer,J) - SoilMoisture(I,NumSoilLayer,J)) * &
             !                     ThicknessSnowSoilLayer(I,NumSoilLayer,J) / SoilTimeStep
             ! SoilMoisture(I,NumSoilLayer,J) = SoilMoistureEqui(I,NumSoilLayer,J)

             ! adjust water table depth in the ficticious layer below
             SoilMoistDeep       = SoilMoistureSat(I,NumSoilLayer,J) * (-SoilMatPotentialSat(I,NumSoilLayer,J) /          &
                                   (-SoilMatPotentialSat(I,NumSoilLayer,J) - ThicknessSnowSoilLayer(I,NumSoilLayer,J)))** &
                                   (1.0/SoilExpCoeffB(I,NumSoilLayer,J))
             WaterTableDepth(I,J)     = min((SoilMoistureToWT(I,J) * ThicknessSnowSoilLayer(I,NumSoilLayer,J) -                 &
                                        SoilMoistDeep * DepthSoilLayer0(NumSoilLayer) +                           &
                                        SoilMoistureSat(I,NumSoilLayer,J) * (DepthSoilLayer0(NumSoilLayer) -          &
                                        ThicknessSnowSoilLayer(I,NumSoilLayer,J))) /                                  &
                                       (SoilMoistureSat(I,NumSoilLayer,J)-SoilMoistDeep), DepthSoilLayer0(NumSoilLayer))
             RechargeGwShallowWT = RechargeGwShallowWT - (WatTblDepthOld-WaterTableDepth(I,J)) *                       &
                                   (SoilMoistureSat(I,NumSoilLayer,J) - SoilMoistDeep)
          endif
       endif
    else if ( WaterTableDepth(I,J) >= (DepthSoilLayer0(NumSoilLayer)-ThicknessSnowSoilLayer(I,NumSoilLayer,J)) ) then
    ! if water table depth was already below the bottom of the resolved soil crust
       WatTblDepthOld = WaterTableDepth(I,J)
       SoilMoistDeep  = SoilMoistureSat(I,NumSoilLayer,J) * (-SoilMatPotentialSat(I,NumSoilLayer,J) /                     &
                        (-SoilMatPotentialSat(I,NumSoilLayer,J) - ThicknessSnowSoilLayer(I,NumSoilLayer,J)))**            &
                        (1.0/SoilExpCoeffB(I,NumSoilLayer,J))
       if ( SoilMoistureToWT(I,J) > SoilMoistDeep ) then
          WaterTableDepth(I,J) = min((SoilMoistureToWT(I,J) * ThicknessSnowSoilLayer(I,NumSoilLayer,J) -                        &
                                 SoilMoistDeep * DepthSoilLayer0(NumSoilLayer) +                                  &
                                 SoilMoistureSat(I,NumSoilLayer,J) * (DepthSoilLayer0(NumSoilLayer) -                 &
                                 ThicknessSnowSoilLayer(I,NumSoilLayer,J))) /                                         &
                                (SoilMoistureSat(I,NumSoilLayer,J)-SoilMoistDeep), DepthSoilLayer0(NumSoilLayer))
          RechargeGwShallowWT = -(WatTblDepthOld-WaterTableDepth(I,J)) * (SoilMoistureSat(I,NumSoilLayer,J)-SoilMoistDeep)
       else
          RechargeGwShallowWT = -(WatTblDepthOld - (DepthSoilLayer0(NumSoilLayer)-ThicknessSnowSoilLayer(I,NumSoilLayer,J))) * &
                                 (SoilMoistureSat(I,NumSoilLayer,J) - SoilMoistDeep)
          WatTblDepthOld      = DepthSoilLayer0(NumSoilLayer) - ThicknessSnowSoilLayer(I,NumSoilLayer,J)
          ! and now even further down
          ThicknessUpLy       = (SoilMoistDeep - SoilMoistureToWT(I,J)) * ThicknessSnowSoilLayer(I,NumSoilLayer,J) /       &
                                (SoilMoistureSat(I,NumSoilLayer,J) - SoilMoistDeep)
          WaterTableDepth(I,J)     = WatTblDepthOld - ThicknessUpLy
          RechargeGwShallowWT = RechargeGwShallowWT - (SoilMoistureSat(I,NumSoilLayer,J)-SoilMoistDeep) * ThicknessUpLy
          SoilMoistureToWT(I,J)    = SoilMoistDeep
       endif
    endif

    if ( (IndAbvWatTbl < NumSoilLayer) .and. (IndAbvWatTbl > 0) ) then
       SoilMoistureToWT(I,J) = SoilMoistureSat(I,IndAbvWatTbl,J)
    else if ( (IndAbvWatTbl < NumSoilLayer) .and. (IndAbvWatTbl <= 0) ) then
       SoilMoistureToWT(I,J) = SoilMoistureSat(I,1,J)
    endif


   enddo
enddo


    end associate

  end subroutine ShallowWaterTableMMF

end module ShallowWaterTableMmfMod
