module IrrigationTriggerMod

!!! Trigger irrigation if soil moisture less than the management allowable deficit (MAD)
!!! and estimate irrigation water depth [m] using current rootzone soil moisture and field 
!!! capacity. There are two options here to trigger the irrigation scheme based on MAD
!!! OptIrrigation = 1 -> if irrigated fraction > threshold fraction
!!! OptIrrigation = 2 -> if irrigated fraction > threshold fraction and within crop season
!!! OptIrrigation = 3 -> if irrigated fraction > threshold fraction and LeafAreaIndex > threshold LeafAreaIndex

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine IrrigationTrigger(noahmp, I, J)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: TRIGGER_IRRIGATION
! Original code: P. Valayamkunnath (NCAR) <prasanth@ucar.edu> (08/06/2020)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------
!$acc routine seq

    implicit none

    type(noahmp_type), intent(inout) :: noahmp
    integer, intent(in)              :: I, J          ! grid indices

! local variable
    logical                          :: FlagIrri           ! flag for irrigation activation
    integer                          :: LoopInd            ! loop index
    real(kind=kind_noahmp)           :: SoilMoistAvail     ! available soil moisture [m] at timestep
    real(kind=kind_noahmp)           :: SoilMoistAvailMax  ! maximum available moisture [m]
    real(kind=kind_noahmp)           :: IrrigationWater    ! irrigation water amount [m]

    FlagIrri = .true.

    ! check if irrigation is can be activated or not
    if ( noahmp%config%nmlist%OptIrrigation == 2 ) then ! activate irrigation if within crop season
       if ( (noahmp%config%domain%DayJulianInYear < noahmp%biochem%param%DatePlanting(I,J)) .or. &
            (noahmp%config%domain%DayJulianInYear > (noahmp%biochem%param%DateHarvest(I,J)-noahmp%water%param%IrriStopDayBfHarvest(I,J))) ) &
          FlagIrri = .false.
    elseif ( noahmp%config%nmlist%OptIrrigation == 3) then ! activate if LeafAreaIndex > threshold LeafAreaIndex
       if ( noahmp%energy%state%LeafAreaIndex(I,J) < noahmp%water%param%IrriTriggerLaiMin(I,J)) FlagIrri = .false.
    elseif ( (noahmp%config%nmlist%OptIrrigation > 3) .or. (noahmp%config%nmlist%OptIrrigation < 1) ) then
       FlagIrri = .false.
    endif

    if ( FlagIrri .eqv. .true. ) then
       ! estimate available water and field capacity for the root zone
       SoilMoistAvail      = 0.0
       SoilMoistAvailMax   = 0.0
       SoilMoistAvail      = (noahmp%water%state%SoilLiqWater(I,1,J) - noahmp%water%param%SoilMoistureWilt(I,1,J)) * &
                              (-1.0) * noahmp%config%domain%DepthSoilLayer(I,1,J)          ! current soil water (m)
       SoilMoistAvailMax   = (noahmp%water%param%SoilMoistureFieldCap(I,1,J) - noahmp%water%param%SoilMoistureWilt(I,1,J)) * &
                              (-1.0) * noahmp%config%domain%DepthSoilLayer(I,1,J)  ! available water (m)
       do LoopInd = 2, noahmp%water%param%NumSoilLayerRoot(I,J)
         SoilMoistAvail    = SoilMoistAvail + &
                              (noahmp%water%state%SoilLiqWater(I,LoopInd,J) - noahmp%water%param%SoilMoistureWilt(I,LoopInd,J)) * &
                              (noahmp%config%domain%DepthSoilLayer(I,LoopInd-1,J) - noahmp%config%domain%DepthSoilLayer(I,LoopInd,J))
         SoilMoistAvailMax = SoilMoistAvailMax + &
                              (noahmp%water%param%SoilMoistureFieldCap(I,LoopInd,J) - noahmp%water%param%SoilMoistureWilt(I,LoopInd,J)) * &
                              (noahmp%config%domain%DepthSoilLayer(I,LoopInd-1,J) - noahmp%config%domain%DepthSoilLayer(I,LoopInd,J))
       enddo

      ! check if root zone soil moisture < SoilWatDeficitAllow (calibratable)
      if ( (SoilMoistAvail/SoilMoistAvailMax) <= noahmp%water%param%SoilWatDeficitAllow(I,J) ) then
         ! amount of water need to be added to bring soil moisture back to
         ! field capacity, i.e., irrigation water amount (m)
         IrrigationWater = (SoilMoistAvailMax - SoilMoistAvail) * noahmp%water%state%IrrigationFracGrid(I,J) * &
                           noahmp%energy%state%VegFrac(I,J)

         ! sprinkler irrigation amount (m) based on 2D IrrigationFracSprinkler
         if ( (noahmp%water%state%IrrigationAmtSprinkler(I,J) == 0.0) .and. &
              (noahmp%water%state%IrrigationFracSprinkler(I,J) > 0.0) .and. &
              (noahmp%config%nmlist%OptIrrigationMethod == 0) ) then
            noahmp%water%state%IrrigationAmtSprinkler(I,J) = noahmp%water%state%IrrigationFracSprinkler(I,J) * IrrigationWater
            noahmp%water%state%IrrigationCntSprinkler(I,J) = noahmp%water%state%IrrigationCntSprinkler(I,J) + 1
         ! sprinkler irrigation amount (m) based on namelist choice
         elseif ( (noahmp%water%state%IrrigationAmtSprinkler(I,J) == 0.0) .and. &
                  (noahmp%config%nmlist%OptIrrigationMethod == 1) ) then
            noahmp%water%state%IrrigationAmtSprinkler(I,J) = IrrigationWater
            noahmp%water%state%IrrigationCntSprinkler(I,J) = noahmp%water%state%IrrigationCntSprinkler(I,J) + 1
         endif

         ! micro irrigation amount (m) based on 2D IrrigationFracMicro
         if ( (noahmp%water%state%IrrigationAmtMicro(I,J) == 0.0) .and. &
              (noahmp%water%state%IrrigationFracMicro(I,J) > 0.0) .and. &
              (noahmp%config%nmlist%OptIrrigationMethod == 0) ) then
            noahmp%water%state%IrrigationAmtMicro(I,J) = noahmp%water%state%IrrigationFracMicro(I,J) * IrrigationWater
            noahmp%water%state%IrrigationCntMicro(I,J) = noahmp%water%state%IrrigationCntMicro(I,J) + 1
         ! micro irrigation amount (m) based on namelist choice
         elseif ( (noahmp%water%state%IrrigationAmtMicro(I,J) == 0.0) .and. &
                  (noahmp%config%nmlist%OptIrrigationMethod == 2) ) then
            noahmp%water%state%IrrigationAmtMicro(I,J) = IrrigationWater
            noahmp%water%state%IrrigationCntMicro(I,J) = noahmp%water%state%IrrigationCntMicro(I,J) + 1
         endif

         ! flood irrigation amount (m): Assumed to saturate top two layers and
         ! third layer to FC. As water moves from one end of the field to
         ! another, surface layers will be saturated.
         ! flood irrigation amount (m) based on 2D IrrigationFracFlood
         if ( (noahmp%water%state%IrrigationAmtFlood(I,J) == 0.0) .and. &
              (noahmp%water%state%IrrigationFracFlood(I,J) > 0.0) .and. &
              (noahmp%config%nmlist%OptIrrigationMethod == 0) ) then
            noahmp%water%state%IrrigationAmtFlood(I,J) = noahmp%water%state%IrrigationFracFlood(I,J) * IrrigationWater * &
                                                         (1.0/(1.0 - noahmp%water%param%IrriFloodLossFrac(I,J)))
            noahmp%water%state%IrrigationCntFlood(I,J) = noahmp%water%state%IrrigationCntFlood(I,J) + 1
         !flood irrigation amount (m) based on namelist choice
         elseif ( (noahmp%water%state%IrrigationAmtFlood(I,J) == 0.0) .and. &
                  (noahmp%config%nmlist%OptIrrigationMethod == 3) ) then
            noahmp%water%state%IrrigationAmtFlood(I,J) = IrrigationWater * &
                                                         (1.0/(1.0 - noahmp%water%param%IrriFloodLossFrac(I,J)))
            noahmp%water%state%IrrigationCntFlood(I,J) = noahmp%water%state%IrrigationCntFlood(I,J) + 1
         endif
      else
         IrrigationWater                                = 0.0
         noahmp%water%state%IrrigationAmtSprinkler(I,J) = 0.0
         noahmp%water%state%IrrigationAmtMicro(I,J)     = 0.0
         noahmp%water%state%IrrigationAmtFlood(I,J)     = 0.0
      endif

    endif

  end subroutine IrrigationTrigger

end module IrrigationTriggerMod
