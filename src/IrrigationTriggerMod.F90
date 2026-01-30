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

! --------------------------------------------------------------------
    associate(                                                                            &
              DepthSoilLayer          => noahmp%config%domain%DepthSoilLayer             ,& ! in,    depth [m] of layer-bottom from soil surface
              DayJulianInYear         => noahmp%config%domain%DayJulianInYear            ,& ! in,    Julian day of the year
              OptIrrigation           => noahmp%config%nmlist%OptIrrigation              ,& ! in,    irrigation option
              OptIrrigationMethod     => noahmp%config%nmlist%OptIrrigationMethod        ,& ! in,    irrigation method option
              DatePlanting            => noahmp%biochem%param%DatePlanting(I,J)          ,& ! in,    Planting day (day of year)
              DateHarvest             => noahmp%biochem%param%DateHarvest(I,J)           ,& ! in,    Harvest date (day of year)
              SoilMoistureWilt        => noahmp%water%param%SoilMoistureWilt             ,& ! in,    wilting point soil moisture [m3/m3]
              SoilMoistureFieldCap    => noahmp%water%param%SoilMoistureFieldCap         ,& ! in,    reference soil moisture (field capacity) (m3/m3)
              NumSoilLayerRoot        => noahmp%water%param%NumSoilLayerRoot(I,J)        ,& ! in,    number of soil layers with root present
              IrriStopDayBfHarvest    => noahmp%water%param%IrriStopDayBfHarvest(I,J)    ,& ! in,    number of days before harvest date to stop irrigation
              IrriTriggerLaiMin       => noahmp%water%param%IrriTriggerLaiMin(I,J)       ,& ! in,    minimum lai to trigger irrigation
              SoilWatDeficitAllow     => noahmp%water%param%SoilWatDeficitAllow(I,J)     ,& ! in,    management allowable deficit (0-1)
              IrriFloodLossFrac       => noahmp%water%param%IrriFloodLossFrac(I,J)       ,& ! in,    factor of flood irrigation loss
              VegFrac                 => noahmp%energy%state%VegFrac(I,J)                ,& ! in,    greeness vegetation fraction
              LeafAreaIndex           => noahmp%energy%state%LeafAreaIndex(I,J)          ,& ! in,    leaf area index [m2/m2]
              IrrigationFracGrid      => noahmp%water%state%IrrigationFracGrid(I,J)      ,& ! in,    irrigated area fraction of a grid
              SoilLiqWater            => noahmp%water%state%SoilLiqWater                 ,& ! in,    soil water content [m3/m3]
              IrrigationFracMicro     => noahmp%water%state%IrrigationFracMicro(I,J)     ,& ! in,    fraction of grid under micro irrigation (0 to 1)
              IrrigationFracFlood     => noahmp%water%state%IrrigationFracFlood(I,J)     ,& ! in,    fraction of grid under flood irrigation (0 to 1)
              IrrigationFracSprinkler => noahmp%water%state%IrrigationFracSprinkler(I,J) ,& ! in,    sprinkler irrigation fraction (0 to 1)
              IrrigationAmtMicro      => noahmp%water%state%IrrigationAmtMicro(I,J)      ,& ! inout, irrigation water amount [m] to be applied, Micro
              IrrigationAmtFlood      => noahmp%water%state%IrrigationAmtFlood(I,J)      ,& ! inout, irrigation water amount [m] to be applied, Flood
              IrrigationAmtSprinkler  => noahmp%water%state%IrrigationAmtSprinkler(I,J)  ,& ! inout, irrigation water amount [m] to be applied, Sprinkler
              IrrigationCntSprinkler  => noahmp%water%state%IrrigationCntSprinkler(I,J)  ,& ! inout, irrigation event number, Sprinkler
              IrrigationCntMicro      => noahmp%water%state%IrrigationCntMicro(I,J)      ,& ! inout, irrigation event number, Micro
              IrrigationCntFlood      => noahmp%water%state%IrrigationCntFlood(I,J)       & ! inout, irrigation event number, Flood
             )
! ----------------------------------------------------------------------

    FlagIrri = .true.

    ! check if irrigation is can be activated or not
    if ( OptIrrigation == 2 ) then ! activate irrigation if within crop season
       if ( (DayJulianInYear < DatePlanting) .or. (DayJulianInYear > (DateHarvest-IrriStopDayBfHarvest)) ) &
          FlagIrri = .false.
    elseif ( OptIrrigation == 3) then ! activate if LeafAreaIndex > threshold LeafAreaIndex
       if ( LeafAreaIndex < IrriTriggerLaiMin) FlagIrri = .false.
    elseif ( (OptIrrigation > 3) .or. (OptIrrigation < 1) ) then
       FlagIrri = .false.
    endif

    if ( FlagIrri .eqv. .true. ) then
       ! estimate available water and field capacity for the root zone
       SoilMoistAvail      = 0.0
       SoilMoistAvailMax   = 0.0
       SoilMoistAvail      = (SoilLiqWater(I,1,J) - SoilMoistureWilt(I,1,J)) * (-1.0) * DepthSoilLayer(I,1,J)          ! current soil water (m) 
       SoilMoistAvailMax   = (SoilMoistureFieldCap(I,1,J) - SoilMoistureWilt(I,1,J)) * (-1.0) * DepthSoilLayer(I,1,J)  ! available water (m)
       do LoopInd = 2, NumSoilLayerRoot
         SoilMoistAvail    = SoilMoistAvail + (SoilLiqWater(I,LoopInd,J) - SoilMoistureWilt(I,LoopInd,J)) * &
                                              (DepthSoilLayer(I,LoopInd-1,J) - DepthSoilLayer(I,LoopInd,J))
         SoilMoistAvailMax = SoilMoistAvailMax + (SoilMoistureFieldCap(I,LoopInd,J) - SoilMoistureWilt(I,LoopInd,J)) * &
                                                 (DepthSoilLayer(I,LoopInd-1,J) - DepthSoilLayer(I,LoopInd,J))
       enddo

      ! check if root zone soil moisture < SoilWatDeficitAllow (calibratable)
      if ( (SoilMoistAvail/SoilMoistAvailMax) <= SoilWatDeficitAllow ) then
         ! amount of water need to be added to bring soil moisture back to 
         ! field capacity, i.e., irrigation water amount (m)
         IrrigationWater = (SoilMoistAvailMax - SoilMoistAvail) * IrrigationFracGrid * VegFrac

         ! sprinkler irrigation amount (m) based on 2D IrrigationFracSprinkler
         if ( (IrrigationAmtSprinkler == 0.0) .and. (IrrigationFracSprinkler > 0.0) .and. (OptIrrigationMethod == 0) ) then
            IrrigationAmtSprinkler = IrrigationFracSprinkler * IrrigationWater
            IrrigationCntSprinkler = IrrigationCntSprinkler + 1
         ! sprinkler irrigation amount (m) based on namelist choice
         elseif ( (IrrigationAmtSprinkler == 0.0) .and. (OptIrrigationMethod == 1) ) then
            IrrigationAmtSprinkler = IrrigationWater
            IrrigationCntSprinkler = IrrigationCntSprinkler + 1
         endif

         ! micro irrigation amount (m) based on 2D IrrigationFracMicro
         if ( (IrrigationAmtMicro == 0.0) .and. (IrrigationFracMicro > 0.0) .and. (OptIrrigationMethod == 0) ) then
            IrrigationAmtMicro = IrrigationFracMicro * IrrigationWater
            IrrigationCntMicro = IrrigationCntMicro + 1
         ! micro irrigation amount (m) based on namelist choice
         elseif ( (IrrigationAmtMicro == 0.0) .and. (OptIrrigationMethod == 2) ) then
            IrrigationAmtMicro = IrrigationWater
            IrrigationCntMicro = IrrigationCntMicro + 1
         endif

         ! flood irrigation amount (m): Assumed to saturate top two layers and 
         ! third layer to FC. As water moves from one end of the field to
         ! another, surface layers will be saturated. 
         ! flood irrigation amount (m) based on 2D IrrigationFracFlood
         if ( (IrrigationAmtFlood == 0.0) .and. (IrrigationFracFlood > 0.0) .and. (OptIrrigationMethod == 0) ) then
            IrrigationAmtFlood = IrrigationFracFlood * IrrigationWater * (1.0/(1.0 - IrriFloodLossFrac))
            IrrigationCntFlood = IrrigationCntFlood + 1
         !flood irrigation amount (m) based on namelist choice
         elseif ( (IrrigationAmtFlood == 0.0) .and. (OptIrrigationMethod == 3) ) then
            IrrigationAmtFlood = IrrigationWater * (1.0/(1.0 - IrriFloodLossFrac))
            IrrigationCntFlood = IrrigationCntFlood + 1
         endif
      else
         IrrigationWater        = 0.0
         IrrigationAmtSprinkler = 0.0
         IrrigationAmtMicro     = 0.0
         IrrigationAmtFlood     = 0.0
      endif

    endif

    end associate

  end subroutine IrrigationTrigger

end module IrrigationTriggerMod
