module SnowLayerDivideMod

!!! Snowpack layer division process
!!! Update snow ice, snow water, snow thickness, snow temperature

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowLayerWaterComboMod, only: SnowLayerWaterCombo

  implicit none

contains

  subroutine SnowLayerDivide(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: DIVIDE
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd                              ! snow layer loop index
    integer                          :: NumSnowLayerTmp                      ! number of snow layer top to bottom
    real(kind=kind_noahmp)           :: SnowThickCombTmp                     ! thickness of the combined [m]
    real(kind=kind_noahmp)           :: SnowIceExtra                         ! extra snow ice to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: SnowLiqExtra                         ! extra snow liquid water to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: SnowFracExtra                        ! fraction of extra snow to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: SnowTempGrad                         ! temperature gradient between two snow layers
    real(kind=kind_noahmp)           :: MassBChydrophoExtra                  ! extra mass of hydrophobic BC in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassBChydrophiExtra                  ! extra mass of hydrophillic BC in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassOChydrophoExtra                  ! extra mass of hydrophobic OC in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassOChydrophiExtra                  ! extra mass of hydrophillic OC in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassDust1Extra                       ! extra mass of dust species 1 in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassDust2Extra                       ! extra mass of dust species 2 in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassDust3Extra                       ! extra mass of dust species 3 in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassDust4Extra                       ! extra mass of dust species 4 in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassDust5Extra                       ! extra mass of dust species 5 in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowThickTmp        ! snow layer thickness [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowIceTmp          ! partial volume of ice [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowLiqTmp          ! partial volume of liquid water [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: TemperatureSnowTmp  ! node temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassBChydrophoTmp   ! mass of hydrophobic Black Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassBChydrophiTmp   ! mass of hydrophillic Black Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassOChydrophoTmp   ! mass of hydrophobic Organic Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassOChydrophiTmp   ! mass of hydrophillic Organic Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassDust1Tmp        ! mass of dust species 1 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassDust2Tmp        ! mass of dust species 2 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassDust3Tmp        ! mass of dust species 3 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassDust4Tmp        ! mass of dust species 4 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MassDust5Tmp        ! mass of dust species 5 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SnowRadiusTmp       ! effective grain radius [microns, m-6]
    integer                          :: I, J                                 ! grid indices

    associate(                                                                       &
              OptSnowAlbedo          => noahmp%config%nmlist%OptSnowAlbedo          ,& ! in,    options for ground snow surface albedo
              NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax        ,& ! in,    maximum number of snow layers
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg        ,& ! inout, actual number of snow layers (negative)
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m]
              TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow     ,& ! inout, snow and soil layer temperature [K]
              SnowIce                => noahmp%water%state%SnowIce                  ,& ! inout, snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater             ,& ! inout, snow layer liquid water [mm]
              MassBChydropho         => noahmp%water%state%MassBChydropho           ,& ! inout, mass of hydrophobic Black Carbon in snow [kg m-2]
              MassBChydrophi         => noahmp%water%state%MassBChydrophi           ,& ! inout, mass of hydrophillic Black Carbon in snow [kg m-2]
              MassOChydropho         => noahmp%water%state%MassOChydropho           ,& ! inout, mass of hydrophobic Organic Carbon in snow [kg m-2]
              MassOChydrophi         => noahmp%water%state%MassOChydrophi           ,& ! inout, mass of hydrophillic Organic Carbon in snow [kg m-2]
              MassDust1              => noahmp%water%state%MassDust1                ,& ! inout, mass of dust species 1 in snow [kg m-2]
              MassDust2              => noahmp%water%state%MassDust2                ,& ! inout, mass of dust species 2 in snow [kg m-2]
              MassDust3              => noahmp%water%state%MassDust3                ,& ! inout, mass of dust species 3 in snow [kg m-2]
              MassDust4              => noahmp%water%state%MassDust4                ,& ! inout, mass of dust species 4 in snow [kg m-2]
              MassDust5              => noahmp%water%state%MassDust5                ,& ! inout, mass of dust species 5 in snow [kg m-2]
              SnowRadius             => noahmp%water%state%SnowRadius                & ! inout, effective grain radius [microns, m-6]
             )

    allocate(SnowThickTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(SnowIceTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(SnowLiqTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(TemperatureSnowTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassBChydrophoTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassBChydrophiTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassOChydrophoTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassOChydrophiTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassDust1Tmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassDust2Tmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassDust3Tmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassDust4Tmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MassDust5Tmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(SnowRadiusTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSnowLayerMax, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    !$acc data create(SnowThickTmp, SnowIceTmp, SnowLiqTmp, TemperatureSnowTmp, &
    !$acc             MassBChydrophoTmp, MassBChydrophiTmp, MassOChydrophoTmp, MassOChydrophiTmp, &
    !$acc             MassDust1Tmp, MassDust2Tmp, MassDust3Tmp, MassDust4Tmp, MassDust5Tmp, SnowRadiusTmp)

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd, NumSnowLayerTmp, SnowThickCombTmp, SnowIceExtra, SnowLiqExtra, SnowFracExtra, SnowTempGrad) &
    !$acc private(MassBChydrophoExtra, MassBChydrophiExtra, MassOChydrophoExtra, MassOChydrophiExtra, MassDust1Extra, MassDust2Extra, MassDust3Extra, MassDust4Extra, MassDust5Extra)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
      if ( NumSnowLayerNeg(I,J) >= 0 ) cycle  ! no snow layers
      
    ! initialization

    !$acc loop seq
    do LoopInd = 1, NumSnowLayerMax
      SnowIceTmp        (I,LoopInd,J) = 0.0
      SnowLiqTmp        (I,LoopInd,J) = 0.0
      TemperatureSnowTmp(I,LoopInd,J) = 0.0
      SnowThickTmp      (I,LoopInd,J) = 0.0

      if ( OptSnowAlbedo == 3 ) then
         MassBChydrophoTmp(I,LoopInd,J) = 0.0
         MassBChydrophiTmp(I,LoopInd,J) = 0.0
         MassOChydrophoTmp(I,LoopInd,J) = 0.0
         MassOChydrophiTmp(I,LoopInd,J) = 0.0
         MassDust1Tmp     (I,LoopInd,J) = 0.0
         MassDust2Tmp     (I,LoopInd,J) = 0.0
         MassDust3Tmp     (I,LoopInd,J) = 0.0
         MassDust4Tmp     (I,LoopInd,J) = 0.0
         MassDust5Tmp     (I,LoopInd,J) = 0.0
         SnowRadiusTmp    (I,LoopInd,J) = 0.0
      endif

       if ( LoopInd <= abs(NumSnowLayerNeg(I,J)) ) then
          SnowThickTmp(I,LoopInd,J)       = ThicknessSnowSoilLayer(I,LoopInd+NumSnowLayerNeg(I,J),J)
          SnowIceTmp(I,LoopInd,J)         = SnowIce(I,LoopInd+NumSnowLayerNeg(I,J),J)
          SnowLiqTmp(I,LoopInd,J)         = SnowLiqWater(I,LoopInd+NumSnowLayerNeg(I,J),J)
          TemperatureSnowTmp(I,LoopInd,J) = TemperatureSoilSnow(I,LoopInd+NumSnowLayerNeg(I,J),J)

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(I,LoopInd,J) = MassBChydropho(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassBChydrophiTmp(I,LoopInd,J) = MassBChydrophi(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassOChydrophoTmp(I,LoopInd,J) = MassOChydropho(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassOChydrophiTmp(I,LoopInd,J) = MassOChydrophi(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassDust1Tmp(I,LoopInd,J)      = MassDust1(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassDust2Tmp(I,LoopInd,J)      = MassDust2(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassDust3Tmp(I,LoopInd,J)      = MassDust3(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassDust4Tmp(I,LoopInd,J)      = MassDust4(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassDust5Tmp(I,LoopInd,J)      = MassDust5(I,LoopInd+NumSnowLayerNeg(I,J),J)
             SnowRadiusTmp(I,LoopInd,J)     = SnowRadius(I,LoopInd+NumSnowLayerNeg(I,J),J)
          endif
       endif
    enddo

    ! start snow layer division
    NumSnowLayerTmp = abs(NumSnowLayerNeg(I,J))

    if ( NumSnowLayerTmp == 1 ) then
       ! Specify a new snow layer
       if ( SnowThickTmp(I,1,J) > 0.05 ) then
          NumSnowLayerTmp       = 2
          SnowThickTmp(I,1,J)       = SnowThickTmp(I,1,J)/2.0
          SnowIceTmp(I,1,J)         = SnowIceTmp(I,1,J)/2.0
          SnowLiqTmp(I,1,J)         = SnowLiqTmp(I,1,J)/2.0
          SnowThickTmp(I,2,J)       = SnowThickTmp(I,1,J)
          SnowIceTmp(I,2,J)         = SnowIceTmp(I,1,J)
          SnowLiqTmp(I,2,J)         = SnowLiqTmp(I,1,J)
          TemperatureSnowTmp(I,2,J) = TemperatureSnowTmp(I,1,J)
         
          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(I,1,J) = MassBChydrophoTmp(I,1,J)/2.0
             MassBChydrophoTmp(I,2,J) = MassBChydrophoTmp(I,1,J)
             MassBChydrophiTmp(I,1,J) = MassBChydrophiTmp(I,1,J)/2.0
             MassBChydrophiTmp(I,2,J) = MassBChydrophiTmp(I,1,J)
             MassOChydrophoTmp(I,1,J) = MassOChydrophoTmp(I,1,J)/2.0
             MassOChydrophoTmp(I,2,J) = MassOChydrophoTmp(I,1,J)
             MassOChydrophiTmp(I,1,J) = MassOChydrophiTmp(I,1,J)/2.0
             MassOChydrophiTmp(I,2,J) = MassOChydrophiTmp(I,1,J)
             MassDust1Tmp(I,1,J)      = MassDust1Tmp(I,1,J)/2.0
             MassDust1Tmp(I,2,J)      = MassDust1Tmp(I,1,J)
             MassDust2Tmp(I,1,J)      = MassDust2Tmp(I,1,J)/2.0
             MassDust2Tmp(I,2,J)      = MassDust2Tmp(I,1,J)
             MassDust3Tmp(I,1,J)      = MassDust3Tmp(I,1,J)/2.0
             MassDust3Tmp(I,2,J)      = MassDust3Tmp(I,1,J)
             MassDust4Tmp(I,1,J)      = MassDust4Tmp(I,1,J)/2.0
             MassDust4Tmp(I,2,J)      = MassDust4Tmp(I,1,J)
             MassDust5Tmp(I,1,J)      = MassDust5Tmp(I,1,J)/2.0
             MassDust5Tmp(I,2,J)      = MassDust5Tmp(I,1,J)
             SnowRadiusTmp(I,2,J)     = SnowRadiusTmp(I,1,J)
          endif
       endif
    endif

    if ( NumSnowLayerTmp > 1 ) then
       if ( SnowThickTmp(I,1,J) > 0.05 ) then     ! maximum allowed thickness (5cm) for top snow layer
          SnowThickCombTmp     = SnowThickTmp(I,1,J) - 0.05
          SnowFracExtra        = SnowThickCombTmp / SnowThickTmp(I,1,J)
          SnowIceExtra         = SnowFracExtra * SnowIceTmp(I,1,J)
          SnowLiqExtra         = SnowFracExtra * SnowLiqTmp(I,1,J)

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoExtra = SnowFracExtra * MassBChydrophoTmp(I,1,J)
             MassBChydrophiExtra = SnowFracExtra * MassBChydrophiTmp(I,1,J)
             MassOChydrophoExtra = SnowFracExtra * MassOChydrophoTmp(I,1,J)
             MassOChydrophiExtra = SnowFracExtra * MassOChydrophiTmp(I,1,J)
             MassDust1Extra      = SnowFracExtra * MassDust1Tmp(I,1,J)
             MassDust2Extra      = SnowFracExtra * MassDust2Tmp(I,1,J)
             MassDust3Extra      = SnowFracExtra * MassDust3Tmp(I,1,J)
             MassDust4Extra      = SnowFracExtra * MassDust4Tmp(I,1,J)
             MassDust5Extra      = SnowFracExtra * MassDust5Tmp(I,1,J)
          endif

          SnowFracExtra        = 0.05 / SnowThickTmp(I,1,J)
          SnowIceTmp(I,1,J)        = SnowFracExtra * SnowIceTmp(I,1,J)
          SnowLiqTmp(I,1,J)        = SnowFracExtra * SnowLiqTmp(I,1,J)
          SnowThickTmp(I,1,J)      = 0.05

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(I,1,J) = SnowFracExtra * MassBChydrophoTmp(I,1,J)
             MassBChydrophiTmp(I,1,J) = SnowFracExtra * MassBChydrophiTmp(I,1,J)
             MassOChydrophoTmp(I,1,J) = SnowFracExtra * MassOChydrophoTmp(I,1,J)
             MassOChydrophiTmp(I,1,J) = SnowFracExtra * MassOChydrophiTmp(I,1,J)
             MassDust1Tmp(I,1,J)      = SnowFracExtra * MassDust1Tmp(I,1,J)
             MassDust2Tmp(I,1,J)      = SnowFracExtra * MassDust2Tmp(I,1,J)
             MassDust3Tmp(I,1,J)      = SnowFracExtra * MassDust3Tmp(I,1,J)
             MassDust4Tmp(I,1,J)      = SnowFracExtra * MassDust4Tmp(I,1,J)
             MassDust5Tmp(I,1,J)      = SnowFracExtra * MassDust5Tmp(I,1,J)

             MassBChydrophoTmp(I,2,J) = MassBChydrophoTmp(I,2,J) + MassBChydrophoExtra
             MassBChydrophiTmp(I,2,J) = MassBChydrophiTmp(I,2,J) + MassBChydrophiExtra
             MassOChydrophoTmp(I,2,J) = MassOChydrophoTmp(I,2,J) + MassOChydrophoExtra
             MassOChydrophiTmp(I,2,J) = MassOChydrophiTmp(I,2,J) + MassOChydrophiExtra
             MassDust1Tmp(I,2,J)      = MassDust1Tmp(I,2,J) + MassDust1Extra 
             MassDust2Tmp(I,2,J)      = MassDust2Tmp(I,2,J) + MassDust2Extra
             MassDust3Tmp(I,2,J)      = MassDust3Tmp(I,2,J) + MassDust3Extra
             MassDust4Tmp(I,2,J)      = MassDust4Tmp(I,2,J) + MassDust4Extra
             MassDust5Tmp(I,2,J)      = MassDust5Tmp(I,2,J) + MassDust5Extra
             SnowRadiusTmp(I,2,J)     = (SnowRadiusTmp(I,2,J)*(SnowLiqTmp(I,2,J)+SnowIceTmp(I,2,J))+SnowRadiusTmp(I,1,J)*(SnowLiqExtra+SnowIceExtra)) / &
                                    (SnowLiqTmp(I,2,J) + SnowIceTmp(I,2,J) + SnowLiqExtra + SnowIceExtra) 
          endif

          ! update combined snow water & temperature
          call SnowLayerWaterCombo(SnowThickTmp(I,2,J), SnowLiqTmp(I,2,J), SnowIceTmp(I,2,J), TemperatureSnowTmp(I,2,J), &
                                   SnowThickCombTmp, SnowLiqExtra, SnowIceExtra, TemperatureSnowTmp(I,1,J))

          ! subdivide a new layer, maximum allowed thickness (20cm) for second snow layer
          if ( (NumSnowLayerTmp <= 2) .and. (SnowThickTmp(I,2,J) > 0.20) ) then  ! MB: change limit
         !if ( (NumSnowLayerTmp <= 2) .and. (SnowThickTmp(I,2,J) > 0.10) ) then
             NumSnowLayerTmp       = 3
             SnowTempGrad          = (TemperatureSnowTmp(I,1,J) - TemperatureSnowTmp(I,2,J)) / &
                                     ((SnowThickTmp(I,1,J)+SnowThickTmp(I,2,J)) / 2.0)
             SnowThickTmp(I,2,J)       = SnowThickTmp(I,2,J) / 2.0
             SnowIceTmp(I,2,J)         = SnowIceTmp(I,2,J) / 2.0
             SnowLiqTmp(I,2,J)         = SnowLiqTmp(I,2,J) / 2.0
             SnowThickTmp(I,3,J)       = SnowThickTmp(I,2,J)
             SnowIceTmp(I,3,J)         = SnowIceTmp(I,2,J)
             SnowLiqTmp(I,3,J)         = SnowLiqTmp(I,2,J)
             TemperatureSnowTmp(I,3,J) = TemperatureSnowTmp(I,2,J) - SnowTempGrad * SnowThickTmp(I,2,J) / 2.0
             if ( TemperatureSnowTmp(I,3,J) >= ConstFreezePoint ) then
                TemperatureSnowTmp(I,3,J) = TemperatureSnowTmp(I,2,J)
             else
                TemperatureSnowTmp(I,2,J) = TemperatureSnowTmp(I,2,J) + SnowTempGrad * SnowThickTmp(I,2,J) / 2.0
             endif

             if ( OptSnowAlbedo == 3 ) then
                MassBChydrophoTmp(I,2,J) = MassBChydrophoTmp(I,2,J) / 2.0
                MassBChydrophoTmp(I,3,J) = MassBChydrophoTmp(I,2,J)
                MassBChydrophiTmp(I,2,J) = MassBChydrophiTmp(I,2,J) / 2.0
                MassBChydrophiTmp(I,3,J) = MassBChydrophiTmp(I,2,J)
                MassOChydrophoTmp(I,2,J) = MassOChydrophoTmp(I,2,J) / 2.0
                MassOChydrophoTmp(I,3,J) = MassOChydrophoTmp(I,2,J)
                MassOChydrophiTmp(I,2,J) = MassOChydrophiTmp(I,2,J) / 2.0
                MassOChydrophiTmp(I,3,J) = MassOChydrophiTmp(I,2,J)
                MassDust1Tmp(I,2,J)      = MassDust1Tmp(I,2,J) / 2.0
                MassDust1Tmp(I,3,J)      = MassDust1Tmp(I,2,J)
                MassDust2Tmp(I,2,J)      = MassDust2Tmp(I,2,J) / 2.0
                MassDust2Tmp(I,3,J)      = MassDust2Tmp(I,2,J)
                MassDust3Tmp(I,2,J)      = MassDust3Tmp(I,2,J) / 2.0
                MassDust3Tmp(I,3,J)      = MassDust3Tmp(I,2,J)
                MassDust4Tmp(I,2,J)      = MassDust4Tmp(I,2,J) / 2.0
                MassDust4Tmp(I,3,J)      = MassDust4Tmp(I,2,J)
                MassDust5Tmp(I,2,J)      = MassDust5Tmp(I,2,J) / 2.0
                MassDust5Tmp(I,3,J)      = MassDust5Tmp(I,2,J)
                SnowRadiusTmp(I,3,J)     = SnowRadiusTmp(I,2,J)
             endif

          endif
       endif ! if(SnowThickTmp(I,1,J) > 0.05)
    endif  ! if (NumSnowLayerTmp > 1)

    if ( NumSnowLayerTmp > 2 ) then
       if ( SnowThickTmp(I,2,J) > 0.2 ) then
          SnowThickCombTmp = SnowThickTmp(I,2,J) - 0.2
          SnowFracExtra    = SnowThickCombTmp / SnowThickTmp(I,2,J)
          SnowIceExtra     = SnowFracExtra * SnowIceTmp(I,2,J)
          SnowLiqExtra     = SnowFracExtra * SnowLiqTmp(I,2,J)

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoExtra = SnowFracExtra * MassBChydrophoTmp(I,2,J)
             MassBChydrophiExtra = SnowFracExtra * MassBChydrophiTmp(I,2,J)
             MassOChydrophoExtra = SnowFracExtra * MassOChydrophoTmp(I,2,J)
             MassOChydrophiExtra = SnowFracExtra * MassOChydrophiTmp(I,2,J)
             MassDust1Extra      = SnowFracExtra * MassDust1Tmp(I,2,J)
             MassDust2Extra      = SnowFracExtra * MassDust2Tmp(I,2,J)
             MassDust3Extra      = SnowFracExtra * MassDust3Tmp(I,2,J)
             MassDust4Extra      = SnowFracExtra * MassDust4Tmp(I,2,J)
             MassDust5Extra      = SnowFracExtra * MassDust5Tmp(I,2,J)
          endif

          SnowFracExtra    = 0.2 / SnowThickTmp(I,2,J)
          SnowIceTmp(I,2,J)    = SnowFracExtra * SnowIceTmp(I,2,J)
          SnowLiqTmp(I,2,J)    = SnowFracExtra * SnowLiqTmp(I,2,J)
          SnowThickTmp(I,2,J)  = 0.2

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(I,2,J) = SnowFracExtra * MassBChydrophoTmp(I,2,J)
             MassBChydrophiTmp(I,2,J) = SnowFracExtra * MassBChydrophiTmp(I,2,J)
             MassOChydrophoTmp(I,2,J) = SnowFracExtra * MassOChydrophoTmp(I,2,J)
             MassOChydrophiTmp(I,2,J) = SnowFracExtra * MassOChydrophiTmp(I,2,J)
             MassDust1Tmp(I,2,J)      = SnowFracExtra * MassDust1Tmp(I,2,J)
             MassDust2Tmp(I,2,J)      = SnowFracExtra * MassDust2Tmp(I,2,J)
             MassDust3Tmp(I,2,J)      = SnowFracExtra * MassDust3Tmp(I,2,J)
             MassDust4Tmp(I,2,J)      = SnowFracExtra * MassDust4Tmp(I,2,J)
             MassDust5Tmp(I,2,J)      = SnowFracExtra * MassDust5Tmp(I,2,J)

             MassBChydrophoTmp(I,3,J) = MassBChydrophoTmp(I,3,J) + MassBChydrophoExtra
             MassBChydrophiTmp(I,3,J) = MassBChydrophiTmp(I,3,J) + MassBChydrophiExtra
             MassOChydrophoTmp(I,3,J) = MassOChydrophoTmp(I,3,J) + MassOChydrophoExtra
             MassOChydrophiTmp(I,3,J) = MassOChydrophiTmp(I,3,J) + MassOChydrophiExtra
             MassDust1Tmp(I,3,J)      = MassDust1Tmp(I,3,J) + MassDust1Extra
             MassDust2Tmp(I,3,J)      = MassDust2Tmp(I,3,J) + MassDust2Extra
             MassDust3Tmp(I,3,J)      = MassDust3Tmp(I,3,J) + MassDust3Extra
             MassDust4Tmp(I,3,J)      = MassDust4Tmp(I,3,J) + MassDust4Extra
             MassDust5Tmp(I,3,J)      = MassDust5Tmp(I,3,J) + MassDust5Extra
             SnowRadiusTmp(I,3,J)     = (SnowRadiusTmp(I,3,J)*(SnowLiqTmp(I,3,J)+SnowIceTmp(I,3,J))+SnowRadiusTmp(I,2,J)*(SnowLiqExtra+SnowIceExtra)) / &
                                    (SnowLiqTmp(I,3,J) + SnowIceTmp(I,3,J) + SnowLiqExtra + SnowIceExtra)
          endif

          ! update combined snow water & temperature
          call SnowLayerWaterCombo(SnowThickTmp(I,3,J), SnowLiqTmp(I,3,J), SnowIceTmp(I,3,J), TemperatureSnowTmp(I,3,J), &
                                   SnowThickCombTmp, SnowLiqExtra, SnowIceExtra, TemperatureSnowTmp(I,2,J))
       endif
    endif

    NumSnowLayerNeg(I,J) = -NumSnowLayerTmp

    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, 0
       ThicknessSnowSoilLayer(I,LoopInd,J) = SnowThickTmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
       SnowIce(I,LoopInd,J)                = SnowIceTmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
       SnowLiqWater(I,LoopInd,J)           = SnowLiqTmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
       TemperatureSoilSnow(I,LoopInd,J)    = TemperatureSnowTmp(I,LoopInd-NumSnowLayerNeg(I,J),J)

       if ( OptSnowAlbedo == 3 ) then
          MassBChydropho(I,LoopInd,J)      = MassBChydrophoTmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
          MassBChydrophi(I,LoopInd,J)      = MassBChydrophiTmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
          MassOChydropho(I,LoopInd,J)      = MassOChydrophoTmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
          MassOChydrophi(I,LoopInd,J)      = MassOChydrophiTmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
          MassDust1(I,LoopInd,J)           = MassDust1Tmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
          MassDust2(I,LoopInd,J)           = MassDust2Tmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
          MassDust3(I,LoopInd,J)           = MassDust3Tmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
          MassDust4(I,LoopInd,J)           = MassDust4Tmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
          MassDust5(I,LoopInd,J)           = MassDust5Tmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
          SnowRadius(I,LoopInd,J)          = SnowRadiusTmp(I,LoopInd-NumSnowLayerNeg(I,J),J)
       endif
    enddo

   end do
   end do

    !$acc end data
    deallocate(SnowThickTmp)
    deallocate(SnowIceTmp)
    deallocate(SnowLiqTmp)
    deallocate(TemperatureSnowTmp)
    deallocate(MassBChydrophoTmp)
    deallocate(MassBChydrophiTmp)
    deallocate(MassOChydrophoTmp)
    deallocate(MassOChydrophiTmp)
    deallocate(MassDust1Tmp)
    deallocate(MassDust2Tmp)
    deallocate(MassDust3Tmp)
    deallocate(MassDust4Tmp)
    deallocate(MassDust5Tmp)
    deallocate(SnowRadiusTmp)

    end associate

  end subroutine SnowLayerDivide

end module SnowLayerDivideMod
