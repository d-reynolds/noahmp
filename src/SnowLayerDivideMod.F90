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
    real(kind=kind_noahmp)           :: SnowThickTmp(1:noahmp%config%domain%NumSnowLayerMax)        ! snow layer thickness [m]
    real(kind=kind_noahmp)           :: SnowIceTmp(1:noahmp%config%domain%NumSnowLayerMax)          ! partial volume of ice [m3/m3]
    real(kind=kind_noahmp)           :: SnowLiqTmp(1:noahmp%config%domain%NumSnowLayerMax)          ! partial volume of liquid water [m3/m3]
    real(kind=kind_noahmp)           :: TemperatureSnowTmp(1:noahmp%config%domain%NumSnowLayerMax)  ! node temperature [K]
    real(kind=kind_noahmp)           :: MassBChydrophoTmp(1:noahmp%config%domain%NumSnowLayerMax)   ! mass of hydrophobic Black Carbon in snow [kg m-2]
    real(kind=kind_noahmp)           :: MassBChydrophiTmp(1:noahmp%config%domain%NumSnowLayerMax)   ! mass of hydrophillic Black Carbon in snow [kg m-2]
    real(kind=kind_noahmp)           :: MassOChydrophoTmp(1:noahmp%config%domain%NumSnowLayerMax)   ! mass of hydrophobic Organic Carbon in snow [kg m-2]
    real(kind=kind_noahmp)           :: MassOChydrophiTmp(1:noahmp%config%domain%NumSnowLayerMax)   ! mass of hydrophillic Organic Carbon in snow [kg m-2]
    real(kind=kind_noahmp)           :: MassDust1Tmp(1:noahmp%config%domain%NumSnowLayerMax)        ! mass of dust species 1 in snow [kg m-2]
    real(kind=kind_noahmp)           :: MassDust2Tmp(1:noahmp%config%domain%NumSnowLayerMax)        ! mass of dust species 2 in snow [kg m-2]
    real(kind=kind_noahmp)           :: MassDust3Tmp(1:noahmp%config%domain%NumSnowLayerMax)        ! mass of dust species 3 in snow [kg m-2]
    real(kind=kind_noahmp)           :: MassDust4Tmp(1:noahmp%config%domain%NumSnowLayerMax)        ! mass of dust species 4 in snow [kg m-2]
    real(kind=kind_noahmp)           :: MassDust5Tmp(1:noahmp%config%domain%NumSnowLayerMax)        ! mass of dust species 5 in snow [kg m-2]
    real(kind=kind_noahmp)           :: SnowRadiusTmp(1:noahmp%config%domain%NumSnowLayerMax)       ! effective grain radius [microns, m-6]
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

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd, NumSnowLayerTmp, SnowThickCombTmp, SnowIceExtra, SnowLiqExtra, SnowFracExtra, SnowTempGrad) &
    !$acc private(MassBChydrophoExtra, MassBChydrophiExtra, MassOChydrophoExtra, MassOChydrophiExtra, MassDust1Extra, MassDust2Extra, MassDust3Extra, MassDust4Extra, MassDust5Extra) &
    !$acc private(SnowThickTmp, SnowIceTmp, SnowLiqTmp, TemperatureSnowTmp, MassBChydrophoTmp, MassBChydrophiTmp, MassOChydrophoTmp, MassOChydrophiTmp, MassDust1Tmp, MassDust2Tmp, MassDust3Tmp, MassDust4Tmp, MassDust5Tmp, SnowRadiusTmp, I, J)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
      if ( NumSnowLayerNeg(I,J) >= 0 ) cycle  ! no snow layers
      
    ! initialization

    !$acc loop seq
    do LoopInd = 1, NumSnowLayerMax
      SnowIceTmp        (LoopInd) = 0.0
      SnowLiqTmp        (LoopInd) = 0.0
      TemperatureSnowTmp(LoopInd) = 0.0
      SnowThickTmp      (LoopInd) = 0.0

      if ( OptSnowAlbedo == 3 ) then
         MassBChydrophoTmp(LoopInd) = 0.0
         MassBChydrophiTmp(LoopInd) = 0.0
         MassOChydrophoTmp(LoopInd) = 0.0
         MassOChydrophiTmp(LoopInd) = 0.0
         MassDust1Tmp     (LoopInd) = 0.0
         MassDust2Tmp     (LoopInd) = 0.0
         MassDust3Tmp     (LoopInd) = 0.0
         MassDust4Tmp     (LoopInd) = 0.0
         MassDust5Tmp     (LoopInd) = 0.0
         SnowRadiusTmp    (LoopInd) = 0.0
      endif

       if ( LoopInd <= abs(NumSnowLayerNeg(I,J)) ) then
          SnowThickTmp(LoopInd)       = ThicknessSnowSoilLayer(I,LoopInd+NumSnowLayerNeg(I,J),J)
          SnowIceTmp(LoopInd)         = SnowIce(I,LoopInd+NumSnowLayerNeg(I,J),J)
          SnowLiqTmp(LoopInd)         = SnowLiqWater(I,LoopInd+NumSnowLayerNeg(I,J),J)
          TemperatureSnowTmp(LoopInd) = TemperatureSoilSnow(I,LoopInd+NumSnowLayerNeg(I,J),J)

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(LoopInd) = MassBChydropho(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassBChydrophiTmp(LoopInd) = MassBChydrophi(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassOChydrophoTmp(LoopInd) = MassOChydropho(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassOChydrophiTmp(LoopInd) = MassOChydrophi(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassDust1Tmp(LoopInd)      = MassDust1(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassDust2Tmp(LoopInd)      = MassDust2(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassDust3Tmp(LoopInd)      = MassDust3(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassDust4Tmp(LoopInd)      = MassDust4(I,LoopInd+NumSnowLayerNeg(I,J),J)
             MassDust5Tmp(LoopInd)      = MassDust5(I,LoopInd+NumSnowLayerNeg(I,J),J)
             SnowRadiusTmp(LoopInd)     = SnowRadius(I,LoopInd+NumSnowLayerNeg(I,J),J)
          endif
       endif
    enddo

    ! start snow layer division
    NumSnowLayerTmp = abs(NumSnowLayerNeg(I,J))

    if ( NumSnowLayerTmp == 1 ) then
       ! Specify a new snow layer
       if ( SnowThickTmp(1) > 0.05 ) then
          NumSnowLayerTmp       = 2
          SnowThickTmp(1)       = SnowThickTmp(1)/2.0
          SnowIceTmp(1)         = SnowIceTmp(1)/2.0
          SnowLiqTmp(1)         = SnowLiqTmp(1)/2.0
          SnowThickTmp(2)       = SnowThickTmp(1)
          SnowIceTmp(2)         = SnowIceTmp(1)
          SnowLiqTmp(2)         = SnowLiqTmp(1)
          TemperatureSnowTmp(2) = TemperatureSnowTmp(1)
         
          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(1) = MassBChydrophoTmp(1)/2.0
             MassBChydrophoTmp(2) = MassBChydrophoTmp(1)
             MassBChydrophiTmp(1) = MassBChydrophiTmp(1)/2.0
             MassBChydrophiTmp(2) = MassBChydrophiTmp(1)
             MassOChydrophoTmp(1) = MassOChydrophoTmp(1)/2.0
             MassOChydrophoTmp(2) = MassOChydrophoTmp(1)
             MassOChydrophiTmp(1) = MassOChydrophiTmp(1)/2.0
             MassOChydrophiTmp(2) = MassOChydrophiTmp(1)
             MassDust1Tmp(1)      = MassDust1Tmp(1)/2.0
             MassDust1Tmp(2)      = MassDust1Tmp(1)
             MassDust2Tmp(1)      = MassDust2Tmp(1)/2.0
             MassDust2Tmp(2)      = MassDust2Tmp(1)
             MassDust3Tmp(1)      = MassDust3Tmp(1)/2.0
             MassDust3Tmp(2)      = MassDust3Tmp(1)
             MassDust4Tmp(1)      = MassDust4Tmp(1)/2.0
             MassDust4Tmp(2)      = MassDust4Tmp(1)
             MassDust5Tmp(1)      = MassDust5Tmp(1)/2.0
             MassDust5Tmp(2)      = MassDust5Tmp(1)
             SnowRadiusTmp(2)     = SnowRadiusTmp(1)
          endif
       endif
    endif

    if ( NumSnowLayerTmp > 1 ) then
       if ( SnowThickTmp(1) > 0.05 ) then     ! maximum allowed thickness (5cm) for top snow layer
          SnowThickCombTmp     = SnowThickTmp(1) - 0.05
          SnowFracExtra        = SnowThickCombTmp / SnowThickTmp(1)
          SnowIceExtra         = SnowFracExtra * SnowIceTmp(1)
          SnowLiqExtra         = SnowFracExtra * SnowLiqTmp(1)

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoExtra = SnowFracExtra * MassBChydrophoTmp(1)
             MassBChydrophiExtra = SnowFracExtra * MassBChydrophiTmp(1)
             MassOChydrophoExtra = SnowFracExtra * MassOChydrophoTmp(1)
             MassOChydrophiExtra = SnowFracExtra * MassOChydrophiTmp(1)
             MassDust1Extra      = SnowFracExtra * MassDust1Tmp(1)
             MassDust2Extra      = SnowFracExtra * MassDust2Tmp(1)
             MassDust3Extra      = SnowFracExtra * MassDust3Tmp(1)
             MassDust4Extra      = SnowFracExtra * MassDust4Tmp(1)
             MassDust5Extra      = SnowFracExtra * MassDust5Tmp(1)
          endif

          SnowFracExtra        = 0.05 / SnowThickTmp(1)
          SnowIceTmp(1)        = SnowFracExtra * SnowIceTmp(1)
          SnowLiqTmp(1)        = SnowFracExtra * SnowLiqTmp(1)
          SnowThickTmp(1)      = 0.05

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(1) = SnowFracExtra * MassBChydrophoTmp(1)
             MassBChydrophiTmp(1) = SnowFracExtra * MassBChydrophiTmp(1)
             MassOChydrophoTmp(1) = SnowFracExtra * MassOChydrophoTmp(1)
             MassOChydrophiTmp(1) = SnowFracExtra * MassOChydrophiTmp(1)
             MassDust1Tmp(1)      = SnowFracExtra * MassDust1Tmp(1)
             MassDust2Tmp(1)      = SnowFracExtra * MassDust2Tmp(1)
             MassDust3Tmp(1)      = SnowFracExtra * MassDust3Tmp(1)
             MassDust4Tmp(1)      = SnowFracExtra * MassDust4Tmp(1)
             MassDust5Tmp(1)      = SnowFracExtra * MassDust5Tmp(1)

             MassBChydrophoTmp(2) = MassBChydrophoTmp(2) + MassBChydrophoExtra
             MassBChydrophiTmp(2) = MassBChydrophiTmp(2) + MassBChydrophiExtra
             MassOChydrophoTmp(2) = MassOChydrophoTmp(2) + MassOChydrophoExtra
             MassOChydrophiTmp(2) = MassOChydrophiTmp(2) + MassOChydrophiExtra
             MassDust1Tmp(2)      = MassDust1Tmp(2) + MassDust1Extra 
             MassDust2Tmp(2)      = MassDust2Tmp(2) + MassDust2Extra
             MassDust3Tmp(2)      = MassDust3Tmp(2) + MassDust3Extra
             MassDust4Tmp(2)      = MassDust4Tmp(2) + MassDust4Extra
             MassDust5Tmp(2)      = MassDust5Tmp(2) + MassDust5Extra
             SnowRadiusTmp(2)     = (SnowRadiusTmp(2)*(SnowLiqTmp(2)+SnowIceTmp(2))+SnowRadiusTmp(1)*(SnowLiqExtra+SnowIceExtra)) / &
                                    (SnowLiqTmp(2) + SnowIceTmp(2) + SnowLiqExtra + SnowIceExtra) 
          endif

          ! update combined snow water & temperature
          call SnowLayerWaterCombo(SnowThickTmp(2), SnowLiqTmp(2), SnowIceTmp(2), TemperatureSnowTmp(2), &
                                   SnowThickCombTmp, SnowLiqExtra, SnowIceExtra, TemperatureSnowTmp(1))

          ! subdivide a new layer, maximum allowed thickness (20cm) for second snow layer
          if ( (NumSnowLayerTmp <= 2) .and. (SnowThickTmp(2) > 0.20) ) then  ! MB: change limit
         !if ( (NumSnowLayerTmp <= 2) .and. (SnowThickTmp(2) > 0.10) ) then
             NumSnowLayerTmp       = 3
             SnowTempGrad          = (TemperatureSnowTmp(1) - TemperatureSnowTmp(2)) / &
                                     ((SnowThickTmp(1)+SnowThickTmp(2)) / 2.0)
             SnowThickTmp(2)       = SnowThickTmp(2) / 2.0
             SnowIceTmp(2)         = SnowIceTmp(2) / 2.0
             SnowLiqTmp(2)         = SnowLiqTmp(2) / 2.0
             SnowThickTmp(3)       = SnowThickTmp(2)
             SnowIceTmp(3)         = SnowIceTmp(2)
             SnowLiqTmp(3)         = SnowLiqTmp(2)
             TemperatureSnowTmp(3) = TemperatureSnowTmp(2) - SnowTempGrad * SnowThickTmp(2) / 2.0
             if ( TemperatureSnowTmp(3) >= ConstFreezePoint ) then
                TemperatureSnowTmp(3) = TemperatureSnowTmp(2)
             else
                TemperatureSnowTmp(2) = TemperatureSnowTmp(2) + SnowTempGrad * SnowThickTmp(2) / 2.0
             endif

             if ( OptSnowAlbedo == 3 ) then
                MassBChydrophoTmp(2) = MassBChydrophoTmp(2) / 2.0
                MassBChydrophoTmp(3) = MassBChydrophoTmp(2)
                MassBChydrophiTmp(2) = MassBChydrophiTmp(2) / 2.0
                MassBChydrophiTmp(3) = MassBChydrophiTmp(2)
                MassOChydrophoTmp(2) = MassOChydrophoTmp(2) / 2.0
                MassOChydrophoTmp(3) = MassOChydrophoTmp(2)
                MassOChydrophiTmp(2) = MassOChydrophiTmp(2) / 2.0
                MassOChydrophiTmp(3) = MassOChydrophiTmp(2)
                MassDust1Tmp(2)      = MassDust1Tmp(2) / 2.0
                MassDust1Tmp(3)      = MassDust1Tmp(2)
                MassDust2Tmp(2)      = MassDust2Tmp(2) / 2.0
                MassDust2Tmp(3)      = MassDust2Tmp(2)
                MassDust3Tmp(2)      = MassDust3Tmp(2) / 2.0
                MassDust3Tmp(3)      = MassDust3Tmp(2)
                MassDust4Tmp(2)      = MassDust4Tmp(2) / 2.0
                MassDust4Tmp(3)      = MassDust4Tmp(2)
                MassDust5Tmp(2)      = MassDust5Tmp(2) / 2.0
                MassDust5Tmp(3)      = MassDust5Tmp(2)
                SnowRadiusTmp(3)     = SnowRadiusTmp(2)
             endif

          endif
       endif ! if(SnowThickTmp(1) > 0.05)
    endif  ! if (NumSnowLayerTmp > 1)

    if ( NumSnowLayerTmp > 2 ) then
       if ( SnowThickTmp(2) > 0.2 ) then
          SnowThickCombTmp = SnowThickTmp(2) - 0.2
          SnowFracExtra    = SnowThickCombTmp / SnowThickTmp(2)
          SnowIceExtra     = SnowFracExtra * SnowIceTmp(2)
          SnowLiqExtra     = SnowFracExtra * SnowLiqTmp(2)

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoExtra = SnowFracExtra * MassBChydrophoTmp(2)
             MassBChydrophiExtra = SnowFracExtra * MassBChydrophiTmp(2)
             MassOChydrophoExtra = SnowFracExtra * MassOChydrophoTmp(2)
             MassOChydrophiExtra = SnowFracExtra * MassOChydrophiTmp(2)
             MassDust1Extra      = SnowFracExtra * MassDust1Tmp(2)
             MassDust2Extra      = SnowFracExtra * MassDust2Tmp(2)
             MassDust3Extra      = SnowFracExtra * MassDust3Tmp(2)
             MassDust4Extra      = SnowFracExtra * MassDust4Tmp(2)
             MassDust5Extra      = SnowFracExtra * MassDust5Tmp(2)
          endif

          SnowFracExtra    = 0.2 / SnowThickTmp(2)
          SnowIceTmp(2)    = SnowFracExtra * SnowIceTmp(2)
          SnowLiqTmp(2)    = SnowFracExtra * SnowLiqTmp(2)
          SnowThickTmp(2)  = 0.2

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(2) = SnowFracExtra * MassBChydrophoTmp(2)
             MassBChydrophiTmp(2) = SnowFracExtra * MassBChydrophiTmp(2)
             MassOChydrophoTmp(2) = SnowFracExtra * MassOChydrophoTmp(2)
             MassOChydrophiTmp(2) = SnowFracExtra * MassOChydrophiTmp(2)
             MassDust1Tmp(2)      = SnowFracExtra * MassDust1Tmp(2)
             MassDust2Tmp(2)      = SnowFracExtra * MassDust2Tmp(2)
             MassDust3Tmp(2)      = SnowFracExtra * MassDust3Tmp(2)
             MassDust4Tmp(2)      = SnowFracExtra * MassDust4Tmp(2)
             MassDust5Tmp(2)      = SnowFracExtra * MassDust5Tmp(2)

             MassBChydrophoTmp(3) = MassBChydrophoTmp(3) + MassBChydrophoExtra
             MassBChydrophiTmp(3) = MassBChydrophiTmp(3) + MassBChydrophiExtra
             MassOChydrophoTmp(3) = MassOChydrophoTmp(3) + MassOChydrophoExtra
             MassOChydrophiTmp(3) = MassOChydrophiTmp(3) + MassOChydrophiExtra
             MassDust1Tmp(3)      = MassDust1Tmp(3) + MassDust1Extra
             MassDust2Tmp(3)      = MassDust2Tmp(3) + MassDust2Extra
             MassDust3Tmp(3)      = MassDust3Tmp(3) + MassDust3Extra
             MassDust4Tmp(3)      = MassDust4Tmp(3) + MassDust4Extra
             MassDust5Tmp(3)      = MassDust5Tmp(3) + MassDust5Extra
             SnowRadiusTmp(3)     = (SnowRadiusTmp(3)*(SnowLiqTmp(3)+SnowIceTmp(3))+SnowRadiusTmp(2)*(SnowLiqExtra+SnowIceExtra)) / &
                                    (SnowLiqTmp(3) + SnowIceTmp(3) + SnowLiqExtra + SnowIceExtra)
          endif

          ! update combined snow water & temperature
          call SnowLayerWaterCombo(SnowThickTmp(3), SnowLiqTmp(3), SnowIceTmp(3), TemperatureSnowTmp(3), &
                                   SnowThickCombTmp, SnowLiqExtra, SnowIceExtra, TemperatureSnowTmp(2))
       endif
    endif

    NumSnowLayerNeg(I,J) = -NumSnowLayerTmp

    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, 0
       ThicknessSnowSoilLayer(I,LoopInd,J) = SnowThickTmp(LoopInd-NumSnowLayerNeg(I,J))
       SnowIce(I,LoopInd,J)                = SnowIceTmp(LoopInd-NumSnowLayerNeg(I,J))
       SnowLiqWater(I,LoopInd,J)           = SnowLiqTmp(LoopInd-NumSnowLayerNeg(I,J))
       TemperatureSoilSnow(I,LoopInd,J)    = TemperatureSnowTmp(LoopInd-NumSnowLayerNeg(I,J))

       if ( OptSnowAlbedo == 3 ) then
          MassBChydropho(I,LoopInd,J)      = MassBChydrophoTmp(LoopInd-NumSnowLayerNeg(I,J))
          MassBChydrophi(I,LoopInd,J)      = MassBChydrophiTmp(LoopInd-NumSnowLayerNeg(I,J))
          MassOChydropho(I,LoopInd,J)      = MassOChydrophoTmp(LoopInd-NumSnowLayerNeg(I,J))
          MassOChydrophi(I,LoopInd,J)      = MassOChydrophiTmp(LoopInd-NumSnowLayerNeg(I,J))
          MassDust1(I,LoopInd,J)           = MassDust1Tmp(LoopInd-NumSnowLayerNeg(I,J))
          MassDust2(I,LoopInd,J)           = MassDust2Tmp(LoopInd-NumSnowLayerNeg(I,J))
          MassDust3(I,LoopInd,J)           = MassDust3Tmp(LoopInd-NumSnowLayerNeg(I,J))
          MassDust4(I,LoopInd,J)           = MassDust4Tmp(LoopInd-NumSnowLayerNeg(I,J))
          MassDust5(I,LoopInd,J)           = MassDust5Tmp(LoopInd-NumSnowLayerNeg(I,J))
          SnowRadius(I,LoopInd,J)          = SnowRadiusTmp(LoopInd-NumSnowLayerNeg(I,J))
       endif
    enddo

   end do
   end do



    end associate

  end subroutine SnowLayerDivide

end module SnowLayerDivideMod
