module GroundWaterTopModelMod

!!! Compute groundwater flow and subsurface runoff based on TOPMODEL (Niu et al., 2007)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine GroundWaterTopModel(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: GROUNDWATER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                       ! grid indices
    integer                          :: LoopInd                    ! loop index
    integer                          :: IndUnsatSoil               ! layer index of the first unsaturated layer
    real(kind=8)                     :: SatDegUnsatSoil            ! degree of saturation of IndUnsatSoil layer
    real(kind=kind_noahmp)           :: SoilMatPotFrz              ! soil matric potential (frozen effects) [mm]
    real(kind=kind_noahmp)           :: AquiferWatConduct          ! aquifer hydraulic conductivity [mm/s]
    real(kind=kind_noahmp)           :: WaterHeadTbl               ! water head at water table [mm]
    real(kind=kind_noahmp)           :: WaterHead                  ! water head at layer above water table [mm]
    real(kind=kind_noahmp)           :: WaterFillPore              ! water used to fill air pore [mm]
    real(kind=kind_noahmp)           :: WatConductAcc              ! sum of SoilWatConductTmp*ThicknessSoil
    real(kind=kind_noahmp)           :: SoilMoistureMin            ! minimum soil moisture [m3/m3]
    real(kind=kind_noahmp)           :: WaterExcessSat             ! excessive water above saturation [mm]
    real(kind=kind_noahmp)           :: ThicknessSoil(noahmp%config%domain%NumSoilLayer)          ! layer thickness [mm]
    real(kind=kind_noahmp)           :: DepthSoilMid(noahmp%config%domain%NumSoilLayer)           ! node depth [m]
    real(kind=kind_noahmp)           :: SoilLiqTmp(noahmp%config%domain%NumSoilLayer)             ! liquid water mass [kg/m2 or mm]
    real(kind=kind_noahmp)           :: SoilEffPorosity(noahmp%config%domain%NumSoilLayer)        ! soil effective porosity
    real(kind=kind_noahmp)           :: SoilWatConductTmp(noahmp%config%domain%NumSoilLayer)      ! hydraulic conductivity [mm/s]
    real(kind=kind_noahmp)           :: SoilMoisture(noahmp%config%domain%NumSoilLayer)           ! total soil water content [m3/m3]

   !$acc parallel loop collapse(2) gang vector present(noahmp) &
   !$acc private(LoopInd, IndUnsatSoil, SatDegUnsatSoil, SoilMatPotFrz, AquiferWatConduct) &
   !$acc private(WaterHeadTbl, WaterHead, WaterFillPore, WatConductAcc, SoilMoistureMin, WaterExcessSat) &
   !$acc private(ThicknessSoil, DepthSoilMid, SoilLiqTmp, SoilEffPorosity, SoilWatConductTmp, SoilMoisture)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points

! --------------------------------------------------------------------
    associate(                                                                     &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer    ,& ! in,    number of soil layers
              SoilTimeStep           => noahmp%config%domain%SoilTimeStep    ,& ! in,    noahmp soil timestep [s]
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer       ,& ! in,    depth of soil layer-bottom [m]
              SoilImpervFracMax      => noahmp%water%state%SoilImpervFracMax(I,J) ,& ! in,    maximum soil imperviousness fraction
              SoilIce                => noahmp%water%state%SoilIce                ,& ! in,    soil ice content [m3/m3]
              SoilWatConductivity    => noahmp%water%state%SoilWatConductivity    ,& ! in,    soil hydraulic conductivity [m/s]
              SoilMoistureSat        => noahmp%water%param%SoilMoistureSat        ,& ! in,    saturated value of soil moisture [m3/m3]
              GridTopoIndex          => noahmp%water%param%GridTopoIndex(I,J)     ,& ! in,    gridcell mean topgraphic index (global mean)
              SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat    ,& ! in,    saturated soil matric potential
              SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB          ,& ! in,    soil B parameter
              SpecYieldGw            => noahmp%water%param%SpecYieldGw(I,J)       ,& ! in,    specific yield [-], default:0.2
              MicroPoreContent       => noahmp%water%param%MicroPoreContent(I,J)  ,& ! in,    microprore content (0.0-1.0), default:0.2
              SoilWatConductivitySat => noahmp%water%param%SoilWatConductivitySat ,& ! in,    saturated soil hydraulic conductivity [m/s]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater           ,& ! inout, soil water content [m3/m3]
              WaterTableDepth        => noahmp%water%state%WaterTableDepth(I,J)   ,& ! inout, water table depth [m]
              WaterStorageAquifer    => noahmp%water%state%WaterStorageAquifer(I,J),& ! inout, water storage in aquifer [mm]
              WaterStorageSoilAqf    => noahmp%water%state%WaterStorageSoilAqf(I,J),& ! inout, water storage in aquifer + saturated soil [mm]
              RunoffDecayFac         => noahmp%water%param%RunoffDecayFac(I,J)    ,& ! inout, runoff decay factor (1/m)
              BaseflowCoeff          => noahmp%water%param%BaseflowCoeff(I,J)     ,& ! inout, baseflow coefficient [mm/s]
              RechargeGw             => noahmp%water%flux%RechargeGw(I,J)         ,& ! out,   groundwater recharge rate [mm/s]
              DischargeGw            => noahmp%water%flux%DischargeGw(I,J)         & ! out,   groundwater discharge rate [mm/s]
             )
! ----------------------------------------------------------------------

    ! initialization
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
      DepthSoilMid(LoopInd)      = 0.0
      ThicknessSoil(LoopInd)     = 0.0
      SoilLiqTmp(LoopInd)        = 0.0
      SoilEffPorosity(LoopInd)   = 0.0
      SoilWatConductTmp(LoopInd) = 0.0
      SoilMoisture(LoopInd)      = 0.0
    enddo

    DischargeGw       = 0.0
    RechargeGw        = 0.0

    ! Derive layer-bottom depth in [mm]; KWM:Derive layer thickness in mm
    ThicknessSoil(1) = -DepthSoilLayer(I,1,J) * 1.0e3
    !$acc loop seq
    do LoopInd = 2, NumSoilLayer
       ThicknessSoil(LoopInd) = 1.0e3 * (DepthSoilLayer(I,LoopInd-1,J) - DepthSoilLayer(I,LoopInd,J))
    enddo

    ! Derive node (middle) depth in [m]; KWM: Positive number, depth below ground surface in m
    DepthSoilMid(1) = -DepthSoilLayer(I,1,J) / 2.0
    !$acc loop seq
    do LoopInd = 2, NumSoilLayer
       DepthSoilMid(LoopInd) = -DepthSoilLayer(I,LoopInd-1,J) + &
                               0.5 * (DepthSoilLayer(I,LoopInd-1,J) - DepthSoilLayer(I,LoopInd,J))
    enddo

    ! Convert volumetric soil moisture to mass
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       SoilMoisture(LoopInd)   = SoilLiqWater(I,LoopInd,J) + SoilIce(I,LoopInd,J)
       SoilLiqTmp(LoopInd)        = SoilLiqWater(I,LoopInd,J) * ThicknessSoil(LoopInd)
       SoilEffPorosity(LoopInd)   = max(0.01, SoilMoistureSat(I,LoopInd,J)-SoilIce(I,LoopInd,J))
       SoilWatConductTmp(LoopInd) = 1.0e3 * SoilWatConductivity(I,LoopInd,J)
    enddo

    ! The layer index of the first unsaturated layer (the layer right above the water table)
    IndUnsatSoil = NumSoilLayer
    !$acc loop seq
    do LoopInd = 2, NumSoilLayer
       if ( WaterTableDepth <= -DepthSoilLayer(I,LoopInd,J) ) then
          IndUnsatSoil = LoopInd - 1
          exit
       endif
    enddo

    ! Groundwater discharge [mm/s]
    !RunoffDecayFac    = 6.0
    !BaseflowCoeff     = 5.0
    !DischargeGw       = (1.0 - SoilImpervFracMax) * BaseflowCoeff * &
    !                    exp(-GridTopoIndex) * exp(-RunoffDecayFac * (WaterTableDepth-2.0))
    ! Update from GY Niu 2022
    RunoffDecayFac    = SoilExpCoeffB(I,IndUnsatSoil,J) / 3.0
    BaseflowCoeff     = SoilWatConductTmp(IndUnsatSoil) * 1.0e3 * exp(3.0)  ! [mm/s]
    DischargeGw       = (1.0 - SoilImpervFracMax) * BaseflowCoeff * exp(-GridTopoIndex) * &
                        exp(-RunoffDecayFac * WaterTableDepth)

    ! Matric potential at the layer above the water table
    SatDegUnsatSoil   = min(1.0, SoilMoisture(IndUnsatSoil)/SoilMoistureSat(I,IndUnsatSoil,J))
    SatDegUnsatSoil   = max(SatDegUnsatSoil, real(0.01,kind=8))
    if (SatDegUnsatSoil < 0.01) SatDegUnsatSoil = 0.01
    SoilMatPotFrz     = -SoilMatPotentialSat(I,IndUnsatSoil,J) * 1000.0 * &
                        SatDegUnsatSoil**(-SoilExpCoeffB(I,IndUnsatSoil,J))
    SoilMatPotFrz     = max(-120000.0, MicroPoreContent*SoilMatPotFrz)

    ! Recharge rate qin to groundwater
    AquiferWatConduct = 2.0 * (SoilWatConductTmp(IndUnsatSoil) * SoilWatConductivitySat(I,IndUnsatSoil,J)*1.0e3) / &
                        (SoilWatConductTmp(IndUnsatSoil) + SoilWatConductivitySat(I,IndUnsatSoil,J)*1.0e3)  ! harmonic average, GY Niu's update 2022
    WaterHeadTbl      = -WaterTableDepth * 1.0e3                 !(mm)
    WaterHead         = SoilMatPotFrz - DepthSoilMid(IndUnsatSoil) * 1.0e3   !(mm)
    RechargeGw        = -AquiferWatConduct * (WaterHeadTbl - WaterHead) / &
                        ((WaterTableDepth-DepthSoilMid(IndUnsatSoil)) * 1.0e3)
    RechargeGw        = max(-10.0/SoilTimeStep, min(10.0/SoilTimeStep, RechargeGw))

    ! Water storage in the aquifer + saturated soil
    WaterStorageSoilAqf = WaterStorageSoilAqf + (RechargeGw - DischargeGw) * SoilTimeStep     !(mm)
    if ( IndUnsatSoil == NumSoilLayer ) then
       WaterStorageAquifer      = WaterStorageAquifer + (RechargeGw - DischargeGw) * SoilTimeStep     !(mm)
       WaterStorageSoilAqf      = WaterStorageAquifer
       WaterTableDepth          = (-DepthSoilLayer(I,NumSoilLayer,J) + 25.0) - &
                                  WaterStorageAquifer / 1000.0 / SpecYieldGw      !(m)
       SoilLiqTmp(NumSoilLayer) = SoilLiqTmp(NumSoilLayer) - RechargeGw * SoilTimeStep        ! [mm]
       SoilLiqTmp(NumSoilLayer) = SoilLiqTmp(NumSoilLayer) + max(0.0, (WaterStorageAquifer-5000.0))
       WaterStorageAquifer      = min(WaterStorageAquifer, 5000.0)
    else
       if ( IndUnsatSoil == NumSoilLayer-1 ) then
          WaterTableDepth = -DepthSoilLayer(I,NumSoilLayer,J) - (WaterStorageSoilAqf - SpecYieldGw*1000.0*25.0) / &
                                                            (SoilEffPorosity(NumSoilLayer)) / 1000.0
       else
          WaterFillPore   = 0.0   ! water used to fill soil air pores
          !$acc loop seq
          do LoopInd = IndUnsatSoil+2, NumSoilLayer
             WaterFillPore = WaterFillPore + SoilEffPorosity(LoopInd) * ThicknessSoil(LoopInd)
          enddo
          WaterTableDepth  = -DepthSoilLayer(I,IndUnsatSoil+1,J) - (WaterStorageSoilAqf - SpecYieldGw*1000.0*25.0 - &
                                                                WaterFillPore) / (SoilEffPorosity(IndUnsatSoil+1)) / 1000.0
       endif
       WatConductAcc = 0.0
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          WatConductAcc = WatConductAcc + SoilWatConductTmp(LoopInd) * ThicknessSoil(LoopInd)
       enddo
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer           ! Removing subsurface runoff
          SoilLiqTmp(LoopInd) = SoilLiqTmp(LoopInd) - DischargeGw * SoilTimeStep * &
                                                      SoilWatConductTmp(LoopInd) * ThicknessSoil(LoopInd) / WatConductAcc
       enddo
    endif
    WaterTableDepth = max(1.5, WaterTableDepth)

    ! Limit SoilLiqTmp to be greater than or equal to SoilMoistureMin
    ! Get water needed to bring SoilLiqTmp equal SoilMoistureMin from lower layer.
    SoilMoistureMin = 0.01
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer-1
       if ( SoilLiqTmp(LoopInd) < 0.0 ) then
          WaterExcessSat = SoilMoistureMin - SoilLiqTmp(LoopInd)
       else
          WaterExcessSat = 0.0
       endif
       SoilLiqTmp(LoopInd  ) = SoilLiqTmp(LoopInd  ) + WaterExcessSat
       SoilLiqTmp(LoopInd+1) = SoilLiqTmp(LoopInd+1) - WaterExcessSat
    enddo
    LoopInd = NumSoilLayer
    if ( SoilLiqTmp(LoopInd) < SoilMoistureMin ) then
       WaterExcessSat   = SoilMoistureMin - SoilLiqTmp(LoopInd)
    else
       WaterExcessSat   = 0.0
    endif
    SoilLiqTmp(LoopInd) = SoilLiqTmp(LoopInd) + WaterExcessSat
    WaterStorageAquifer = WaterStorageAquifer - WaterExcessSat
    WaterStorageSoilAqf = WaterStorageSoilAqf - WaterExcessSat

    ! update soil moisture
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
        SoilLiqWater(I,LoopInd,J) = SoilLiqTmp(LoopInd) / ThicknessSoil(LoopInd)
    enddo

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine GroundWaterTopModel

end module GroundWaterTopModelMod
