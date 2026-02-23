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
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ThicknessSoil          ! layer thickness [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: DepthSoilMid           ! node depth [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilLiqTmp             ! liquid water mass [kg/m2 or mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilEffPorosity        ! soil effective porosity
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilWatConductTmp      ! hydraulic conductivity [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilMoisture           ! total soil water content [m3/m3]

    associate(                                                                     &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer    ,& ! in,    number of soil layers
              SoilTimeStep           => noahmp%config%domain%SoilTimeStep    ,& ! in,    noahmp soil timestep [s]
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer       ,& ! in,    depth of soil layer-bottom [m]
              SoilImpervFracMax      => noahmp%water%state%SoilImpervFracMax ,& ! in,    maximum soil imperviousness fraction
              SoilIce                => noahmp%water%state%SoilIce                ,& ! in,    soil ice content [m3/m3]
              SoilWatConductivity    => noahmp%water%state%SoilWatConductivity    ,& ! in,    soil hydraulic conductivity [m/s]
              SoilMoistureSat        => noahmp%water%param%SoilMoistureSat        ,& ! in,    saturated value of soil moisture [m3/m3]
              GridTopoIndex          => noahmp%water%param%GridTopoIndex     ,& ! in,    gridcell mean topgraphic index (global mean)
              SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat    ,& ! in,    saturated soil matric potential
              SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB          ,& ! in,    soil B parameter
              SpecYieldGw            => noahmp%water%param%SpecYieldGw       ,& ! in,    specific yield [-], default:0.2
              MicroPoreContent       => noahmp%water%param%MicroPoreContent  ,& ! in,    microprore content (0.0-1.0), default:0.2
              SoilWatConductivitySat => noahmp%water%param%SoilWatConductivitySat ,& ! in,    saturated soil hydraulic conductivity [m/s]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater           ,& ! inout, soil water content [m3/m3]
              WaterTableDepth        => noahmp%water%state%WaterTableDepth   ,& ! inout, water table depth [m]
              WaterStorageAquifer    => noahmp%water%state%WaterStorageAquifer,& ! inout, water storage in aquifer [mm]
              WaterStorageSoilAqf    => noahmp%water%state%WaterStorageSoilAqf,& ! inout, water storage in aquifer + saturated soil [mm]
              RunoffDecayFac         => noahmp%water%param%RunoffDecayFac    ,& ! inout, runoff decay factor (1/m)
              BaseflowCoeff          => noahmp%water%param%BaseflowCoeff     ,& ! inout, baseflow coefficient [mm/s]
              RechargeGw             => noahmp%water%flux%RechargeGw         ,& ! out,   groundwater recharge rate [mm/s]
              DischargeGw            => noahmp%water%flux%DischargeGw         & ! out,   groundwater discharge rate [mm/s]
             )

    allocate(ThicknessSoil(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(DepthSoilMid(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(SoilLiqTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(SoilEffPorosity(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(SoilWatConductTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(SoilMoisture(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    !$acc data create(ThicknessSoil, DepthSoilMid, SoilLiqTmp, SoilEffPorosity, SoilWatConductTmp, SoilMoisture)

   !$acc parallel loop collapse(2) gang vector default(present) &
   !$acc private(LoopInd, IndUnsatSoil, SatDegUnsatSoil, SoilMatPotFrz, AquiferWatConduct) &
   !$acc private(WaterHeadTbl, WaterHead, WaterFillPore, WatConductAcc, SoilMoistureMin, WaterExcessSat)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points


    ! initialization
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
      DepthSoilMid(I,LoopInd,J)      = 0.0
      ThicknessSoil(I,LoopInd,J)     = 0.0
      SoilLiqTmp(I,LoopInd,J)        = 0.0
      SoilEffPorosity(I,LoopInd,J)   = 0.0
      SoilWatConductTmp(I,LoopInd,J) = 0.0
      SoilMoisture(I,LoopInd,J)      = 0.0
    enddo

    DischargeGw(I,J)       = 0.0
    RechargeGw(I,J)        = 0.0

    ! Derive layer-bottom depth in [mm]; KWM:Derive layer thickness in mm
    ThicknessSoil(I,1,J) = -DepthSoilLayer(I,1,J) * 1.0e3
    !$acc loop seq
    do LoopInd = 2, NumSoilLayer
       ThicknessSoil(I,LoopInd,J) = 1.0e3 * (DepthSoilLayer(I,LoopInd-1,J) - DepthSoilLayer(I,LoopInd,J))
    enddo

    ! Derive node (middle) depth in [m]; KWM: Positive number, depth below ground surface in m
    DepthSoilMid(I,1,J) = -DepthSoilLayer(I,1,J) / 2.0
    !$acc loop seq
    do LoopInd = 2, NumSoilLayer
       DepthSoilMid(I,LoopInd,J) = -DepthSoilLayer(I,LoopInd-1,J) + &
                               0.5 * (DepthSoilLayer(I,LoopInd-1,J) - DepthSoilLayer(I,LoopInd,J))
    enddo

    ! Convert volumetric soil moisture to mass
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       SoilMoisture(I,LoopInd,J)   = SoilLiqWater(I,LoopInd,J) + SoilIce(I,LoopInd,J)
       SoilLiqTmp(I,LoopInd,J)        = SoilLiqWater(I,LoopInd,J) * ThicknessSoil(I,LoopInd,J)
       SoilEffPorosity(I,LoopInd,J)   = max(0.01, SoilMoistureSat(I,LoopInd,J)-SoilIce(I,LoopInd,J))
       SoilWatConductTmp(I,LoopInd,J) = 1.0e3 * SoilWatConductivity(I,LoopInd,J)
    enddo

    ! The layer index of the first unsaturated layer (the layer right above the water table)
    IndUnsatSoil = NumSoilLayer
    !$acc loop seq
    do LoopInd = 2, NumSoilLayer
       if ( WaterTableDepth(I,J) <= -DepthSoilLayer(I,LoopInd,J) ) then
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
    RunoffDecayFac(I,J)    = SoilExpCoeffB(I,IndUnsatSoil,J) / 3.0
    BaseflowCoeff(I,J)     = SoilWatConductTmp(I,IndUnsatSoil,J) * 1.0e3 * exp(3.0)  ! [mm/s]
    DischargeGw(I,J)       = (1.0 - SoilImpervFracMax(I,J)) * BaseflowCoeff(I,J) * exp(-GridTopoIndex(I,J)) * &
                        exp(-RunoffDecayFac(I,J) * WaterTableDepth(I,J))

    ! Matric potential at the layer above the water table
    SatDegUnsatSoil   = min(1.0, SoilMoisture(I,IndUnsatSoil,J)/SoilMoistureSat(I,IndUnsatSoil,J))
    SatDegUnsatSoil   = max(SatDegUnsatSoil, real(0.01,kind=8))
    if (SatDegUnsatSoil < 0.01) SatDegUnsatSoil = 0.01
    SoilMatPotFrz     = -SoilMatPotentialSat(I,IndUnsatSoil,J) * 1000.0 * &
                        SatDegUnsatSoil**(-SoilExpCoeffB(I,IndUnsatSoil,J))
    SoilMatPotFrz     = max(-120000.0, MicroPoreContent(I,J)*SoilMatPotFrz)

    ! Recharge rate qin to groundwater
    AquiferWatConduct = 2.0 * (SoilWatConductTmp(I,IndUnsatSoil,J) * SoilWatConductivitySat(I,IndUnsatSoil,J)*1.0e3) / &
                        (SoilWatConductTmp(I,IndUnsatSoil,J) + SoilWatConductivitySat(I,IndUnsatSoil,J)*1.0e3)  ! harmonic average, GY Niu's update 2022
    WaterHeadTbl      = -WaterTableDepth(I,J) * 1.0e3                 !(mm)
    WaterHead         = SoilMatPotFrz - DepthSoilMid(I,IndUnsatSoil,J) * 1.0e3   !(mm)
    RechargeGw(I,J)        = -AquiferWatConduct * (WaterHeadTbl - WaterHead) / &
                        ((WaterTableDepth(I,J)-DepthSoilMid(I,IndUnsatSoil,J)) * 1.0e3)
    RechargeGw(I,J)        = max(-10.0/SoilTimeStep, min(10.0/SoilTimeStep, RechargeGw(I,J)))

    ! Water storage in the aquifer + saturated soil
    WaterStorageSoilAqf(I,J) = WaterStorageSoilAqf(I,J) + (RechargeGw(I,J) - DischargeGw(I,J)) * SoilTimeStep     !(mm)
    if ( IndUnsatSoil == NumSoilLayer ) then
       WaterStorageAquifer(I,J)      = WaterStorageAquifer(I,J) + (RechargeGw(I,J) - DischargeGw(I,J)) * SoilTimeStep     !(mm)
       WaterStorageSoilAqf(I,J)      = WaterStorageAquifer(I,J)
       WaterTableDepth(I,J)          = (-DepthSoilLayer(I,NumSoilLayer,J) + 25.0) - &
                                  WaterStorageAquifer(I,J) / 1000.0 / SpecYieldGw(I,J)      !(m)
       SoilLiqTmp(I,NumSoilLayer,J) = SoilLiqTmp(I,NumSoilLayer,J) - RechargeGw(I,J) * SoilTimeStep        ! [mm]
       SoilLiqTmp(I,NumSoilLayer,J) = SoilLiqTmp(I,NumSoilLayer,J) + max(0.0, (WaterStorageAquifer(I,J)-5000.0))
       WaterStorageAquifer(I,J)      = min(WaterStorageAquifer(I,J), 5000.0)
    else
       if ( IndUnsatSoil == NumSoilLayer-1 ) then
          WaterTableDepth(I,J) = -DepthSoilLayer(I,NumSoilLayer,J) - (WaterStorageSoilAqf(I,J) - SpecYieldGw(I,J)*1000.0*25.0) / &
                                                            (SoilEffPorosity(I,NumSoilLayer,J)) / 1000.0
       else
          WaterFillPore   = 0.0   ! water used to fill soil air pores
          !$acc loop seq
          do LoopInd = IndUnsatSoil+2, NumSoilLayer
             WaterFillPore = WaterFillPore + SoilEffPorosity(I,LoopInd,J) * ThicknessSoil(I,LoopInd,J)
          enddo
          WaterTableDepth(I,J)  = -DepthSoilLayer(I,IndUnsatSoil+1,J) - (WaterStorageSoilAqf(I,J) - SpecYieldGw(I,J)*1000.0*25.0 - &
                                                                WaterFillPore) / (SoilEffPorosity(I,IndUnsatSoil+1,J)) / 1000.0
       endif
       WatConductAcc = 0.0
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          WatConductAcc = WatConductAcc + SoilWatConductTmp(I,LoopInd,J) * ThicknessSoil(I,LoopInd,J)
       enddo
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer           ! Removing subsurface runoff
          SoilLiqTmp(I,LoopInd,J) = SoilLiqTmp(I,LoopInd,J) - DischargeGw(I,J) * SoilTimeStep * &
                                                      SoilWatConductTmp(I,LoopInd,J) * ThicknessSoil(I,LoopInd,J) / WatConductAcc
       enddo
    endif
    WaterTableDepth(I,J) = max(1.5, WaterTableDepth(I,J))

    ! Limit SoilLiqTmp to be greater than or equal to SoilMoistureMin
    ! Get water needed to bring SoilLiqTmp equal SoilMoistureMin from lower layer.
    SoilMoistureMin = 0.01
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer-1
       if ( SoilLiqTmp(I,LoopInd,J) < 0.0 ) then
          WaterExcessSat = SoilMoistureMin - SoilLiqTmp(I,LoopInd,J)
       else
          WaterExcessSat = 0.0
       endif
       SoilLiqTmp(I,LoopInd,J) = SoilLiqTmp(I,LoopInd,J) + WaterExcessSat
       SoilLiqTmp(I,LoopInd+1,J) = SoilLiqTmp(I,LoopInd+1,J) - WaterExcessSat
    enddo
    LoopInd = NumSoilLayer
    if ( SoilLiqTmp(I,LoopInd,J) < SoilMoistureMin ) then
       WaterExcessSat   = SoilMoistureMin - SoilLiqTmp(I,LoopInd,J)
    else
       WaterExcessSat   = 0.0
    endif
    SoilLiqTmp(I,LoopInd,J) = SoilLiqTmp(I,LoopInd,J) + WaterExcessSat
    WaterStorageAquifer(I,J) = WaterStorageAquifer(I,J) - WaterExcessSat
    WaterStorageSoilAqf(I,J) = WaterStorageSoilAqf(I,J) - WaterExcessSat

    ! update soil moisture
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
        SoilLiqWater(I,LoopInd,J) = SoilLiqTmp(I,LoopInd,J) / ThicknessSoil(I,LoopInd,J)
    enddo


      end do
    end do
   !$acc end parallel loop
    !$acc end data
    deallocate(ThicknessSoil, DepthSoilMid, SoilLiqTmp, SoilEffPorosity, SoilWatConductTmp, SoilMoisture)


    end associate

  end subroutine GroundWaterTopModel

end module GroundWaterTopModelMod
