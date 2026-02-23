module TileDrainageHooghoudtMod

!!! Calculate tile drainage discharge [mm] based on Hooghoudt's equation

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use TileDrainageEquiDepthMod, only : TileDrainageEquiDepth
  use WaterTableDepthSearchMod, only : WaterTableDepthSearch
  use WaterTableEquilibriumMod, only : WaterTableEquilibrium

  implicit none

contains

  subroutine TileDrainageHooghoudt(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: TILE_HOOGHOUDT
! Original code: P. Valayamkunnath (NCAR)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                                   ! grid indices
    integer                          :: IndSoil                                ! soil layer loop index
    integer                          :: NumDrain                               ! number of drains
    real(kind=kind_noahmp)           :: ThickSatZoneTot                        ! total thickness of saturated zone
    real(kind=kind_noahmp)           :: LateralFlow                            ! lateral flow
    real(kind=kind_noahmp)           :: DepthToLayerTop                        ! depth to top of the layer
    real(kind=kind_noahmp)           :: WatTblTmp1                             ! temporary water table variable
    real(kind=kind_noahmp)           :: WatTblTmp2                             ! temporary water table variable
    real(kind=kind_noahmp)           :: LateralWatCondAve                      ! average lateral hydruaic conductivity
    real(kind=kind_noahmp)           :: DrainWatHgtAbvImp                      ! Height of water table in the drain Above Impermeable Layer
    real(kind=kind_noahmp)           :: DepthSfcToImp                          ! Effective Depth to impermeable layer from soil surface
    real(kind=kind_noahmp)           :: HgtDrnToWatTbl                         ! Effective Height between water level in drains to water table MiDpoint
    real(kind=kind_noahmp)           :: DrainCoeffTmp                          ! Drainage Coefficient
    real(kind=kind_noahmp)           :: TileDrainTmp                           ! temporary drainage discharge
    real(kind=kind_noahmp)           :: DrainDepthToImpTmp                     ! drain depth to impermeable layer
    real(kind=kind_noahmp)           :: WatExcFieldCapTot                      ! amount of water over field capacity
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: ThickSatZone
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: LateralWatCondTmp
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: WatExcFieldCapTmp
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilLiqWaterAftDrain

! ----------------------------------------------------------------------------

#ifndef WRF_HYDRO
    ! WaterTableDepthSearch has its own internal parallel loop
    call WaterTableDepthSearch(noahmp)
#endif

        associate(                                                                 &
                  NumSoilLayer         => noahmp%config%domain%NumSoilLayer       ,& ! in,    number of soil layers
                  DepthSoilLayer       => noahmp%config%domain%DepthSoilLayer     ,& ! in,    depth [m] of layer-bottom from soil surface
                  SoilTimeStep         => noahmp%config%domain%SoilTimeStep       ,& ! in,    noahmp soil timestep [s]
                  GridSize             => noahmp%config%domain%GridSize           ,& ! in,    noahmp model grid spacing [m]
                  ThicknessSoilLayer   => noahmp%config%domain%ThicknessSoilLayer ,& ! in,    soil layer thickness [m]
                  SoilMoistureFieldCap => noahmp%water%param%SoilMoistureFieldCap ,& ! in,    reference soil moisture (field capacity) [m3/m3]
                  TileDrainCoeff       => noahmp%water%param%TileDrainCoeff  ,& ! in,    drainage coefficent [m/day]
                  DrainDepthToImperv   => noahmp%water%param%DrainDepthToImperv,& ! in,    Actual depth to impermeable layer from surface [m]
                  LateralWatCondFac    => noahmp%water%param%LateralWatCondFac,& ! in,    multiplication factor to determine lateral hydraulic conductivity
                  TileDrainDepth       => noahmp%water%param%TileDrainDepth  ,& ! in,    Depth of drain [m]
                  DrainTubeDist        => noahmp%water%param%DrainTubeDist   ,& ! in,    distance between two drain tubes or tiles [m]
                  DrainTubeRadius      => noahmp%water%param%DrainTubeRadius ,& ! in,    effective radius of drains [m]
                  SoilWatConductivity  => noahmp%water%state%SoilWatConductivity  ,& ! in,    soil hydraulic conductivity [m/s]
                  SoilIce              => noahmp%water%state%SoilIce              ,& ! in,    soil ice content [m3/m3]
                  WaterTableHydro      => noahmp%water%state%WaterTableHydro ,& ! in,    water table depth estimated in WRF-Hydro fine grids [m]
                  SoilLiqWater         => noahmp%water%state%SoilLiqWater         ,& ! inout, soil water content [m3/m3]
                  SoilMoisture         => noahmp%water%state%SoilMoisture         ,& ! inout, total soil moisture [m3/m3]
                  WaterTableDepth      => noahmp%water%state%WaterTableDepth ,& ! inout, water table depth [m]
                  TileDrain            => noahmp%water%flux%TileDrain         & ! inout, tile drainage [mm/s]
                 )

    allocate(ThickSatZone(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(LateralWatCondTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(WatExcFieldCapTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(SoilLiqWaterAftDrain(noahmp%config%domain%ITS:noahmp%config%domain%ITE, 1:NumSoilLayer, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    !$acc data create(ThickSatZone, LateralWatCondTmp, WatExcFieldCapTmp, SoilLiqWaterAftDrain)

    !$acc parallel loop collapse(2) gang vector default(present) &
    !$acc private(IndSoil, NumDrain, ThickSatZoneTot, LateralFlow, DepthToLayerTop) &
    !$acc private(WatTblTmp1, WatTblTmp2, LateralWatCondAve, DrainWatHgtAbvImp) &
    !$acc private(DepthSfcToImp, HgtDrnToWatTbl, DrainCoeffTmp, TileDrainTmp) &
    !$acc private(DrainDepthToImpTmp, WatExcFieldCapTot)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         !cycle condition copied from SoilWaterMainMod before GPU port
         if ( noahmp%water%state%TileDrainFrac (I,J) <= 0.1 ) cycle


        ! initialization
        !$acc loop seq
         do IndSoil = 1, NumSoilLayer
            ThicknessSoilLayer(I,IndSoil,J)  = 0.0
         enddo
        !$acc loop seq
        do IndSoil = 1, NumSoilLayer
           ThickSatZone(I,IndSoil,J)         = 0.0
           LateralWatCondTmp(I,IndSoil,J)    = 0.0
           WatExcFieldCapTmp(I,IndSoil,J)    = 0.0
           SoilLiqWaterAftDrain(I,IndSoil,J) = 0.0
        enddo
        DepthToLayerTop      = 0.0
        LateralFlow          = 0.0
        ThickSatZoneTot      = 0.0
        DrainCoeffTmp        = TileDrainCoeff(I,J) * 1000.0 * SoilTimeStep / (24.0 * 3600.0)  ! m per day to mm per timestep

        ! Thickness of soil layers
        !$acc loop seq
        do IndSoil = 1, NumSoilLayer
           if ( IndSoil == 1 ) then
              ThicknessSoilLayer(I,IndSoil,J) = -1.0 * DepthSoilLayer(I,IndSoil,J)
           else
              ThicknessSoilLayer(I,IndSoil,J) = (DepthSoilLayer(I,IndSoil-1,J) - DepthSoilLayer(I,IndSoil,J))
           endif
        enddo

#ifdef WRF_HYDRO
        ! Depth to water table from WRF-HYDRO, m
        WatTblTmp2 = WaterTableHydro(I,J)
#else
        WatTblTmp2 = WaterTableDepth(I,J)
#endif

        if ( WatTblTmp2 > DrainDepthToImperv(I,J)) WatTblTmp2 = DrainDepthToImperv(I,J)

        ! Depth of saturated zone
        !$acc loop seq
        do IndSoil = 1, NumSoilLayer
           if ( WatTblTmp2 > (-1.0*DepthSoilLayer(I,IndSoil,J)) ) then
              ThickSatZone(I,IndSoil,J) = 0.0
           else
              ThickSatZone(I,IndSoil,J) = (-1.0 * DepthSoilLayer(I,IndSoil,J)) - WatTblTmp2
              WatTblTmp1            = (-1.0 * DepthSoilLayer(I,IndSoil,J)) - DepthToLayerTop
              if ( ThickSatZone(I,IndSoil,J) > WatTblTmp1 ) ThickSatZone(I,IndSoil,J) = WatTblTmp1
           endif
           DepthToLayerTop = -1.0 * DepthSoilLayer(I,IndSoil,J)
        enddo

        ! amount of water over field capacity
        WatExcFieldCapTot = 0.0
        !$acc loop seq
        do IndSoil = 1, NumSoilLayer
           WatExcFieldCapTmp(I,IndSoil,J) = (SoilLiqWater(I,IndSoil,J) - (SoilMoistureFieldCap(I,IndSoil,J)-SoilIce(I,IndSoil,J))) * &
                                        ThicknessSoilLayer(I,IndSoil,J) * 1000.0
           if ( WatExcFieldCapTmp(I,IndSoil,J) < 0.0 ) WatExcFieldCapTmp(I,IndSoil,J) = 0.0
           WatExcFieldCapTot = WatExcFieldCapTot + WatExcFieldCapTmp(I,IndSoil,J)
        enddo

        ! lateral hydraulic conductivity and total lateral flow
        !$acc loop seq
        do IndSoil = 1, NumSoilLayer
           LateralWatCondTmp(I,IndSoil,J) = SoilWatConductivity(I,IndSoil,J) * LateralWatCondFac(I,J) * SoilTimeStep  ! m/s to m/timestep
           LateralFlow                = LateralFlow + (ThickSatZone(I,IndSoil,J) * LateralWatCondTmp(I,IndSoil,J))
           ThickSatZoneTot            = ThickSatZoneTot + ThickSatZone(I,IndSoil,J)
        enddo
        if ( ThickSatZoneTot < 0.001 ) ThickSatZoneTot = 0.001                                               ! unit is m
        if ( LateralFlow < 0.001 )     LateralFlow     = 0.0                                                 ! unit is m
        LateralWatCondAve  = LateralFlow / ThickSatZoneTot                                                   ! lateral hydraulic conductivity per timestep
        DrainDepthToImpTmp = DrainDepthToImperv(I,J) - TileDrainDepth(I,J)

        call TileDrainageEquiDepth(DrainDepthToImpTmp, DrainTubeDist(I,J), DrainTubeRadius(I,J), DrainWatHgtAbvImp)

        DepthSfcToImp  = DrainWatHgtAbvImp + TileDrainDepth(I,J)
        HgtDrnToWatTbl = TileDrainDepth(I,J) - WatTblTmp2
        if ( HgtDrnToWatTbl <= 0.0 ) then
           TileDrain(I,J) = 0.0
        else
           TileDrain(I,J) = ((8.0*LateralWatCondAve*DrainWatHgtAbvImp*HgtDrnToWatTbl) + &
                       (4.0*LateralWatCondAve*HgtDrnToWatTbl*HgtDrnToWatTbl)) / (DrainTubeDist(I,J)*DrainTubeDist(I,J))
        endif
        TileDrain(I,J)    = TileDrain(I,J) * 1000.0                                                                     ! m per timestep to mm/timestep /one tile
        if ( TileDrain(I,J) <= 0.0 ) TileDrain(I,J) = 0.0
        if ( TileDrain(I,J) > DrainCoeffTmp ) TileDrain(I,J) = DrainCoeffTmp
        NumDrain  = int(GridSize / DrainTubeDist(I,J))
        TileDrain(I,J) = TileDrain(I,J) * NumDrain
        if ( TileDrain(I,J) > WatExcFieldCapTot ) TileDrain(I,J) = WatExcFieldCapTot

        ! update soil moisture after drainage: moisture drains from top to bottom
        TileDrainTmp = TileDrain(I,J)
        !$acc loop seq
        do IndSoil = 1, NumSoilLayer
           if ( TileDrainTmp > 0.0) then
              if ( (ThickSatZone(I,IndSoil,J) > 0.0) .and. (WatExcFieldCapTmp(I,IndSoil,J) > 0.0) ) then
                 SoilLiqWaterAftDrain(I,IndSoil,J) = WatExcFieldCapTmp(I,IndSoil,J) - TileDrainTmp                    ! remaining water after tile drain
                 if ( SoilLiqWaterAftDrain(I,IndSoil,J) > 0.0 ) then
                    SoilLiqWater(I,IndSoil,J) = (SoilMoistureFieldCap(I,IndSoil,J) - SoilIce(I,IndSoil,J)) + &
                                                 SoilLiqWaterAftDrain(I,IndSoil,J) / (ThicknessSoilLayer(I,IndSoil,J) * 1000.0)
                    SoilMoisture(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) + SoilIce(I,IndSoil,J)
                    exit
                 else
                    SoilLiqWater(I,IndSoil,J) = SoilMoistureFieldCap(I,IndSoil,J) - SoilIce(I,IndSoil,J)
                    SoilMoisture(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) + SoilIce(I,IndSoil,J)
                    TileDrainTmp              = TileDrainTmp - WatExcFieldCapTmp(I,IndSoil,J)
                 endif
              endif
           endif
        enddo

        TileDrain(I,J) = TileDrain(I,J) / SoilTimeStep            ! mm/s


      end do
    end do
    !$acc end parallel loop
    !$acc end data
    deallocate(ThickSatZone, LateralWatCondTmp, WatExcFieldCapTmp, SoilLiqWaterAftDrain)


        end associate

  end subroutine TileDrainageHooghoudt

end module TileDrainageHooghoudtMod
