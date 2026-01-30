module BalanceErrorCheckMod

!!! Check water and energy balance and report error (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

!!!! Water balance check initialization
  subroutine BalanceWaterInit(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in NOAHMP_SFLX)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J         ! grid indices
    integer                          :: LoopInd      ! loop index

! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                       &
                  NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,  number of soil layers
                  SurfaceType            => noahmp%config%domain%SurfaceType(I,J)       ,& ! in,  surface type 1-soil; 2-lake
                  ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,  thickness of snow/soil layers [m]
                  CanopyLiqWater         => noahmp%water%state%CanopyLiqWater(I,J)      ,& ! in,  canopy intercepted liquid water [mm]
                  CanopyIce              => noahmp%water%state%CanopyIce(I,J)           ,& ! in,  canopy intercepted ice [mm]
                  SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv(I,J)      ,& ! in,  snow water equivalent [mm]
                  SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! in,  total soil moisture [m3/m3]
                  WaterStorageAquifer    => noahmp%water%state%WaterStorageAquifer(I,J) ,& ! in,  water storage in aquifer [mm]
                  WaterStorageWetland    => noahmp%water%state%WaterStorageWetland(I,J) ,& ! in,  water storage in wetland [mm]
                  WaterStorageTotBeg     => noahmp%water%state%WaterStorageTotBeg(I,J)   & ! out, total water storage [mm] at the beginning
                 )
! ----------------------------------------------------------------------

    ! compute total water storage before NoahMP processes
    if ( SurfaceType == 1 ) then  ! soil
       WaterStorageTotBeg = CanopyLiqWater + CanopyIce + SnowWaterEquiv + WaterStorageAquifer + WaterStorageWetland
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          WaterStorageTotBeg = WaterStorageTotBeg + SoilMoisture(I,LoopInd,J) * &
                                                     ThicknessSnowSoilLayer(I,LoopInd,J) * 1000.0
       enddo
    endif

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine BalanceWaterInit


!!!! Water balance check and report error
  subroutine BalanceWaterCheck(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ERROR
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! Note: write() statements disabled on GPU - diagnostics should be post-processed
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J         ! grid indices
    integer                          :: LoopInd      ! loop index

! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                        &
                  NumSoilLayer            => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
                  SurfaceType             => noahmp%config%domain%SurfaceType(I,J)       ,& ! in,    surface type 1-soil; 2-lake
                  ThicknessSnowSoilLayer  => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
                  MainTimeStep            => noahmp%config%domain%MainTimeStep           ,& ! in,    main noahmp timestep [s]
                  FlagCropland            => noahmp%config%domain%FlagCropland(I,J)      ,& ! in,    flag to identify croplands
                  FlagSoilProcess         => noahmp%config%domain%FlagSoilProcess        ,& ! in,    flag to calculate soil process
                  IrriFracThreshold       => noahmp%water%param%IrriFracThreshold(I,J)        ,& ! in,    irrigation fraction parameter
                  IrrigationFracGrid      => noahmp%water%state%IrrigationFracGrid(I,J)  ,& ! in,    total input irrigation fraction
                  WaterTableDepth         => noahmp%water%state%WaterTableDepth(I,J)     ,& ! in,    water table depth [m]
                  CanopyLiqWater          => noahmp%water%state%CanopyLiqWater(I,J)      ,& ! in,    canopy intercepted liquid water [mm]
                  CanopyIce               => noahmp%water%state%CanopyIce(I,J)           ,& ! in,    canopy intercepted ice [mm]
                  SnowWaterEquiv          => noahmp%water%state%SnowWaterEquiv(I,J)      ,& ! in,    snow water equivalent [mm]
                  SoilMoisture            => noahmp%water%state%SoilMoisture             ,& ! in,    total soil moisture [m3/m3]
                  WaterStorageAquifer     => noahmp%water%state%WaterStorageAquifer(I,J) ,& ! in,    water storage in aquifer [mm]
                  WaterStorageWetland     => noahmp%water%state%WaterStorageWetland(I,J) ,& ! in,    water storage in wetland [mm]
                  WaterStorageTotBeg      => noahmp%water%state%WaterStorageTotBeg(I,J)  ,& ! in,    total water storage [mm] at the beginning
                  PrecipTotRefHeight      => noahmp%water%flux%PrecipTotRefHeight(I,J)   ,& ! in,    total precipitation [mm/s] at reference height
                  EvapCanopyNet           => noahmp%water%flux%EvapCanopyNet(I,J)        ,& ! in,    evaporation of intercepted water [mm/s]
                  Transpiration           => noahmp%water%flux%Transpiration(I,J)        ,& ! in,    transpiration rate [mm/s]
                  EvapGroundNet           => noahmp%water%flux%EvapGroundNet(I,J)        ,& ! in,    net ground (soil/snow) evaporation [mm/s]
                  RunoffSurface           => noahmp%water%flux%RunoffSurface(I,J)        ,& ! in,    surface runoff [mm/dt_soil] per soil timestep
                  RunoffSubsurface        => noahmp%water%flux%RunoffSubsurface(I,J)     ,& ! in,    subsurface runoff [mm/dt_soil] per soil timestep
                  TileDrain               => noahmp%water%flux%TileDrain(I,J)            ,& ! in,    tile drainage [mm/dt_soil] per soil timestep
                  IrrigationRateSprinkler => noahmp%water%flux%IrrigationRateSprinkler(I,J) ,& ! in,    rate of irrigation by sprinkler [m/timestep]
                  IrrigationRateMicro     => noahmp%water%flux%IrrigationRateMicro(I,J)  ,& ! in,    micro irrigation water rate [m/timestep]
                  IrrigationRateFlood     => noahmp%water%flux%IrrigationRateFlood(I,J)  ,& ! in,    flood irrigation water rate [m/timestep]
                  SfcWaterTotChgAcc       => noahmp%water%flux%SfcWaterTotChgAcc(I,J)    ,& ! inout, accumulated snow,soil,canopy water change per soil timestep [mm]
                  PrecipTotAcc            => noahmp%water%flux%PrecipTotAcc(I,J)         ,& ! inout, accumulated precipitation per soil timestep [mm]
                  EvapCanopyNetAcc        => noahmp%water%flux%EvapCanopyNetAcc(I,J)     ,& ! inout, accumulated net canopy evaporation per soil timestep [mm]
                  TranspirationAcc        => noahmp%water%flux%TranspirationAcc(I,J)     ,& ! inout, accumulated transpiration per soil timestep [mm]
                  EvapGroundNetAcc        => noahmp%water%flux%EvapGroundNetAcc(I,J)     ,& ! inout, accumulated net ground evaporation per soil timestep [mm]
                  WaterStorageTotEnd      => noahmp%water%state%WaterStorageTotEnd(I,J)  ,& ! out,   total water storage [mm] at the end
                  WaterBalanceError       => noahmp%water%state%WaterBalanceError(I,J)    & ! out,   water balance error [mm] per time step
                 )
! ----------------------------------------------------------------------

    ! before water balance check, add irrigation water to precipitation
    if ( (FlagCropland .eqv. .true.) .and. (IrrigationFracGrid >= IrriFracThreshold) ) then
       PrecipTotRefHeight = PrecipTotRefHeight + IrrigationRateSprinkler * 1000.0 / MainTimeStep  ! irrigation
    endif

    ! only water balance check for every soil timestep
    ! Error in water balance should be < 0.1 mm
    if ( SurfaceType == 1 ) then   ! soil
       WaterStorageTotEnd = CanopyLiqWater + CanopyIce + SnowWaterEquiv + WaterStorageAquifer + WaterStorageWetland
      !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          WaterStorageTotEnd = WaterStorageTotEnd + SoilMoisture(I,LoopInd,J) * &
                                                     ThicknessSnowSoilLayer(I,LoopInd,J) * 1000.0
       enddo
       ! accumualted water change (only for canopy and snow during non-soil timestep)
       SfcWaterTotChgAcc = SfcWaterTotChgAcc + (WaterStorageTotEnd - WaterStorageTotBeg)  ! snow, canopy, and soil water change
       PrecipTotAcc      = PrecipTotAcc      + PrecipTotRefHeight * MainTimeStep          ! accumulated precip 
       EvapCanopyNetAcc  = EvapCanopyNetAcc  + EvapCanopyNet      * MainTimeStep          ! accumulated canopy evapo
       TranspirationAcc  = TranspirationAcc  + Transpiration      * MainTimeStep          ! accumulated transpiration
       EvapGroundNetAcc  = EvapGroundNetAcc  + EvapGroundNet      * MainTimeStep          ! accumulated soil evapo

       ! check water balance at soil timestep
       if ( FlagSoilProcess .eqv. .true. ) then
          WaterBalanceError = SfcWaterTotChgAcc - (PrecipTotAcc + IrrigationRateMicro*1000.0 + IrrigationRateFlood*1000.0 - &
                              EvapCanopyNetAcc - TranspirationAcc - EvapGroundNetAcc - RunoffSurface - RunoffSubsurface -   &
                              TileDrain)
#if !defined(WRF_HYDRO) && !defined(_OPENACC)
          ! Note: Diagnostic writes disabled for GPU execution
          if ( abs(WaterBalanceError) > 0.1 ) then
             if ( WaterBalanceError > 0 ) then
                write(*,*) "The model is gaining water (WaterBalanceError is positive)"
             else
                write(*,*) "The model is losing water (WaterBalanceError is negative)"
             endif
             write(*,*) 'WaterBalanceError = ',WaterBalanceError, "kg m{-2} timestep{-1}"
             write(*, &
                  '("  GridIndexI  GridIndexJ  SfcWaterTotChgAcc  PrecipTotRefHeightAcc  IrrigationRateMicro       &
                       IrrigationRateFlood  EvapCanopyNetAcc  EvapGroundNetAcc  TranspirationAcc  RunoffSurface    &
                       RunoffSubsurface  WaterTableDepth  TileDrain WaterStorageWetland  ")')
             write(*,'(i6,i6,f10.3,11f10.5)') I, J, SfcWaterTotChgAcc, PrecipTotAcc,                               &
                                              IrrigationRateMicro*1000.0, IrrigationRateFlood*1000.0,              &
                                              EvapCanopyNetAcc, EvapGroundNetAcc, TranspirationAcc, RunoffSurface, &
                                              RunoffSubsurface, WaterTableDepth, TileDrain, WaterStorageWetland
             stop "Error: Water budget problem in NoahMP LSM"
          endif
#endif
       endif ! FlagSoilProcess

    else ! water point
       WaterBalanceError = 0.0
    endif

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine BalanceWaterCheck


!!!! Energy balance check and error report
  subroutine BalanceEnergyCheck(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ERROR
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! Note: write() statements disabled on GPU - diagnostics should be post-processed
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J         ! grid indices

! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                 &
                  OptSnowAlbedo        => noahmp%config%nmlist%OptSnowAlbedo      ,& ! in,  options for ground snow surface albedo
                  RadSwDownRefHeight   => noahmp%forcing%RadSwDownRefHeight(I,J)  ,& ! in,  downward shortwave radiation [W/m2] at reference height
                  VegFrac              => noahmp%energy%state%VegFrac(I,J)        ,& ! in,  greeness vegetation fraction
                  RadSwAbsSfc          => noahmp%energy%flux%RadSwAbsSfc(I,J)     ,& ! in,  total absorbed solar radiation [W/m2]
                  RadSwAbsSnowSoilLayer=> noahmp%energy%flux%RadSwAbsSnowSoilLayer,& ! in,  total absorbed solar radiation by snow/soil for each layer [W/m2]
                  RadSwReflSfc         => noahmp%energy%flux%RadSwReflSfc(I,J)    ,& ! in,  total reflected solar radiation [W/m2]
                  RadSwReflVeg         => noahmp%energy%flux%RadSwReflVeg(I,J)    ,& ! in,  reflected solar radiation by vegetation [W/m2]
                  RadSwReflGrd         => noahmp%energy%flux%RadSwReflGrd(I,J)    ,& ! in,  reflected solar radiation by ground [W/m2]
                  RadLwNetSfc          => noahmp%energy%flux%RadLwNetSfc(I,J)     ,& ! in,  total net longwave rad [W/m2] (+ to atm)
                  HeatSensibleSfc      => noahmp%energy%flux%HeatSensibleSfc(I,J) ,& ! in,  total sensible heat [W/m2] (+ to atm)
                  HeatLatentCanopy     => noahmp%energy%flux%HeatLatentCanopy(I,J),& ! in,  canopy latent heat flux [W/m2] (+ to atm)
                  HeatLatentGrd        => noahmp%energy%flux%HeatLatentGrd(I,J)   ,& ! in,  total ground latent heat [W/m2] (+ to atm)
                  HeatLatentTransp     => noahmp%energy%flux%HeatLatentTransp(I,J),& ! in,  latent heat flux from transpiration [W/m2] (+ to atm)
                  HeatGroundTot        => noahmp%energy%flux%HeatGroundTot(I,J)   ,& ! in,  total ground heat flux [W/m2] (+ to soil/snow)
                  SnowCoverFrac        => noahmp%water%state%SnowCoverFrac(I,J)   ,& ! in,  snow cover fraction
                  RadSwAbsVeg          => noahmp%energy%flux%RadSwAbsVeg(I,J)     ,& ! in,  solar radiation absorbed by vegetation [W/m2]
                  RadSwAbsGrd          => noahmp%energy%flux%RadSwAbsGrd(I,J)     ,& ! in,  solar radiation absorbed by ground [W/m2]
                  HeatPrecipAdvSfc     => noahmp%energy%flux%HeatPrecipAdvSfc(I,J),& ! in,  precipitation advected heat - total [W/m2]
                  HeatPrecipAdvBareGrd => noahmp%energy%flux%HeatPrecipAdvBareGrd(I,J) ,& ! in,  precipitation advected heat - bare ground net [W/m2]
                  HeatPrecipAdvVegGrd  => noahmp%energy%flux%HeatPrecipAdvVegGrd(I,J)  ,& ! in,  precipitation advected heat - under canopy net [W/m2]
                  HeatPrecipAdvCanopy  => noahmp%energy%flux%HeatPrecipAdvCanopy(I,J)  ,& ! in,  precipitation advected heat - vegetation net [W/m2]
                  HeatLatentIrriEvap   => noahmp%energy%flux%HeatLatentIrriEvap(I,J)   ,& ! in,  latent heating due to sprinkler evaporation [W/m2]
                  HeatCanStorageChg    => noahmp%energy%flux%HeatCanStorageChg(I,J)    ,& ! in,  canopy heat storage change [W/m2]
                  RadSwPenetrateGrd    => noahmp%energy%flux%RadSwPenetrateGrd    ,& ! in,  light penetrating through soil and snowpack [W/m2]
                  EnergyBalanceError   => noahmp%energy%state%EnergyBalanceError(I,J)  ,& ! out, error in surface energy balance [W/m2]
                  RadSwBalanceError    => noahmp%energy%state%RadSwBalanceError(I,J)    & ! out, error in shortwave radiation balance [W/m2]
                 )
! ----------------------------------------------------------------------

    ! error in shortwave radiation balance should be <0.01 W/m2
    RadSwBalanceError = RadSwDownRefHeight - (RadSwAbsSfc + RadSwReflSfc)

#ifndef _OPENACC
    ! Note: Diagnostic writes disabled for GPU execution
    if ( abs(RadSwBalanceError) > 0.01 ) then
       write(*,*) "GridIndexI, GridIndexJ              = ", I, J
       write(*,*) "RadSwBalanceError                   = ", RadSwBalanceError
       write(*,*) "VEGETATION ---------"
       write(*,*) "RadSwDownRefHeight * VegFrac        = ", RadSwDownRefHeight*VegFrac
       write(*,*) "VegFrac*RadSwAbsVeg + RadSwAbsGrd   = ", VegFrac*RadSwAbsVeg+RadSwAbsGrd
       write(*,*) "VegFrac*RadSwReflVeg + RadSwReflGrd = ", VegFrac*RadSwReflVeg+RadSwReflGrd
       write(*,*) "GROUND -------"
       write(*,*) "(1 - VegFrac) * RadSwDownRefHeight  = ", (1.0-VegFrac)*RadSwDownRefHeight
       write(*,*) "(1 - VegFrac) * RadSwAbsGrd         = ", (1.0-VegFrac)*RadSwAbsGrd
       write(*,*) "(1 - VegFrac) * RadSwReflGrd        = ", (1.0-VegFrac)*RadSwReflGrd
       write(*,*) "RadSwReflVeg                        = ", RadSwReflVeg
       write(*,*) "RadSwReflGrd                        = ", RadSwReflGrd
       write(*,*) "RadSwReflSfc                        = ", RadSwReflSfc
       write(*,*) "RadSwAbsVeg                         = ", RadSwAbsVeg
       write(*,*) "RadSwAbsGrd                         = ", RadSwAbsGrd
       write(*,*) "RadSwAbsSfc                         = ", RadSwAbsSfc
       stop "Error: Solar radiation budget problem in NoahMP LSM"
    endif

    ! SNICAR
    if ( OptSnowAlbedo == 3 ) then
       if ( abs(RadSwAbsGrd-sum(RadSwAbsSnowSoilLayer(I,:,J)))>0.001 ) then ! original check is 0.0001, precision issue
          write(*,*) "RadSwAbsGrd gridmean                            = ", RadSwAbsGrd
          write(*,*) "sum(RadSwAbsSnowSoilLayer) gridmean             = ", sum(RadSwAbsSnowSoilLayer(I,:,J))
          write(*,*) "RadSwAbsSnowSoilLayer gridmean                  = ", RadSwAbsSnowSoilLayer(I,:,J)
          write(*,*) "RadSwAbsGrd-sum(RadSwAbsSnowSoilLayer) gridmean = ", RadSwAbsGrd-sum(RadSwAbsSnowSoilLayer(I,:,J))
          stop "Error: SNICAR snow albedo radiation budget problem in NoahMP LSM"
       endif
    endif
#endif

    ! error in surface energy balance should be <0.01 W/m2
    EnergyBalanceError = RadSwAbsVeg + RadSwAbsGrd + HeatPrecipAdvSfc -                     &
                        (RadLwNetSfc + HeatSensibleSfc + HeatLatentCanopy + HeatLatentGrd + &
                         HeatLatentTransp + HeatGroundTot + HeatLatentIrriEvap + HeatCanStorageChg)

#ifndef _OPENACC
    ! Note: Diagnostic writes disabled for GPU execution
    if ( abs(EnergyBalanceError) > 0.01 ) then
       write(*,*) 'EnergyBalanceError = ', EnergyBalanceError, ' at GridIndexI,GridIndexJ: ', I, J
       write(*,'(a17,F10.4)' ) "Net solar:        ", RadSwAbsSfc
       write(*,'(a17,F10.4)' ) "Net longwave:     ", RadLwNetSfc
       write(*,'(a17,F10.4)' ) "Total sensible:   ", HeatSensibleSfc
       write(*,'(a17,F10.4)' ) "Canopy evap:      ", HeatLatentCanopy
       write(*,'(a17,F10.4)' ) "Ground evap:      ", HeatLatentGrd
       write(*,'(a17,F10.4)' ) "Transpiration:    ", HeatLatentTransp
       write(*,'(a17,F10.4)' ) "Total ground:     ", HeatGroundTot
       write(*,'(a17,F10.4)' ) "Sprinkler:        ", HeatLatentIrriEvap
       write(*,'(a17,F10.4)' ) "Canopy heat storage change: ", HeatCanStorageChg
       write(*,'(a17,4F10.4)') "Precip advected:  ", HeatPrecipAdvSfc,HeatPrecipAdvCanopy,HeatPrecipAdvVegGrd,HeatPrecipAdvBareGrd
       write(*,'(a17,F10.4)' ) "Veg fraction:     ", VegFrac
       write(*,'(a17,F10.4)' ) "Light through soil/snow layer total:  ", sum(RadSwPenetrateGrd(I,:,J))
       write(*,'(a17,4F10.4)') "Light through soil/snow layer:  ", RadSwPenetrateGrd(I,:,J)
       stop "Error: Energy budget problem in NoahMP LSM"
    endif
#endif

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine BalanceEnergyCheck

end module BalanceErrorCheckMod
