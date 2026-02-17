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
    associate(                                                                              &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,  number of soil layers
              SurfaceType            => noahmp%config%domain%SurfaceType            ,& ! in,  surface type 1-soil; 2-lake
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,  thickness of snow/soil layers [m]
              CanopyLiqWater         => noahmp%water%state%CanopyLiqWater           ,& ! in,  canopy intercepted liquid water [mm]
              CanopyIce              => noahmp%water%state%CanopyIce                ,& ! in,  canopy intercepted ice [mm]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv           ,& ! in,  snow water equivalent [mm]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! in,  total soil moisture [m3/m3]
              WaterStorageAquifer    => noahmp%water%state%WaterStorageAquifer      ,& ! in,  water storage in aquifer [mm]
              WaterStorageWetland    => noahmp%water%state%WaterStorageWetland      ,& ! in,  water storage in wetland [mm]
              WaterStorageTotBeg     => noahmp%water%state%WaterStorageTotBeg        & ! out, total water storage [mm] at the beginning
             )

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! compute total water storage before NoahMP processes
    if ( (noahmp%config%domain%IndicatorIceSfc(I,J) == 0) .and. (SurfaceType(I,J) == 1) ) then  ! soil
       WaterStorageTotBeg(I,J) = CanopyLiqWater(I,J) + CanopyIce(I,J) + SnowWaterEquiv(I,J) + WaterStorageAquifer(I,J) + WaterStorageWetland(I,J)
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          WaterStorageTotBeg(I,J) = WaterStorageTotBeg(I,J) + SoilMoisture(I,LoopInd,J) * &
                                                     ThicknessSnowSoilLayer(I,LoopInd,J) * 1000.0
       enddo
    else if (noahmp%config%domain%IndicatorIceSfc(I,J) == -1) then ! ice point
          ! compute total glacier water storage before NoahMP processes
          ! need more work on including glacier ice mass underneath snow
          WaterStorageTotBeg(I,J) = SnowWaterEquiv(I,J)
    endif


      end do
    end do
    !$acc end parallel loop

    end associate

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
    associate(                                                                              &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,  number of soil layers
              SurfaceType            => noahmp%config%domain%SurfaceType            ,& ! in,  surface type 1-soil; 2-lake
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,  thickness of snow/soil layers [m]
              MainTimeStep           => noahmp%config%domain%MainTimeStep           ,& ! in,  main noahmp timestep [s]
              FlagCropland           => noahmp%config%domain%FlagCropland           ,& ! in,  flag to identify croplands
              FlagSoilProcess        => noahmp%config%domain%FlagSoilProcess        ,& ! in,  flag to calculate soil process
              IrriFracThreshold      => noahmp%water%param%IrriFracThreshold        ,& ! in,  irrigation fraction parameter
              CanopyLiqWater         => noahmp%water%state%CanopyLiqWater           ,& ! in,  canopy intercepted liquid water [mm]
              CanopyIce              => noahmp%water%state%CanopyIce                ,& ! in,  canopy intercepted ice [mm]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv           ,& ! in,  snow water equivalent [mm]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! in,  total soil moisture [m3/m3]
              WaterStorageAquifer    => noahmp%water%state%WaterStorageAquifer      ,& ! in,  water storage in aquifer [mm]
              WaterStorageWetland    => noahmp%water%state%WaterStorageWetland      ,& ! in,  water storage in wetland [mm]
              IrrigationFracGrid     => noahmp%water%state%IrrigationFracGrid      ,& ! in,  total input irrigation fraction
              WaterTableDepth        => noahmp%water%state%WaterTableDepth          ,& ! in,  water table depth [m]
              WaterStorageTotBeg     => noahmp%water%state%WaterStorageTotBeg       ,& ! in,  total water storage [mm] at the beginning
              WaterStorageTotEnd     => noahmp%water%state%WaterStorageTotEnd       ,& ! out, total water storage [mm] at the end
              WaterBalanceError      => noahmp%water%state%WaterBalanceError        ,& ! out, water balance error [mm] per time step
              PrecipTotRefHeight     => noahmp%water%flux%PrecipTotRefHeight        ,& ! in,  total precipitation [mm/s] at reference height
              EvapCanopyNet          => noahmp%water%flux%EvapCanopyNet             ,& ! in,  evaporation of intercepted water [mm/s]
              Transpiration          => noahmp%water%flux%Transpiration             ,& ! in,  transpiration rate [mm/s]
              EvapGroundNet          => noahmp%water%flux%EvapGroundNet             ,& ! in,  net ground (soil/snow) evaporation [mm/s]
              RunoffSurface          => noahmp%water%flux%RunoffSurface             ,& ! in,  surface runoff [mm/dt_soil] per soil timestep
              RunoffSubsurface       => noahmp%water%flux%RunoffSubsurface         ,& ! in,  subsurface runoff [mm/dt_soil] per soil timestep
              TileDrain              => noahmp%water%flux%TileDrain                 ,& ! in,  tile drainage [mm/dt_soil] per soil timestep
              IrrigationRateSprinkler => noahmp%water%flux%IrrigationRateSprinkler ,& ! in,  rate of irrigation by sprinkler [m/timestep]
              IrrigationRateMicro    => noahmp%water%flux%IrrigationRateMicro      ,& ! in,  micro irrigation water rate [m/timestep]
              IrrigationRateFlood    => noahmp%water%flux%IrrigationRateFlood      ,& ! in,  flood irrigation water rate [m/timestep]
              SfcWaterTotChgAcc      => noahmp%water%flux%SfcWaterTotChgAcc        ,& ! inout, accumulated snow,soil,canopy water change per soil timestep [mm]
              PrecipTotAcc           => noahmp%water%flux%PrecipTotAcc             ,& ! inout, accumulated precipitation per soil timestep [mm]
              EvapCanopyNetAcc       => noahmp%water%flux%EvapCanopyNetAcc         ,& ! inout, accumulated net canopy evaporation per soil timestep [mm]
              TranspirationAcc       => noahmp%water%flux%TranspirationAcc         ,& ! inout, accumulated transpiration per soil timestep [mm]
              EvapGroundNetAcc       => noahmp%water%flux%EvapGroundNetAcc          & ! inout, accumulated net ground evaporation per soil timestep [mm]
             )

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! before water balance check, add irrigation water to precipitation
    if ( (FlagCropland(I,J) .eqv. .true.) .and. (IrrigationFracGrid(I,J) >= IrriFracThreshold(I,J)) ) then
       PrecipTotRefHeight(I,J) = PrecipTotRefHeight(I,J) + IrrigationRateSprinkler(I,J) * 1000.0 / MainTimeStep  ! irrigation
    endif

    ! only water balance check for every soil timestep
    ! Error in water balance should be < 0.1 mm
    if ( (noahmp%config%domain%IndicatorIceSfc(I,J) == 0) .and. (SurfaceType(I,J) == 1) ) then   ! soil
       WaterStorageTotEnd(I,J) = CanopyLiqWater(I,J) + CanopyIce(I,J) + SnowWaterEquiv(I,J) + WaterStorageAquifer(I,J) + WaterStorageWetland(I,J)
      !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          WaterStorageTotEnd(I,J) = WaterStorageTotEnd(I,J) + SoilMoisture(I,LoopInd,J) * &
                                                     ThicknessSnowSoilLayer(I,LoopInd,J) * 1000.0
       enddo
       ! accumualted water change (only for canopy and snow during non-soil timestep)
       SfcWaterTotChgAcc(I,J) = SfcWaterTotChgAcc(I,J) + (WaterStorageTotEnd(I,J) - WaterStorageTotBeg(I,J))  ! snow, canopy, and soil water change
       PrecipTotAcc(I,J)      = PrecipTotAcc(I,J)      + PrecipTotRefHeight(I,J) * MainTimeStep          ! accumulated precip 
       EvapCanopyNetAcc(I,J)  = EvapCanopyNetAcc(I,J)  + EvapCanopyNet(I,J)      * MainTimeStep          ! accumulated canopy evapo
       TranspirationAcc(I,J)  = TranspirationAcc(I,J)  + Transpiration(I,J)      * MainTimeStep          ! accumulated transpiration
       EvapGroundNetAcc(I,J)  = EvapGroundNetAcc(I,J)  + EvapGroundNet(I,J)      * MainTimeStep          ! accumulated soil evapo

       ! check water balance at soil timestep
       if ( FlagSoilProcess .eqv. .true. ) then
          WaterBalanceError(I,J) = SfcWaterTotChgAcc(I,J) - (PrecipTotAcc(I,J) + IrrigationRateMicro(I,J)*1000.0 + IrrigationRateFlood(I,J)*1000.0 - &
                              EvapCanopyNetAcc(I,J) - TranspirationAcc(I,J) - EvapGroundNetAcc(I,J) - RunoffSurface(I,J) - RunoffSubsurface(I,J) -   &
                              TileDrain(I,J))
#if !defined(WRF_HYDRO) && !defined(_OPENACC)
          ! Note: Diagnostic writes disabled for GPU execution
          if ( abs(WaterBalanceError(I,J)) > 0.1 ) then
             if ( WaterBalanceError(I,J) > 0 ) then
                write(*,*) "The model is gaining water (WaterBalanceError(I,J) is positive)"
             else
                write(*,*) "The model is losing water (WaterBalanceError(I,J) is negative)"
             endif
             write(*,*) 'WaterBalanceError(I,J) = ',WaterBalanceError(I,J), "kg m{-2} timestep{-1}"
             write(*, &
                  '("  GridIndexI  GridIndexJ  SfcWaterTotChgAcc(I,J)  PrecipTotRefHeightAcc  IrrigationRateMicro(I,J)       &
                       IrrigationRateFlood(I,J)  EvapCanopyNetAcc(I,J)  EvapGroundNetAcc(I,J)  TranspirationAcc(I,J)  RunoffSurface(I,J)    &
                       RunoffSubsurface(I,J)  WaterTableDepth(I,J)  TileDrain(I,J) WaterStorageWetland(I,J)  ")')
             write(*,'(i6,i6,f10.3,11f10.5)') I, J, SfcWaterTotChgAcc(I,J), PrecipTotAcc(I,J),                               &
                                              IrrigationRateMicro(I,J)*1000.0, IrrigationRateFlood(I,J)*1000.0,              &
                                              EvapCanopyNetAcc(I,J), EvapGroundNetAcc(I,J), TranspirationAcc(I,J), RunoffSurface(I,J), &
                                              RunoffSubsurface(I,J), WaterTableDepth(I,J), TileDrain(I,J), WaterStorageWetland(I,J)
             stop "Error: Water budget problem in NoahMP LSM"
          endif
#endif
       endif ! FlagSoilProcess

    else if (noahmp%config%domain%IndicatorIceSfc(I,J) == -1) then ! ice point
      ! Error in water balance should be < 0.1 mm
      ! compute total glacier water storage before NoahMP processes
      ! need more work on including glacier ice mass underneath snow
      WaterStorageTotEnd(I,J) = SnowWaterEquiv(I,J)
      WaterBalanceError(I,J)  = WaterStorageTotEnd(I,J) - WaterStorageTotBeg(I,J) - &
                           (PrecipTotRefHeight(I,J) - EvapGroundNet(I,J) - RunoffSurface(I,J) - RunoffSubsurface(I,J)) * MainTimeStep

#if !defined(WRF_HYDRO) && !defined(_OPENACC)
      if ( abs(WaterBalanceError(I,J)) > 0.1 ) then
         if ( WaterBalanceError(I,J) > 0) then
            write(*,*) "The model is gaining water (WaterBalanceError(I,J) is positive)"
         else
            write(*,*) "The model is losing water (WaterBalanceError(I,J) is negative)"
         endif
         write(*,*) "WaterBalanceError(I,J) = ",WaterBalanceError(I,J), "kg m{-2} timestep{-1}"
         write(*, &
            '("  GridIndexI   GridIndexJ     WaterStorageTotEnd(I,J)  WaterStorageTotBeg(I,J)  PrecipTotRefHeight(I,J)  &
                  EvapGroundNet(I,J)  RunoffSurface(I,J)  RunoffSubsurface(I,J)")')
         write(*,'(i6,1x,i6,1x,2f15.3,9f11.5)') I, J, WaterStorageTotEnd(I,J), WaterStorageTotBeg(I,J), &
                                                PrecipTotRefHeight(I,J)*MainTimeStep, EvapGroundNet(I,J)*MainTimeStep,    &
                                                RunoffSurface(I,J)*MainTimeStep, RunoffSubsurface(I,J)*MainTimeStep
         stop "Error: Water budget problem in NoahMP LSM"
      endif
#endif
    else if (SurfaceType(I,J) == 2) then ! water point
         WaterBalanceError(I,J) = 0.0
    endif


      end do
    end do
    !$acc end parallel loop

    end associate

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
    associate(                                                                              &
              OptSnowAlbedo          => noahmp%config%nmlist%OptSnowAlbedo          ,& ! in,  options for ground snow surface albedo
              RadSwDownRefHeight     => noahmp%forcing%RadSwDownRefHeight           ,& ! in,  downward shortwave radiation [W/m2] at reference height
              VegFrac                => noahmp%energy%state%VegFrac                 ,& ! in,  greeness vegetation fraction
              RadSwAbsSfc            => noahmp%energy%flux%RadSwAbsSfc              ,& ! in,  total absorbed solar radiation [W/m2]
              RadSwAbsSnowSoilLayer  => noahmp%energy%flux%RadSwAbsSnowSoilLayer   ,& ! in,  total absorbed solar radiation by snow/soil for each layer [W/m2]
              RadSwReflSfc           => noahmp%energy%flux%RadSwReflSfc             ,& ! in,  total reflected solar radiation [W/m2]
              RadSwReflVeg           => noahmp%energy%flux%RadSwReflVeg             ,& ! in,  reflected solar radiation by vegetation [W/m2]
              RadSwReflGrd           => noahmp%energy%flux%RadSwReflGrd             ,& ! in,  reflected solar radiation by ground [W/m2]
              RadLwNetSfc            => noahmp%energy%flux%RadLwNetSfc              ,& ! in,  total net longwave rad [W/m2] (+ to atm)
              HeatSensibleSfc        => noahmp%energy%flux%HeatSensibleSfc          ,& ! in,  total sensible heat [W/m2] (+ to atm)
              HeatLatentCanopy       => noahmp%energy%flux%HeatLatentCanopy         ,& ! in,  canopy latent heat flux [W/m2] (+ to atm)
              HeatLatentGrd          => noahmp%energy%flux%HeatLatentGrd            ,& ! in,  total ground latent heat [W/m2] (+ to atm)
              HeatLatentTransp       => noahmp%energy%flux%HeatLatentTransp         ,& ! in,  latent heat flux from transpiration [W/m2] (+ to atm)
              HeatGroundTot          => noahmp%energy%flux%HeatGroundTot            ,& ! in,  total ground heat flux [W/m2] (+ to soil/snow)
              SnowCoverFrac          => noahmp%water%state%SnowCoverFrac            ,& ! in,  snow cover fraction
              RadSwAbsVeg            => noahmp%energy%flux%RadSwAbsVeg              ,& ! in,  solar radiation absorbed by vegetation [W/m2]
              RadSwAbsGrd            => noahmp%energy%flux%RadSwAbsGrd              ,& ! in,  solar radiation absorbed by ground [W/m2]
              HeatPrecipAdvSfc       => noahmp%energy%flux%HeatPrecipAdvSfc         ,& ! in,  precipitation advected heat - total [W/m2]
              HeatPrecipAdvBareGrd   => noahmp%energy%flux%HeatPrecipAdvBareGrd     ,& ! in,  precipitation advected heat - bare ground net [W/m2]
              HeatPrecipAdvVegGrd    => noahmp%energy%flux%HeatPrecipAdvVegGrd      ,& ! in,  precipitation advected heat - under canopy net [W/m2]
              HeatPrecipAdvCanopy    => noahmp%energy%flux%HeatPrecipAdvCanopy      ,& ! in,  precipitation advected heat - vegetation net [W/m2]
              HeatLatentIrriEvap     => noahmp%energy%flux%HeatLatentIrriEvap       ,& ! in,  latent heating due to sprinkler evaporation [W/m2]
              HeatCanStorageChg      => noahmp%energy%flux%HeatCanStorageChg        ,& ! in,  canopy heat storage change [W/m2]
              RadSwPenetrateGrd      => noahmp%energy%flux%RadSwPenetrateGrd        ,& ! in,  light penetrating through soil and snowpack [W/m2]
              EnergyBalanceError     => noahmp%energy%state%EnergyBalanceError      ,& ! out, error in surface energy balance [W/m2]
              RadSwBalanceError      => noahmp%energy%state%RadSwBalanceError        & ! out, error in shortwave radiation balance [W/m2]
             )

    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! error in shortwave radiation balance should be <0.01 W/m2
    RadSwBalanceError(I,J) = RadSwDownRefHeight(I,J) - (RadSwAbsSfc(I,J) + RadSwReflSfc(I,J))

  if (noahmp%config%domain%IndicatorIceSfc(I,J) == 0) then
#ifndef _OPENACC
    ! Note: Diagnostic writes disabled for GPU execution
    if ( abs(RadSwBalanceError(I,J)) > 0.01 ) then
       write(*,*) "GridIndexI, GridIndexJ              = ", I, J
       write(*,*) "RadSwBalanceError(I,J)                   = ", RadSwBalanceError(I,J)
       write(*,*) "VEGETATION ---------"
       write(*,*) "RadSwDownRefHeight(I,J) * VegFrac(I,J)        = ", RadSwDownRefHeight(I,J)*VegFrac(I,J)
       write(*,*) "VegFrac(I,J)*RadSwAbsVeg(I,J) + RadSwAbsGrd(I,J)   = ", VegFrac(I,J)*RadSwAbsVeg(I,J)+RadSwAbsGrd(I,J)
       write(*,*) "VegFrac(I,J)*RadSwReflVeg(I,J) + RadSwReflGrd(I,J) = ", VegFrac(I,J)*RadSwReflVeg(I,J)+RadSwReflGrd(I,J)
       write(*,*) "GROUND -------"
       write(*,*) "(1 - VegFrac(I,J)) * RadSwDownRefHeight(I,J)  = ", (1.0-VegFrac(I,J))*RadSwDownRefHeight(I,J)
       write(*,*) "(1 - VegFrac(I,J)) * RadSwAbsGrd(I,J)         = ", (1.0-VegFrac(I,J))*RadSwAbsGrd(I,J)
       write(*,*) "(1 - VegFrac(I,J)) * RadSwReflGrd(I,J)        = ", (1.0-VegFrac(I,J))*RadSwReflGrd(I,J)
       write(*,*) "RadSwReflVeg(I,J)                        = ", RadSwReflVeg(I,J)
       write(*,*) "RadSwReflGrd(I,J)                        = ", RadSwReflGrd(I,J)
       write(*,*) "RadSwReflSfc(I,J)                        = ", RadSwReflSfc(I,J)
       write(*,*) "RadSwAbsVeg(I,J)                         = ", RadSwAbsVeg(I,J)
       write(*,*) "RadSwAbsGrd(I,J)                         = ", RadSwAbsGrd(I,J)
       write(*,*) "RadSwAbsSfc(I,J)                         = ", RadSwAbsSfc(I,J)
       stop "Error: Solar radiation budget problem in NoahMP LSM"
    endif

    ! SNICAR
    if ( OptSnowAlbedo == 3 ) then
       if ( abs(RadSwAbsGrd(I,J)-sum(RadSwAbsSnowSoilLayer(I,:,J)))>0.001 ) then ! original check is 0.0001, precision issue
          write(*,*) "RadSwAbsGrd(I,J) gridmean                            = ", RadSwAbsGrd(I,J)
          write(*,*) "sum(RadSwAbsSnowSoilLayer) gridmean             = ", sum(RadSwAbsSnowSoilLayer(I,:,J))
          write(*,*) "RadSwAbsSnowSoilLayer gridmean                  = ", RadSwAbsSnowSoilLayer(I,:,J)
          write(*,*) "RadSwAbsGrd(I,J)-sum(RadSwAbsSnowSoilLayer) gridmean = ", RadSwAbsGrd(I,J)-sum(RadSwAbsSnowSoilLayer(I,:,J))
          stop "Error: SNICAR snow albedo radiation budget problem in NoahMP LSM"
       endif
    endif
#endif

    ! error in surface energy balance should be <0.01 W/m2
    EnergyBalanceError(I,J) = RadSwAbsVeg(I,J) + RadSwAbsGrd(I,J) + HeatPrecipAdvSfc(I,J) -                     &
                        (RadLwNetSfc(I,J) + HeatSensibleSfc(I,J) + HeatLatentCanopy(I,J) + HeatLatentGrd(I,J) + &
                         HeatLatentTransp(I,J) + HeatGroundTot(I,J) + HeatLatentIrriEvap(I,J) + HeatCanStorageChg(I,J))

#ifndef _OPENACC
    ! Note: Diagnostic writes disabled for GPU execution
    if ( abs(EnergyBalanceError(I,J)) > 0.01 ) then
       write(*,*) 'EnergyBalanceError(I,J) = ', EnergyBalanceError(I,J), ' at GridIndexI,GridIndexJ: ', I, J
       write(*,'(a17,F10.4)' ) "Net solar:        ", RadSwAbsSfc(I,J)
       write(*,'(a17,F10.4)' ) "Net longwave:     ", RadLwNetSfc(I,J)
       write(*,'(a17,F10.4)' ) "Total sensible:   ", HeatSensibleSfc(I,J)
       write(*,'(a17,F10.4)' ) "Canopy evap:      ", HeatLatentCanopy(I,J)
       write(*,'(a17,F10.4)' ) "Ground evap:      ", HeatLatentGrd(I,J)
       write(*,'(a17,F10.4)' ) "Transpiration(I,J):    ", HeatLatentTransp(I,J)
       write(*,'(a17,F10.4)' ) "Total ground:     ", HeatGroundTot(I,J)
       write(*,'(a17,F10.4)' ) "Sprinkler:        ", HeatLatentIrriEvap(I,J)
       write(*,'(a17,F10.4)' ) "Canopy heat storage change: ", HeatCanStorageChg(I,J)
       write(*,'(a17,4F10.4)') "Precip advected:  ", HeatPrecipAdvSfc(I,J),HeatPrecipAdvCanopy(I,J),HeatPrecipAdvVegGrd(I,J),HeatPrecipAdvBareGrd(I,J)
       write(*,'(a17,F10.4)' ) "Veg fraction:     ", VegFrac(I,J)
       write(*,'(a17,F10.4)' ) "Light through soil/snow layer total:  ", sum(RadSwPenetrateGrd(I,:,J))
       write(*,'(a17,4F10.4)') "Light through soil/snow layer:  ", RadSwPenetrateGrd(I,:,J)
       stop "Error: Energy budget problem in NoahMP LSM"
    endif
#endif
  else if (noahmp%config%domain%IndicatorIceSfc(I,J) == -1) then ! ice point
    ! print out diagnostics when error is large
#ifdef _OPENACC
    ! Skip error checking on GPU
#else
    if ( abs(RadSwBalanceError(I,J)) > 0.01 ) then
       write(*,*) "GridIndexI, GridIndexJ = ", I, J
       write(*,*) "RadSwBalanceError(I,J)      = ", RadSwBalanceError(I,J)
       write(*,*) "RadSwDownRefHeight(I,J)     = ", RadSwDownRefHeight(I,J)
       write(*,*) "RadSwReflSfc(I,J)           = ", RadSwReflSfc(I,J)
       write(*,*) "RadSwAbsGrd(I,J)            = ", RadSwAbsGrd(I,J)
       write(*,*) "RadSwAbsSfc(I,J)            = ", RadSwAbsSfc(I,J)
       stop "Error: Solar radiation budget problem in NoahMP LSM"
    endif

    ! SNICAR
    if ( OptSnowAlbedo == 3 ) then
       if ( abs(RadSwAbsGrd(I,J)-sum(RadSwAbsSnowSoilLayer(I,:,J)))>0.001 ) then ! original check is 0.0001, precision issue
          write(*,*) "RadSwAbsGrd(I,J) gridmean                            = ", RadSwAbsGrd(I,J)
          write(*,*) "sum(RadSwAbsSnowSoilLayer) gridmean             = ", sum(RadSwAbsSnowSoilLayer(I,:,J))
          write(*,*) "RadSwAbsSnowSoilLayer gridmean                  = ", RadSwAbsSnowSoilLayer(I,:,J)
          write(*,*) "RadSwAbsGrd(I,J)-sum(RadSwAbsSnowSoilLayer) gridmean = ", RadSwAbsGrd(I,J)-sum(RadSwAbsSnowSoilLayer(I,:,J))
          stop "Error: SNICAR snow albedo radiation budget problem in NoahMP LSM"
       endif
    endif

    ! error in surface energy balance should be <0.01 W/m2
    EnergyBalanceError(I,J) = RadSwAbsGrd(I,J) + HeatPrecipAdvSfc(I,J) - (RadLwNetSfc(I,J) + HeatSensibleSfc(I,J) + HeatLatentGrd(I,J) + HeatGroundTot(I,J))
    ! print out diagnostics when error is large
    if ( abs(EnergyBalanceError(I,J)) > 0.01 ) then
       write(*,*) 'EnergyBalanceError(I,J) = ', EnergyBalanceError(I,J), ' at GridIndexI,GridIndexJ: ', I, J
       write(*,'(a17,F10.4)' ) "Net longwave:       ", RadLwNetSfc(I,J)
       write(*,'(a17,F10.4)' ) "Total sensible:     ", HeatSensibleSfc(I,J)
       write(*,'(a17,F10.4)' ) "Ground evap:        ", HeatLatentGrd(I,J)
       write(*,'(a17,F10.4)' ) "Total ground:       ", HeatGroundTot(I,J)
       write(*,'(a17,4F10.4)') "Precip advected:    ", HeatPrecipAdvSfc(I,J)
       write(*,'(a17,F10.4)' ) "absorbed shortwave: ", RadSwAbsGrd(I,J)
       stop "Error: Surface energy budget problem in NoahMP LSM"
    endif
#endif

  endif

      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine BalanceEnergyCheck

end module BalanceErrorCheckMod
