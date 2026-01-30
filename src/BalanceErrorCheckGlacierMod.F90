module BalanceErrorCheckGlacierMod

!!! Check glacier water and energy balance and report error

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

!!!! Water balance check initialization
  subroutine BalanceWaterInitGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in NOAHMP_GLACIER)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: I, J                        ! grid indices

! --------------------------------------------------------------------
   !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                             &
              SnowWaterEquiv     => noahmp%water%state%SnowWaterEquiv(I,J)     ,& ! in,  snow water equivalent [mm]
              WaterStorageTotBeg => noahmp%water%state%WaterStorageTotBeg(I,J)  & ! out, total water storage [mm] at the beginning
             )
! ----------------------------------------------------------------------

    ! compute total glacier water storage before NoahMP processes
    ! need more work on including glacier ice mass underneath snow
    WaterStorageTotBeg = SnowWaterEquiv

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine BalanceWaterInitGlacier


!!!! Water balance check and report error
  subroutine BalanceWaterCheckGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ERROR_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: I, J                        ! grid indices

! --------------------------------------------------------------------
   !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                             &
              MainTimeStep       => noahmp%config%domain%MainTimeStep     ,& ! in,  main noahmp timestep [s]
              SnowWaterEquiv     => noahmp%water%state%SnowWaterEquiv(I,J)     ,& ! in,  snow water equivalent [mm]
              WaterStorageTotBeg => noahmp%water%state%WaterStorageTotBeg(I,J) ,& ! in,  total water storage [mm] at the beginning
              PrecipTotRefHeight => noahmp%water%flux%PrecipTotRefHeight(I,J)  ,& ! in,  total precipitation [mm/s] at reference height
              EvapGroundNet      => noahmp%water%flux%EvapGroundNet(I,J)       ,& ! in,  net ground evaporation [mm/s]
              RunoffSurface      => noahmp%water%flux%RunoffSurface(I,J)       ,& ! in,  surface runoff [mm/s]
              RunoffSubsurface   => noahmp%water%flux%RunoffSubsurface(I,J)    ,& ! in,  subsurface runoff [mm/s]
              WaterStorageTotEnd => noahmp%water%state%WaterStorageTotEnd(I,J) ,& ! out, total water storage [mm] at the end
              WaterBalanceError  => noahmp%water%state%WaterBalanceError(I,J)   & ! out, water balance error [mm] per time step
             )
! ----------------------------------------------------------------------

    ! Error in water balance should be < 0.1 mm
    ! compute total glacier water storage before NoahMP processes
    ! need more work on including glacier ice mass underneath snow
    WaterStorageTotEnd = SnowWaterEquiv
    WaterBalanceError  = WaterStorageTotEnd - WaterStorageTotBeg - &
                         (PrecipTotRefHeight - EvapGroundNet - RunoffSurface - RunoffSubsurface) * MainTimeStep

#if !defined(WRF_HYDRO) && !defined(_OPENACC)
    if ( abs(WaterBalanceError) > 0.1 ) then
       if ( WaterBalanceError > 0) then
          write(*,*) "The model is gaining water (WaterBalanceError is positive)"
       else
          write(*,*) "The model is losing water (WaterBalanceError is negative)"
       endif
       write(*,*) "WaterBalanceError = ",WaterBalanceError, "kg m{-2} timestep{-1}"
       write(*, &
           '("  GridIndexI   GridIndexJ     WaterStorageTotEnd  WaterStorageTotBeg  PrecipTotRefHeight  &
                EvapGroundNet  RunoffSurface  RunoffSubsurface")')
       write(*,'(i6,1x,i6,1x,2f15.3,9f11.5)') I, J, WaterStorageTotEnd, WaterStorageTotBeg, &
                                              PrecipTotRefHeight*MainTimeStep, EvapGroundNet*MainTimeStep,    &
                                              RunoffSurface*MainTimeStep, RunoffSubsurface*MainTimeStep
       stop "Error: Water budget problem in NoahMP LSM"
    endif
#endif

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine BalanceWaterCheckGlacier


!!!! Energy balance check and error report
  subroutine BalanceEnergyCheckGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ERROR_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: I, J                        ! grid indices

! --------------------------------------------------------------------
   !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                                 &
              OptSnowAlbedo        => noahmp%config%nmlist%OptSnowAlbedo      ,& ! in,  options for ground snow surface albedo
              RadSwDownRefHeight   => noahmp%forcing%RadSwDownRefHeight(I,J)       ,& ! in,  downward shortwave radiation [W/m2] at reference height
              RadSwAbsSfc          => noahmp%energy%flux%RadSwAbsSfc(I,J)          ,& ! in,  total absorbed solar radiation [W/m2]
              RadSwAbsSnowSoilLayer=> noahmp%energy%flux%RadSwAbsSnowSoilLayer     ,& ! in,  total absorbed solar radiation by snow/soil for each layer [W/m2]
              RadSwReflSfc         => noahmp%energy%flux%RadSwReflSfc(I,J)         ,& ! in,  total reflected solar radiation [W/m2]
              RadLwNetSfc          => noahmp%energy%flux%RadLwNetSfc(I,J)          ,& ! in,  total net longwave rad [W/m2] (+ to atm)
              HeatSensibleSfc      => noahmp%energy%flux%HeatSensibleSfc(I,J)      ,& ! in,  total sensible heat [W/m2] (+ to atm)
              HeatLatentGrd        => noahmp%energy%flux%HeatLatentGrd(I,J)        ,& ! in,  total ground latent heat [W/m2] (+ to atm)
              HeatGroundTot        => noahmp%energy%flux%HeatGroundTot(I,J)        ,& ! in,  total ground heat flux [W/m2] (+ to soil/snow)
              RadSwAbsGrd          => noahmp%energy%flux%RadSwAbsGrd(I,J)          ,& ! in,  solar radiation absorbed by ground [W/m2]
              HeatPrecipAdvSfc     => noahmp%energy%flux%HeatPrecipAdvSfc(I,J)     ,& ! in,  precipitation advected heat - total [W/m2]
              EnergyBalanceError   => noahmp%energy%state%EnergyBalanceError(I,J)  ,& ! out, error in surface energy balance [W/m2]
              RadSwBalanceError    => noahmp%energy%state%RadSwBalanceError(I,J)    & ! out, error in shortwave radiation balance [W/m2]
             )
! ----------------------------------------------------------------------

    ! error in shortwave radiation balance should be <0.01 W/m2
    RadSwBalanceError = RadSwDownRefHeight - (RadSwAbsSfc + RadSwReflSfc)

    ! print out diagnostics when error is large
#ifdef _OPENACC
    ! Skip error checking on GPU
#else
    if ( abs(RadSwBalanceError) > 0.01 ) then
       write(*,*) "GridIndexI, GridIndexJ = ", I, J
       write(*,*) "RadSwBalanceError      = ", RadSwBalanceError
       write(*,*) "RadSwDownRefHeight     = ", RadSwDownRefHeight
       write(*,*) "RadSwReflSfc           = ", RadSwReflSfc
       write(*,*) "RadSwAbsGrd            = ", RadSwAbsGrd
       write(*,*) "RadSwAbsSfc            = ", RadSwAbsSfc
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

    ! error in surface energy balance should be <0.01 W/m2
    EnergyBalanceError = RadSwAbsGrd + HeatPrecipAdvSfc - (RadLwNetSfc + HeatSensibleSfc + HeatLatentGrd + HeatGroundTot)
    ! print out diagnostics when error is large
    if ( abs(EnergyBalanceError) > 0.01 ) then
       write(*,*) 'EnergyBalanceError = ', EnergyBalanceError, ' at GridIndexI,GridIndexJ: ', I, J
       write(*,'(a17,F10.4)' ) "Net longwave:       ", RadLwNetSfc
       write(*,'(a17,F10.4)' ) "Total sensible:     ", HeatSensibleSfc
       write(*,'(a17,F10.4)' ) "Ground evap:        ", HeatLatentGrd
       write(*,'(a17,F10.4)' ) "Total ground:       ", HeatGroundTot
       write(*,'(a17,4F10.4)') "Precip advected:    ", HeatPrecipAdvSfc
       write(*,'(a17,F10.4)' ) "absorbed shortwave: ", RadSwAbsGrd
       stop "Error: Surface energy budget problem in NoahMP LSM"
    endif
#endif

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine BalanceEnergyCheckGlacier

end module BalanceErrorCheckGlacierMod
