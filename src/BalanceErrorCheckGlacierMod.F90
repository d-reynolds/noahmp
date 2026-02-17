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
    associate(&
              SnowWaterEquiv     => noahmp%water%state%SnowWaterEquiv ,& ! in,  snow water equivalent [mm]
              WaterStorageTotBeg => noahmp%water%state%WaterStorageTotBeg  & ! out, total water storage [mm] at the beginning
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! compute total glacier water storage before NoahMP processes
    ! need more work on including glacier ice mass underneath snow
    WaterStorageTotBeg(I,J) = SnowWaterEquiv(I,J)


      end do
    end do
   !$acc end parallel loop

    end associate

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
    associate(&
              SnowWaterEquiv     => noahmp%water%state%SnowWaterEquiv ,& ! in,  snow water equivalent [mm]
              WaterStorageTotBeg => noahmp%water%state%WaterStorageTotBeg ,& ! in,  total water storage [mm] at the beginning
              WaterStorageTotEnd => noahmp%water%state%WaterStorageTotEnd ,& ! out, total water storage [mm] at the end
              WaterBalanceError  => noahmp%water%state%WaterBalanceError ,& ! out, water balance error [mm] per time step
              MainTimeStep       => noahmp%config%domain%MainTimeStep ,& ! in,  main noahmp timestep [s]
              PrecipTotRefHeight => noahmp%water%flux%PrecipTotRefHeight ,& ! in,  total precipitation [mm/s] at reference height
              EvapGroundNet      => noahmp%water%flux%EvapGroundNet ,& ! in,  net ground evaporation [mm/s]
              RunoffSurface      => noahmp%water%flux%RunoffSurface ,& ! in,  surface runoff [mm/s]
              RunoffSubsurface   => noahmp%water%flux%RunoffSubsurface  & ! in,  subsurface runoff [mm/s]
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


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


      end do
    end do
   !$acc end parallel loop

    end associate

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
    associate(&
              OptSnowAlbedo        => noahmp%config%nmlist%OptSnowAlbedo ,& ! in,  options for ground snow surface albedo
              RadSwDownRefHeight   => noahmp%forcing%RadSwDownRefHeight ,& ! in,  downward shortwave radiation [W/m2] at reference height
              RadSwAbsSfc          => noahmp%energy%flux%RadSwAbsSfc ,& ! in,  total absorbed solar radiation [W/m2]
              RadSwAbsSnowSoilLayer=> noahmp%energy%flux%RadSwAbsSnowSoilLayer ,& ! in,  total absorbed solar radiation by snow/soil for each layer [W/m2]
              RadSwReflSfc         => noahmp%energy%flux%RadSwReflSfc ,& ! in,  total reflected solar radiation [W/m2]
              RadLwNetSfc          => noahmp%energy%flux%RadLwNetSfc ,& ! in,  total net longwave rad [W/m2] (+ to atm)
              HeatSensibleSfc      => noahmp%energy%flux%HeatSensibleSfc ,& ! in,  total sensible heat [W/m2] (+ to atm)
              HeatLatentGrd        => noahmp%energy%flux%HeatLatentGrd ,& ! in,  total ground latent heat [W/m2] (+ to atm)
              HeatGroundTot        => noahmp%energy%flux%HeatGroundTot ,& ! in,  total ground heat flux [W/m2] (+ to soil/snow)
              RadSwAbsGrd          => noahmp%energy%flux%RadSwAbsGrd ,& ! in,  solar radiation absorbed by ground [W/m2]
              HeatPrecipAdvSfc     => noahmp%energy%flux%HeatPrecipAdvSfc ,& ! in,  precipitation advected heat - total [W/m2]
              EnergyBalanceError   => noahmp%energy%state%EnergyBalanceError ,& ! out, error in surface energy balance [W/m2]
              RadSwBalanceError    => noahmp%energy%state%RadSwBalanceError  & ! out, error in shortwave radiation balance [W/m2]
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! error in shortwave radiation balance should be <0.01 W/m2
    RadSwBalanceError(I,J) = RadSwDownRefHeight(I,J) - (RadSwAbsSfc(I,J) + RadSwReflSfc(I,J))

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


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine BalanceEnergyCheckGlacier

end module BalanceErrorCheckGlacierMod
