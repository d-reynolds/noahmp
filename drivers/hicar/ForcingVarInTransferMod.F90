module ForcingVarInTransferMod

!!! Transfer input 2-D NoahmpIO Forcing variables to 1-D column variable
!!! 1-D variables should be first defined in /src/ForcingVarType.F90
!!! 2-D variables should be first defined in NoahmpIOVarType.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type
  use NoahmpVarType

  implicit none

contains

!=== initialize with input data or table values

  subroutine ForcingVarInTransfer(noahmp, NoahmpIO)

    implicit none

    type(NoahmpIO_type), intent(inout) :: NoahmpIO
    type(noahmp_type),   intent(inout) :: noahmp
    
    ! local variables
    integer                             :: I, J
    real(kind=kind_noahmp)              :: PrecipOtherRefHeight  ! other precipitation, e.g. fog [mm/s] at reference height
    real(kind=kind_noahmp)              :: PrecipTotalRefHeight  ! total precipitation [mm/s] at reference height

    !$acc parallel loop collapse(2) present(noahmp, NoahmpIO)
    do J = noahmp%config%domain%JTS, noahmp%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
    noahmp%forcing%TemperatureAirRefHeight(I,J) = NoahmpIO%T_PHY(I,1,J)
    noahmp%forcing%WindEastwardRefHeight(I,J)   = NoahmpIO%U_PHY(I,1,J)
    noahmp%forcing%WindNorthwardRefHeight(I,J)  = NoahmpIO%V_PHY(I,1,J)
    noahmp%forcing%SpecHumidityRefHeight(I,J)   = NoahmpIO%QV_CURR(I,1,J)/(1.0+NoahmpIO%QV_CURR(I,1,J))  ! convert from mixing ratio to specific humidity
    noahmp%forcing%PressureAirRefHeight(I,J)    = (NoahmpIO%P8W(I,NoahmpIO%KTS,J) + NoahmpIO%P8W(I,NoahmpIO%KTS+1,J)) * 0.5 ! air pressure at middle point of lowest atmos model layer
    noahmp%forcing%PressureAirSurface(I,J)      = NoahmpIO%P8W      (I,1,J)
    noahmp%forcing%RadLwDownRefHeight(I,J)      = NoahmpIO%GLW      (I,J)
    noahmp%forcing%RadSwDownRefHeight(I,J)      = NoahmpIO%SWDOWN   (I,J)
    noahmp%forcing%TemperatureSoilBottom(I,J)   = NoahmpIO%TMN      (I,J)

    ! treat different precipitation types
    PrecipTotalRefHeight                   = NoahmpIO%RAINBL   (I,J) / NoahmpIO%DTBL                ! convert precip unit from mm/timestep to mm/s
    noahmp%forcing%PrecipConvRefHeight(I,J)     = NoahmpIO%MP_RAINC (I,J) / NoahmpIO%DTBL
    noahmp%forcing%PrecipNonConvRefHeight(I,J)  = NoahmpIO%MP_RAINNC(I,J) / NoahmpIO%DTBL
    noahmp%forcing%PrecipShConvRefHeight(I,J)   = NoahmpIO%MP_SHCV  (I,J) / NoahmpIO%DTBL
    noahmp%forcing%PrecipSnowRefHeight(I,J)     = NoahmpIO%MP_SNOW  (I,J) / NoahmpIO%DTBL
    noahmp%forcing%PrecipGraupelRefHeight(I,J)  = NoahmpIO%MP_GRAUP (I,J) / NoahmpIO%DTBL
    noahmp%forcing%PrecipHailRefHeight(I,J)     = NoahmpIO%MP_HAIL  (I,J) / NoahmpIO%DTBL
    ! treat other precipitation (e.g. fog) contained in total precipitation
    PrecipOtherRefHeight                   = PrecipTotalRefHeight - noahmp%forcing%PrecipConvRefHeight(I,J) - &
                                             noahmp%forcing%PrecipNonConvRefHeight(I,J) - noahmp%forcing%PrecipShConvRefHeight(I,J)
    PrecipOtherRefHeight                   = max(0.0, PrecipOtherRefHeight)
    noahmp%forcing%PrecipNonConvRefHeight(I,J)  = noahmp%forcing%PrecipNonConvRefHeight(I,J) + PrecipOtherRefHeight
    noahmp%forcing%PrecipSnowRefHeight(I,J)     = noahmp%forcing%PrecipSnowRefHeight(I,J) + PrecipOtherRefHeight * NoahmpIO%SR(I,J)

    ! downward solar radiation direct/diffuse and visible/NIR partition
    noahmp%forcing%RadSwDirFrac(I,J)            = NoahmpIO%RadSwDirFrac(I,J)
    noahmp%forcing%RadSwVisFrac(I,J)            = NoahmpIO%RadSwVisFrac(I,J)

    ! SNICAR aerosol deposition flux forcing
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( noahmp%config%nmlist%FlagSnicarAerosolReadTable .eqv. .true. ) then 
          noahmp%forcing%DepBChydropho(I,J)     = NoahmpIO%DepBChydropho_TABLE
          noahmp%forcing%DepBChydrophi(I,J)     = NoahmpIO%DepBChydrophi_TABLE
          noahmp%forcing%DepOChydropho(I,J)     = NoahmpIO%DepOChydropho_TABLE
          noahmp%forcing%DepOChydrophi(I,J)     = NoahmpIO%DepOChydrophi_TABLE
          noahmp%forcing%DepDust1(I,J)          = NoahmpIO%DepDust1_TABLE
          noahmp%forcing%DepDust2(I,J)          = NoahmpIO%DepDust2_TABLE
          noahmp%forcing%DepDust3(I,J)          = NoahmpIO%DepDust3_TABLE
          noahmp%forcing%DepDust4(I,J)          = NoahmpIO%DepDust4_TABLE
          noahmp%forcing%DepDust5(I,J)          = NoahmpIO%DepDust5_TABLE
       else
          noahmp%forcing%DepBChydropho(I,J)     = NoahmpIO%DepBChydrophoXY(I,J)
          noahmp%forcing%DepBChydrophi(I,J)     = NoahmpIO%DepBChydrophiXY(I,J)
          noahmp%forcing%DepOChydropho(I,J)     = NoahmpIO%DepOChydrophoXY(I,J)
          noahmp%forcing%DepOChydrophi(I,J)     = NoahmpIO%DepOChydrophiXY(I,J)
          noahmp%forcing%DepDust1(I,J)          = NoahmpIO%DepDust1XY(I,J)
          noahmp%forcing%DepDust2(I,J)          = NoahmpIO%DepDust2XY(I,J)
          noahmp%forcing%DepDust3(I,J)          = NoahmpIO%DepDust3XY(I,J)
          noahmp%forcing%DepDust4(I,J)          = NoahmpIO%DepDust4XY(I,J)
          noahmp%forcing%DepDust5(I,J)          = NoahmpIO%DepDust5XY(I,J)
       endif
    endif

  end do
  end do

 
  end subroutine ForcingVarInTransfer

end module ForcingVarInTransferMod
