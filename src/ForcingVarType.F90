module ForcingVarType

!!! Define 2D Noah-MP forcing variables
!!! Forcing variable initialization is done in ForcingVarInitMod.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

  use Machine

  implicit none
  save
  private

  type, public :: forcing_type

    ! All forcing fields are now 2D arrays (I,J)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SpecHumidityRefHeight    ! Specific humidity [kg water vapor / kg moist air] forcing at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureAirRefHeight  ! Air temperature [K] forcing at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WindEastwardRefHeight    ! wind speed [m/s] in eastward dir at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: WindNorthwardRefHeight   ! wind speed [m/s] in northward dir at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadSwDownRefHeight       ! downward shortwave radiation [W/m2] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadLwDownRefHeight       ! downward longwave radiation [W/m2] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PressureAirRefHeight     ! air pressure [Pa] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PressureAirSurface       ! air pressure [Pa] at surface-atmosphere interface (lowest atmos model boundary)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PrecipConvRefHeight      ! convective precipitation rate [mm/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PrecipNonConvRefHeight   ! non-convective precipitation rate [mm/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PrecipShConvRefHeight    ! shallow convective precipitation rate [mm/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PrecipSnowRefHeight      ! snowfall rate [mm/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PrecipGraupelRefHeight   ! graupel rate [mm/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: PrecipHailRefHeight      ! hail rate [mm/s] at reference height
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TemperatureSoilBottom    ! bottom boundary condition for soil temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DepBChydropho            ! hydrophobic Black Carbon deposition flux [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DepBChydrophi            ! hydrophillic Black Carbon deposition flux [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DepOChydropho            ! hydrophobic Organic Carbon deposition flux [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DepOChydrophi            ! hydrophillic Organic Carbon deposition flux [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DepDust1                 ! dust species 1 deposition flux [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DepDust2                 ! dust species 2 deposition flux [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DepDust3                 ! dust species 3 deposition flux [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DepDust4                 ! dust species 4 deposition flux [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DepDust5                 ! dust species 5 deposition flux [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadSwVisFrac             ! fraction of visible band radiation
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RadSwDirFrac             ! fraction of direct raidation

  end type forcing_type

end module ForcingVarType
