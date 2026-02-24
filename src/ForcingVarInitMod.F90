module ForcingVarInitMod

!!! Initialize column (1-D) Noah-MP forcing variables
!!! Forcing variables should be first defined in ForcingVarType.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpVarType

  implicit none

contains

!=== initialize with default values
  subroutine ForcingVarInitDefault(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer :: I, J      ! grid indices

    ! Allocate 3D crop parameter arrays and transfer to GPU
    associate(                                                                      &
              NumCropGrowStage => noahmp%config%domain%NumCropGrowStage ,&
               ITS => noahmp%config%domain%ITS, ITE => noahmp%config%domain%ITE ,&
               JTS => noahmp%config%domain%JTS, JTE => noahmp%config%domain%JTE  &
             )

    ! Phase 1: Host allocations
    if (.not.(allocated(noahmp%forcing%SpecHumidityRefHeight))) &
      allocate(noahmp%forcing%SpecHumidityRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%TemperatureAirRefHeight))) &
      allocate(noahmp%forcing%TemperatureAirRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%WindEastwardRefHeight))) &
      allocate(noahmp%forcing%WindEastwardRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%WindNorthwardRefHeight))) &
      allocate(noahmp%forcing%WindNorthwardRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%RadLwDownRefHeight))) &
      allocate(noahmp%forcing%RadLwDownRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%RadSwDownRefHeight))) &
      allocate(noahmp%forcing%RadSwDownRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%PrecipConvRefHeight))) &
      allocate(noahmp%forcing%PrecipConvRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%PrecipNonConvRefHeight))) &
      allocate(noahmp%forcing%PrecipNonConvRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%PrecipShConvRefHeight))) &
      allocate(noahmp%forcing%PrecipShConvRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%PrecipSnowRefHeight))) &
      allocate(noahmp%forcing%PrecipSnowRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%PrecipGraupelRefHeight))) &
      allocate(noahmp%forcing%PrecipGraupelRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%PrecipHailRefHeight))) &
      allocate(noahmp%forcing%PrecipHailRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%PressureAirSurface))) &
      allocate(noahmp%forcing%PressureAirSurface(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%PressureAirRefHeight))) &
      allocate(noahmp%forcing%PressureAirRefHeight(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%TemperatureSoilBottom))) &
      allocate(noahmp%forcing%TemperatureSoilBottom(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%DepBChydropho))) &
      allocate(noahmp%forcing%DepBChydropho(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%DepBChydrophi))) &
      allocate(noahmp%forcing%DepBChydrophi(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%DepOChydropho))) &
      allocate(noahmp%forcing%DepOChydropho(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%DepOChydrophi))) &
      allocate(noahmp%forcing%DepOChydrophi(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%DepDust1))) &
      allocate(noahmp%forcing%DepDust1(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%DepDust2))) &
      allocate(noahmp%forcing%DepDust2(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%DepDust3))) &
      allocate(noahmp%forcing%DepDust3(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%DepDust4))) &
      allocate(noahmp%forcing%DepDust4(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%DepDust5))) &
      allocate(noahmp%forcing%DepDust5(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%RadSwVisFrac))) &
      allocate(noahmp%forcing%RadSwVisFrac(ITS:ITE,JTS:JTE))
    if (.not.(allocated(noahmp%forcing%RadSwDirFrac))) &
      allocate(noahmp%forcing%RadSwDirFrac(ITS:ITE,JTS:JTE))

    ! Phase 2: Batched async device create
    !$acc enter data create(                              &
    !$acc   noahmp%forcing%SpecHumidityRefHeight,         &
    !$acc   noahmp%forcing%TemperatureAirRefHeight,       &
    !$acc   noahmp%forcing%WindEastwardRefHeight,         &
    !$acc   noahmp%forcing%WindNorthwardRefHeight,        &
    !$acc   noahmp%forcing%RadLwDownRefHeight,            &
    !$acc   noahmp%forcing%RadSwDownRefHeight,            &
    !$acc   noahmp%forcing%PrecipConvRefHeight,           &
    !$acc   noahmp%forcing%PrecipNonConvRefHeight,        &
    !$acc   noahmp%forcing%PrecipShConvRefHeight,         &
    !$acc   noahmp%forcing%PrecipSnowRefHeight,           &
    !$acc   noahmp%forcing%PrecipGraupelRefHeight,        &
    !$acc   noahmp%forcing%PrecipHailRefHeight,           &
    !$acc   noahmp%forcing%PressureAirSurface,            &
    !$acc   noahmp%forcing%PressureAirRefHeight,          &
    !$acc   noahmp%forcing%TemperatureSoilBottom,         &
    !$acc   noahmp%forcing%DepBChydropho,                 &
    !$acc   noahmp%forcing%DepBChydrophi,                 &
    !$acc   noahmp%forcing%DepOChydropho,                 &
    !$acc   noahmp%forcing%DepOChydrophi,                 &
    !$acc   noahmp%forcing%DepDust1,                      &
    !$acc   noahmp%forcing%DepDust2,                      &
    !$acc   noahmp%forcing%DepDust3,                      &
    !$acc   noahmp%forcing%DepDust4,                      &
    !$acc   noahmp%forcing%DepDust5,                      &
    !$acc   noahmp%forcing%RadSwVisFrac,                  &
    !$acc   noahmp%forcing%RadSwDirFrac                   &
    !$acc   ) async(NOAHMP_ACC_QUEUE)

    end associate

    ! Phase 3: Wait for device memory, then initialize
    !$acc wait(NOAHMP_ACC_QUEUE)

    associate(                                                                      &
              NumCropGrowStage => noahmp%config%domain%NumCropGrowStage ,&
               ITS => noahmp%config%domain%ITS, ITE => noahmp%config%domain%ITE ,&
               JTS => noahmp%config%domain%JTS, JTE => noahmp%config%domain%JTE ,&
                  SpecHumidityRefHeight   => noahmp%forcing%SpecHumidityRefHeight ,&
                  TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight ,&
                  WindEastwardRefHeight   => noahmp%forcing%WindEastwardRefHeight ,&
                  WindNorthwardRefHeight  => noahmp%forcing%WindNorthwardRefHeight ,&
                  RadLwDownRefHeight      => noahmp%forcing%RadLwDownRefHeight ,&
                  RadSwDownRefHeight      => noahmp%forcing%RadSwDownRefHeight ,&
                  PrecipConvRefHeight     => noahmp%forcing%PrecipConvRefHeight ,&
                  PrecipNonConvRefHeight  => noahmp%forcing%PrecipNonConvRefHeight ,&
                  PrecipShConvRefHeight   => noahmp%forcing%PrecipShConvRefHeight ,&
                  PrecipSnowRefHeight     => noahmp%forcing%PrecipSnowRefHeight ,&
                  PrecipGraupelRefHeight  => noahmp%forcing%PrecipGraupelRefHeight ,&
                  PrecipHailRefHeight     => noahmp%forcing%PrecipHailRefHeight ,&
                  PressureAirSurface      => noahmp%forcing%PressureAirSurface ,&
                  PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight ,&
                  TemperatureSoilBottom   => noahmp%forcing%TemperatureSoilBottom ,&
                  DepBChydropho           => noahmp%forcing%DepBChydropho ,&
                  DepBChydrophi           => noahmp%forcing%DepBChydrophi ,&
                  DepOChydropho           => noahmp%forcing%DepOChydropho ,&
                  DepOChydrophi           => noahmp%forcing%DepOChydrophi ,&
                  DepDust1                => noahmp%forcing%DepDust1 ,&
                  DepDust2                => noahmp%forcing%DepDust2 ,&
                  DepDust3                => noahmp%forcing%DepDust3 ,&
                  DepDust4                => noahmp%forcing%DepDust4 ,&
                  DepDust5                => noahmp%forcing%DepDust5 ,&
                  RadSwVisFrac            => noahmp%forcing%RadSwVisFrac ,&
                  RadSwDirFrac            => noahmp%forcing%RadSwDirFrac  &
             )

    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


        SpecHumidityRefHeight(I,J)   = undefined_real
        TemperatureAirRefHeight(I,J) = undefined_real
        WindEastwardRefHeight(I,J)   = undefined_real
        WindNorthwardRefHeight(I,J)  = undefined_real
        RadLwDownRefHeight(I,J)      = undefined_real
        RadSwDownRefHeight(I,J)      = undefined_real
        PrecipConvRefHeight(I,J)     = undefined_real
        PrecipNonConvRefHeight(I,J)  = undefined_real
        PrecipShConvRefHeight(I,J)   = undefined_real
        PrecipSnowRefHeight(I,J)     = undefined_real
        PrecipGraupelRefHeight(I,J)  = undefined_real
        PrecipHailRefHeight(I,J)     = undefined_real
        PressureAirSurface(I,J)      = undefined_real
        PressureAirRefHeight(I,J)    = undefined_real
        TemperatureSoilBottom(I,J)   = undefined_real

        DepBChydropho(I,J)           = undefined_real
        DepBChydrophi(I,J)           = undefined_real
        DepOChydropho(I,J)           = undefined_real
        DepOChydrophi(I,J)           = undefined_real
        DepDust1(I,J)                = undefined_real
        DepDust2(I,J)                = undefined_real
        DepDust3(I,J)                = undefined_real
        DepDust4(I,J)                = undefined_real
        DepDust5(I,J)                = undefined_real
        RadSwVisFrac(I,J)            = undefined_real
        RadSwDirFrac(I,J)            = undefined_real


      end do
    end do
    !$acc end parallel loop


    end associate

  end subroutine ForcingVarInitDefault


  subroutine ForcingVarExitDevice(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    !$acc exit data delete(               &
    !$acc   noahmp%forcing%SpecHumidityRefHeight,   &
    !$acc   noahmp%forcing%TemperatureAirRefHeight,   &
    !$acc   noahmp%forcing%WindEastwardRefHeight,   &
    !$acc   noahmp%forcing%WindNorthwardRefHeight,   &
    !$acc   noahmp%forcing%RadLwDownRefHeight,   &
    !$acc   noahmp%forcing%RadSwDownRefHeight,   &
    !$acc   noahmp%forcing%PrecipConvRefHeight,   &
    !$acc   noahmp%forcing%PrecipNonConvRefHeight,   &
    !$acc   noahmp%forcing%PrecipShConvRefHeight,   &
    !$acc   noahmp%forcing%PrecipSnowRefHeight,   &
    !$acc   noahmp%forcing%PrecipGraupelRefHeight,   &
    !$acc   noahmp%forcing%PrecipHailRefHeight,   &
    !$acc   noahmp%forcing%PressureAirSurface,   &
    !$acc   noahmp%forcing%PressureAirRefHeight,   &
    !$acc   noahmp%forcing%TemperatureSoilBottom,   &
    !$acc   noahmp%forcing%DepBChydropho,   &
    !$acc   noahmp%forcing%DepBChydrophi,   &
    !$acc   noahmp%forcing%DepOChydropho,   &
    !$acc   noahmp%forcing%DepOChydrophi,   &
    !$acc   noahmp%forcing%DepDust1,   &
    !$acc   noahmp%forcing%DepDust2,   &
    !$acc   noahmp%forcing%DepDust3,   &
    !$acc   noahmp%forcing%DepDust4,   &
    !$acc   noahmp%forcing%DepDust5,   &
    !$acc   noahmp%forcing%RadSwVisFrac,   &
    !$acc   noahmp%forcing%RadSwDirFrac    &
    !$acc   )

  end subroutine ForcingVarExitDevice
end module ForcingVarInitMod
