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
    associate( NumCropGrowStage => noahmp%config%domain%NumCropGrowStage ,&
               ITS => noahmp%config%domain%ITS, ITE => noahmp%config%domain%ITE ,&
               JTS => noahmp%config%domain%JTS, JTE => noahmp%config%domain%JTE )

    if (.not.(allocated(noahmp%forcing%SpecHumidityRefHeight))) then
      allocate(noahmp%forcing%SpecHumidityRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%SpecHumidityRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%TemperatureAirRefHeight))) then
      allocate(noahmp%forcing%TemperatureAirRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%TemperatureAirRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%WindEastwardRefHeight))) then
      allocate(noahmp%forcing%WindEastwardRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%WindEastwardRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%WindNorthwardRefHeight))) then
      allocate(noahmp%forcing%WindNorthwardRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%WindNorthwardRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%RadLwDownRefHeight))) then
      allocate(noahmp%forcing%RadLwDownRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%RadLwDownRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%RadSwDownRefHeight))) then
      allocate(noahmp%forcing%RadSwDownRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%RadSwDownRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%PrecipConvRefHeight))) then
      allocate(noahmp%forcing%PrecipConvRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%PrecipConvRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%PrecipNonConvRefHeight))) then
      allocate(noahmp%forcing%PrecipNonConvRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%PrecipNonConvRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%PrecipShConvRefHeight))) then
      allocate(noahmp%forcing%PrecipShConvRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%PrecipShConvRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%PrecipSnowRefHeight))) then
      allocate(noahmp%forcing%PrecipSnowRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%PrecipSnowRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%PrecipGraupelRefHeight))) then
      allocate(noahmp%forcing%PrecipGraupelRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%PrecipGraupelRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%PrecipHailRefHeight))) then
      allocate(noahmp%forcing%PrecipHailRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%PrecipHailRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%PressureAirSurface))) then
      allocate(noahmp%forcing%PressureAirSurface(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%PressureAirSurface)
    end if

    if (.not.(allocated(noahmp%forcing%PressureAirRefHeight))) then
      allocate(noahmp%forcing%PressureAirRefHeight(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%PressureAirRefHeight)
    end if

    if (.not.(allocated(noahmp%forcing%TemperatureSoilBottom))) then
      allocate(noahmp%forcing%TemperatureSoilBottom(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%TemperatureSoilBottom)
    end if

    if (.not.(allocated(noahmp%forcing%DepBChydropho))) then
      allocate(noahmp%forcing%DepBChydropho(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%DepBChydropho)
    end if

    if (.not.(allocated(noahmp%forcing%DepBChydrophi))) then
      allocate(noahmp%forcing%DepBChydrophi(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%DepBChydrophi)
    end if

    if (.not.(allocated(noahmp%forcing%DepOChydropho))) then
      allocate(noahmp%forcing%DepOChydropho(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%DepOChydropho)
    end if

    if (.not.(allocated(noahmp%forcing%DepOChydrophi))) then
      allocate(noahmp%forcing%DepOChydrophi(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%DepOChydrophi)
    end if

    if (.not.(allocated(noahmp%forcing%DepDust1))) then
      allocate(noahmp%forcing%DepDust1(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%DepDust1)
    end if

    if (.not.(allocated(noahmp%forcing%DepDust2))) then
      allocate(noahmp%forcing%DepDust2(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%DepDust2)
    end if

    if (.not.(allocated(noahmp%forcing%DepDust3))) then
      allocate(noahmp%forcing%DepDust3(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%DepDust3)
    end if

    if (.not.(allocated(noahmp%forcing%DepDust4))) then
      allocate(noahmp%forcing%DepDust4(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%DepDust4)
    end if

    if (.not.(allocated(noahmp%forcing%DepDust5))) then
      allocate(noahmp%forcing%DepDust5(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%DepDust5)
    end if

    if (.not.(allocated(noahmp%forcing%RadSwVisFrac))) then
      allocate(noahmp%forcing%RadSwVisFrac(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%RadSwVisFrac)
    end if

    if (.not.(allocated(noahmp%forcing%RadSwDirFrac))) then
      allocate(noahmp%forcing%RadSwDirFrac(ITS:ITE,JTS:JTE))
      !$acc enter data create(noahmp%forcing%RadSwDirFrac)
    end if

    end associate

    !$acc parallel loop collapse(2) gang vector present(noahmp%forcing)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                            &
                  SpecHumidityRefHeight   => noahmp%forcing%SpecHumidityRefHeight(I,J)       ,&
                  TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight(I,J)     ,&
                  WindEastwardRefHeight   => noahmp%forcing%WindEastwardRefHeight(I,J)       ,&
                  WindNorthwardRefHeight  => noahmp%forcing%WindNorthwardRefHeight(I,J)      ,&
                  RadLwDownRefHeight      => noahmp%forcing%RadLwDownRefHeight(I,J)          ,&
                  RadSwDownRefHeight      => noahmp%forcing%RadSwDownRefHeight(I,J)          ,&
                  PrecipConvRefHeight     => noahmp%forcing%PrecipConvRefHeight(I,J)         ,&
                  PrecipNonConvRefHeight  => noahmp%forcing%PrecipNonConvRefHeight(I,J)      ,&
                  PrecipShConvRefHeight   => noahmp%forcing%PrecipShConvRefHeight(I,J)       ,&
                  PrecipSnowRefHeight     => noahmp%forcing%PrecipSnowRefHeight(I,J)         ,&
                  PrecipGraupelRefHeight  => noahmp%forcing%PrecipGraupelRefHeight(I,J)      ,&
                  PrecipHailRefHeight     => noahmp%forcing%PrecipHailRefHeight(I,J)         ,&
                  PressureAirSurface      => noahmp%forcing%PressureAirSurface(I,J)          ,&
                  PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight(I,J)        ,&
                  TemperatureSoilBottom   => noahmp%forcing%TemperatureSoilBottom(I,J)       ,&
                  DepBChydropho           => noahmp%forcing%DepBChydropho(I,J)               ,&
                  DepBChydrophi           => noahmp%forcing%DepBChydrophi(I,J)               ,&
                  DepOChydropho           => noahmp%forcing%DepOChydropho(I,J)               ,&
                  DepOChydrophi           => noahmp%forcing%DepOChydrophi(I,J)               ,&
                  DepDust1                => noahmp%forcing%DepDust1(I,J)                    ,&
                  DepDust2                => noahmp%forcing%DepDust2(I,J)                    ,&
                  DepDust3                => noahmp%forcing%DepDust3(I,J)                    ,&
                  DepDust4                => noahmp%forcing%DepDust4(I,J)                    ,&
                  DepDust5                => noahmp%forcing%DepDust5(I,J)                    ,&
                  RadSwVisFrac            => noahmp%forcing%RadSwVisFrac(I,J)                ,&
                  RadSwDirFrac            => noahmp%forcing%RadSwDirFrac(I,J)                 &
                 )

        SpecHumidityRefHeight   = undefined_real
        TemperatureAirRefHeight = undefined_real
        WindEastwardRefHeight   = undefined_real
        WindNorthwardRefHeight  = undefined_real
        RadLwDownRefHeight      = undefined_real
        RadSwDownRefHeight      = undefined_real
        PrecipConvRefHeight     = undefined_real
        PrecipNonConvRefHeight  = undefined_real
        PrecipShConvRefHeight   = undefined_real
        PrecipSnowRefHeight     = undefined_real
        PrecipGraupelRefHeight  = undefined_real
        PrecipHailRefHeight     = undefined_real
        PressureAirSurface      = undefined_real
        PressureAirRefHeight    = undefined_real
        TemperatureSoilBottom   = undefined_real

        DepBChydropho           = undefined_real
        DepBChydrophi           = undefined_real
        DepOChydropho           = undefined_real
        DepOChydrophi           = undefined_real
        DepDust1                = undefined_real
        DepDust2                = undefined_real
        DepDust3                = undefined_real
        DepDust4                = undefined_real
        DepDust5                = undefined_real
        RadSwVisFrac            = undefined_real
        RadSwDirFrac            = undefined_real

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine ForcingVarInitDefault

end module ForcingVarInitMod
