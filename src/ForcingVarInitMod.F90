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

    !$acc parallel loop collapse(2) gang vector present(noahmp)
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
