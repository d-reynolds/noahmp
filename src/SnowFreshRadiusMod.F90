module SnowFreshRadiusMod

!!! Compute fresh fallen snow grain size for SNICAR albedo calculation
!!! Returns fresh snow grain radius, which is linearly dependent on temperature.
!!! This is implemented to remedy an outstanding bias that SNICAR has in initial
!!! grain size. See e.g. Sandells et al, 2017 for a discussion (10.5194/tc-11-229-2017).
!!! Yang et al. (2017), 10.1016/j.jqsrt.2016.03.033
!!! discusses grain size observations, which suggest a temperature dependence.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowFreshRadius(noahmp)

! ------------------------ Code history -----------------------------------
! Implementation: T.-S. Lin, C. He, et al. (2025, JHM)
! Adapted from CTSM function: FreshSnowRadius
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J   ! grid indices
    real(kind=kind_noahmp)           :: Tmin   ! start of linear ramp
    real(kind=kind_noahmp)           :: Tmax   ! end of linear ramp

! --------------------------------------------------------------------
    associate(                                                                      &
              TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight,& ! in,  air temperature [K] at reference height
              SnowRadiusMin           => noahmp%water%param%SnowRadiusMin      ,& ! in,  minimum allowed snow effective radius (also cold "fresh snow" value) [microns]
              FreshSnowRadiusMax      => noahmp%water%param%FreshSnowRadiusMax ,& ! in,  maximum warm fresh snow effective radius [microns]
              SnowRadiusFresh         => noahmp%water%state%SnowRadiusFresh     & ! out, fresh snow radius [microns]
             )

   !$acc parallel loop collapse(2) gang vector default(present) private(Tmin,Tmax)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        ! solar radiation process is only done if there is light
        if ( noahmp%config%domain%CosSolarZenithAngle(I,J) <= 0 ) cycle


    Tmin = ConstFreezePoint - 30.0
    Tmax = ConstFreezePoint

    if ( FreshSnowRadiusMax(I,J) <= SnowRadiusMin(I,J) )then
        SnowRadiusFresh(I,J) = SnowRadiusMin(I,J)
    else
        if (TemperatureAirRefHeight(I,J) < Tmin) then
            SnowRadiusFresh(I,J) = SnowRadiusMin(I,J)
        else if (TemperatureAirRefHeight(I,J) > Tmax) then
            SnowRadiusFresh(I,J) = FreshSnowRadiusMax(I,J)
        else
            SnowRadiusFresh(I,J) = (Tmax - TemperatureAirRefHeight(I,J)) / (Tmax - Tmin) * SnowRadiusMin(I,J) + &
                              (TemperatureAirRefHeight(I,J) - Tmin) / (Tmax - Tmin) * FreshSnowRadiusMax(I,J)
        end if
    end if


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine SnowFreshRadius

end module SnowFreshRadiusMod
