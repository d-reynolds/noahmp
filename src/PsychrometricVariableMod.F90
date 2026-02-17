module PsychrometricVariableMod

!!! Compute psychrometric variables for canopy and ground (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine PsychrometricVariable(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in ENERGY subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer :: I, J  ! grid indices

! --------------------------------------------------------------------
        associate(                                                                     &
                  PressureAirRefHeight => noahmp%forcing%PressureAirRefHeight   ,& ! in,  air pressure [Pa] at reference height
                  TemperatureCanopy    => noahmp%energy%state%TemperatureCanopy ,& ! in,  vegetation temperature [K]
                  TemperatureGrd       => noahmp%energy%state%TemperatureGrd    ,& ! in,  ground temperature [K]
                  LatHeatVapCanopy     => noahmp%energy%state%LatHeatVapCanopy  ,& ! out, latent heat of vaporization/subli [J/kg], canopy
                  LatHeatVapGrd        => noahmp%energy%state%LatHeatVapGrd     ,& ! out, latent heat of vaporization/subli [J/kg], ground
                  FlagFrozenCanopy     => noahmp%energy%state%FlagFrozenCanopy  ,& ! out, used to define latent heat pathway
                  FlagFrozenGround     => noahmp%energy%state%FlagFrozenGround  ,& ! out, frozen ground (logical) to define latent heat pathway
                  PsychConstCanopy     => noahmp%energy%state%PsychConstCanopy  ,& ! out, psychrometric constant [Pa/K], canopy
                  PsychConstGrd        => noahmp%energy%state%PsychConstGrd      & ! out, psychrometric constant [Pa/K], ground
                 )

    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    if (noahmp%config%domain%IndicatorIceSfc(I,J) == 0) then

      ! for canopy  ! Barlage: add distinction between ground and vegetation in v3.6
      if ( TemperatureCanopy(I,J) > ConstFreezePoint ) then
        LatHeatVapCanopy(I,J) = ConstLatHeatEvap
        FlagFrozenCanopy(I,J) = .false.
      else
        LatHeatVapCanopy(I,J) = ConstLatHeatSublim
        FlagFrozenCanopy(I,J) = .true.
      endif
      PsychConstCanopy(I,J)    = ConstHeatCapacAir * PressureAirRefHeight(I,J) / (0.622*LatHeatVapCanopy(I,J))

      ! for ground
      if ( TemperatureGrd(I,J) > ConstFreezePoint ) then
        LatHeatVapGrd(I,J)    = ConstLatHeatEvap
        FlagFrozenGround(I,J) = .false.
      else
        LatHeatVapGrd(I,J)    = ConstLatHeatSublim
        FlagFrozenGround(I,J) = .true.
      endif
      PsychConstGrd(I,J)       = ConstHeatCapacAir * PressureAirRefHeight(I,J) / (0.622*LatHeatVapGrd(I,J))
    else if (noahmp%config%domain%IndicatorIceSfc(I,J) == -1) then ! glacier ice surface
      LatHeatVapGrd(I,J) = ConstLatHeatSublim
      PsychConstGrd(I,J) = ConstHeatCapacAir * PressureAirRefHeight(I,J) / (0.622 * LatHeatVapGrd(I,J))
    endif

      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine PsychrometricVariable

end module PsychrometricVariableMod
