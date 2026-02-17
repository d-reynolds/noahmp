module PrecipitationHeatAdvectGlacierMod

!!! Estimate heat flux advected from precipitation to glacier ground

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine PrecipitationHeatAdvectGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: none (adapted from PRECIP_HEAT)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                ! grid indices
    real(kind=kind_noahmp)           :: HeatPrcpAirToGrd    ! precipitation advected heat - air to ground [W/m2]

! --------------------------------------------------------------------
    associate(                                                                       &
              TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight ,& ! in,  air temperature [K] at reference height
              TemperatureGrd          => noahmp%energy%state%TemperatureGrd     ,& ! in,  ground temperature [K]
              RainfallRefHeight       => noahmp%water%flux%RainfallRefHeight    ,& ! in,  total liquid rainfall [mm/s] before interception
              SnowfallRefHeight       => noahmp%water%flux%SnowfallRefHeight    ,& ! in,  total snowfall [mm/s] before interception
              SnowfallGround          => noahmp%water%flux%SnowfallGround       ,& ! out, snowfall at ground surface [mm/s]
              RainfallGround          => noahmp%water%flux%RainfallGround       ,& ! out, rainfall at ground surface [mm/s]
              HeatPrecipAdvBareGrd    => noahmp%energy%flux%HeatPrecipAdvBareGrd & ! out, precipitation advected heat - bare ground net [W/m2]
             )

   !$acc parallel loop collapse(2) gang vector default(present) private(HeatPrcpAirToGrd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! initialization
    HeatPrcpAirToGrd     = 0.0
    HeatPrecipAdvBareGrd(I,J) = 0.0
    RainfallGround(I,J)       = RainfallRefHeight(I,J)
    SnowfallGround(I,J)       = SnowfallRefHeight(I,J)

    ! Heat advection for liquid rainfall
    HeatPrcpAirToGrd     = RainfallGround(I,J) * (ConstHeatCapacWater/1000.0) * (TemperatureAirRefHeight(I,J) - TemperatureGrd(I,J))

    ! Heat advection for snowfall
    HeatPrcpAirToGrd     = HeatPrcpAirToGrd + &
                           SnowfallGround(I,J) * (ConstHeatCapacIce/1000.0) * (TemperatureAirRefHeight(I,J) - TemperatureGrd(I,J))

    ! net heat advection
    HeatPrecipAdvBareGrd(I,J) = HeatPrcpAirToGrd

    ! Put some artificial limits here for stability
    HeatPrecipAdvBareGrd(I,J) = max(HeatPrecipAdvBareGrd(I,J), -20.0)
    HeatPrecipAdvBareGrd(I,J) = min(HeatPrecipAdvBareGrd(I,J),  20.0)


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine PrecipitationHeatAdvectGlacier

end module PrecipitationHeatAdvectGlacierMod
