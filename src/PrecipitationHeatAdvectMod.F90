module PrecipitationHeatAdvectMod

!!! Estimate heat flux advected from precipitation to vegetation and ground (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine PrecipitationHeatAdvect(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: PRECIP_HEAT
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! The water and heat portions of PRECIP_HEAT are separated in refactored code
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: I, J                ! grid indices
    real(kind=kind_noahmp)           :: HeatPrcpAirToCan    ! precipitation advected heat - air to canopy [W/m2]
    real(kind=kind_noahmp)           :: HeatPrcpCanToGrd    ! precipitation advected heat - canopy to ground [W/m2]
    real(kind=kind_noahmp)           :: HeatPrcpAirToGrd    ! precipitation advected heat - air to ground [W/m2]

! --------------------------------------------------------------------
        associate(                                                                          &
                  TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight  ,& ! in,  air temperature [K] at reference height
                  TemperatureCanopy       => noahmp%energy%state%TemperatureCanopy   ,& ! in,  vegetation temperature [K]
                  TemperatureGrd          => noahmp%energy%state%TemperatureGrd      ,& ! in,  ground temperature [K]
                  VegFrac                 => noahmp%energy%state%VegFrac             ,& ! in,  greeness vegetation fraction
                  RainfallRefHeight       => noahmp%water%flux%RainfallRefHeight     ,& ! in,  total liquid rainfall [mm/s] before interception
                  SnowfallRefHeight       => noahmp%water%flux%SnowfallRefHeight     ,& ! in,  total snowfall [mm/s] before interception
                  DripCanopyRain          => noahmp%water%flux%DripCanopyRain        ,& ! in,  drip rate for intercepted rain [mm/s]
                  ThroughfallRain         => noahmp%water%flux%ThroughfallRain       ,& ! in,  throughfall for rain [mm/s]
                  DripCanopySnow          => noahmp%water%flux%DripCanopySnow        ,& ! in,  drip (unloading) rate for intercepted snow [mm/s]
                  ThroughfallSnow         => noahmp%water%flux%ThroughfallSnow       ,& ! in,  throughfall of snowfall [mm/s]
                  SnowfallGround          => noahmp%water%flux%SnowfallGround       ,& ! out, snowfall at ground surface [mm/s]
                  RainfallGround          => noahmp%water%flux%RainfallGround       ,& ! out, rainfall at ground surface [mm/s]
                  HeatPrecipAdvCanopy     => noahmp%energy%flux%HeatPrecipAdvCanopy  ,& ! out, precipitation advected heat - vegetation net [W/m2]
                  HeatPrecipAdvVegGrd     => noahmp%energy%flux%HeatPrecipAdvVegGrd  ,& ! out, precipitation advected heat - under canopy net [W/m2]
                  HeatPrecipAdvBareGrd    => noahmp%energy%flux%HeatPrecipAdvBareGrd  & ! out, precipitation advected heat - bare ground net [W/m2]
                 )

    !$acc parallel loop collapse(2) gang vector default(present) &
    !$acc private(HeatPrcpAirToCan, HeatPrcpCanToGrd, HeatPrcpAirToGrd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

  if (noahmp%config%domain%IndicatorIceSfc(I,J) == 0) then

    ! initialization
    HeatPrcpAirToCan     = 0.0
    HeatPrcpCanToGrd     = 0.0
    HeatPrcpAirToGrd     = 0.0
    HeatPrecipAdvCanopy(I,J)  = 0.0
    HeatPrecipAdvVegGrd(I,J)  = 0.0
    HeatPrecipAdvBareGrd(I,J) = 0.0

    ! Heat advection for liquid rainfall
    HeatPrcpAirToCan = VegFrac(I,J) * RainfallRefHeight(I,J) * (ConstHeatCapacWater/1000.0) * (TemperatureAirRefHeight(I,J)-TemperatureCanopy(I,J))
    HeatPrcpCanToGrd = DripCanopyRain(I,J) * (ConstHeatCapacWater/1000.0) * (TemperatureCanopy(I,J)-TemperatureGrd(I,J))
    HeatPrcpAirToGrd = ThroughfallRain(I,J) * (ConstHeatCapacWater/1000.0) * (TemperatureAirRefHeight(I,J)-TemperatureGrd(I,J))

    ! Heat advection for snowfall
    HeatPrcpAirToCan = HeatPrcpAirToCan + &
                       VegFrac(I,J) * SnowfallRefHeight(I,J) * (ConstHeatCapacIce/1000.0) * (TemperatureAirRefHeight(I,J)-TemperatureCanopy(I,J))
    HeatPrcpCanToGrd = HeatPrcpCanToGrd + &
                       DripCanopySnow(I,J) * (ConstHeatCapacIce/1000.0) * (TemperatureCanopy(I,J)-TemperatureGrd(I,J))
    HeatPrcpAirToGrd = HeatPrcpAirToGrd + &
                       ThroughfallSnow(I,J) * (ConstHeatCapacIce/1000.0) * (TemperatureAirRefHeight(I,J)-TemperatureGrd(I,J))

    ! net heat advection
    HeatPrecipAdvCanopy(I,J)  = HeatPrcpAirToCan - HeatPrcpCanToGrd
    HeatPrecipAdvVegGrd(I,J)  = HeatPrcpCanToGrd
    HeatPrecipAdvBareGrd(I,J) = HeatPrcpAirToGrd

    ! adjust for VegFrac
    if ( (VegFrac(I,J) > 0.0) .and. (VegFrac(I,J) < 1.0) ) then
       HeatPrecipAdvVegGrd(I,J)  = HeatPrecipAdvVegGrd(I,J) / VegFrac(I,J)                  ! these will be multiplied by fraction later
       HeatPrecipAdvBareGrd(I,J) = HeatPrecipAdvBareGrd(I,J) / (1.0-VegFrac(I,J))
    elseif ( VegFrac(I,J) <= 0.0 ) then
       HeatPrecipAdvBareGrd(I,J) = HeatPrecipAdvVegGrd(I,J) + HeatPrecipAdvBareGrd(I,J)     ! for case of canopy getting buried
       HeatPrecipAdvVegGrd(I,J)  = 0.0
       HeatPrecipAdvCanopy(I,J)  = 0.0
    elseif ( VegFrac(I,J) >= 1.0 ) then
       HeatPrecipAdvBareGrd(I,J) = 0.0
    endif

    ! Put some artificial limits here for stability
    HeatPrecipAdvCanopy(I,J)  = max(HeatPrecipAdvCanopy(I,J) , -20.0)
    HeatPrecipAdvCanopy(I,J)  = min(HeatPrecipAdvCanopy(I,J) ,  20.0)
    HeatPrecipAdvVegGrd(I,J)  = max(HeatPrecipAdvVegGrd(I,J) , -20.0)
    HeatPrecipAdvVegGrd(I,J)  = min(HeatPrecipAdvVegGrd(I,J) ,  20.0)
    HeatPrecipAdvBareGrd(I,J) = max(HeatPrecipAdvBareGrd(I,J), -20.0)
    HeatPrecipAdvBareGrd(I,J) = min(HeatPrecipAdvBareGrd(I,J),  20.0)

  else if (noahmp%config%domain%IndicatorIceSfc(I,J) == -1) then

    ! initialization for glacier points
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

  endif

      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine PrecipitationHeatAdvect

end module PrecipitationHeatAdvectMod
