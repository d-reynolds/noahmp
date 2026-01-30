module CanopyWaterInterceptMod

!!! Canopy water processes for snow and rain interception (2D GPU-optimized)
!!! Subsequent hydrological process for intercepted water is done in CanopyHydrologyMod.F90

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine CanopyWaterIntercept(noahmp)

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
    integer                          :: I, J                   ! grid indices
    real(kind=kind_noahmp)           :: IceDripFacTemp         ! temperature factor for unloading rate
    real(kind=kind_noahmp)           :: IceDripFacWind         ! wind factor for unloading rate
    real(kind=kind_noahmp)           :: CanopySnowDrip         ! canopy snow/ice unloading 

! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp) &
    !$acc private(IceDripFacTemp, IceDripFacWind, CanopySnowDrip)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                       &
                  SurfaceType            => noahmp%config%domain%SurfaceType(I,J)      ,& ! in,    surface type 1-soil; 2-lake
                  MainTimeStep           => noahmp%config%domain%MainTimeStep          ,& ! in,    noahmp main time step [s]
                  WindEastwardRefHeight  => noahmp%forcing%WindEastwardRefHeight(I,J)  ,& ! in,    wind speed [m/s] in eastward direction at reference height
                  WindNorthwardRefHeight => noahmp%forcing%WindNorthwardRefHeight(I,J) ,& ! in,    wind speed [m/s] in northward direction at reference height
                  LeafAreaIndEff         => noahmp%energy%state%LeafAreaIndEff(I,J)    ,& ! in,    leaf area index, after burying by snow
                  StemAreaIndEff         => noahmp%energy%state%StemAreaIndEff(I,J)    ,& ! in,    stem area index, after burying by snow
                  VegFrac                => noahmp%energy%state%VegFrac(I,J)           ,& ! in,    greeness vegetation fraction
                  TemperatureCanopy      => noahmp%energy%state%TemperatureCanopy(I,J) ,& ! in,    vegetation temperature [K]
                  TemperatureGrd         => noahmp%energy%state%TemperatureGrd(I,J)    ,& ! in,    ground temperature [K]
                  CanopyLiqHoldCap       => noahmp%water%param%CanopyLiqHoldCap(I,J)   ,& ! in,    maximum intercepted liquid water per unit veg area index [mm]
                  RainfallRefHeight      => noahmp%water%flux%RainfallRefHeight(I,J)   ,& ! in,    total liquid rainfall [mm/s] before interception
                  SnowfallRefHeight      => noahmp%water%flux%SnowfallRefHeight(I,J)   ,& ! in,    total snowfall [mm/s] before interception
                  SnowfallDensity        => noahmp%water%state%SnowfallDensity(I,J)    ,& ! in,    bulk density of snowfall [kg/m3]
                  PrecipAreaFrac         => noahmp%water%state%PrecipAreaFrac(I,J)     ,& ! in,    fraction of the gridcell that receives precipitation
                  CanopyLiqWater         => noahmp%water%state%CanopyLiqWater(I,J)     ,& ! inout, intercepted canopy liquid water [mm]
                  CanopyIce              => noahmp%water%state%CanopyIce(I,J)          ,& ! inout, intercepted canopy ice [mm]
                  CanopyWetFrac          => noahmp%water%state%CanopyWetFrac(I,J)      ,& ! out,   wetted or snowed fraction of the canopy
                  CanopyTotalWater       => noahmp%water%state%CanopyTotalWater(I,J)   ,& ! out,   total canopy intercepted water [mm]
                  CanopyIceMax           => noahmp%water%state%CanopyIceMax(I,J)       ,& ! out,   canopy capacity for snow interception [mm]
                  CanopyLiqWaterMax      => noahmp%water%state%CanopyLiqWaterMax(I,J)  ,& ! out,   canopy capacity for rain interception [mm]
                  InterceptCanopyRain    => noahmp%water%flux%InterceptCanopyRain(I,J) ,& ! out,   interception rate for rain [mm/s]
                  DripCanopyRain         => noahmp%water%flux%DripCanopyRain(I,J)      ,& ! out,   drip rate for intercepted rain [mm/s]
                  ThroughfallRain        => noahmp%water%flux%ThroughfallRain(I,J)     ,& ! out,   throughfall for rain [mm/s]
                  InterceptCanopySnow    => noahmp%water%flux%InterceptCanopySnow(I,J) ,& ! out,   interception (loading) rate for snowfall [mm/s]
                  DripCanopySnow         => noahmp%water%flux%DripCanopySnow(I,J)      ,& ! out,   drip (unloading) rate for intercepted snow [mm/s]
                  ThroughfallSnow        => noahmp%water%flux%ThroughfallSnow(I,J)     ,& ! out,   throughfall of snowfall [mm/s]
                  RainfallGround         => noahmp%water%flux%RainfallGround(I,J)      ,& ! out,   rainfall at ground surface [mm/s]
                  SnowfallGround         => noahmp%water%flux%SnowfallGround(I,J)      ,& ! out,   snowfall at ground surface [mm/s]
                  SnowDepthIncr          => noahmp%water%flux%SnowDepthIncr(I,J)        & ! out,   snow depth increasing rate [m/s] due to snowfall
                 )
! ----------------------------------------------------------------------

    ! initialization
    InterceptCanopyRain = 0.0
    DripCanopyRain      = 0.0
    ThroughfallRain     = 0.0
    InterceptCanopySnow = 0.0
    DripCanopySnow      = 0.0
    ThroughfallSnow     = 0.0
    RainfallGround      = 0.0
    SnowfallGround      = 0.0
    SnowDepthIncr       = 0.0
    CanopySnowDrip      = 0.0
    IceDripFacTemp      = 0.0
    IceDripFacWind      = 0.0

    ! ----------------------- canopy liquid water ------------------------------
    ! maximum canopy water
    CanopyLiqWaterMax =  VegFrac * CanopyLiqHoldCap * (LeafAreaIndEff + StemAreaIndEff)

    ! average rain interception and throughfall
    if ( (LeafAreaIndEff+StemAreaIndEff) > 0.0 ) then
       InterceptCanopyRain = VegFrac * RainfallRefHeight * PrecipAreaFrac  ! max interception capability
       InterceptCanopyRain = min( InterceptCanopyRain, (CanopyLiqWaterMax-CanopyLiqWater)/MainTimeStep * &
                                  (1.0-exp(-RainfallRefHeight*MainTimeStep/CanopyLiqWaterMax)) )
       InterceptCanopyRain = max( InterceptCanopyRain, 0.0 )
       DripCanopyRain      = VegFrac * RainfallRefHeight - InterceptCanopyRain
       ThroughfallRain     = (1.0 - VegFrac) * RainfallRefHeight
       CanopyLiqWater      = max( 0.0, CanopyLiqWater + InterceptCanopyRain*MainTimeStep )
    else
       InterceptCanopyRain = 0.0
       DripCanopyRain      = 0.0
       ThroughfallRain     = RainfallRefHeight
       if ( CanopyLiqWater > 0.0 ) then   ! canopy gets buried by rain
          DripCanopyRain   = DripCanopyRain + CanopyLiqWater / MainTimeStep
          CanopyLiqWater   = 0.0
       endif
    endif

    ! ----------------------- canopy ice ------------------------------
    ! maximum canopy ice
    CanopyIceMax = VegFrac * 6.6 * (0.27 + 46.0/SnowfallDensity) * (LeafAreaIndEff + StemAreaIndEff)

    ! average snow interception and throughfall
    if ( (LeafAreaIndEff+StemAreaIndEff) > 0.0 ) then
       InterceptCanopySnow = VegFrac * SnowfallRefHeight * PrecipAreaFrac
       InterceptCanopySnow = min( InterceptCanopySnow, (CanopyIceMax-CanopyIce)/MainTimeStep * &
                                  (1.0-exp(-SnowfallRefHeight*MainTimeStep/CanopyIceMax)) )
       InterceptCanopySnow = max( InterceptCanopySnow, 0.0 )
       IceDripFacTemp      = max( 0.0, (TemperatureCanopy - 270.15) / 1.87e5 )
       IceDripFacWind      = sqrt(WindEastwardRefHeight**2.0 + WindNorthwardRefHeight**2.0) / 1.56e5
       ! MB: changed below to reflect the rain assumption that all precip gets intercepted 
       CanopySnowDrip      = max( 0.0, CanopyIce ) * (IceDripFacWind + IceDripFacTemp)
       CanopySnowDrip      = min( CanopyIce/MainTimeStep + InterceptCanopySnow, CanopySnowDrip) ! add constraint to keep water balance
       DripCanopySnow      = (VegFrac * SnowfallRefHeight - InterceptCanopySnow) + CanopySnowDrip
       ThroughfallSnow     = (1.0 - VegFrac) * SnowfallRefHeight
       CanopyIce           = max( 0.0, CanopyIce + (InterceptCanopySnow-CanopySnowDrip)*MainTimeStep )
    else
       InterceptCanopySnow = 0.0
       DripCanopySnow      = 0.0
       ThroughfallSnow     = SnowfallRefHeight
       if ( CanopyIce > 0.0 ) then   ! canopy gets buried by snow
          DripCanopySnow   = DripCanopySnow + CanopyIce / MainTimeStep
          CanopyIce        = 0.0
       endif
    endif

    ! wetted fraction of canopy
    if ( CanopyIce > 0.0 ) then
       CanopyWetFrac  = max( 0.0, CanopyIce ) / max( CanopyIceMax, 1.0e-06 )
    else
       CanopyWetFrac  = max( 0.0, CanopyLiqWater ) / max( CanopyLiqWaterMax, 1.0e-06 )
    endif
    CanopyWetFrac     = min( CanopyWetFrac, 1.0 ) ** 0.667

    ! total canopy water
    CanopyTotalWater  = CanopyLiqWater + CanopyIce

    ! rain or snow on the ground
    RainfallGround    = DripCanopyRain + ThroughfallRain
    SnowfallGround    = DripCanopySnow + ThroughfallSnow
    SnowDepthIncr     = SnowfallGround / SnowfallDensity
    if ( (SurfaceType == 2) .and. (TemperatureGrd > ConstFreezePoint) ) then
       SnowfallGround = 0.0
       SnowDepthIncr  = 0.0
    endif

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine CanopyWaterIntercept

end module CanopyWaterInterceptMod
