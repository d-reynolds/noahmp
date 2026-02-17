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
        associate(                                                                       &
                  SurfaceType            => noahmp%config%domain%SurfaceType      ,& ! in,    surface type 1-soil; 2-lake
                  MainTimeStep           => noahmp%config%domain%MainTimeStep          ,& ! in,    noahmp main time step [s]
                  WindEastwardRefHeight  => noahmp%forcing%WindEastwardRefHeight  ,& ! in,    wind speed [m/s] in eastward direction at reference height
                  WindNorthwardRefHeight => noahmp%forcing%WindNorthwardRefHeight ,& ! in,    wind speed [m/s] in northward direction at reference height
                  LeafAreaIndEff         => noahmp%energy%state%LeafAreaIndEff    ,& ! in,    leaf area index, after burying by snow
                  StemAreaIndEff         => noahmp%energy%state%StemAreaIndEff    ,& ! in,    stem area index, after burying by snow
                  VegFrac                => noahmp%energy%state%VegFrac           ,& ! in,    greeness vegetation fraction
                  TemperatureCanopy      => noahmp%energy%state%TemperatureCanopy ,& ! in,    vegetation temperature [K]
                  TemperatureGrd         => noahmp%energy%state%TemperatureGrd    ,& ! in,    ground temperature [K]
                  CanopyLiqHoldCap       => noahmp%water%param%CanopyLiqHoldCap   ,& ! in,    maximum intercepted liquid water per unit veg area index [mm]
                  RainfallRefHeight      => noahmp%water%flux%RainfallRefHeight   ,& ! in,    total liquid rainfall [mm/s] before interception
                  SnowfallRefHeight      => noahmp%water%flux%SnowfallRefHeight   ,& ! in,    total snowfall [mm/s] before interception
                  SnowfallDensity        => noahmp%water%state%SnowfallDensity    ,& ! in,    bulk density of snowfall [kg/m3]
                  PrecipAreaFrac         => noahmp%water%state%PrecipAreaFrac     ,& ! in,    fraction of the gridcell that receives precipitation
                  CanopyLiqWater         => noahmp%water%state%CanopyLiqWater     ,& ! inout, intercepted canopy liquid water [mm]
                  CanopyIce              => noahmp%water%state%CanopyIce          ,& ! inout, intercepted canopy ice [mm]
                  CanopyWetFrac          => noahmp%water%state%CanopyWetFrac      ,& ! out,   wetted or snowed fraction of the canopy
                  CanopyTotalWater       => noahmp%water%state%CanopyTotalWater   ,& ! out,   total canopy intercepted water [mm]
                  CanopyIceMax           => noahmp%water%state%CanopyIceMax       ,& ! out,   canopy capacity for snow interception [mm]
                  CanopyLiqWaterMax      => noahmp%water%state%CanopyLiqWaterMax  ,& ! out,   canopy capacity for rain interception [mm]
                  InterceptCanopyRain    => noahmp%water%flux%InterceptCanopyRain ,& ! out,   interception rate for rain [mm/s]
                  DripCanopyRain         => noahmp%water%flux%DripCanopyRain      ,& ! out,   drip rate for intercepted rain [mm/s]
                  ThroughfallRain        => noahmp%water%flux%ThroughfallRain     ,& ! out,   throughfall for rain [mm/s]
                  InterceptCanopySnow    => noahmp%water%flux%InterceptCanopySnow ,& ! out,   interception (loading) rate for snowfall [mm/s]
                  DripCanopySnow         => noahmp%water%flux%DripCanopySnow      ,& ! out,   drip (unloading) rate for intercepted snow [mm/s]
                  ThroughfallSnow        => noahmp%water%flux%ThroughfallSnow     ,& ! out,   throughfall of snowfall [mm/s]
                  RainfallGround         => noahmp%water%flux%RainfallGround      ,& ! out,   rainfall at ground surface [mm/s]
                  SnowfallGround         => noahmp%water%flux%SnowfallGround      ,& ! out,   snowfall at ground surface [mm/s]
                  SnowDepthIncr          => noahmp%water%flux%SnowDepthIncr        & ! out,   snow depth increasing rate [m/s] due to snowfall
                 )

    !$acc parallel loop collapse(2) gang vector default(present) &
    !$acc private(IceDripFacTemp, IceDripFacWind, CanopySnowDrip)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! initialization
    InterceptCanopyRain(I,J) = 0.0
    DripCanopyRain(I,J)      = 0.0
    ThroughfallRain(I,J)     = 0.0
    InterceptCanopySnow(I,J) = 0.0
    DripCanopySnow(I,J)      = 0.0
    ThroughfallSnow(I,J)     = 0.0
    RainfallGround(I,J)      = 0.0
    SnowfallGround(I,J)      = 0.0
    SnowDepthIncr(I,J)       = 0.0
    CanopySnowDrip      = 0.0
    IceDripFacTemp      = 0.0
    IceDripFacWind      = 0.0

    ! ----------------------- canopy liquid water ------------------------------
    ! maximum canopy water
    CanopyLiqWaterMax(I,J) =  VegFrac(I,J) * CanopyLiqHoldCap(I,J) * (LeafAreaIndEff(I,J) + StemAreaIndEff(I,J))

    ! average rain interception and throughfall
    if ( (LeafAreaIndEff(I,J)+StemAreaIndEff(I,J)) > 0.0 ) then
       InterceptCanopyRain(I,J) = VegFrac(I,J) * RainfallRefHeight(I,J) * PrecipAreaFrac(I,J)  ! max interception capability
       InterceptCanopyRain(I,J) = min( InterceptCanopyRain(I,J), (CanopyLiqWaterMax(I,J)-CanopyLiqWater(I,J))/MainTimeStep * &
                                  (1.0-exp(-RainfallRefHeight(I,J)*MainTimeStep/CanopyLiqWaterMax(I,J))) )
       InterceptCanopyRain(I,J) = max( InterceptCanopyRain(I,J), 0.0 )
       DripCanopyRain(I,J)      = VegFrac(I,J) * RainfallRefHeight(I,J) - InterceptCanopyRain(I,J)
       ThroughfallRain(I,J)     = (1.0 - VegFrac(I,J)) * RainfallRefHeight(I,J)
       CanopyLiqWater(I,J)      = max( 0.0, CanopyLiqWater(I,J) + InterceptCanopyRain(I,J)*MainTimeStep )
    else
       InterceptCanopyRain(I,J) = 0.0
       DripCanopyRain(I,J)      = 0.0
       ThroughfallRain(I,J)     = RainfallRefHeight(I,J)
       if ( CanopyLiqWater(I,J) > 0.0 ) then   ! canopy gets buried by rain
          DripCanopyRain(I,J)   = DripCanopyRain(I,J) + CanopyLiqWater(I,J) / MainTimeStep
          CanopyLiqWater(I,J)   = 0.0
       endif
    endif

    ! ----------------------- canopy ice ------------------------------
    ! maximum canopy ice
    CanopyIceMax(I,J) = VegFrac(I,J) * 6.6 * (0.27 + 46.0/SnowfallDensity(I,J)) * (LeafAreaIndEff(I,J) + StemAreaIndEff(I,J))

    ! average snow interception and throughfall
    if ( (LeafAreaIndEff(I,J)+StemAreaIndEff(I,J)) > 0.0 ) then
       InterceptCanopySnow(I,J) = VegFrac(I,J) * SnowfallRefHeight(I,J) * PrecipAreaFrac(I,J)
       InterceptCanopySnow(I,J) = min( InterceptCanopySnow(I,J), (CanopyIceMax(I,J)-CanopyIce(I,J))/MainTimeStep * &
                                  (1.0-exp(-SnowfallRefHeight(I,J)*MainTimeStep/CanopyIceMax(I,J))) )
       InterceptCanopySnow(I,J) = max( InterceptCanopySnow(I,J), 0.0 )
       IceDripFacTemp      = max( 0.0, (TemperatureCanopy(I,J) - 270.15) / 1.87e5 )
       IceDripFacWind      = sqrt(WindEastwardRefHeight(I,J)**2.0 + WindNorthwardRefHeight(I,J)**2.0) / 1.56e5
       ! MB: changed below to reflect the rain assumption that all precip gets intercepted 
       CanopySnowDrip      = max( 0.0, CanopyIce(I,J) ) * (IceDripFacWind + IceDripFacTemp)
       CanopySnowDrip      = min( CanopyIce(I,J)/MainTimeStep + InterceptCanopySnow(I,J), CanopySnowDrip) ! add constraint to keep water balance
       DripCanopySnow(I,J)      = (VegFrac(I,J) * SnowfallRefHeight(I,J) - InterceptCanopySnow(I,J)) + CanopySnowDrip
       ThroughfallSnow(I,J)     = (1.0 - VegFrac(I,J)) * SnowfallRefHeight(I,J)
       CanopyIce(I,J)           = max( 0.0, CanopyIce(I,J) + (InterceptCanopySnow(I,J)-CanopySnowDrip)*MainTimeStep )
    else
       InterceptCanopySnow(I,J) = 0.0
       DripCanopySnow(I,J)      = 0.0
       ThroughfallSnow(I,J)     = SnowfallRefHeight(I,J)
       if ( CanopyIce(I,J) > 0.0 ) then   ! canopy gets buried by snow
          DripCanopySnow(I,J)   = DripCanopySnow(I,J) + CanopyIce(I,J) / MainTimeStep
          CanopyIce(I,J)        = 0.0
       endif
    endif

    ! wetted fraction of canopy
    if ( CanopyIce(I,J) > 0.0 ) then
       CanopyWetFrac(I,J)  = max( 0.0, CanopyIce(I,J) ) / max( CanopyIceMax(I,J), 1.0e-06 )
    else
       CanopyWetFrac(I,J)  = max( 0.0, CanopyLiqWater(I,J) ) / max( CanopyLiqWaterMax(I,J), 1.0e-06 )
    endif
    CanopyWetFrac(I,J)     = min( CanopyWetFrac(I,J), 1.0 ) ** 0.667

    ! total canopy water
    CanopyTotalWater(I,J)  = CanopyLiqWater(I,J) + CanopyIce(I,J)

    ! rain or snow on the ground
    RainfallGround(I,J)    = DripCanopyRain(I,J) + ThroughfallRain(I,J)
    SnowfallGround(I,J)    = DripCanopySnow(I,J) + ThroughfallSnow(I,J)
    SnowDepthIncr(I,J)     = SnowfallGround(I,J) / SnowfallDensity(I,J)
    if ( (SurfaceType(I,J) == 2) .and. (TemperatureGrd(I,J) > ConstFreezePoint) ) then
       SnowfallGround(I,J) = 0.0
       SnowDepthIncr(I,J)  = 0.0
    endif


      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine CanopyWaterIntercept

end module CanopyWaterInterceptMod
