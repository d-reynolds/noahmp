module CanopyHydrologyMod

!!! Canopy Hydrology processes for intercepted rain and snow water
!!! Canopy liquid water evaporation and dew; canopy ice water sublimation and frost
  
  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine CanopyHydrology(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: CANWATER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: I, J         ! grid indices

! --------------------------------------------------------------------
    associate(                                                            &
              MainTimeStep      => noahmp%config%domain%MainTimeStep          ,& ! in,    noahmp main time step [s]
              HeatLatentCanopy  => noahmp%energy%flux%HeatLatentCanopy   ,& ! in,    canopy latent heat flux [W/m2] (+ to atm)
              HeatLatentTransp  => noahmp%energy%flux%HeatLatentTransp   ,& ! in,    latent heat flux from transpiration [W/m2] (+ to atm)
              LeafAreaIndEff    => noahmp%energy%state%LeafAreaIndEff    ,& ! in,    leaf area index, after burying by snow
              StemAreaIndEff    => noahmp%energy%state%StemAreaIndEff    ,& ! in,    stem area index, after burying by snow
              FlagFrozenCanopy  => noahmp%energy%state%FlagFrozenCanopy  ,& ! in,    used to define latent heat pathway
              VegFrac           => noahmp%energy%state%VegFrac           ,& ! in,    greeness vegetation fraction
              SnowfallDensity   => noahmp%water%state%SnowfallDensity    ,& ! in,    bulk density of snowfall [kg/m3]
              CanopyLiqHoldCap  => noahmp%water%param%CanopyLiqHoldCap   ,& ! in,    maximum intercepted liquid water per unit veg area index [mm]
              CanopyLiqWater    => noahmp%water%state%CanopyLiqWater     ,& ! inout, intercepted canopy liquid water [mm]
              CanopyIce         => noahmp%water%state%CanopyIce          ,& ! inout, intercepted canopy ice [mm]
              TemperatureCanopy => noahmp%energy%state%TemperatureCanopy ,& ! inout, vegetation temperature [K]
              CanopyTotalWater  => noahmp%water%state%CanopyTotalWater   ,& ! out,   total canopy intercepted water [mm]
              CanopyWetFrac     => noahmp%water%state%CanopyWetFrac      ,& ! out,   wetted or snowed fraction of the canopy
              CanopyIceMax      => noahmp%water%state%CanopyIceMax       ,& ! out,   canopy capacity for snow interception [mm]
              CanopyLiqWaterMax => noahmp%water%state%CanopyLiqWaterMax  ,& ! out,   canopy capacity for rain interception [mm]
              EvapCanopyNet     => noahmp%water%flux%EvapCanopyNet       ,& ! out,   evaporation of intercepted total water [mm/s]
              Transpiration     => noahmp%water%flux%Transpiration       ,& ! out,   transpiration rate [mm/s]
              EvapCanopyLiq     => noahmp%water%flux%EvapCanopyLiq       ,& ! out,   canopy liquid water evaporation rate [mm/s]
              DewCanopyLiq      => noahmp%water%flux%DewCanopyLiq        ,& ! out,   canopy liquid water dew rate [mm/s]
              FrostCanopyIce    => noahmp%water%flux%FrostCanopyIce      ,& ! out,   canopy ice frost rate [mm/s]
              SublimCanopyIce   => noahmp%water%flux%SublimCanopyIce     ,& ! out,   canopy ice sublimation rate [mm/s]
              MeltCanopyIce     => noahmp%water%flux%MeltCanopyIce       ,& ! out,   canopy ice melting rate [mm/s]
              FreezeCanopyLiq   => noahmp%water%flux%FreezeCanopyLiq      & ! out,   canopy water freezing rate [mm/s]
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! initialization for out-only variables
    EvapCanopyNet(I,J)     = 0.0
    Transpiration(I,J)     = 0.0
    EvapCanopyLiq(I,J)     = 0.0
    DewCanopyLiq(I,J)      = 0.0
    FrostCanopyIce(I,J)    = 0.0
    SublimCanopyIce(I,J)   = 0.0
    MeltCanopyIce(I,J)     = 0.0
    FreezeCanopyLiq(I,J)   = 0.0
    CanopyLiqWaterMax(I,J) = 0.0
    CanopyIceMax(I,J)      = 0.0
    CanopyWetFrac(I,J)     = 0.0
    CanopyTotalWater(I,J)  = 0.0

    ! canopy liquid water
    ! maximum canopy intercepted water
    CanopyLiqWaterMax(I,J) =  VegFrac(I,J) * CanopyLiqHoldCap(I,J) * (LeafAreaIndEff(I,J) + StemAreaIndEff(I,J))

    ! canopy evaporation, transpiration, and dew
    if ( FlagFrozenCanopy(I,J) .eqv. .false. ) then    ! Barlage: change to FlagFrozenCanopy(I,J)
       Transpiration(I,J)   = max( HeatLatentTransp(I,J)/ConstLatHeatEvap, 0.0 )
       EvapCanopyLiq(I,J)   = max( HeatLatentCanopy(I,J)/ConstLatHeatEvap, 0.0 )
       DewCanopyLiq(I,J)    = abs( min( HeatLatentCanopy(I,J)/ConstLatHeatEvap, 0.0 ) )
       SublimCanopyIce(I,J) = 0.0
       FrostCanopyIce(I,J)  = 0.0
    else
       Transpiration(I,J)   = max( HeatLatentTransp(I,J)/ConstLatHeatSublim, 0.0 )
       EvapCanopyLiq(I,J)   = 0.0
       DewCanopyLiq(I,J)    = 0.0
       SublimCanopyIce(I,J) = max( HeatLatentCanopy(I,J)/ConstLatHeatSublim, 0.0 )
       FrostCanopyIce(I,J)  = abs( min( HeatLatentCanopy(I,J)/ConstLatHeatSublim, 0.0 ) )
    endif

    ! canopy water balance. for convenience allow dew to bring CanopyLiqWater above
    ! maxh2o or else would have to re-adjust drip
    EvapCanopyLiq(I,J)  = min( CanopyLiqWater(I,J)/MainTimeStep, EvapCanopyLiq(I,J) )
    CanopyLiqWater(I,J) = max( 0.0, CanopyLiqWater(I,J)+(DewCanopyLiq(I,J)-EvapCanopyLiq(I,J))*MainTimeStep )
    if ( CanopyLiqWater(I,J) <= 1.0e-06 ) CanopyLiqWater(I,J) = 0.0

    ! canopy ice 
    ! maximum canopy intercepted ice
    CanopyIceMax(I,J) = VegFrac(I,J) * 6.6 * (0.27 + 46.0/SnowfallDensity(I,J)) * (LeafAreaIndEff(I,J) + StemAreaIndEff(I,J))

    ! canopy sublimation and frost
    SublimCanopyIce(I,J) = min( CanopyIce(I,J)/MainTimeStep, SublimCanopyIce(I,J) )
    CanopyIce(I,J)       = max( 0.0, CanopyIce(I,J)+(FrostCanopyIce(I,J)-SublimCanopyIce(I,J))*MainTimeStep )
    if ( CanopyIce(I,J) <= 1.0e-6 ) CanopyIce(I,J) = 0.0

    ! wetted fraction of canopy
    if ( (CanopyIce(I,J) > 0.0) .and. (CanopyIce(I,J) >= CanopyLiqWater(I,J)) ) then
       CanopyWetFrac(I,J) = max(0.0,CanopyIce(I,J)) / max(CanopyIceMax(I,J),1.0e-06)
    else
       CanopyWetFrac(I,J) = max(0.0,CanopyLiqWater(I,J)) / max(CanopyLiqWaterMax(I,J),1.0e-06)
    endif
    CanopyWetFrac(I,J)    = min(CanopyWetFrac(I,J), 1.0) ** 0.667
    CanopyTotalWater(I,J) = CanopyLiqWater(I,J) + CanopyIce(I,J)

    ! phase change
    ! canopy ice melting
    if ( (CanopyIce(I,J) > 1.0e-6) .and. (TemperatureCanopy(I,J) > ConstFreezePoint) ) then
       MeltCanopyIce(I,J)     = min( CanopyIce(I,J)/MainTimeStep, (TemperatureCanopy(I,J)-ConstFreezePoint) * ConstHeatCapacIce * &
                                CanopyIce(I,J) / ConstDensityIce / (MainTimeStep*ConstLatHeatFusion) )
       CanopyIce(I,J)         = max( 0.0, CanopyIce(I,J) - MeltCanopyIce(I,J)*MainTimeStep )
       CanopyLiqWater(I,J)    = max( 0.0, CanopyTotalWater(I,J) - CanopyIce(I,J) )
       TemperatureCanopy(I,J) = CanopyWetFrac(I,J)*ConstFreezePoint + (1.0 - CanopyWetFrac(I,J))*TemperatureCanopy(I,J)
    endif

    ! canopy water refreeezing
    if ( (CanopyLiqWater(I,J) > 1.0e-6) .and. (TemperatureCanopy(I,J) < ConstFreezePoint) ) then
       FreezeCanopyLiq(I,J)   = min( CanopyLiqWater(I,J)/MainTimeStep, (ConstFreezePoint-TemperatureCanopy(I,J)) * ConstHeatCapacWater * &
                                CanopyLiqWater(I,J) / ConstDensityWater / (MainTimeStep*ConstLatHeatFusion) )
       CanopyLiqWater(I,J)    = max( 0.0, CanopyLiqWater(I,J) - FreezeCanopyLiq(I,J)*MainTimeStep )
       CanopyIce(I,J)         = max( 0.0, CanopyTotalWater(I,J) - CanopyLiqWater(I,J) )
       TemperatureCanopy(I,J) = CanopyWetFrac(I,J)*ConstFreezePoint + (1.0 - CanopyWetFrac(I,J))*TemperatureCanopy(I,J)
    endif

    ! update total canopy water
    CanopyTotalWater(I,J) = CanopyLiqWater(I,J) + CanopyIce(I,J)

    ! total canopy net evaporation
    EvapCanopyNet(I,J)    = EvapCanopyLiq(I,J) + SublimCanopyIce(I,J) - DewCanopyLiq(I,J) - FrostCanopyIce(I,J)


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine CanopyHydrology

end module CanopyHydrologyMod
