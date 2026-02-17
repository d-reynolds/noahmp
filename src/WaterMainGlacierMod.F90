module WaterMainGlacierMod

!!! Main glacier water module including all water relevant processes
!!! snowpack water -> ice water -> runoff

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowWaterMainGlacierMod, only : SnowWaterMainGlacier

  implicit none

contains

  subroutine WaterMainGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: WATER_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J        ! grid indices
    integer                          :: LoopInd     ! loop index
    real(kind=kind_noahmp)           :: WatReplaceSublim  ! replacement water due to sublimation of glacier
    real(kind=kind_noahmp)  :: SoilIceTmp(1:noahmp%config%domain%NumSoilLayer)       ! temporary glacier ice content [m3/m3]
    real(kind=kind_noahmp)  :: SoilLiqWaterTmp(1:noahmp%config%domain%NumSoilLayer)  ! temporary glacier liquid water content [m3/m3]

    ! glacier ice water processing with OpenACC
    !$acc data create(SoilIceTmp, SoilLiqWaterTmp)
    associate(                                                                      &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer ,& ! in,    number of soil layers
              VaporizeGrd            => noahmp%water%flux%VaporizeGrd ,& ! inout, ground vaporize rate total (evap+sublim) [mm/s]
              CondenseVapGrd         => noahmp%water%flux%CondenseVapGrd ,& ! inout, ground vapor condense rate total (dew+frost) [mm/s]
              SnowfallGround         => noahmp%water%flux%SnowfallGround ,& ! in,    snowfall on the ground [mm/s]
              SnowfallDensity        => noahmp%water%state%SnowfallDensity ,& ! in,    bulk density of snowfall [kg/m3]
              LatHeatVapGrd          => noahmp%energy%state%LatHeatVapGrd ,& ! in,    latent heat of vaporization/subli [J/kg], ground
              HeatLatentGrd          => noahmp%energy%flux%HeatLatentGrd ,& ! inout, total ground latent heat [W/m2] (+ to atm)
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv ,& ! inout, snow water equivalent [mm]
              SnowWaterEquivPrev     => noahmp%water%state%SnowWaterEquivPrev ,& ! inout, snow water equivalent at last time step [mm]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater ,& ! inout, glacier water content [m3/m3]
              SoilIce                => noahmp%water%state%SoilIce ,& ! inout, glacier ice moisture [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture ,& ! inout, total glacier water [m3/m3]
              FrostSnowSfcIce        => noahmp%water%flux%FrostSnowSfcIce ,& ! inout, snow surface frost rate [mm/s]
              SublimSnowSfcIce       => noahmp%water%flux%SublimSnowSfcIce ,& ! inout, snow surface sublimation rate [mm/s]
              GlacierExcessFlow      => noahmp%water%flux%GlacierExcessFlow ,& ! inout, glacier snow excess flow [mm/s]
              SnowDepthIncr          => noahmp%water%flux%SnowDepthIncr ,& ! out,   snow depth increasing rate [m/s] due to snowfall
              EvapGroundNet          => noahmp%water%flux%EvapGroundNet ,& ! out,   net direct ground evaporation [mm/s]
              RunoffSurface          => noahmp%water%flux%RunoffSurface ,& ! out,   surface runoff [mm/s]
              RunoffSubsurface       => noahmp%water%flux%RunoffSubsurface ,& ! out,   subsurface runoff [mm/s]
              OptGlacierTreatment    => noahmp%config%nmlist%OptGlacierTreatment ,& ! in,    option for glacier treatment
              MainTimeStep           => noahmp%config%domain%MainTimeStep ,& ! in,    noahmp main time step [s]
              RainfallGround         => noahmp%water%flux%RainfallGround ,& ! in,    ground surface rain rate [mm/s]
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg ,& ! inout, actual number of snow layers (negative)
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/glacier layers [m]
              PondSfcThinSnwMelt     => noahmp%water%state%PondSfcThinSnwMelt ,& ! inout, surface ponding [mm] from snowmelt when thin snow has no layer
              WaterHeadSfc           => noahmp%water%state%WaterHeadSfc ,& ! inout, surface water head [mm)]
              SoilSfcInflow          => noahmp%water%flux%SoilSfcInflow ,& ! inout, water input on glacier/soil surface [m/s]
              SnowBotOutflow         => noahmp%water%flux%SnowBotOutflow ,& ! out,   total water (snowmelt + rain through pack) out of snowpack bottom [mm/s]
              PondSfcThinSnwComb     => noahmp%water%state%PondSfcThinSnwComb ,& ! out,   surface ponding [mm] from liquid in thin snow layer combination
              PondSfcThinSnwTrans    => noahmp%water%state%PondSfcThinSnwTrans  & ! out,   surface ponding [mm] from thin snow liquid during transition from multilayer to no layer
             )

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS,noahmp%config%domain%ITE
    SoilIceTmp         = 0.0
    SoilLiqWaterTmp    = 0.0
    GlacierExcessFlow(I,J)  = 0.0
    RunoffSubsurface(I,J)   = 0.0
    RunoffSurface(I,J)      = 0.0
    SnowDepthIncr(I,J)      = 0.0

    ! prepare for water process
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
      SoilIce(I,LoopInd,J)         = max(0.0, SoilMoisture(I,LoopInd,J)-SoilLiqWater(I,LoopInd,J))
      SoilIceTmp(LoopInd)      = SoilIce(I,LoopInd,J)
      SoilLiqWaterTmp(LoopInd) = SoilLiqWater(I,LoopInd,J)
    enddo
    SnowWaterEquivPrev(I,J) = SnowWaterEquiv(I,J)

    ! compute soil/snow surface evap/dew rate based on energy flux
    VaporizeGrd(I,J)        = max(HeatLatentGrd(I,J)/LatHeatVapGrd(I,J), 0.0)       ! positive part of ground latent heat; Barlage change to ground v3.6
    CondenseVapGrd(I,J)     = abs(min(HeatLatentGrd(I,J)/LatHeatVapGrd(I,J), 0.0))  ! negative part of ground latent heat
    EvapGroundNet(I,J)      = VaporizeGrd(I,J) - CondenseVapGrd(I,J)

    ! snow height increase
    SnowDepthIncr(I,J)      = SnowfallGround(I,J) / SnowfallDensity(I,J)

    ! ground sublimation and evaporation
    SublimSnowSfcIce(I,J)   = VaporizeGrd(I,J)

    ! ground frost and dew
    FrostSnowSfcIce(I,J)    = CondenseVapGrd(I,J)

    enddo
   enddo

    ! snowpack water processs
    call SnowWaterMainGlacier(noahmp)

    ! glacier ice water processing with OpenACC
    !$acc parallel loop collapse(2) gang vector default(present) private(WatReplaceSublim) firstprivate(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS,noahmp%config%domain%ITE

    ! total surface input water to glacier ice
    SoilSfcInflow(I,J) = (PondSfcThinSnwMelt(I,J) + PondSfcThinSnwComb(I,J) + PondSfcThinSnwTrans(I,J)) / MainTimeStep * 0.001  ! convert units (mm/s -> m/s)
    if ( NumSnowLayerNeg(I,J) == 0 ) then
       SoilSfcInflow(I,J) = SoilSfcInflow(I,J) + (SnowBotOutflow(I,J) + RainfallGround(I,J)) * 0.001
    else
       SoilSfcInflow(I,J) = SoilSfcInflow(I,J) + SnowBotOutflow(I,J) * 0.001
    endif
#ifdef WRF_HYDRO
    SoilSfcInflow(I,J) = SoilSfcInflow(I,J) + WaterHeadSfc(I,J) / MainTimeStep * 0.001
#endif

    ! surface runoff
    RunoffSurface(I,J) = SoilSfcInflow(I,J) * 1000.0   ! mm/s

    ! glacier ice water
    if ( OptGlacierTreatment == 1 ) then
       WatReplaceSublim = 0.0
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          WatReplaceSublim = WatReplaceSublim + ThicknessSnowSoilLayer(I,LoopInd,J)*(SoilIce(I,LoopInd,J) - &
                             SoilIceTmp(LoopInd) + SoilLiqWater(I,LoopInd,J) - SoilLiqWaterTmp(LoopInd))
       enddo
       WatReplaceSublim = WatReplaceSublim * 1000.0 / MainTimeStep     ! convert to [mm/s]
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          SoilIce(I,LoopInd,J) = min(1.0, SoilIceTmp(LoopInd))
       enddo
    elseif ( OptGlacierTreatment == 2 ) then
       WatReplaceSublim = 0.0
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          SoilIce(I,LoopInd,J) = 1.0
       enddo
    endif

    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       SoilLiqWater(I,LoopInd,J) = 1.0 - SoilIce(I,LoopInd,J)
    enddo

    ! use RunoffSubsurface as a water balancer, GlacierExcessFlow is snow that disappears, WatReplaceSublim is
    ! water from below that replaces glacier loss
    if ( OptGlacierTreatment == 1 ) then
       RunoffSubsurface(I,J) = GlacierExcessFlow(I,J) + WatReplaceSublim
    elseif ( OptGlacierTreatment == 2 ) then
       RunoffSubsurface(I,J) = GlacierExcessFlow(I,J)
       VaporizeGrd(I,J)      = SublimSnowSfcIce(I,J)
       CondenseVapGrd(I,J)   = FrostSnowSfcIce(I,J)
    endif

    if ( OptGlacierTreatment == 2 ) then
       EvapGroundNet(I,J) = VaporizeGrd(I,J) - CondenseVapGrd(I,J)
       HeatLatentGrd(I,J) = EvapGroundNet(I,J) * LatHeatVapGrd(I,J)
    endif

#ifndef _OPENACC
    if ( maxval(SoilIce) < 0.0001 ) then
       write(*,*) "GLACIER HAS MELTED AT: ", I, J, " ARE YOU SURE THIS SHOULD BE A GLACIER POINT?"
    endif
#endif

      end do
    end do
    !$acc end parallel loop

    !$acc end data


    end associate

  end subroutine WaterMainGlacier

end module WaterMainGlacierMod
