module SnowpackCompactionMod

!!! Snowpack compaction process (2D GPU-optimized)
!!! Update snow depth via compaction due to destructive metamorphism, overburden, & melt

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowpackCompaction(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: COMPACT
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                   ! grid indices
    integer                          :: LoopInd                ! snow layer loop index
    real(kind=kind_noahmp)           :: SnowBurden             ! pressure of overlying snow [kg/m2]
    real(kind=kind_noahmp)           :: SnowCompactAgeExpFac   ! EXPF=exp(-c4*(273.15-TemperatureSoilSnow))
    real(kind=kind_noahmp)           :: TempDiff               ! ConstFreezePoint - TemperatureSoilSnow[K]
    real(kind=kind_noahmp)           :: SnowVoid               ! void (1 - SnowIce - SnowLiqWater)
    real(kind=kind_noahmp)           :: SnowWatTotTmp          ! water mass (ice + liquid) [kg/m2]
    real(kind=kind_noahmp)           :: SnowIceDens            ! partial density of ice [kg/m3]

! --------------------------------------------------------------------
        associate(                                                                       &
                  MainTimeStep           => noahmp%config%domain%MainTimeStep           ,& ! in,    noahmp main time step [s]
                  TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow     ,& ! in,    snow and soil layer temperature [K]
                  SnowIce                => noahmp%water%state%SnowIce                  ,& ! in,    snow layer ice [mm]
                  SnowLiqWater           => noahmp%water%state%SnowLiqWater             ,& ! in,    snow layer liquid water [mm]
                  IndexPhaseChange       => noahmp%water%state%IndexPhaseChange         ,& ! in,    phase change index [0-none;1-melt;2-refreeze]
                  SnowIceFracPrev        => noahmp%water%state%SnowIceFracPrev          ,& ! in,    ice fraction in snow layers at previous timestep
                  SnowCompactBurdenFac   => noahmp%water%param%SnowCompactBurdenFac ,& ! in,    snow overburden compaction parameter [m3/kg]
                  SnowCompactAgingFac1   => noahmp%water%param%SnowCompactAgingFac1 ,& ! in,    snow desctructive metamorphism compaction factor1 [1/s]
                  SnowCompactAgingFac2   => noahmp%water%param%SnowCompactAgingFac2 ,& ! in,    snow desctructive metamorphism compaction factor2 [1/k]
                  SnowCompactAgingFac3   => noahmp%water%param%SnowCompactAgingFac3 ,& ! in,    snow desctructive metamorphism compaction factor3
                  SnowCompactAgingMax    => noahmp%water%param%SnowCompactAgingMax ,& ! in,    maximum destructive metamorphism compaction [kg/m3]
                  SnowViscosityCoeff     => noahmp%water%param%SnowViscosityCoeff       ,& ! in,    snow viscosity coeff [kg s/m2],Anderson1979:0.52e6~1.38e6
                  NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg        ,& ! inout, actual number of snow layers (negative)
                  ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m]
                  CompactionSnowAging    => noahmp%water%flux%CompactionSnowAging       ,& ! out,   rate of compaction due to destructive metamorphism [1/s]
                  CompactionSnowBurden   => noahmp%water%flux%CompactionSnowBurden      ,& ! out,   rate of compaction of snowpack due to overburden [1/s]
                  CompactionSnowMelt     => noahmp%water%flux%CompactionSnowMelt        ,& ! out,   rate of compaction of snowpack due to melt [1/s]
                  CompactionSnowTot      => noahmp%water%flux%CompactionSnowTot         ,& ! out,   change in fractional-thickness due to compaction [1/s]
                  SnowIceFrac            => noahmp%water%state%SnowIceFrac               & ! out,   fraction of ice in snow layers at current time step
                 )

    !$acc parallel loop collapse(2) gang vector default(present) &
    !$acc private(LoopInd, SnowBurden, SnowCompactAgeExpFac, TempDiff, SnowVoid, SnowWatTotTmp, SnowIceDens)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

      if ( NumSnowLayerNeg(I,J) >= 0 ) cycle  ! no snow layers

        ! initialization for out-only variables
        !$acc loop seq
        do LoopInd = NumSnowLayerNeg(I,J)+1, 0
          CompactionSnowAging(I,LoopInd,J)  = 0.0
          CompactionSnowBurden(I,LoopInd,J) = 0.0
          CompactionSnowMelt(I,LoopInd,J)   = 0.0
          CompactionSnowTot(I,LoopInd,J)    = 0.0
          SnowIceFrac(I,LoopInd,J)          = 0.0
        enddo
        
        ! start snow compaction
        SnowBurden = 0.0
        !$acc loop seq
        do LoopInd = NumSnowLayerNeg(I,J)+1, 0

          SnowWatTotTmp        = SnowIce(I,LoopInd,J) + SnowLiqWater(I,LoopInd,J)
          SnowIceFrac(I,LoopInd,J) = SnowIce(I,LoopInd,J) / SnowWatTotTmp
          SnowVoid             = 1.0 - (SnowIce(I,LoopInd,J)/ConstDensityIce + SnowLiqWater(I,LoopInd,J)/ConstDensityWater) / &
                                       ThicknessSnowSoilLayer(I,LoopInd,J)

          ! Allow compaction only for non-saturated node and higher ice lens node.
          if ( (SnowVoid > 0.001) .and. (SnowIce(I,LoopInd,J) > 0.1) ) then
             SnowIceDens = SnowIce(I,LoopInd,J) / ThicknessSnowSoilLayer(I,LoopInd,J)
             TempDiff    = max(0.0, ConstFreezePoint-TemperatureSoilSnow(I,LoopInd,J))

             ! Settling/compaction as a result of destructive metamorphism
             SnowCompactAgeExpFac         = exp(-SnowCompactAgingFac2(I,J) * TempDiff)
             CompactionSnowAging(I,LoopInd,J) = -SnowCompactAgingFac1(I,J) * SnowCompactAgeExpFac
             if ( SnowIceDens > SnowCompactAgingMax(I,J) ) &
                CompactionSnowAging(I,LoopInd,J) = CompactionSnowAging(I,LoopInd,J) * exp(-46.0e-3*(SnowIceDens-SnowCompactAgingMax(I,J)))
             if ( SnowLiqWater(I,LoopInd,J) > (0.01*ThicknessSnowSoilLayer(I,LoopInd,J)) ) &
                CompactionSnowAging(I,LoopInd,J) = CompactionSnowAging(I,LoopInd,J) * SnowCompactAgingFac3(I,J)                ! Liquid water term

             ! Compaction due to overburden
             CompactionSnowBurden(I,LoopInd,J) = -(SnowBurden + 0.5*SnowWatTotTmp) * &
                                       exp(-0.08*TempDiff-SnowCompactBurdenFac(I,J)*SnowIceDens) / SnowViscosityCoeff(I,J)  ! 0.5*SnowWatTotTmp -> self-burden

             ! Compaction occurring during melt
             if ( IndexPhaseChange(I,LoopInd,J) == 1 ) then
                CompactionSnowMelt(I,LoopInd,J) = max(0.0, (SnowIceFracPrev(I,LoopInd,J)-SnowIceFrac(I,LoopInd,J)) / &
                                                       max(1.0e-6, SnowIceFracPrev(I,LoopInd,J)))
                CompactionSnowMelt(I,LoopInd,J) = -CompactionSnowMelt(I,LoopInd,J) / MainTimeStep   ! sometimes too large
             else
                CompactionSnowMelt(I,LoopInd,J) = 0.0
             endif

             ! Time rate of fractional change in snow thickness (units of s-1)
             CompactionSnowTot(I,LoopInd,J) = (CompactionSnowAging(I,LoopInd,J) + CompactionSnowBurden(I,LoopInd,J) + &
                                           CompactionSnowMelt(I,LoopInd,J) ) * MainTimeStep
             CompactionSnowTot(I,LoopInd,J) = max(-0.5, CompactionSnowTot(I,LoopInd,J))

             ! The change in DZ due to compaction
             ThicknessSnowSoilLayer(I,LoopInd,J) = ThicknessSnowSoilLayer(I,LoopInd,J) * (1.0 + CompactionSnowTot(I,LoopInd,J))
             ThicknessSnowSoilLayer(I,LoopInd,J) = max(ThicknessSnowSoilLayer(I,LoopInd,J), &
                                               SnowIce(I,LoopInd,J)/ConstDensityIce + SnowLiqWater(I,LoopInd,J)/ConstDensityWater)

             ! Constrain snow density to a reasonable range (50~500 kg/m3)
             ThicknessSnowSoilLayer(I,LoopInd,J) = min( max( ThicknessSnowSoilLayer(I,LoopInd,J),&
                                                        (SnowIce(I,LoopInd,J)+SnowLiqWater(I,LoopInd,J))/500.0 ), &
                                                   (SnowIce(I,LoopInd,J)+SnowLiqWater(I,LoopInd,J))/50.0 )
          endif

          ! Pressure of overlying snow
          SnowBurden = SnowBurden + SnowWatTotTmp

        enddo


      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine SnowpackCompaction

end module SnowpackCompactionMod
