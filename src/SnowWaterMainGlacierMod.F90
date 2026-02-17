module SnowWaterMainGlacierMod

!!! Main glacier snow water module including all snowpack processes
!!! Snowfall -> Snowpack compaction -> Snow layer combination -> Snow layer division -> Snow Hydrology

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowfallBelowCanopyMod,      only : SnowfallAfterCanopyIntercept
  use SnowpackCompactionMod,       only : SnowpackCompaction
  use SnowpackCompactionAR24Mod,   only : SnowpackCompactionAR24
  use SnowLayerCombineMod,         only : SnowLayerCombine
  use SnowLayerDivideMod,          only : SnowLayerDivide
  use SnowpackHydrologyGlacierMod, only : SnowpackHydrologyGlacier
  use SnowAerosolSnicarMod,        only : SnowAerosolSnicar

  implicit none

contains

  subroutine SnowWaterMainGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SNOWWATER_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: LoopInd                 ! do loop/array indices
    real(kind=kind_noahmp)           :: SnowDensBulk            ! bulk density of snow [kg/m3]
    real(kind=kind_noahmp)           :: GlacierExcessRemainFrac ! fraction of mass remaining after glacier excess flow
    integer                          :: I, J                    ! grid indices

    associate(                                                                      &
              GlacierExcessFlow      => noahmp%water%flux%GlacierExcessFlow ,& ! out,   glacier excess flow [mm/s]
              PondSfcThinSnwComb     => noahmp%water%state%PondSfcThinSnwComb ,& ! out,   surface ponding [mm] from liquid in thin snow layer combination
              PondSfcThinSnwTrans    => noahmp%water%state%PondSfcThinSnwTrans ,& ! out,   surface ponding [mm] from thin snow when changing from multilayer to no layer
              NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax ,& ! in,    maximum number of snow layers
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer ,& ! in,    number of soil layers
              MainTimeStep           => noahmp%config%domain%MainTimeStep ,& ! in,    noahmp main time step [s]
              OptSnowAlbedo          => noahmp%config%nmlist%OptSnowAlbedo ,& ! in,    options for ground snow surface albedo
              SnoWatEqvMaxGlacier    => noahmp%water%param%SnoWatEqvMaxGlacier ,& ! in,    Maximum SWE allowed at glaciers [mm]
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m]
              DepthSnowSoilLayer     => noahmp%config%domain%DepthSnowSoilLayer ,& ! inout, depth of snow/soil layer-bottom [m]
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg ,& ! inout, actual number of snow layers (negative)
              TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow ,& ! inout, snow and soil layer temperature [K]
              SnowDepth              => noahmp%water%state%SnowDepth ,& ! inout, snow depth [m]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv ,& ! inout, snow water equivalent [mm]
              SnowIce                => noahmp%water%state%SnowIce ,& ! inout, snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater ,& ! inout, snow layer liquid water [mm]
              MassBChydropho         => noahmp%water%state%MassBChydropho ,& ! inout, mass of hydrophobic Black Carbon in snow [kg m-2]
              MassBChydrophi         => noahmp%water%state%MassBChydrophi ,& ! inout, mass of hydrophillic Black Carbon in snow [kg m-2]
              MassOChydropho         => noahmp%water%state%MassOChydropho ,& ! inout, mass of hydrophobic Organic Carbon in snow [kg m-2]
              MassOChydrophi         => noahmp%water%state%MassOChydrophi ,& ! inout, mass of hydrophillic Organic Carbon in snow [kg m-2]
              MassDust1              => noahmp%water%state%MassDust1 ,& ! inout, mass of dust species 1 in snow [kg m-2]
              MassDust2              => noahmp%water%state%MassDust2 ,& ! inout, mass of dust species 2 in snow [kg m-2]
              MassDust3              => noahmp%water%state%MassDust3 ,& ! inout, mass of dust species 3 in snow [kg m-2]
              MassDust4              => noahmp%water%state%MassDust4 ,& ! inout, mass of dust species 4 in snow [kg m-2]
              MassDust5              => noahmp%water%state%MassDust5 ,& ! inout, mass of dust species 5 in snow [kg m-2]
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer  & ! inout, depth [m] of layer-bottom from soil surface
             )

    !$acc parallel loop gang vector collapse(2) default(present) private(LoopInd, SnowDensBulk, GlacierExcessRemainFrac)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


         ! initialize out-only variables
         GlacierExcessFlow(I,J)       = 0.0
         PondSfcThinSnwComb(I,J)      = 0.0
         PondSfcThinSnwTrans(I,J)     = 0.0

         end do
      end do

    GlacierExcessRemainFrac = 1.0

    ! snowfall
    call SnowfallAfterCanopyIntercept(noahmp)

    ! do following snow layer compaction, combination, and division only for multi-layer snowpack

    ! snowpack compaction (option: 1->original,Anderson1976; 2->new,Abolafia-Rosenzweig2024)
    if (noahmp%config%nmlist%OptSnowCompaction == 1) call SnowpackCompaction(noahmp)
    if (noahmp%config%nmlist%OptSnowCompaction == 2) call SnowpackCompactionAR24(noahmp)

    ! snow layer combination
    !$acc parallel loop gang vector collapse(2) default(present) private(LoopInd, SnowDensBulk, GlacierExcessRemainFrac)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
        call SnowLayerCombine(noahmp, I, J)
      end do
    end do

    ! snow layer division
    call SnowLayerDivide(noahmp)

    ! snow hydrology for all snow cases
    call SnowpackHydrologyGlacier(noahmp)

    !$acc parallel loop gang vector collapse(2) default(present) private(LoopInd, SnowDensBulk, GlacierExcessRemainFrac)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! set empty snow layer properties to zero
    !$acc loop seq
    do LoopInd = -NumSnowLayerMax+1, NumSnowLayerNeg(I,J)
       SnowIce(I,LoopInd,J)                = 0.0
       SnowLiqWater(I,LoopInd,J)           = 0.0
       TemperatureSoilSnow(I,LoopInd,J)    = 0.0
       ThicknessSnowSoilLayer(I,LoopInd,J) = 0.0
       DepthSnowSoilLayer(I,LoopInd,J)     = 0.0

       if ( (OptSnowAlbedo == 3) .and. (NumSnowLayerNeg(I,J) < 0) ) then
          MassBChydropho(I,LoopInd,J)      = 0.0
          MassBChydrophi(I,LoopInd,J)      = 0.0
          MassOChydropho(I,LoopInd,J)      = 0.0
          MassOChydrophi(I,LoopInd,J)      = 0.0
          MassDust1(I,LoopInd,J)           = 0.0
          MassDust2(I,LoopInd,J)           = 0.0
          MassDust3(I,LoopInd,J)           = 0.0
          MassDust4(I,LoopInd,J)           = 0.0
          MassDust5(I,LoopInd,J)           = 0.0
       endif
    enddo

    ! to obtain equilibrium state of snow in glacier region
    if ( SnowWaterEquiv(I,J) > SnoWatEqvMaxGlacier(I,J) ) then 
       SnowDensBulk              = SnowIce(I,0,J) / ThicknessSnowSoilLayer(I,0,J)
       GlacierExcessFlow(I,J)         = SnowWaterEquiv(I,J) - SnoWatEqvMaxGlacier(I,J)
       SnowIce(I,0,J)                = SnowIce(I,0,J)  - GlacierExcessFlow(I,J)
       ThicknessSnowSoilLayer(I,0,J) = ThicknessSnowSoilLayer(I,0,J) - GlacierExcessFlow(I,J) / SnowDensBulk
       GlacierExcessFlow(I,J)         = GlacierExcessFlow(I,J) / MainTimeStep

       if ( OptSnowAlbedo == 3 ) then
          GlacierExcessRemainFrac = SnowIce(I,0,J) / (SnowIce(I,0,J) + GlacierExcessFlow(I,J))
          MassBChydropho(I,0,J)       = MassBChydropho(I,0,J) * GlacierExcessRemainFrac
          MassBChydrophi(I,0,J)       = MassBChydrophi(I,0,J) * GlacierExcessRemainFrac
          MassOChydropho(I,0,J)       = MassOChydropho(I,0,J) * GlacierExcessRemainFrac
          MassOChydrophi(I,0,J)       = MassOChydrophi(I,0,J) * GlacierExcessRemainFrac
          MassDust1(I,0,J)            = MassDust1(I,0,J) * GlacierExcessRemainFrac
          MassDust2(I,0,J)            = MassDust2(I,0,J) * GlacierExcessRemainFrac
          MassDust3(I,0,J)            = MassDust3(I,0,J) * GlacierExcessRemainFrac
          MassDust4(I,0,J)            = MassDust4(I,0,J) * GlacierExcessRemainFrac
          MassDust5(I,0,J)            = MassDust5(I,0,J) * GlacierExcessRemainFrac
       endif
    endif
   enddo
enddo

    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) call SnowAerosolSnicar(noahmp)

    !$acc parallel loop gang vector collapse(2) default(present) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
    ! sum up snow mass for layered snow
    if ( NumSnowLayerNeg(I,J) < 0 ) then  ! MB: only do for multi-layer
       SnowWaterEquiv(I,J) = 0.0
       !$acc loop seq
       do LoopInd = NumSnowLayerNeg(I,J)+1, 0
          SnowWaterEquiv(I,J) = SnowWaterEquiv(I,J) + SnowIce(I,LoopInd,J) + SnowLiqWater(I,LoopInd,J)
       enddo
    endif

    ! Reset DepthSnowSoilLayer and ThicknessSnowSoilLayer
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, 0
       ThicknessSnowSoilLayer(I,LoopInd,J) = -ThicknessSnowSoilLayer(I,LoopInd,J)
    enddo

    ThicknessSnowSoilLayer(I,1,J) = DepthSoilLayer(I,1,J)
    !$acc loop seq
    do LoopInd = 2, NumSoilLayer
       ThicknessSnowSoilLayer(I,LoopInd,J) = DepthSoilLayer(I,LoopInd,J) - DepthSoilLayer(I,LoopInd-1,J)
    enddo

    DepthSnowSoilLayer(I,NumSnowLayerNeg(I,J)+1,J) = ThicknessSnowSoilLayer(I,NumSnowLayerNeg(I,J)+1,J)
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+2, NumSoilLayer
       DepthSnowSoilLayer(I,LoopInd,J) = DepthSnowSoilLayer(I,LoopInd-1,J) + ThicknessSnowSoilLayer(I,LoopInd,J)
    enddo

    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       ThicknessSnowSoilLayer(I,LoopInd,J) = -ThicknessSnowSoilLayer(I,LoopInd,J)
    enddo

    ! Update SnowDepth for multi-layer snow
    if ( NumSnowLayerNeg(I,J) < 0 ) then
       SnowDepth(I,J) = 0.0
       !$acc loop seq
       do LoopInd = NumSnowLayerNeg(I,J)+1, 0
          SnowDepth(I,J) = SnowDepth(I,J) + ThicknessSnowSoilLayer(I,LoopInd,J)
       enddo
    endif

    ! update snow quantity
    if ( (SnowDepth(I,J) <= 1.0e-6) .or. (SnowWaterEquiv(I,J) <= 1.0e-6) ) then
       SnowDepth(I,J)      = 0.0
       SnowWaterEquiv(I,J) = 0.0
    endif


      end do
   end do


    end associate

  end subroutine SnowWaterMainGlacier

end module SnowWaterMainGlacierMod
