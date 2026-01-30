module TileDrainageSimpleMod

!!! Calculate tile drainage discharge [mm] based on simple model

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine TileDrainageSimple(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: TILE_DRAIN
! Original code: P. Valayamkunnath (NCAR)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                              ! grid indices
    integer                          :: IndSoil                           ! soil layer loop index
    real(kind=kind_noahmp)           :: DrainWatVolTot                    ! temporary variable for drainage volume [mm]
    real(kind=kind_noahmp)           :: DrainCoeffTmp                     ! temporary variable for drainage
    real(kind=kind_noahmp)           :: DrainWatTmp                       ! temporary variable for drainage
    real(kind=kind_noahmp)           :: WatExcFieldCap(1:noahmp%config%domain%NumSoilLayer)     ! temp variable for volume of water above field capacity
    real(kind=kind_noahmp)           :: SoilFieldCapLiq(1:noahmp%config%domain%NumSoilLayer)    ! Available field capacity = field capacity - SoilIce [m3/m3]
    real(kind=kind_noahmp)           :: DrainFracTmp(1:noahmp%config%domain%NumSoilLayer)       ! tile drainage fraction

! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp) &
    !$acc private(IndSoil, DrainWatVolTot, DrainCoeffTmp, DrainWatTmp) &
    !$acc private(WatExcFieldCap, SoilFieldCapLiq, DrainFracTmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         !cycle condition copied from SoilWaterMainMod before GPU port
         if ( noahmp%water%state%TileDrainFrac (I,J) <= 0.3 ) cycle
        associate(                                                                 &
                  NumSoilLayer         => noahmp%config%domain%NumSoilLayer       ,& ! in,    number of soil layers
                  DepthSoilLayer       => noahmp%config%domain%DepthSoilLayer     ,& ! in,    depth [m] of layer-bottom from soil surface
                  SoilTimeStep         => noahmp%config%domain%SoilTimeStep       ,& ! in,    noahmp soil timestep [s]
                  ThicknessSoilLayer   => noahmp%config%domain%ThicknessSoilLayer ,& ! in,    soil layer thickness [m]
                  TileDrainCoeffSp     => noahmp%water%param%TileDrainCoeffSp(I,J),& ! in,    drainage coefficient [mm/d]
                  DrainSoilLayerInd    => noahmp%water%param%DrainSoilLayerInd(I,J),& ! in,    starting soil layer for drainage
                  TileDrainTubeDepth   => noahmp%water%param%TileDrainTubeDepth(I,J),& ! in,    depth of drain tube from the soil surface
                  DrainFacSoilWat      => noahmp%water%param%DrainFacSoilWat(I,J) ,& ! in,    drainage factor for soil moisture
                  SoilMoistureFieldCap => noahmp%water%param%SoilMoistureFieldCap ,& ! in,    reference soil moisture (field capacity) [m3/m3]
                  SoilIce              => noahmp%water%state%SoilIce              ,& ! in,    soil ice content [m3/m3]
                  SoilLiqWater         => noahmp%water%state%SoilLiqWater         ,& ! inout, soil water content [m3/m3]
                  SoilMoisture         => noahmp%water%state%SoilMoisture         ,& ! inout, total soil moisture [m3/m3]
                  TileDrain            => noahmp%water%flux%TileDrain(I,J)         & ! out,   tile drainage [mm/s]
                 )
! ----------------------------------------------------------------------

        ! initialization
        !$acc loop seq
         do IndSoil = 1, NumSoilLayer
            ThicknessSoilLayer(I,IndSoil,J)  = 0.0
         enddo
        DrainFracTmp          = 0.0
        SoilFieldCapLiq       = 0.0
        DrainWatVolTot        = 0.0
        WatExcFieldCap        = 0.0
        TileDrain             = 0.0
        DrainWatTmp           = 0.0
        DrainCoeffTmp         = TileDrainCoeffSp * SoilTimeStep / (24.0 * 3600.0)

        !$acc loop seq
        do IndSoil = 1, NumSoilLayer
           if ( IndSoil == 1 ) then
              ThicknessSoilLayer(I,IndSoil,J) = -1.0 * DepthSoilLayer(I,IndSoil,J)
           else
              ThicknessSoilLayer(I,IndSoil,J) = DepthSoilLayer(I,IndSoil-1,J) - DepthSoilLayer(I,IndSoil,J)
           endif
        enddo

        if ( DrainSoilLayerInd == 0 ) then ! drainage from one specified layer in NoahmpTable.TBL
           IndSoil                  = int(TileDrainTubeDepth)
           SoilFieldCapLiq(IndSoil) = SoilMoistureFieldCap(I,IndSoil,J) - SoilIce(I,IndSoil,J)
           WatExcFieldCap(IndSoil)  = (SoilLiqWater(I,IndSoil,J) - (DrainFacSoilWat*SoilFieldCapLiq(IndSoil))) * &
                                      ThicknessSoilLayer(I,IndSoil,J) * 1000.0 ! mm
           if ( WatExcFieldCap(IndSoil) > 0.0 ) then
              if ( WatExcFieldCap(IndSoil) > DrainCoeffTmp ) WatExcFieldCap(IndSoil) = DrainCoeffTmp
              DrainWatVolTot            = DrainWatVolTot  + WatExcFieldCap(IndSoil)
              SoilLiqWater(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) - &
                                          (WatExcFieldCap(IndSoil) / (ThicknessSoilLayer(I,IndSoil,J) * 1000.0))
              SoilMoisture(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) + SoilIce(I,IndSoil,J)
           endif
        else if ( DrainSoilLayerInd == 1 ) then
           !$acc loop seq
           do IndSoil = 1, 2
              SoilFieldCapLiq(IndSoil) = SoilMoistureFieldCap(I,IndSoil,J) - SoilIce(I,IndSoil,J)
              WatExcFieldCap(IndSoil)  = (SoilLiqWater(I,IndSoil,J) - (DrainFacSoilWat*SoilFieldCapLiq(IndSoil))) * &
                                         ThicknessSoilLayer(I,IndSoil,J) * 1000.0 ! mm
              if ( WatExcFieldCap(IndSoil) < 0.0 ) WatExcFieldCap(IndSoil) = 0.0
              DrainWatTmp = DrainWatTmp + WatExcFieldCap(IndSoil)
           enddo
           !$acc loop seq
           do IndSoil = 1, 2
              if ( WatExcFieldCap(IndSoil) /= 0.0 ) then
                 DrainFracTmp(IndSoil) = WatExcFieldCap(IndSoil) / DrainWatTmp
              endif
           enddo
           if ( DrainWatTmp > 0.0 ) then
              if ( DrainWatTmp > DrainCoeffTmp ) DrainWatTmp = DrainCoeffTmp
              DrainWatVolTot = DrainWatVolTot + DrainWatTmp
              !$acc loop seq
              do IndSoil = 1, 2
                 WatExcFieldCap(IndSoil)   = DrainFracTmp(IndSoil) * DrainWatTmp
                 SoilLiqWater(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) - &
                                             (WatExcFieldCap(IndSoil) / (ThicknessSoilLayer(I,IndSoil,J) * 1000.0))
                 SoilMoisture(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) + SoilIce(I,IndSoil,J)
              enddo
           endif
        else if ( DrainSoilLayerInd == 2 ) then
           !$acc loop seq
           do IndSoil = 1, 3
              SoilFieldCapLiq(IndSoil) = SoilMoistureFieldCap(I,IndSoil,J) - SoilIce(I,IndSoil,J)
              WatExcFieldCap(IndSoil)  = (SoilLiqWater(I,IndSoil,J) - (DrainFacSoilWat*SoilFieldCapLiq(IndSoil))) * &
                                         ThicknessSoilLayer(I,IndSoil,J) * 1000.0
              if ( WatExcFieldCap(IndSoil) < 0.0 ) WatExcFieldCap(IndSoil) = 0.0
              DrainWatTmp = DrainWatTmp + WatExcFieldCap(IndSoil)
           enddo
           !$acc loop seq
           do IndSoil = 1, 3
              if ( WatExcFieldCap(IndSoil) /= 0.0 ) then
                 DrainFracTmp(IndSoil) = WatExcFieldCap(IndSoil) / DrainWatTmp
              endif
           enddo
           if ( DrainWatTmp > 0.0 ) then
              if ( DrainWatTmp > DrainCoeffTmp ) DrainWatTmp = DrainCoeffTmp
              DrainWatVolTot = DrainWatVolTot + DrainWatTmp
              !$acc loop seq
              do IndSoil = 1, 3
                 WatExcFieldCap(IndSoil)   = DrainFracTmp(IndSoil) * DrainWatTmp
                 SoilLiqWater(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) - &
                                             (WatExcFieldCap(IndSoil) / (ThicknessSoilLayer(I,IndSoil,J) * 1000.0))
                 SoilMoisture(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) + SoilIce(I,IndSoil,J)
              enddo
           endif
        else if ( DrainSoilLayerInd == 3 ) then
           !$acc loop seq
           do IndSoil = 2, 3
              SoilFieldCapLiq(IndSoil) = SoilMoistureFieldCap(I,IndSoil,J) - SoilIce(I,IndSoil,J)
              WatExcFieldCap(IndSoil)  = (SoilLiqWater(I,IndSoil,J) - (DrainFacSoilWat*SoilFieldCapLiq(IndSoil))) * &
                                         ThicknessSoilLayer(I,IndSoil,J) * 1000.0
              if ( WatExcFieldCap(IndSoil) < 0.0 ) WatExcFieldCap(IndSoil) = 0.0
              DrainWatTmp = DrainWatTmp + WatExcFieldCap(IndSoil)
           enddo
           !$acc loop seq
           do IndSoil = 2, 3
              if ( WatExcFieldCap(IndSoil) /= 0.0 ) then
                 DrainFracTmp(IndSoil) = WatExcFieldCap(IndSoil) / DrainWatTmp
              endif
           enddo
           if ( DrainWatTmp > 0.0 ) then
              if ( DrainWatTmp > DrainCoeffTmp ) DrainWatTmp = DrainCoeffTmp
              DrainWatVolTot = DrainWatVolTot + DrainWatTmp
              !$acc loop seq
              do IndSoil = 2, 3
                 WatExcFieldCap(IndSoil)   = DrainFracTmp(IndSoil) * DrainWatTmp
                 SoilLiqWater(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) - &
                                             (WatExcFieldCap(IndSoil) / (ThicknessSoilLayer(I,IndSoil,J) * 1000.0))
                 SoilMoisture(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) + SoilIce(I,IndSoil,J)
              enddo
           endif
        else if ( DrainSoilLayerInd == 4 ) then
           !$acc loop seq
           do IndSoil = 3, 4
              SoilFieldCapLiq(IndSoil) = SoilMoistureFieldCap(I,IndSoil,J) - SoilIce(I,IndSoil,J)
              WatExcFieldCap(IndSoil)  = (SoilLiqWater(I,IndSoil,J) - (DrainFacSoilWat*SoilFieldCapLiq(IndSoil))) * &
                                         ThicknessSoilLayer(I,IndSoil,J) * 1000.0
              if ( WatExcFieldCap(IndSoil) < 0.0 ) WatExcFieldCap(IndSoil) = 0.0
              DrainWatTmp = DrainWatTmp + WatExcFieldCap(IndSoil)
           enddo
           !$acc loop seq
           do IndSoil = 3, 4
              if ( WatExcFieldCap(IndSoil) /= 0.0 ) then
                 DrainFracTmp(IndSoil) = WatExcFieldCap(IndSoil) / DrainWatTmp
              endif
           enddo
           if ( DrainWatTmp > 0.0 ) then
              if ( DrainWatTmp > DrainCoeffTmp ) DrainWatTmp = DrainCoeffTmp
              DrainWatVolTot = DrainWatVolTot + DrainWatTmp
              !$acc loop seq
              do IndSoil = 3, 4
                 WatExcFieldCap(IndSoil)   = DrainFracTmp(IndSoil) * DrainWatTmp
                 SoilLiqWater(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) - (WatExcFieldCap(IndSoil) / &
                                             (ThicknessSoilLayer(I,IndSoil,J) * 1000.0))
                 SoilMoisture(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) + SoilIce(I,IndSoil,J)
              enddo
           endif
        else if ( DrainSoilLayerInd == 5 ) then ! from all the four layers
           !$acc loop seq
           do IndSoil = 1, 4
              SoilFieldCapLiq(IndSoil) = SoilMoistureFieldCap(I,IndSoil,J) - SoilIce(I,IndSoil,J)
              WatExcFieldCap(IndSoil)  = (SoilLiqWater(I,IndSoil,J) - (DrainFacSoilWat*SoilFieldCapLiq(IndSoil))) * &
                                         ThicknessSoilLayer(I,IndSoil,J) * 1000.0
              if ( WatExcFieldCap(IndSoil) < 0.0 ) WatExcFieldCap(IndSoil) = 0.0
              DrainWatTmp = DrainWatTmp + WatExcFieldCap(IndSoil)
           enddo
           !$acc loop seq
           do IndSoil = 1, 4
              if ( WatExcFieldCap(IndSoil) /= 0.0 ) then
                 DrainFracTmp(IndSoil) = WatExcFieldCap(IndSoil) / DrainWatTmp
              endif
           enddo
           if ( DrainWatTmp > 0.0 ) then
              if ( DrainWatTmp > DrainCoeffTmp ) DrainWatTmp = DrainCoeffTmp
              DrainWatVolTot = DrainWatVolTot + DrainWatTmp
              !$acc loop seq
              do IndSoil = 1, 4
                 WatExcFieldCap(IndSoil)   = DrainFracTmp(IndSoil) * DrainWatTmp
                 SoilLiqWater(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) - (WatExcFieldCap(IndSoil) / &
                                             (ThicknessSoilLayer(I,IndSoil,J) * 1000.0))
                 SoilMoisture(I,IndSoil,J) = SoilLiqWater(I,IndSoil,J) + SoilIce(I,IndSoil,J)
              enddo
           endif
        endif

        TileDrain = DrainWatVolTot / SoilTimeStep

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine TileDrainageSimple

end module TileDrainageSimpleMod
