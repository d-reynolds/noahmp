module GroundThermalPropertyGlacierMod

!!! Compute snow and glacier ice thermal conductivity and heat capacity

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowThermalPropertyMod,       only : SnowThermalProperty
  use GlacierIceThermalPropertyMod, only : GlacierIceThermalProperty

  implicit none

contains

  subroutine GroundThermalPropertyGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: THERMOPROP_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                :: I, J           ! grid indices
    integer                :: LoopInd        ! loop index

! --------------------------------------------------------------------

    ! compute snow thermal conductivity and heat capacity
    call SnowThermalProperty(noahmp)

    ! compute glacier ice thermal properties (using Noah glacial ice approximations)
    call GlacierIceThermalProperty(noahmp)

   !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                                             &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer                  ,& ! in,  number of soil layers
              MainTimeStep           => noahmp%config%domain%MainTimeStep                  ,& ! in,  main noahmp timestep [s]
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer        ,& ! in,  thickness of snow/soil layers [m]
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg(I,J)          ,& ! in,  actual number of snow layers (negative)
              SnowDepth              => noahmp%water%state%SnowDepth(I,J)                  ,& ! in,  snow depth [m]
              ThermConductSoilSnow   => noahmp%energy%state%ThermConductSoilSnow           ,& ! out, thermal conductivity [W/m/K] for all soil & snow
              HeatCapacSoilSnow      => noahmp%energy%state%HeatCapacSoilSnow              ,& ! out, heat capacity [J/m3/K] for all soil & snow
              PhaseChgFacSoilSnow    => noahmp%energy%state%PhaseChgFacSoilSnow            ,& ! out, energy factor for soil & snow phase change
              HeatCapacVolSnow       => noahmp%energy%state%HeatCapacVolSnow               ,& ! out, snow layer volumetric specific heat [J/m3/K]
              ThermConductSnow       => noahmp%energy%state%ThermConductSnow               ,& ! out, snow layer thermal conductivity [W/m/K]
              HeatCapacGlaIce        => noahmp%energy%state%HeatCapacGlaIce                ,& ! out, glacier ice layer volumetric specific heat [J/m3/K]
              ThermConductGlaIce     => noahmp%energy%state%ThermConductGlaIce              & ! out, glacier ice layer thermal conductivity [W/m/K]
             )
! ----------------------------------------------------------------------

    ! initialize
    !$acc loop seq
    do LoopInd = -2, NumSoilLayer
       HeatCapacSoilSnow(I,LoopInd,J)    = 0.0
       ThermConductSoilSnow(I,LoopInd,J) = 0.0
    enddo

    ! copy snow thermal properties
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg+1, 0
       ThermConductSoilSnow(I,LoopInd,J) = ThermConductSnow(I,LoopInd,J)
       HeatCapacSoilSnow(I,LoopInd,J)    = HeatCapacVolSnow(I,LoopInd,J)
    enddo

    ! copy glacier ice thermal properties
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       ThermConductSoilSnow(I,LoopInd,J) = ThermConductGlaIce(I,LoopInd,J)
       HeatCapacSoilSnow(I,LoopInd,J)    = HeatCapacGlaIce(I,LoopInd,J)
    enddo

    ! combine a temporary variable used for melting/freezing of snow and glacier ice
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       PhaseChgFacSoilSnow(I,LoopInd,J) = MainTimeStep / (HeatCapacSoilSnow(I,LoopInd,J)*ThicknessSnowSoilLayer(I,LoopInd,J))
    enddo

    ! snow/glacier ice interface
    if ( NumSnowLayerNeg == 0 ) then
       ThermConductSoilSnow(I,1,J) = (ThermConductSoilSnow(I,1,J)*ThicknessSnowSoilLayer(I,1,J) + 0.35*SnowDepth) / &
                                     (SnowDepth + ThicknessSnowSoilLayer(I,1,J))
    else
       ThermConductSoilSnow(I,1,J) = (ThermConductSoilSnow(I,1,J)*ThicknessSnowSoilLayer(I,1,J) + &
                                      ThermConductSoilSnow(I,0,J)*ThicknessSnowSoilLayer(I,0,J)) / &
                                     (ThicknessSnowSoilLayer(I,0,J) + ThicknessSnowSoilLayer(I,1,J))
    endif

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine GroundThermalPropertyGlacier

end module GroundThermalPropertyGlacierMod
