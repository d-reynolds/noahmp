module GlacierIceThermalPropertyMod

!!! Compute glacier ice thermal conductivity based on Noah scheme

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine GlacierIceThermalProperty(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: none (embedded in ENERGY_GLACIER)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                :: I, J               ! grid indices
    integer                :: LoopInd1, LoopInd2 ! loop index
    real(kind=kind_noahmp) :: DepthIceLayerMid   ! mid-point ice layer depth

! --------------------------------------------------------------------
    associate(                                                                          &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer               ,& ! in,  number of soil layers
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer     ,& ! in,  thickness of snow/soil layers [m]
              HeatCapacGlaIce        => noahmp%energy%state%HeatCapacGlaIce             ,& ! out, glacier ice layer volumetric specific heat [J/m3/K]
              ThermConductGlaIce     => noahmp%energy%state%ThermConductGlaIce           & ! out, glacier ice layer thermal conductivity [W/m/K]
             )

   !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd1,LoopInd2,DepthIceLayerMid)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    !$acc loop seq
    do LoopInd1 = 1, NumSoilLayer
       DepthIceLayerMid = 0.5 * ThicknessSnowSoilLayer(I,LoopInd1,J)
       !$acc loop seq
       do LoopInd2 = 1, LoopInd1-1
          DepthIceLayerMid = DepthIceLayerMid + ThicknessSnowSoilLayer(I,LoopInd2,J)
       enddo
       HeatCapacGlaIce(I,LoopInd1,J)    = 1.0e6 * (0.8194 + 0.1309 * DepthIceLayerMid)
       ThermConductGlaIce(I,LoopInd1,J) = 0.32333 + (0.10073 * DepthIceLayerMid)
    enddo


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine GlacierIceThermalProperty

end module GlacierIceThermalPropertyMod
