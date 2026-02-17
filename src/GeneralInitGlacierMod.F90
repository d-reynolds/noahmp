module GeneralInitGlacierMod

!!! General initialization for glacier variables

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
 
  implicit none

contains

  subroutine GeneralInitGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in NOAHMP_GLACIER)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J      ! grid indices
    integer                          :: LoopInd   ! loop index

    associate(                                                                       &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer            ,& ! in,  number of soil layers
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg    ,& ! in,  actual number of snow layers (negative)
              DepthSnowSoilLayer     => noahmp%config%domain%DepthSnowSoilLayer      ,& ! in,  depth of snow/soil layer-bottom [m]
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer   & ! out, thickness of snow/soil layers [m]
             )

   !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! initialize snow/soil layer thickness
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       if ( LoopInd == (NumSnowLayerNeg(I,J)+1) ) then
          ThicknessSnowSoilLayer(I,LoopInd,J) = - DepthSnowSoilLayer(I,LoopInd,J)
       else
          ThicknessSnowSoilLayer(I,LoopInd,J) = DepthSnowSoilLayer(I,LoopInd-1,J) - DepthSnowSoilLayer(I,LoopInd,J)
       endif
    enddo


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine GeneralInitGlacier

end module GeneralInitGlacierMod
