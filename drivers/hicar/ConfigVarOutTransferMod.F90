module ConfigVarOutTransferMod

!!! To transfer 1D Noah-MP column Config variables to 2D NoahmpIO for output

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type
  use NoahmpVarType

  implicit none

contains

!=== Transfer model states to output=====

  subroutine ConfigVarOutTransfer(noahmp, NoahmpIO)

    implicit none

    type(NoahmpIO_type) , intent(inout) :: NoahmpIO
    type(noahmp_type),    intent(inout) :: noahmp

    integer :: I, J
! ----------------------------------------------------------------------
    associate(                                                         &
              NumSnowLayerMax => noahmp%config%domain%NumSnowLayerMax ,&
              NumSoilLayer    => noahmp%config%domain%NumSoilLayer     &
             )
! ----------------------------------------------------------------------

    !$acc parallel loop collapse(2) default(present) private(NumSnowLayerMax, NumSoilLayer)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if (NoahmpIO%XLAND(I,J) - 1.5 >= 0.0) cycle ! Do out write output for open water points
    ! config domain variables
    NoahmpIO%ISNOWXY(I,J)  = noahmp%config%domain%NumSnowLayerNeg(I,J)
    NoahmpIO%ZSNSOXY(I,-NumSnowLayerMax+1:NumSoilLayer,J) = &
                            noahmp%config%domain%DepthSnowSoilLayer(I,-NumSnowLayerMax+1:NumSoilLayer,J)
    NoahmpIO%FORCZLSM(I,J) = noahmp%config%domain%RefHeightAboveSfc(I,J)

    end do
    end do
    !$acc end parallel loop

    end associate

  end subroutine ConfigVarOutTransfer

end module ConfigVarOutTransferMod
