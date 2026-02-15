module BiochemVarOutTransferMod

!!! Transfer column (1-D) biochemistry variables to 2D NoahmpIO for output

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type
  use NoahmpVarType

  implicit none

contains

!=== Transfer model states to output =====

  subroutine BiochemVarOutTransfer(noahmp, NoahmpIO)

    implicit none

    type(noahmp_type),   intent(inout) :: noahmp
    type(NoahmpIO_type), intent(inout) :: NoahmpIO

    integer :: I, J
! ---------------------------------------------------------------------
    !$acc parallel loop collapse(2) present(noahmp, NoahmpIO)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
! ---------------------------------------------------------------------
        if (NoahmpIO%XLAND(I,J) - 1.5 >= 0.0) cycle ! Do out write output for open water points

    ! biochem state variables
    NoahmpIO%LFMASSXY(I,J) = noahmp%biochem%state%LeafMass(I,J)
    NoahmpIO%RTMASSXY(I,J) = noahmp%biochem%state%RootMass(I,J)
    NoahmpIO%STMASSXY(I,J) = noahmp%biochem%state%StemMass(I,J)
    NoahmpIO%WOODXY  (I,J) = noahmp%biochem%state%WoodMass(I,J)
    NoahmpIO%STBLCPXY(I,J) = noahmp%biochem%state%CarbonMassDeepSoil(I,J)
    NoahmpIO%FASTCPXY(I,J) = noahmp%biochem%state%CarbonMassShallowSoil(I,J)
    NoahmpIO%GDDXY   (I,J) = noahmp%biochem%state%GrowDegreeDay(I,J)
    NoahmpIO%PGSXY   (I,J) = noahmp%biochem%state%PlantGrowStage(I,J)
    NoahmpIO%GRAINXY (I,J) = noahmp%biochem%state%GrainMass(I,J)

    ! biochem flux variables
    NoahmpIO%NEEXY   (I,J) = noahmp%biochem%flux%NetEcoExchange(I,J)
    NoahmpIO%GPPXY   (I,J) = noahmp%biochem%flux%GrossPriProduction(I,J)
    NoahmpIO%NPPXY   (I,J) = noahmp%biochem%flux%NetPriProductionTot(I,J)
    NoahmpIO%PSNXY   (I,J) = noahmp%biochem%flux%PhotosynTotal(I,J)

    end do
  end do

  end subroutine BiochemVarOutTransfer

end module BiochemVarOutTransferMod
