module ForcingVarOutTransferMod

!!! Transfer column (1-D) Noah-MP forcing variables to 2D NoahmpIO for output

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

  subroutine ForcingVarOutTransfer(noahmp, NoahmpIO)

    implicit none

    type(noahmp_type),   intent(inout) :: noahmp
    type(NoahmpIO_type), intent(inout) :: NoahmpIO

    integer :: I, J

    !$acc parallel loop collapse(2) present(noahmp, NoahmpIO)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if (NoahmpIO%XLAND(I,J) - 1.5 >= 0.0) cycle ! Do out write output for open water points

    NoahmpIO%FORCTLSM  (I,J) = noahmp%forcing%TemperatureAirRefHeight(I,J)
    NoahmpIO%FORCQLSM  (I,J) = noahmp%forcing%SpecHumidityRefHeight(I,J)
    NoahmpIO%FORCPLSM  (I,J) = noahmp%forcing%PressureAirRefHeight(I,J)
    NoahmpIO%FORCWLSM  (I,J) = sqrt(noahmp%forcing%WindEastwardRefHeight(I,J)**2 + &
                                    noahmp%forcing%WindNorthwardRefHeight(I,J)**2)
    NoahmpIO%RadSwDirFrac(I,J) = noahmp%forcing%RadSwDirFrac(I,J)
    NoahmpIO%RadSwVisFrac(I,J) = noahmp%forcing%RadSwVisFrac(I,J)

      end do
      end do


  end subroutine ForcingVarOutTransfer

end module ForcingVarOutTransferMod
