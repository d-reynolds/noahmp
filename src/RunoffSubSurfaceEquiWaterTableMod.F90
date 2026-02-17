module RunoffSubSurfaceEquiWaterTableMod

!!! Calculate subsurface runoff using equilibrium water table depth (Niu et al., 2005)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use WaterTableEquilibriumMod, only : WaterTableEquilibrium

  implicit none

contains

  subroutine RunoffSubSurfaceEquiWaterTable(noahmp)

! ------------------------ Code history --------------------------------------------------
! Originally embeded in SOILWATER subroutine instead of as a separate subroutine
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J      ! grid indices

    ! compute equilibrium water table depth
    call WaterTableEquilibrium(noahmp)

    associate(                                                                      &
              SoilImpervFracMax => noahmp%water%state%SoilImpervFracMax       ,& ! in,    maximum soil imperviousness fraction
              GridTopoIndex     => noahmp%water%param%GridTopoIndex           ,& ! in,    gridcell mean topgraphic index (global mean)
              RunoffDecayFac    => noahmp%water%param%RunoffDecayFac          ,& ! in,    runoff decay factor [m-1]
              BaseflowCoeff     => noahmp%water%param%BaseflowCoeff           ,& ! inout, baseflow coefficient [mm/s]
              WaterTableDepth   => noahmp%water%state%WaterTableDepth         ,& ! out,   water table depth [m]
              RunoffSubsurface  => noahmp%water%flux%RunoffSubsurface          & ! out,   subsurface runoff [mm/s]
             )

    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points


    ! set parameter values specific for this scheme
    RunoffDecayFac(I,J) = 2.0
    BaseflowCoeff(I,J)  = 4.0


    ! compuate subsurface runoff mm/s
    RunoffSubsurface(I,J) = (1.0 - SoilImpervFracMax(I,J)) * BaseflowCoeff(I,J) * &
                       exp(-GridTopoIndex(I,J)) * exp(-RunoffDecayFac(I,J) * WaterTableDepth(I,J))


      end do
    end do
    !$acc end parallel loop



    end associate

  end subroutine RunoffSubSurfaceEquiWaterTable

end module RunoffSubSurfaceEquiWaterTableMod
