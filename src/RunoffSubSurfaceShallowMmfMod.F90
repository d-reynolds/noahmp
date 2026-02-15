module RunoffSubSurfaceShallowMmfMod

!!! Calculate subsurface runoff based on MMF groundwater scheme

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use ShallowWaterTableMmfMod, only : ShallowWaterTableMMF

  implicit none

contains

  subroutine RunoffSubSurfaceShallowWaterMMF(noahmp)

! ------------------------ Code history --------------------------------------------------
! Originally embeded in WATER subroutine instead of as a separate subroutine
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    ! local variable
    integer                          :: I, J      ! grid indices

    ! compute shallow water table and moisture
    call ShallowWaterTableMMF(noahmp)

    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points

! --------------------------------------------------------------------
    associate(                                                                    &
              NumSoilLayer        => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              SoilIce             => noahmp%water%state%SoilIce                  ,& ! in,    soil ice content [m3/m3]
              DrainSoilBot        => noahmp%water%flux%DrainSoilBot(I,J)         ,& ! in,    soil bottom drainage [mm/s]
              SoilLiqWater        => noahmp%water%state%SoilLiqWater             ,& ! inout, soil water content [m3/m3]
              SoilMoisture        => noahmp%water%state%SoilMoisture             ,& ! inout, total soil water content [m3/m3]
              WaterStorageAquifer => noahmp%water%state%WaterStorageAquifer(I,J) ,& ! inout, water storage in aquifer [mm]
              RunoffSubsurface    => noahmp%water%flux%RunoffSubsurface(I,J)      & ! out,   subsurface runoff [mm/s] 
             )
! ----------------------------------------------------------------------

    ! update moisture
    SoilLiqWater(I,NumSoilLayer,J) = SoilMoisture(I,NumSoilLayer,J) - SoilIce(I,NumSoilLayer,J)

    ! compute subsurface runoff
    RunoffSubsurface    = RunoffSubsurface + DrainSoilBot 
    WaterStorageAquifer = 0.0

    end associate

      end do
    end do

  end subroutine RunoffSubSurfaceShallowWaterMMF

end module RunoffSubSurfaceShallowMmfMod
