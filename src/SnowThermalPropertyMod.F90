module SnowThermalPropertyMod

!!! Compute snowpack thermal conductivity and volumetric specific heat (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowThermalProperty(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: CSNOW
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd                        ! loop index
    integer                          :: I, J                           ! grid indices
    real(kind=kind_noahmp)           :: SnowDensBulk                   ! bulk density of snow [kg/m3]

! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                         &
                  NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg(I,J)     ,& ! in,  actual number of snow layers (negative)
                  ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer   ,& ! in,  thickness of snow/soil layers [m] (3D)
                  OptSnowThermConduct    => noahmp%config%nmlist%OptSnowThermConduct      ,& ! in,  options for snow thermal conductivity schemes
                  SnowIce                => noahmp%water%state%SnowIce                    ,& ! in,  snow layer ice [mm] (3D)
                  SnowLiqWater           => noahmp%water%state%SnowLiqWater               ,& ! in,  snow layer liquid water [mm] (3D)
                  SnowIceVol             => noahmp%water%state%SnowIceVol                 ,& ! out, partial volume of snow ice [m3/m3] (3D)
                  SnowLiqWaterVol        => noahmp%water%state%SnowLiqWaterVol            ,& ! out, partial volume of snow liquid water [m3/m3] (3D)
                  SnowEffPorosity        => noahmp%water%state%SnowEffPorosity            ,& ! out, snow effective porosity [m3/m3] (3D)
                  HeatCapacVolSnow       => noahmp%energy%state%HeatCapacVolSnow          ,& ! out, snow layer volumetric specific heat [J/m3/K] (3D)
                  ThermConductSnow       => noahmp%energy%state%ThermConductSnow           & ! out, snow layer thermal conductivity [W/m/K] (3D)
                 )
! ----------------------------------------------------------------------

        !  effective porosity of snow
        !$acc loop seq
        do LoopInd = NumSnowLayerNeg+1, 0
           SnowIceVol(I,LoopInd,J)      = min(1.0, SnowIce(I,LoopInd,J)/(ThicknessSnowSoilLayer(I,LoopInd,J)*ConstDensityIce))
           SnowEffPorosity(I,LoopInd,J) = 1.0 - SnowIceVol(I,LoopInd,J)
           SnowLiqWaterVol(I,LoopInd,J) = min(SnowEffPorosity(I,LoopInd,J), &
                                              SnowLiqWater(I,LoopInd,J)/(ThicknessSnowSoilLayer(I,LoopInd,J)*ConstDensityWater))
        enddo

        ! thermal capacity of snow
        !$acc loop seq
        do LoopInd = NumSnowLayerNeg+1, 0
           SnowDensBulk              = (SnowIce(I,LoopInd,J) + SnowLiqWater(I,LoopInd,J)) / ThicknessSnowSoilLayer(I,LoopInd,J)
           HeatCapacVolSnow(I,LoopInd,J) = ConstHeatCapacIce*SnowIceVol(I,LoopInd,J) + ConstHeatCapacWater*SnowLiqWaterVol(I,LoopInd,J)
          !HeatCapacVolSnow(I,LoopInd,J) = 0.525e06  ! constant
        enddo

        ! thermal conductivity of snow
        !$acc loop seq
        do LoopInd = NumSnowLayerNeg+1, 0
           SnowDensBulk = (SnowIce(I,LoopInd,J) + SnowLiqWater(I,LoopInd,J)) / ThicknessSnowSoilLayer(I,LoopInd,J)
           if (OptSnowThermConduct == 1) &
              ThermConductSnow(I,LoopInd,J) = 3.2217e-6 * SnowDensBulk**2.0                      ! Stieglitz(yen,1965)
           if (OptSnowThermConduct == 2) &
              ThermConductSnow(I,LoopInd,J) = 2e-2 + 2.5e-6*SnowDensBulk*SnowDensBulk   ! Anderson, 1976
           if (OptSnowThermConduct == 3) &
              ThermConductSnow(I,LoopInd,J) = 0.35                                                        ! constant
           if (OptSnowThermConduct == 4) &
              ThermConductSnow(I,LoopInd,J) = 2.576e-6 * SnowDensBulk**2.0 + 0.074               ! Verseghy (1991)
           if (OptSnowThermConduct == 5) &
              ThermConductSnow(I,LoopInd,J) = 2.22 * (SnowDensBulk/1000.0)**1.88                 ! Douvill(Yen, 1981)
        enddo

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine SnowThermalProperty

end module SnowThermalPropertyMod
