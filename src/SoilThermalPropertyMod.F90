module SoilThermalPropertyMod

!!! Compute soil thermal conductivity based on Peters-Lidard et al. (1998) (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SoilThermalProperty(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: TDFCND
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! If the soil has any moisture content compute a partial sum/product
! otherwise use a constant value which works well with most soils
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd                       ! loop index
    integer                          :: I, J                          ! grid indices
    real(kind=kind_noahmp)           :: KerstenFac                    ! Kersten number
    real(kind=kind_noahmp)           :: SoilGamFac                    ! temporary soil GAMMD factor
    real(kind=kind_noahmp)           :: ThermConductSoilDry           ! thermal conductivity for dry soil
    real(kind=kind_noahmp)           :: ThermConductSoilSat           ! thermal conductivity for saturated soil
    real(kind=kind_noahmp)           :: ThermConductSolid             ! thermal conductivity for the solids
    real(kind=kind_noahmp)           :: SoilSatRatio                  ! saturation ratio
    real(kind=kind_noahmp)           :: SoilWatFracSat                ! saturated soil water fraction
    real(kind=kind_noahmp)           :: SoilWatFrac                   ! soil water fraction
    real(kind=kind_noahmp)           :: SoilIceTmp                    ! temporal soil ice

! --------------------------------------------------------------------
        associate(                                                          &
                  NumSoilLayer     => noahmp%config%domain%NumSoilLayer    ,& ! in,  number of soil layers
                  SoilMoistureSat  => noahmp%water%param%SoilMoistureSat   ,& ! in,  saturated value of soil moisture [m3/m3] (3D)
                  SoilHeatCapacity => noahmp%energy%param%SoilHeatCapacity ,& ! in,  soil volumetric specific heat [J/m3/K]
                  SoilQuartzFrac   => noahmp%energy%param%SoilQuartzFrac   ,& ! in,  soil quartz content (3D)
                  SoilMoisture     => noahmp%water%state%SoilMoisture      ,& ! in,  total soil moisture [m3/m3] (3D)
                  SoilLiqWater     => noahmp%water%state%SoilLiqWater      ,& ! in,  soil water content [m3/m3] (3D)
                  HeatCapacVolSoil => noahmp%energy%state%HeatCapacVolSoil ,& ! out, soil layer volumetric specific heat [J/m3/K] (3D)
                  ThermConductSoil => noahmp%energy%state%ThermConductSoil  & ! out, soil layer thermal conductivity [W/m/K] (3D)
                 )

    !$acc parallel loop collapse(2) gang vector default(present) private(KerstenFac, LoopInd, SoilGamFac, SoilIceTmp, &
    !$acc SoilSatRatio, SoilWatFrac, SoilWatFracSat, ThermConductSoilDry, ThermConductSoilSat, ThermConductSolid)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        !$acc loop seq
        do LoopInd = 1, NumSoilLayer

           ! ==== soil heat capacity
           SoilIceTmp                        = SoilMoisture(I,LoopInd,J) - SoilLiqWater(I,LoopInd,J)
           HeatCapacVolSoil(I,LoopInd,J) = SoilLiqWater(I,LoopInd,J) * ConstHeatCapacWater +                            &
                                       (1.0 - SoilMoistureSat(I,LoopInd,J)) * SoilHeatCapacity(I,J) +                    &
                                       (SoilMoistureSat(I,LoopInd,J) - SoilMoisture(I,LoopInd,J)) * ConstHeatCapacAir + &
                                       SoilIceTmp * ConstHeatCapacIce

           ! ==== soil thermal conductivity
           SoilSatRatio = SoilMoisture(I,LoopInd,J) / SoilMoistureSat(I,LoopInd,J) ! SATURATION RATIO

           ! UNFROZEN FRACTION (FROM 1., i.e., 100%LIQUID, TO 0. (100% FROZEN))
           ThermConductSolid = (ConstThermConductQuartz ** SoilQuartzFrac(I,LoopInd,J)) * &
                               (ConstThermConductSoilOth ** (1.0 - SoilQuartzFrac(I,LoopInd,J)))

           ! UNFROZEN VOLUME FOR SATURATION (POROSITY*SoilWatFrac)
           SoilWatFrac = 1.0    ! Prevent divide by zero (suggested by D. Mocko)
           if ( SoilMoisture(I,LoopInd,J) > 0.0 ) SoilWatFrac = SoilLiqWater(I,LoopInd,J) / SoilMoisture(I,LoopInd,J)
           SoilWatFracSat = SoilWatFrac * SoilMoistureSat(I,LoopInd,J)

           ! SATURATED THERMAL CONDUCTIVITY
           ThermConductSoilSat = ThermConductSolid ** (1.0-SoilMoistureSat(I,LoopInd,J)) * &
                                 ConstThermConductIce ** (SoilMoistureSat(I,LoopInd,J)-SoilWatFracSat) * &
                                 ConstThermConductWater ** (SoilWatFracSat)

           ! DRY THERMAL CONDUCTIVITY IN W.M-1.K-1
           SoilGamFac          = (1.0 - SoilMoistureSat(I,LoopInd,J)) * 2700.0
           ThermConductSoilDry = (0.135 * SoilGamFac + 64.7) / (2700.0 - 0.947 * SoilGamFac)

           ! THE KERSTEN NUMBER KerstenFac
           if ( (SoilLiqWater(I,LoopInd,J)+0.0005) < SoilMoisture(I,LoopInd,J) ) then ! FROZEN
              KerstenFac = SoilSatRatio
           else  ! UNFROZEN
              ! KERSTEN NUMBER (USING "FINE" FORMULA, VALID FOR SOILS CONTAINING AT
              ! LEAST 5% OF PARTICLES WITH DIAMETER LESS THAN 2.E-6 METERS.)
              ! (FOR "COARSE" FORMULA, SEE PETERS-LIDARD ET AL., 1998).
              if ( SoilSatRatio > 0.1 ) then
                 KerstenFac = log10(SoilSatRatio) + 1.0
              else
                 KerstenFac = 0.0
              endif
           endif

           !  THERMAL CONDUCTIVITY
           ThermConductSoil(I,LoopInd,J) = KerstenFac*(ThermConductSoilSat-ThermConductSoilDry) + ThermConductSoilDry

        enddo ! LoopInd


      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine SoilThermalProperty

end module SoilThermalPropertyMod
