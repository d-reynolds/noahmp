module IrrigationSprinklerMod

!!! Estimate irrigation water depth (m) based on sprinkler method 
!!! Reference: chapter 11 of NRCS, Part 623 National Engineering Handbook. 
!!! Irrigation water will be applied over the canopy, affecting  present soil moisture, 
!!! infiltration rate of the soil, and evaporative loss, which should be executed before canopy process.
 
  use Machine
  use CheckNanMod
  use NoahmpVarType
  use ConstantDefineMod
  use IrrigationInfilPhilipMod, only : IrrigationInfilPhilip

  implicit none

contains

  subroutine IrrigationSprinkler(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: SPRINKLER_IRRIGATION
! Original code: P. Valayamkunnath (NCAR) <prasanth@ucar.edu> (08/06/2020)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J              ! grid indices
    integer                          :: LoopInd           ! loop index
    logical                          :: FlagNan           ! NaN value flag: if NaN, return true
    real(kind=kind_noahmp) :: InfilRateSfc(noahmp%config%domain%ITS:noahmp%config%domain%ITE,noahmp%config%domain%JTS:noahmp%config%domain%JTE)   ! surface infiltration rate [m/s]
    real(kind=kind_noahmp)           :: IrriRateTmp       ! temporary irrigation rate [m/timestep]
    real(kind=kind_noahmp)           :: WindSpdTot        ! total wind speed [m/s]
    real(kind=kind_noahmp)           :: IrriLossTmp       ! temporary irrigation water loss [%]
    real(kind=kind_noahmp)           :: PressureVaporSat  ! satuarated vapor pressure [Pa]

    !$acc data create(InfilRateSfc)
    ! estimate infiltration rate based on Philips Eq.
    call IrrigationInfilPhilip(noahmp, noahmp%config%domain%MainTimeStep, InfilRateSfc)

    associate(                                                                       &
              NumSoilLayer            => noahmp%config%domain%NumSoilLayer          ,& ! in,    number of soil layers
              MainTimeStep            => noahmp%config%domain%MainTimeStep          ,& ! in,    noahmp main time step [s]
              TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight,& ! in,    air temperature [K] at reference height
              WindEastwardRefHeight   => noahmp%forcing%WindEastwardRefHeight  ,& ! in,    wind speed [m/s] in eastward direction at reference height
              WindNorthwardRefHeight  => noahmp%forcing%WindNorthwardRefHeight ,& ! in,    wind speed [m/s] in northward direction at reference height
              PressureVaporRefHeight  => noahmp%energy%state%PressureVaporRefHeight,& ! in,    vapor pressure air [Pa]
              IrriSprinklerRate       => noahmp%water%param%IrriSprinklerRate  ,& ! in,    sprinkler irrigation rate [mm/h]
              IrrigationFracSprinkler => noahmp%water%state%IrrigationFracSprinkler,& ! in,    sprinkler irrigation fraction (0 to 1)
              SoilMoisture            => noahmp%water%state%SoilMoisture            ,& ! in,    total soil moisture [m3/m3]
              SoilLiqWater            => noahmp%water%state%SoilLiqWater            ,& ! in,    soil water content [m3/m3]
              HeatLatentIrriEvap      => noahmp%energy%flux%HeatLatentIrriEvap ,& ! inout, latent heating due to sprinkler evaporation [W/m2]
              EvapIrriSprinkler       => noahmp%water%flux%EvapIrriSprinkler   ,& ! inout, evaporation of irrigation water, sprinkler [mm/s]
              RainfallRefHeight       => noahmp%water%flux%RainfallRefHeight   ,& ! inout, rainfall [mm/s] at reference height
              IrrigationRateSprinkler => noahmp%water%flux%IrrigationRateSprinkler,& ! inout, rate of irrigation by sprinkler [m/timestep]
              IrriEvapLossSprinkler   => noahmp%water%flux%IrriEvapLossSprinkler,& ! inout, loss of irrigation water to evaporation,sprinkler [m/timestep]
              IrrigationAmtSprinkler  => noahmp%water%state%IrrigationAmtSprinkler,& ! inout, irrigation water amount [m] to be applied, Sprinkler
              PrecipAreaFrac          => noahmp%water%state%PrecipAreaFrac     ,& ! inout, fraction of area receiving precipitation
              SoilIce                 => noahmp%water%state%SoilIce                  & ! out,   soil ice content [m3/m3]
             )

   !$acc parallel loop collapse(2) gang vector default(present) &
   !$acc private(LoopInd, FlagNan, IrriRateTmp, WindSpdTot, IrriLossTmp, PressureVaporSat)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         ! skip if not cropland or no irrigation water left
         if ( .not.(noahmp%config%domain%FlagCropland(I,J) .and. (noahmp%water%state%IrrigationAmtSprinkler(I,J) > 0.0)) ) cycle 

    ! initialize
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       SoilIce(I,LoopInd,J) = max(0.0, SoilMoisture(I,LoopInd,J)-SoilLiqWater(I,LoopInd,J))
    enddo

    ! irrigation rate of sprinkler
    IrriRateTmp             = IrriSprinklerRate(I,J) * (1.0/1000.0) * MainTimeStep / 3600.0              ! NRCS rate/time step - calibratable
    IrrigationRateSprinkler(I,J) = min(InfilRateSfc(I,J)*MainTimeStep, IrrigationAmtSprinkler(I,J), IrriRateTmp)   ! Limit irrigation rate to minimum of infiltration rate
                                                                                                    ! and to the NRCS recommended rate
    ! evaporative loss from droplets: Based on Bavi et al., (2009). Evaporation 
    ! losses from sprinkler irrigation systems under various operating 
    ! conditions. Journal of Applied Sciences, 9(3), 597-600.
    WindSpdTot       = sqrt((WindEastwardRefHeight(I,J)**2.0) + (WindNorthwardRefHeight(I,J)**2.0))
    PressureVaporSat = 610.8 * exp((17.27*(TemperatureAirRefHeight(I,J)-273.15)) / (237.3+(TemperatureAirRefHeight(I,J)-273.15)))

    if ( TemperatureAirRefHeight(I,J) > 273.15 ) then ! Equation (3)
       IrriLossTmp   = 4.375 * (exp(0.106*WindSpdTot)) * (((PressureVaporSat-PressureVaporRefHeight(I,J))*0.01)**(-0.092)) * &
                       ((TemperatureAirRefHeight(I,J)-273.15)**(-0.102))
    else ! Equation (4)
       IrriLossTmp   = 4.337 * (exp(0.077*WindSpdTot)) * (((PressureVaporSat-PressureVaporRefHeight(I,J))*0.01)**(-0.098))
    endif
    ! Old PGI Fortran compiler does not support ISNAN function
    call CheckRealNaN(IrriLossTmp, FlagNan)
    if ( FlagNan .eqv. .true. ) IrriLossTmp = 4.0                           ! In case if IrriLossTmp is NaN
    if ( (IrriLossTmp > 100.0) .or. (IrriLossTmp < 0.0) ) IrriLossTmp = 4.0 ! In case if IrriLossTmp is out of range

    ! Sprinkler water [m] for sprinkler fraction 
    IrrigationRateSprinkler(I,J)    = IrrigationRateSprinkler(I,J) * IrrigationFracSprinkler(I,J)
    if ( IrrigationRateSprinkler(I,J) >= IrrigationAmtSprinkler(I,J) ) then
       IrrigationRateSprinkler(I,J) = IrrigationAmtSprinkler(I,J)
       IrrigationAmtSprinkler(I,J)  = 0.0
    else
       IrrigationAmtSprinkler(I,J)  = IrrigationAmtSprinkler(I,J) - IrrigationRateSprinkler(I,J)
    endif

    IrriEvapLossSprinkler(I,J)      = IrrigationRateSprinkler(I,J) * IrriLossTmp * (1.0/100.0)
    IrrigationRateSprinkler(I,J)    = IrrigationRateSprinkler(I,J) - IrriEvapLossSprinkler(I,J)

    ! include sprinkler water to total rain for canopy process later
    RainfallRefHeight(I,J)  = RainfallRefHeight(I,J) + (IrrigationRateSprinkler(I,J) * 1000.0 / MainTimeStep)
    ! assuming sprinkler irrigated water is like large scale precipitation, so PrecipAreaFrac = 1.0
    PrecipAreaFrac(I,J) = 1.0

    ! cooling and humidification due to sprinkler evaporation, per m^2 calculation 
    HeatLatentIrriEvap(I,J) = IrriEvapLossSprinkler(I,J) * 1000.0 * ConstLatHeatEvap / MainTimeStep   ! heat used for evaporation [W/m2]
    EvapIrriSprinkler(I,J)  = IrriEvapLossSprinkler(I,J) * 1000.0 / MainTimeStep                      ! sprinkler evaporation [mm/s]


      end do
    end do
   !$acc end parallel loop
   !$acc end data

    end associate

  end subroutine IrrigationSprinkler

end module IrrigationSprinklerMod
