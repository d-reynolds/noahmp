module CropPhotosynthesisMod

!!! Compute crop photosynthesis

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine CropPhotosynthesis(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: PSN_CROP  
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------
        
    implicit none
       
    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J             ! grid indices
    real(kind=kind_noahmp)           :: PhotosynRad      ! photosynthetically active radiation (w/m2) 1 W m-2 = 0.0864 MJ m-2 day-1
    real(kind=kind_noahmp)           :: Co2AssimMax      ! Maximum CO2 assimulation rate g CO2/m2/s  
    real(kind=kind_noahmp)           :: Co2AssimTot      ! CO2 Assimilation g CO2/m2/s
    real(kind=kind_noahmp)           :: TemperatureAirC  ! air temperature degC
    real(kind=kind_noahmp)           :: L1               ! Three Gaussian method
    real(kind=kind_noahmp)           :: L2               ! Three Gaussian method
    real(kind=kind_noahmp)           :: L3               ! Three Gaussian method
    real(kind=kind_noahmp)           :: I1               ! Three Gaussian method
    real(kind=kind_noahmp)           :: I2               ! Three Gaussian method
    real(kind=kind_noahmp)           :: I3               ! Three Gaussian method
    real(kind=kind_noahmp)           :: A1               ! Three Gaussian method
    real(kind=kind_noahmp)           :: A2               ! Three Gaussian method
    real(kind=kind_noahmp)           :: A3               ! Three Gaussian method

    associate(                                                                     &
              RadSwDownRefHeight    => noahmp%forcing%RadSwDownRefHeight     ,& ! in,  downward shortwave radiation [W/m2] at reference height
              TemperatureAir2m      => noahmp%energy%state%TemperatureAir2m  ,& ! in,  2-m air temperature [K]
              LeafAreaIndex         => noahmp%energy%state%LeafAreaIndex     ,& ! in,  leaf area index, unadjusted for burying by snow
              PhotosynRadFrac       => noahmp%biochem%param%PhotosynRadFrac  ,& ! in,  Fraction of incoming radiation to photosynthetically active radiation
              TempMinCarbonAssim    => noahmp%biochem%param%TempMinCarbonAssim,& ! in,  Minimum temperature for CO2 assimilation [C]
              TempMaxCarbonAssim    => noahmp%biochem%param%TempMaxCarbonAssim,& ! in,  CO2 assim. linearly increasing until reaching this temperature [C]
              TempMaxCarbonAssimMax => noahmp%biochem%param%TempMaxCarbonAssimMax,& ! in,  CO2 assim. remain at CarbonAssimRefMax until reaching this temperature [C]
              CarbonAssimRefMax     => noahmp%biochem%param%CarbonAssimRefMax,& ! in,  reference maximum CO2 assimilation rate
              LightExtCoeff         => noahmp%biochem%param%LightExtCoeff    ,& ! in,  light extinction coefficient
              LightUseEfficiency    => noahmp%biochem%param%LightUseEfficiency,& ! in,  initial light use efficiency
              CarbonAssimReducFac   => noahmp%biochem%param%CarbonAssimReducFac,& ! in,  CO2 assimulation reduction factor(0-1) (caused by e.g.pest,weeds)
              PhotosynCrop          => noahmp%biochem%flux%PhotosynCrop       & ! out, crop photosynthesis [umol co2/m2/s]
             )

   !$acc parallel loop collapse(2) gang vector default(present) &
   !$acc private(PhotosynRad, Co2AssimMax, Co2AssimTot, TemperatureAirC, L1, L2, L3, I1, I2, I3, A1, A2, A3)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! initialize
    TemperatureAirC = TemperatureAir2m(I,J) - 273.15
    PhotosynRad     = PhotosynRadFrac(I,J) * RadSwDownRefHeight(I,J) * 0.0036  !w to MJ m-2

    ! compute Maximum CO2 assimulation rate g/co2/s
    if ( TemperatureAirC < TempMinCarbonAssim(I,J) ) then
       Co2AssimMax = 1.0e-10
    elseif ( (TemperatureAirC >= TempMinCarbonAssim(I,J)) .and. (TemperatureAirC < TempMaxCarbonAssim(I,J)) ) then
       Co2AssimMax = (TemperatureAirC - TempMinCarbonAssim(I,J)) * CarbonAssimRefMax(I,J) / (TempMaxCarbonAssim(I,J) - TempMinCarbonAssim(I,J))
    elseif ( (TemperatureAirC >= TempMaxCarbonAssim(I,J)) .and. (TemperatureAirC < TempMaxCarbonAssimMax(I,J)) ) then
       Co2AssimMax = CarbonAssimRefMax(I,J)
    else
       Co2AssimMax = CarbonAssimRefMax(I,J) - 0.2 * (TemperatureAir2m(I,J) - TempMaxCarbonAssimMax(I,J))
    endif              
    Co2AssimMax    = max(Co2AssimMax, 0.01)

    ! compute coefficients
    if ( LeafAreaIndex(I,J) <= 0.05 ) then
       L1 = 0.1127 * 0.05   ! use initial LeafAreaIndex(0.05), avoid error
       L2 = 0.5    * 0.05
       L3 = 0.8873 * 0.05
    else
       L1 = 0.1127 * LeafAreaIndex(I,J)
       L2 = 0.5    * LeafAreaIndex(I,J)
       L3 = 0.8873 * LeafAreaIndex(I,J)
    endif

    I1 = LightExtCoeff(I,J) * PhotosynRad * exp(-LightExtCoeff(I,J) * L1)
    I2 = LightExtCoeff(I,J) * PhotosynRad * exp(-LightExtCoeff(I,J) * L2)
    I3 = LightExtCoeff(I,J) * PhotosynRad * exp(-LightExtCoeff(I,J) * L3)
    I1 = max(I1, 1.0e-10)
    I2 = max(I2, 1.0e-10)
    I3 = max(I3, 1.0e-10)
    A1 = Co2AssimMax * (1 - exp(-LightUseEfficiency(I,J) * I1 / Co2AssimMax))
    A2 = Co2AssimMax * (1 - exp(-LightUseEfficiency(I,J) * I2 / Co2AssimMax)) * 1.6
    A3 = Co2AssimMax * (1 - exp(-LightUseEfficiency(I,J) * I3 / Co2AssimMax))

    ! compute photosynthesis rate
    if ( LeafAreaIndex(I,J) <= 0.05 ) then
       Co2AssimTot = (A1 + A2 + A3) / 3.6 * 0.05
    elseif ( (LeafAreaIndex(I,J) > 0.05) .and. (LeafAreaIndex(I,J) <= 4.0) ) then
       Co2AssimTot = (A1 + A2 + A3) / 3.6 * LeafAreaIndex(I,J)
    else
       Co2AssimTot = (A1 + A2 + A3) / 3.6 * 4
    endif    
    Co2AssimTot    = Co2AssimTot * CarbonAssimReducFac(I,J) ! Attainable 
    PhotosynCrop(I,J)   = 6.313 * Co2AssimTot               ! (1/44) * 1000000)/3600 = 6.313


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine CropPhotosynthesis

end module CropPhotosynthesisMod
