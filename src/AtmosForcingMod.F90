module AtmosForcingMod

!!! Process input atmospheric forcing variables (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine ProcessAtmosForcing(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ATM
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local varibles
    integer                          :: I, J                   ! grid indices
    integer                          :: LoopInd                ! loop index
    integer, parameter               :: LoopNum = 10           ! iterations for Twet calculation
    real(kind=kind_noahmp)           :: PrecipFrozenTot        ! total frozen precipitation [mm/s] ! MB/AN : v3.7
    real(kind=kind_noahmp)           :: RadDirFrac             ! direct solar radiation fraction
    real(kind=kind_noahmp)           :: RadVisFrac             ! visible band solar radiation fraction
    real(kind=kind_noahmp)           :: VapPresSat             ! saturated vapor pressure of air
    real(kind=kind_noahmp)           :: LatHeatVap             ! latent heat of vapor/sublimation
    real(kind=kind_noahmp)           :: PsychConst             ! (cp*p)/(eps*L), psychrometric coefficient
    real(kind=kind_noahmp)           :: TemperatureDegC        ! air temperature [C]
    real(kind=kind_noahmp)           :: TemperatureWetBulb     ! wetbulb temperature

! ------------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp) &
    !$acc private(LoopInd, PrecipFrozenTot, RadDirFrac, RadVisFrac, VapPresSat, &
    !$acc         LatHeatVap, PsychConst, TemperatureDegC, TemperatureWetBulb)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                        &
                  CosSolarZenithAngle     => noahmp%config%domain%CosSolarZenithAngle(I,J)    ,& ! in,  cosine solar zenith angle [0-1]
                  OptRainSnowPartition    => noahmp%config%nmlist%OptRainSnowPartition        ,& ! in,  rain-snow partition physics option
                  PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight(I,J)         ,& ! in,  air pressure [Pa] at reference height
                  TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight(I,J)      ,& ! in,  air temperature [K] at reference height
                  SpecHumidityRefHeight   => noahmp%forcing%SpecHumidityRefHeight(I,J)        ,& ! in,  specific humidity [kg/kg] forcing at reference height
                  PrecipConvRefHeight     => noahmp%forcing%PrecipConvRefHeight(I,J)          ,& ! in,  convective precipitation rate [mm/s] at reference height
                  PrecipNonConvRefHeight  => noahmp%forcing%PrecipNonConvRefHeight(I,J)       ,& ! in,  non-convective precipitation rate [mm/s] at reference height
                  PrecipShConvRefHeight   => noahmp%forcing%PrecipShConvRefHeight(I,J)        ,& ! in,  shallow convective precipitation rate [mm/s] at reference height
                  PrecipSnowRefHeight     => noahmp%forcing%PrecipSnowRefHeight(I,J)          ,& ! in,  snowfall rate [mm/s] at reference height
                  PrecipGraupelRefHeight  => noahmp%forcing%PrecipGraupelRefHeight(I,J)       ,& ! in,  graupel rate [mm/s] at reference height
                  PrecipHailRefHeight     => noahmp%forcing%PrecipHailRefHeight(I,J)          ,& ! in,  hail rate [mm/s] at reference height
                  RadSwDownRefHeight      => noahmp%forcing%RadSwDownRefHeight(I,J)           ,& ! in,  downward shortwave radiation [W/m2] at reference height
                  WindEastwardRefHeight   => noahmp%forcing%WindEastwardRefHeight(I,J)        ,& ! in,  wind speed [m/s] in eastward direction at reference height
                  WindNorthwardRefHeight  => noahmp%forcing%WindNorthwardRefHeight(I,J)       ,& ! in,  wind speed [m/s] in northward direction at reference height
                  RadSwVisFrac            => noahmp%forcing%RadSwVisFrac(I,J)                 ,& ! in,  downward solar radiation visible band fraction
                  RadSwDirFrac            => noahmp%forcing%RadSwDirFrac(I,J)                 ,& ! in,  downward solar radiation direct beam fraction
                  SnowfallDensityMax      => noahmp%water%param%SnowfallDensityMax(I,J)       ,& ! in,  maximum fresh snowfall density [kg/m3]
                  TemperaturePotRefHeight => noahmp%energy%state%TemperaturePotRefHeight(I,J) ,& ! out, surface potential temperature [K]
                  PressureVaporRefHeight  => noahmp%energy%state%PressureVaporRefHeight(I,J)  ,& ! out, vapor pressure air [Pa] at reference height
                  DensityAirRefHeight     => noahmp%energy%state%DensityAirRefHeight(I,J)     ,& ! out, density air [kg/m3]
                  WindSpdRefHeight        => noahmp%energy%state%WindSpdRefHeight(I,J)        ,& ! out, wind speed [m/s] at reference height
                  RadSwDownDir            => noahmp%energy%flux%RadSwDownDir                  ,& ! out, incoming direct solar radiation [W/m2]
                  RadSwDownDif            => noahmp%energy%flux%RadSwDownDif                  ,& ! out, incoming diffuse solar radiation [W/m2]
                  RainfallRefHeight       => noahmp%water%flux%RainfallRefHeight(I,J)         ,& ! out, rainfall [mm/s] at reference height
                  SnowfallRefHeight       => noahmp%water%flux%SnowfallRefHeight(I,J)         ,& ! out, liquid equivalent snowfall [mm/s] at reference height
                  PrecipTotRefHeight      => noahmp%water%flux%PrecipTotRefHeight(I,J)        ,& ! out, total precipitation [mm/s] at reference height
                  PrecipConvTotRefHeight  => noahmp%water%flux%PrecipConvTotRefHeight(I,J)    ,& ! out, total convective precipitation [mm/s] at reference height
                  PrecipLargeSclRefHeight => noahmp%water%flux%PrecipLargeSclRefHeight(I,J)   ,& ! out, large-scale precipitation [mm/s] at reference height
                  PrecipAreaFrac          => noahmp%water%state%PrecipAreaFrac(I,J)           ,& ! out, fraction of area receiving precipitation
                  FrozenPrecipFrac        => noahmp%water%state%FrozenPrecipFrac(I,J)         ,& ! out, frozen precipitation fraction
                  SnowfallDensity         => noahmp%water%state%SnowfallDensity(I,J)           & ! out, bulk density of snowfall [kg/m3]
                 )
! ----------------------------------------------------------------------

    ! surface air variables
    TemperaturePotRefHeight = TemperatureAirRefHeight * &
                              (PressureAirRefHeight / PressureAirRefHeight) ** (ConstGasDryAir / ConstHeatCapacAir) 
    PressureVaporRefHeight  = SpecHumidityRefHeight * PressureAirRefHeight / (0.622 + 0.378*SpecHumidityRefHeight)
    DensityAirRefHeight     = (PressureAirRefHeight - 0.378*PressureVaporRefHeight) / &
                              (ConstGasDryAir * TemperatureAirRefHeight)

    ! downward solar radiation
    RadDirFrac = 0.7
    RadVisFrac = 0.5

    if ( RadSwDirFrac >= 0.0 .and. RadSwDirFrac <= 1.0 ) then
       RadDirFrac   = RadSwDirFrac
    else
       RadSwDirFrac = RadDirFrac
    endif

    if ( RadSwVisFrac >= 0.0 .and. RadSwVisFrac <= 1.0 ) then 
       RadVisFrac   = RadSwVisFrac
    else
       RadSwVisFrac = RadVisFrac
    endif

    if ( CosSolarZenithAngle <= 0.0 ) RadSwDownRefHeight = 0.0                      ! filter by solar zenith angle
    RadSwDownDir(I,1,J) = RadSwDownRefHeight * RadDirFrac       * RadVisFrac        ! direct  vis
    RadSwDownDir(I,2,J) = RadSwDownRefHeight * RadDirFrac       * (1.0-RadVisFrac)  ! direct  nir
    RadSwDownDif(I,1,J) = RadSwDownRefHeight * (1.0-RadDirFrac) * RadVisFrac        ! diffuse vis
    RadSwDownDif(I,2,J) = RadSwDownRefHeight * (1.0-RadDirFrac) * (1.0-RadVisFrac)  ! diffuse nir

    ! precipitation
    PrecipTotRefHeight = PrecipConvRefHeight + PrecipNonConvRefHeight + PrecipShConvRefHeight
    if ( OptRainSnowPartition == 4 ) then
       PrecipConvTotRefHeight  = PrecipConvRefHeight + PrecipShConvRefHeight
       PrecipLargeSclRefHeight = PrecipNonConvRefHeight
    else
       PrecipConvTotRefHeight  = 0.10 * PrecipTotRefHeight
       PrecipLargeSclRefHeight = 0.90 * PrecipTotRefHeight
    endif

    ! fractional area that receives precipitation (see, Niu et al. 2005)
    PrecipAreaFrac = 0.0
    if ( (PrecipConvTotRefHeight+PrecipLargeSclRefHeight) > 0.0 ) then
       PrecipAreaFrac = (PrecipConvTotRefHeight + PrecipLargeSclRefHeight) / &
                        (10.0*PrecipConvTotRefHeight + PrecipLargeSclRefHeight)
    endif

    ! partition precipitation into rain and snow. Moved from CANWAT MB/AN: v3.7
    ! Jordan (1991)
    if ( OptRainSnowPartition == 1 ) then
       if ( TemperatureAirRefHeight > (ConstFreezePoint+2.5) ) then
          FrozenPrecipFrac = 0.0
       else
          if ( TemperatureAirRefHeight <= (ConstFreezePoint+0.5) ) then
             FrozenPrecipFrac = 1.0
          elseif ( TemperatureAirRefHeight <= (ConstFreezePoint+2.0) ) then
             FrozenPrecipFrac = 1.0 - (-54.632 + 0.2*TemperatureAirRefHeight)
          else
             FrozenPrecipFrac = 0.6
          endif
       endif
    endif

    ! BATS scheme
    if ( OptRainSnowPartition == 2 ) then
       if ( TemperatureAirRefHeight >= (ConstFreezePoint+2.2) ) then
          FrozenPrecipFrac = 0.0
       else
          FrozenPrecipFrac = 1.0
       endif
    endif

    ! Simple temperature scheme
    if ( OptRainSnowPartition == 3 ) then
       if ( TemperatureAirRefHeight >= ConstFreezePoint ) then
          FrozenPrecipFrac = 0.0
       else
          FrozenPrecipFrac = 1.0
       endif
    endif

    ! Use WRF microphysics output
    ! Hedstrom NR and JW Pomeroy (1998), Hydrol. Processes, 12, 1611-1625
    SnowfallDensity = min( SnowfallDensityMax, 67.92+51.25*exp((TemperatureAirRefHeight-ConstFreezePoint)/2.59) )   ! fresh snow density !MB/AN: change to MIN  
    if ( OptRainSnowPartition == 4 ) then
       PrecipFrozenTot = PrecipSnowRefHeight + PrecipGraupelRefHeight + PrecipHailRefHeight
       if ( (PrecipNonConvRefHeight > 0.0) .and. (PrecipFrozenTot > 0.0) ) then
          FrozenPrecipFrac  = min( 1.0, PrecipFrozenTot/PrecipNonConvRefHeight )
          FrozenPrecipFrac  = max( 0.0, FrozenPrecipFrac )
          SnowfallDensity   = SnowfallDensity     * (PrecipSnowRefHeight/PrecipFrozenTot)    + &
                              ConstDensityGraupel * (PrecipGraupelRefHeight/PrecipFrozenTot) + &
                              ConstDensityHail    * (PrecipHailRefHeight/PrecipFrozenTot)
       else
          FrozenPrecipFrac  = 0.0
       endif
    endif

    ! wet-bulb scheme (Wang et al., 2019 GRL), C.He, 12/18/2020, R. Abolafia-Rosnezweig, 02/01/2024
    if ( OptRainSnowPartition == 5 ) then

        if ( TemperatureAirRefHeight >= (ConstFreezePoint+8) ) then !avoid numerical errors when temperature is high
            FrozenPrecipFrac = 0.0
        else
            TemperatureDegC = min( 50.0, max(-50.0,(TemperatureAirRefHeight-ConstFreezePoint)) )    ! Kelvin to degree Celsius with limit -50 to +50
            if ( TemperatureAirRefHeight > ConstFreezePoint ) then
                LatHeatVap = ConstLatHeatEvap
            else
                LatHeatVap = ConstLatHeatSublim
            endif
            PsychConst            = ConstHeatCapacAir * PressureAirRefHeight / (0.622 * LatHeatVap)
            TemperatureWetBulb    = TemperatureDegC - 5.0    ! first guess wetbulb temperature
            !$acc loop seq
            do LoopInd = 1, LoopNum
                VapPresSat         = 610.8 * exp( (17.27*TemperatureWetBulb) / (237.3+TemperatureWetBulb) )
                TemperatureWetBulb = TemperatureWetBulb - (VapPresSat - PressureVaporRefHeight) / PsychConst   ! Wang et al., 2019 GRL Eq.2
            enddo
            FrozenPrecipFrac      = 1.0 / (1.0 + 6.99e-5 * exp(2.0*(TemperatureWetBulb+3.97)))                ! Wang et al., 2019 GRL Eq. 1
        endif
    endif

    ! rain-snow partitioning
    RainfallRefHeight = PrecipTotRefHeight * (1.0 - FrozenPrecipFrac)
    SnowfallRefHeight = PrecipTotRefHeight * FrozenPrecipFrac

    ! wind speed at reference height for turbulence calculation
    WindSpdRefHeight = max(sqrt(WindEastwardRefHeight**2.0 + WindNorthwardRefHeight**2.0), 1.0)

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine ProcessAtmosForcing

end module AtmosForcingMod
