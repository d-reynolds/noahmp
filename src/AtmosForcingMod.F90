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
        associate(                                                                        &
                  CosSolarZenithAngle     => noahmp%config%domain%CosSolarZenithAngle    ,& ! in,  cosine solar zenith angle [0-1]
                  OptRainSnowPartition    => noahmp%config%nmlist%OptRainSnowPartition        ,& ! in,  rain-snow partition physics option
                  PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight         ,& ! in,  air pressure [Pa] at reference height
                  TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight      ,& ! in,  air temperature [K] at reference height
                  SpecHumidityRefHeight   => noahmp%forcing%SpecHumidityRefHeight        ,& ! in,  specific humidity [kg/kg] forcing at reference height
                  PrecipConvRefHeight     => noahmp%forcing%PrecipConvRefHeight          ,& ! in,  convective precipitation rate [mm/s] at reference height
                  PrecipNonConvRefHeight  => noahmp%forcing%PrecipNonConvRefHeight       ,& ! in,  non-convective precipitation rate [mm/s] at reference height
                  PrecipShConvRefHeight   => noahmp%forcing%PrecipShConvRefHeight        ,& ! in,  shallow convective precipitation rate [mm/s] at reference height
                  PrecipSnowRefHeight     => noahmp%forcing%PrecipSnowRefHeight          ,& ! in,  snowfall rate [mm/s] at reference height
                  PrecipGraupelRefHeight  => noahmp%forcing%PrecipGraupelRefHeight       ,& ! in,  graupel rate [mm/s] at reference height
                  PrecipHailRefHeight     => noahmp%forcing%PrecipHailRefHeight          ,& ! in,  hail rate [mm/s] at reference height
                  RadSwDownRefHeight      => noahmp%forcing%RadSwDownRefHeight           ,& ! in,  downward shortwave radiation [W/m2] at reference height
                  WindEastwardRefHeight   => noahmp%forcing%WindEastwardRefHeight        ,& ! in,  wind speed [m/s] in eastward direction at reference height
                  WindNorthwardRefHeight  => noahmp%forcing%WindNorthwardRefHeight       ,& ! in,  wind speed [m/s] in northward direction at reference height
                  RadSwVisFrac            => noahmp%forcing%RadSwVisFrac                 ,& ! in,  downward solar radiation visible band fraction
                  RadSwDirFrac            => noahmp%forcing%RadSwDirFrac                 ,& ! in,  downward solar radiation direct beam fraction
                  SnowfallDensityMax      => noahmp%water%param%SnowfallDensityMax       ,& ! in,  maximum fresh snowfall density [kg/m3]
                  TemperaturePotRefHeight => noahmp%energy%state%TemperaturePotRefHeight ,& ! out, surface potential temperature [K]
                  PressureVaporRefHeight  => noahmp%energy%state%PressureVaporRefHeight  ,& ! out, vapor pressure air [Pa] at reference height
                  DensityAirRefHeight     => noahmp%energy%state%DensityAirRefHeight     ,& ! out, density air [kg/m3]
                  WindSpdRefHeight        => noahmp%energy%state%WindSpdRefHeight        ,& ! out, wind speed [m/s] at reference height
                  RadSwDownDir            => noahmp%energy%flux%RadSwDownDir                  ,& ! out, incoming direct solar radiation [W/m2]
                  RadSwDownDif            => noahmp%energy%flux%RadSwDownDif                  ,& ! out, incoming diffuse solar radiation [W/m2]
                  RainfallRefHeight       => noahmp%water%flux%RainfallRefHeight         ,& ! out, rainfall [mm/s] at reference height
                  SnowfallRefHeight       => noahmp%water%flux%SnowfallRefHeight         ,& ! out, liquid equivalent snowfall [mm/s] at reference height
                  PrecipTotRefHeight      => noahmp%water%flux%PrecipTotRefHeight        ,& ! out, total precipitation [mm/s] at reference height
                  PrecipConvTotRefHeight  => noahmp%water%flux%PrecipConvTotRefHeight    ,& ! out, total convective precipitation [mm/s] at reference height
                  PrecipLargeSclRefHeight => noahmp%water%flux%PrecipLargeSclRefHeight   ,& ! out, large-scale precipitation [mm/s] at reference height
                  PrecipAreaFrac          => noahmp%water%state%PrecipAreaFrac           ,& ! out, fraction of area receiving precipitation
                  FrozenPrecipFrac        => noahmp%water%state%FrozenPrecipFrac         ,& ! out, frozen precipitation fraction
                  SnowfallDensity         => noahmp%water%state%SnowfallDensity           & ! out, bulk density of snowfall [kg/m3]
                 )

    !$acc parallel loop collapse(2) gang vector default(present) &
    !$acc private(LoopInd, PrecipFrozenTot, RadDirFrac, RadVisFrac, VapPresSat, &
    !$acc         LatHeatVap, PsychConst, TemperatureDegC, TemperatureWetBulb)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! surface air variables
    TemperaturePotRefHeight(I,J) = TemperatureAirRefHeight(I,J) * &
                              (PressureAirRefHeight(I,J) / PressureAirRefHeight(I,J)) ** (ConstGasDryAir / ConstHeatCapacAir) 
    PressureVaporRefHeight(I,J)  = SpecHumidityRefHeight(I,J) * PressureAirRefHeight(I,J) / (0.622 + 0.378*SpecHumidityRefHeight(I,J))
    DensityAirRefHeight(I,J)     = (PressureAirRefHeight(I,J) - 0.378*PressureVaporRefHeight(I,J)) / &
                              (ConstGasDryAir * TemperatureAirRefHeight(I,J))

    ! downward solar radiation
    RadDirFrac = 0.7
    RadVisFrac = 0.5

    if ( RadSwDirFrac(I,J) >= 0.0 .and. RadSwDirFrac(I,J) <= 1.0 ) then
       RadDirFrac   = RadSwDirFrac(I,J)
    else
       RadSwDirFrac(I,J) = RadDirFrac
    endif

    if ( RadSwVisFrac(I,J) >= 0.0 .and. RadSwVisFrac(I,J) <= 1.0 ) then 
       RadVisFrac   = RadSwVisFrac(I,J)
    else
       RadSwVisFrac(I,J) = RadVisFrac
    endif

    if ( CosSolarZenithAngle(I,J) <= 0.0 ) RadSwDownRefHeight(I,J) = 0.0                      ! filter by solar zenith angle
    RadSwDownDir(I,1,J) = RadSwDownRefHeight(I,J) * RadDirFrac       * RadVisFrac        ! direct  vis
    RadSwDownDir(I,2,J) = RadSwDownRefHeight(I,J) * RadDirFrac       * (1.0-RadVisFrac)  ! direct  nir
    RadSwDownDif(I,1,J) = RadSwDownRefHeight(I,J) * (1.0-RadDirFrac) * RadVisFrac        ! diffuse vis
    RadSwDownDif(I,2,J) = RadSwDownRefHeight(I,J) * (1.0-RadDirFrac) * (1.0-RadVisFrac)  ! diffuse nir

    ! precipitation
    PrecipTotRefHeight(I,J) = PrecipConvRefHeight(I,J) + PrecipNonConvRefHeight(I,J) + PrecipShConvRefHeight(I,J)
    if ( OptRainSnowPartition == 4 ) then
       PrecipConvTotRefHeight(I,J)  = PrecipConvRefHeight(I,J) + PrecipShConvRefHeight(I,J)
       PrecipLargeSclRefHeight(I,J) = PrecipNonConvRefHeight(I,J)
    else
       PrecipConvTotRefHeight(I,J)  = 0.10 * PrecipTotRefHeight(I,J)
       PrecipLargeSclRefHeight(I,J) = 0.90 * PrecipTotRefHeight(I,J)
    endif

    ! fractional area that receives precipitation (see, Niu et al. 2005)
    PrecipAreaFrac(I,J) = 0.0
    if ( (PrecipConvTotRefHeight(I,J)+PrecipLargeSclRefHeight(I,J)) > 0.0 ) then
       PrecipAreaFrac(I,J) = (PrecipConvTotRefHeight(I,J) + PrecipLargeSclRefHeight(I,J)) / &
                        (10.0*PrecipConvTotRefHeight(I,J) + PrecipLargeSclRefHeight(I,J))
    endif

    ! partition precipitation into rain and snow. Moved from CANWAT MB/AN: v3.7
    ! Jordan (1991)
    if ( OptRainSnowPartition == 1 ) then
       if ( TemperatureAirRefHeight(I,J) > (ConstFreezePoint+2.5) ) then
          FrozenPrecipFrac(I,J) = 0.0
       else
          if ( TemperatureAirRefHeight(I,J) <= (ConstFreezePoint+0.5) ) then
             FrozenPrecipFrac(I,J) = 1.0
          elseif ( TemperatureAirRefHeight(I,J) <= (ConstFreezePoint+2.0) ) then
             FrozenPrecipFrac(I,J) = 1.0 - (-54.632 + 0.2*TemperatureAirRefHeight(I,J))
          else
             FrozenPrecipFrac(I,J) = 0.6
          endif
       endif
    endif

    ! BATS scheme
    if ( OptRainSnowPartition == 2 ) then
       if ( TemperatureAirRefHeight(I,J) >= (ConstFreezePoint+2.2) ) then
          FrozenPrecipFrac(I,J) = 0.0
       else
          FrozenPrecipFrac(I,J) = 1.0
       endif
    endif

    ! Simple temperature scheme
    if ( OptRainSnowPartition == 3 ) then
       if ( TemperatureAirRefHeight(I,J) >= ConstFreezePoint ) then
          FrozenPrecipFrac(I,J) = 0.0
       else
          FrozenPrecipFrac(I,J) = 1.0
       endif
    endif

    ! Use WRF microphysics output
    ! Hedstrom NR and JW Pomeroy (1998), Hydrol. Processes, 12, 1611-1625
    SnowfallDensity(I,J) = min( SnowfallDensityMax(I,J), 67.92+51.25*exp((TemperatureAirRefHeight(I,J)-ConstFreezePoint)/2.59) )   ! fresh snow density !MB/AN: change to MIN  
    if ( OptRainSnowPartition == 4 ) then
       PrecipFrozenTot = PrecipSnowRefHeight(I,J) + PrecipGraupelRefHeight(I,J) + PrecipHailRefHeight(I,J)
       if ( (PrecipNonConvRefHeight(I,J) > 0.0) .and. (PrecipFrozenTot > 0.0) ) then
          FrozenPrecipFrac(I,J)  = min( 1.0, PrecipFrozenTot/PrecipNonConvRefHeight(I,J) )
          FrozenPrecipFrac(I,J)  = max( 0.0, FrozenPrecipFrac(I,J) )
          SnowfallDensity(I,J)   = SnowfallDensity(I,J)     * (PrecipSnowRefHeight(I,J)/PrecipFrozenTot)    + &
                              ConstDensityGraupel * (PrecipGraupelRefHeight(I,J)/PrecipFrozenTot) + &
                              ConstDensityHail    * (PrecipHailRefHeight(I,J)/PrecipFrozenTot)
       else
          FrozenPrecipFrac(I,J)  = 0.0
       endif
    endif

    ! wet-bulb scheme (Wang et al., 2019 GRL), C.He, 12/18/2020, R. Abolafia-Rosnezweig, 02/01/2024
    if ( OptRainSnowPartition == 5 ) then

        if ( TemperatureAirRefHeight(I,J) >= (ConstFreezePoint+8) ) then !avoid numerical errors when temperature is high
            FrozenPrecipFrac(I,J) = 0.0
        else
            TemperatureDegC = min( 50.0, max(-50.0,(TemperatureAirRefHeight(I,J)-ConstFreezePoint)) )    ! Kelvin to degree Celsius with limit -50 to +50
            if ( TemperatureAirRefHeight(I,J) > ConstFreezePoint ) then
                LatHeatVap = ConstLatHeatEvap
            else
                LatHeatVap = ConstLatHeatSublim
            endif
            PsychConst            = ConstHeatCapacAir * PressureAirRefHeight(I,J) / (0.622 * LatHeatVap)
            TemperatureWetBulb    = TemperatureDegC - 5.0    ! first guess wetbulb temperature
            !$acc loop seq
            do LoopInd = 1, LoopNum
                VapPresSat         = 610.8 * exp( (17.27*TemperatureWetBulb) / (237.3+TemperatureWetBulb) )
                TemperatureWetBulb = TemperatureWetBulb - (VapPresSat - PressureVaporRefHeight(I,J)) / PsychConst   ! Wang et al., 2019 GRL Eq.2
            enddo
            FrozenPrecipFrac(I,J)      = 1.0 / (1.0 + 6.99e-5 * exp(2.0*(TemperatureWetBulb+3.97)))                ! Wang et al., 2019 GRL Eq. 1
        endif
    endif

    ! rain-snow partitioning
    RainfallRefHeight(I,J) = PrecipTotRefHeight(I,J) * (1.0 - FrozenPrecipFrac(I,J))
    SnowfallRefHeight(I,J) = PrecipTotRefHeight(I,J) * FrozenPrecipFrac(I,J)

    ! wind speed at reference height for turbulence calculation
    WindSpdRefHeight(I,J) = max(sqrt(WindEastwardRefHeight(I,J)**2.0 + WindNorthwardRefHeight(I,J)**2.0), 1.0)


      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine ProcessAtmosForcing

end module AtmosForcingMod
