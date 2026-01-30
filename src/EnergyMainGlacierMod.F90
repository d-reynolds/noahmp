module EnergyMainGlacierMod

!!! Main energy module for glacier points including all energy relevant processes
!!! snow thermal property -> radiation -> ground heat flux -> snow temperature solver -> snow/ice phase change

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowCoverGlacierMod,                   only : SnowCoverGlacier
  use GroundRoughnessPropertyGlacierMod,     only : GroundRoughnessPropertyGlacier
  use GroundThermalPropertyGlacierMod,       only : GroundThermalPropertyGlacier
  use SurfaceAlbedoGlacierMod,               only : SurfaceAlbedoGlacier
  use SurfaceRadiationGlacierMod,            only : SurfaceRadiationGlacier
  use SurfaceEmissivityGlacierMod,           only : SurfaceEmissivityGlacier
  use ResistanceGroundEvaporationGlacierMod, only : ResistanceGroundEvaporationGlacier
  use PsychrometricVariableGlacierMod,       only : PsychrometricVariableGlacier
  use SurfaceEnergyFluxGlacierMod,           only : SurfaceEnergyFluxGlacier
  use GlacierTemperatureMainMod,             only : GlacierTemperatureMain
  use GlacierPhaseChangeMod,                 only : GlacierPhaseChange

  implicit none

contains

  subroutine EnergyMainGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ENERGY_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer :: I, J  ! grid indices

! --------------------------------------------------------------------

    ! glaicer snow cover fraction
    call SnowCoverGlacier(noahmp)

    ! ground and surface roughness length and reference height
    call GroundRoughnessPropertyGlacier(noahmp)

    ! Thermal properties of snow and glacier ice
    call GroundThermalPropertyGlacier(noahmp)

    ! Glacier surface shortwave abeldo
    call SurfaceAlbedoGlacier(noahmp)

    ! Glacier surface shortwave radiation
    call SurfaceRadiationGlacier(noahmp)

    ! longwave emissivity for glacier surface
    call SurfaceEmissivityGlacier(noahmp)

    ! glacier surface resistance for ground evaporation/sublimation
    call ResistanceGroundEvaporationGlacier(noahmp)

    ! set psychrometric variable/constant
    call PsychrometricVariableGlacier(noahmp)

    ! temperatures and energy fluxes of glacier ground
    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         associate(                                                  &
              TemperatureGrd         => noahmp%energy%state%TemperatureGrd(I,J),         & ! inout, ground temperature [K]
              TemperatureGrdBare    => noahmp%energy%state%TemperatureGrdBare(I,J),     & ! inout, bare ground temperature [K]
              ExchCoeffMomSfc       => noahmp%energy%state%ExchCoeffMomSfc(I,J),       & ! inout, exchange coefficient [m/s] for momentum, surface, grid mean
              ExchCoeffMomBare      => noahmp%energy%state%ExchCoeffMomBare(I,J),      & ! out,   exchange coefficient [m/s] for momentum, bare ground
              ExchCoeffShSfc        => noahmp%energy%state%ExchCoeffShSfc(I,J),         & ! inout, exchange coefficient [m/s] for heat, surface, grid mean
              ExchCoeffShBare       => noahmp%energy%state%ExchCoeffShBare(I,J)         & ! out,   exchange coefficient [m/s] for heat, bare ground
             )
    TemperatureGrdBare = TemperatureGrd
    ExchCoeffMomBare   = ExchCoeffMomSfc
    ExchCoeffShBare    = ExchCoeffShSfc

    end associate
      end do
    end do

    call SurfaceEnergyFluxGlacier(noahmp)

    ! Grid-level computations requiring 2D parallel loop
    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                         &
                  RadLwDownRefHeight     => noahmp%forcing%RadLwDownRefHeight(I,J)        ,& ! in,    downward longwave radiation [W/m2] at reference height
                  HeatPrecipAdvBareGrd   => noahmp%energy%flux%HeatPrecipAdvBareGrd(I,J)  ,& ! in,    precipitation advected heat - bare ground net [W/m2]
                  TemperatureSfc         => noahmp%energy%state%TemperatureSfc(I,J)       ,& ! inout, surface temperature [K]
                  TemperatureGrd         => noahmp%energy%state%TemperatureGrd(I,J)       ,& ! inout, ground temperature [K]
                  SpecHumiditySfc        => noahmp%energy%state%SpecHumiditySfc(I,J)      ,& ! inout, specific humidity at bare surface
                  SpecHumiditySfcMean    => noahmp%energy%state%SpecHumiditySfcMean(I,J)  ,& ! inout, specific humidity at surface grid mean
                  ExchCoeffMomSfc        => noahmp%energy%state%ExchCoeffMomSfc(I,J)      ,& ! inout, exchange coefficient [m/s] for momentum, surface, grid mean
                  ExchCoeffShSfc         => noahmp%energy%state%ExchCoeffShSfc(I,J)       ,& ! inout, exchange coefficient [m/s] for heat, surface, grid mean
                  SnowDepth              => noahmp%water%state%SnowDepth(I,J)             ,& ! inout, snow depth [m]
                  RoughLenMomSfcToAtm    => noahmp%energy%state%RoughLenMomSfcToAtm(I,J)  ,& ! out,   roughness length, momentum, surface, sent to coupled model
                  WindStressEwSfc        => noahmp%energy%state%WindStressEwSfc(I,J)      ,& ! out,   wind stress: east-west [N/m2] grid mean
                  WindStressNsSfc        => noahmp%energy%state%WindStressNsSfc(I,J)      ,& ! out,   wind stress: north-south [N/m2] grid mean
                  TemperatureRadSfc      => noahmp%energy%state%TemperatureRadSfc(I,J)    ,& ! out,   radiative temperature [K]
                  TemperatureAir2m       => noahmp%energy%state%TemperatureAir2m(I,J)     ,& ! out,   grid mean 2-m air temperature [K]
                  TemperatureAir2mBare   => noahmp%energy%state%TemperatureAir2mBare(I,J) ,& ! out,   2 m height air temperature [K] bare ground
                  EmissivitySfc          => noahmp%energy%state%EmissivitySfc(I,J)        ,& ! out,   surface emissivity
                  RoughLenMomGrd         => noahmp%energy%state%RoughLenMomGrd(I,J)       ,& ! out,   roughness length, momentum, ground [m]
                  WindStressEwBare       => noahmp%energy%state%WindStressEwBare(I,J)     ,& ! out,   wind stress: east-west [N/m2] bare ground
                  WindStressNsBare       => noahmp%energy%state%WindStressNsBare(I,J)     ,& ! out,   wind stress: north-south [N/m2] bare ground
                  SpecHumidity2mBare     => noahmp%energy%state%SpecHumidity2mBare(I,J)   ,& ! out,   bare ground 2-m water vapor mixing ratio
                  SpecHumidity2m         => noahmp%energy%state%SpecHumidity2m(I,J)       ,& ! out,   grid mean 2-m water vapor mixing ratio
                  TemperatureGrdBare     => noahmp%energy%state%TemperatureGrdBare(I,J)   ,& ! out,   bare ground temperature [K]
                  ExchCoeffMomBare       => noahmp%energy%state%ExchCoeffMomBare(I,J)     ,& ! out,   exchange coeff [m/s] for momentum, above ZeroPlaneDisp, bare ground
                  ExchCoeffShBare        => noahmp%energy%state%ExchCoeffShBare(I,J)      ,& ! out,   exchange coeff [m/s] for heat, above ZeroPlaneDisp, bare ground
                  RadLwNetSfc            => noahmp%energy%flux%RadLwNetSfc(I,J)           ,& ! out,   total net longwave rad [W/m2] (+ to atm)
                  HeatSensibleSfc        => noahmp%energy%flux%HeatSensibleSfc(I,J)       ,& ! out,   total sensible heat [W/m2] (+ to atm)
                  HeatLatentGrd          => noahmp%energy%flux%HeatLatentGrd(I,J)         ,& ! out,   total ground latent heat [W/m2] (+ to atm)
                  HeatGroundTot          => noahmp%energy%flux%HeatGroundTot(I,J)         ,& ! out,   total ground heat flux [W/m2] (+ to soil/snow)
                  HeatPrecipAdvSfc       => noahmp%energy%flux%HeatPrecipAdvSfc(I,J)      ,& ! out,   precipitation advected heat - total [W/m2]
                  RadLwEmitSfc           => noahmp%energy%flux%RadLwEmitSfc(I,J)          ,& ! out,   emitted outgoing IR [W/m2]
                  RadLwNetBareGrd        => noahmp%energy%flux%RadLwNetBareGrd(I,J)       ,& ! out,   net longwave rad [W/m2] bare ground (+ to atm)
                  HeatSensibleBareGrd    => noahmp%energy%flux%HeatSensibleBareGrd(I,J)   ,& ! out,   sensible heat flux [W/m2] bare ground (+ to atm)
                  HeatLatentBareGrd      => noahmp%energy%flux%HeatLatentBareGrd(I,J)     ,& ! out,   latent heat flux [W/m2] bare ground (+ to atm)
                  HeatGroundBareGrd      => noahmp%energy%flux%HeatGroundBareGrd(I,J)      & ! out,   bare ground heat flux [W/m2] (+ to soil/snow)
                 )
! ----------------------------------------------------------------------

        ! assign glacier bare ground quantity to grid-level quantity
        ! Energy balance at glacier (bare) ground: 
        ! RadSwAbsGrd + HeatPrecipAdvBareGrd = RadLwNetBareGrd + HeatSensibleBareGrd + HeatLatentBareGrd + HeatGroundBareGrd
        WindStressEwSfc     = WindStressEwBare
        WindStressNsSfc     = WindStressNsBare
        RadLwNetSfc         = RadLwNetBareGrd
        HeatSensibleSfc     = HeatSensibleBareGrd
        HeatLatentGrd       = HeatLatentBareGrd
        HeatGroundTot       = HeatGroundBareGrd
        TemperatureGrd      = TemperatureGrdBare
        TemperatureAir2m    = TemperatureAir2mBare
        HeatPrecipAdvSfc    = HeatPrecipAdvBareGrd
        TemperatureSfc      = TemperatureGrd
        ExchCoeffMomSfc     = ExchCoeffMomBare
        ExchCoeffShSfc      = ExchCoeffShBare
        SpecHumiditySfcMean = SpecHumiditySfc
        SpecHumidity2m      = SpecHumidity2mBare
        RoughLenMomSfcToAtm = RoughLenMomGrd

        ! emitted longwave radiation and physical check
        RadLwEmitSfc = RadLwDownRefHeight + RadLwNetSfc
#ifndef _OPENACC
        if ( RadLwEmitSfc <= 0.0 ) then
           write(*,*) "emitted longwave <0; skin T may be wrong due to inconsistent"
           write(*,*) "RadLwDownRefHeight = ", RadLwDownRefHeight, "RadLwNetSfc = ", RadLwNetSfc, "SnowDepth = ", SnowDepth
           stop "Error: Longwave radiation budget problem in NoahMP LSM"
        endif
#endif

        ! radiative temperature: subtract from the emitted IR the
        ! reflected portion of the incoming longwave radiation, so just
        ! considering the IR originating/emitted in the ground system.
        ! Old TemperatureRadSfc calculation not taking into account Emissivity:
        ! TemperatureRadSfc = (RadLwEmitSfc/ConstStefanBoltzmann)**0.25
        TemperatureRadSfc = ((RadLwEmitSfc - (1.0 - EmissivitySfc)*RadLwDownRefHeight) / &
                             (EmissivitySfc * ConstStefanBoltzmann)) ** 0.25

        end associate
      end do
    end do

    ! compute snow and glacier ice temperature
    call GlacierTemperatureMain(noahmp)

    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         associate(                                                  &
              OptSnowSoilTempTime    => noahmp%config%nmlist%OptSnowSoilTempTime      ,& ! in,    options for snow/soil temperature time scheme
              TemperatureGrd         => noahmp%energy%state%TemperatureGrd(I,J),         & ! inout, ground temperature [K]
              TemperatureSfc         => noahmp%energy%state%TemperatureSfc(I,J),         & ! inout, ground temperature [K]
              TemperatureGrdBare    => noahmp%energy%state%TemperatureGrdBare(I,J),     & ! inout, bare ground temperature [K]
              SnowDepth              => noahmp%water%state%SnowDepth(I,J)              & ! inout, snow depth [m]
             )
        ! adjusting suface temperature based on snow condition
        if ( OptSnowSoilTempTime == 2 ) then
           if ( (SnowDepth > 0.05) .and. (TemperatureGrd > ConstFreezePoint) ) then
              TemperatureGrdBare = ConstFreezePoint
              TemperatureGrd     = TemperatureGrdBare
              TemperatureSfc     = TemperatureGrdBare
           endif
        endif

        end associate
      end do
    end do

    ! Phase change and Energy released or consumed by snow & glacier ice
    call GlacierPhaseChange(noahmp)

    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

         associate(                                                  &
                  RadSwDownRefHeight     => noahmp%forcing%RadSwDownRefHeight(I,J)        ,& ! in,    downward shortwave radiation [W/m2] at reference height
                  AlbedoSfc              => noahmp%energy%state%AlbedoSfc(I,J)            ,& ! out,   total shortwave surface albedo
                  RadSwReflSfc           => noahmp%energy%flux%RadSwReflSfc(I,J)           & ! out,   total reflected solar radiation [W/m2]
             )
        ! update total surface albedo
        if ( RadSwDownRefHeight > 0.0 ) then
           AlbedoSfc = RadSwReflSfc / RadSwDownRefHeight
        else
           AlbedoSfc = undefined_real
        endif

        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine EnergyMainGlacier

end module EnergyMainGlacierMod
