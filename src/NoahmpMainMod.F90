module NoahmpMainMod

!!! Main NoahMP module including all column model processes
!!! atmos forcing -> canopy intercept -> precip heat advect -> main energy -> main water -> main biogeochemistry -> balance check

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use AtmosForcingMod,            only : ProcessAtmosForcing
  use GeneralInitMod,             only : GeneralInit
  use PhenologyMainMod,           only : PhenologyMain
  use IrrigationPrepareMod,       only : IrrigationPrepare
  use IrrigationSprinklerMod,     only : IrrigationSprinkler
  use CanopyWaterInterceptMod,    only : CanopyWaterIntercept
  use PrecipitationHeatAdvectMod, only : PrecipitationHeatAdvect
  use EnergyMainMod,              only : EnergyMain
  use WaterMainMod,               only : WaterMain
  use BiochemNatureVegMainMod,    only : BiochemNatureVegMain
  use BiochemCropMainMod,         only : BiochemCropMain
  use BalanceErrorCheckMod,       only : BalanceWaterInit, BalanceWaterCheck, BalanceEnergyCheck 
 
  implicit none

contains

  subroutine NoahmpMain(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: NOAHMP_SFLX
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    ! NoahmpMain is the top-level driver that calls other physics subroutines
    ! Each called subroutine has its own OpenACC parallel regions
    ! No I,J loops needed here as parallelism is handled in the called routines

    !---------------------------------------------------------------------
    ! Atmospheric forcing processing
    !--------------------------------------------------------------------- 

    call ProcessAtmosForcing(noahmp)

    !---------------------------------------------------------------------
    ! General initialization to prepare key variables
    !--------------------------------------------------------------------- 

    call GeneralInit(noahmp)

    !---------------------------------------------------------------------
    ! Prepare for water balance check
    !--------------------------------------------------------------------- 

    call BalanceWaterInit(noahmp)

    !---------------------------------------------------------------------
    ! Phenology
    !--------------------------------------------------------------------- 

    call PhenologyMain(noahmp)

    !---------------------------------------------------------------------
    ! Irrigation prepare including trigger
    !--------------------------------------------------------------------- 

    call IrrigationPrepare(noahmp)

    !---------------------------------------------------------------------
    ! Sprinkler irrigation
    !--------------------------------------------------------------------- 

    call IrrigationSprinkler(noahmp)

    !---------------------------------------------------------------------
    ! Canopy water interception and precip heat advection
    !--------------------------------------------------------------------- 

    call CanopyWaterIntercept(noahmp)
    call PrecipitationHeatAdvect(noahmp)

    !---------------------------------------------------------------------
    ! Energy processes
    !--------------------------------------------------------------------- 

    call EnergyMain(noahmp)

    !---------------------------------------------------------------------
    ! Water processes
    !--------------------------------------------------------------------- 

    call WaterMain(noahmp)

    !---------------------------------------------------------------------
    ! Biochem processes (crop and carbon)
    !--------------------------------------------------------------------- 

    ! for generic vegetation - handled within BiochemNatureVegMain
    call BiochemNatureVegMain(noahmp)

    ! for explicit crop treatment - handled within BiochemCropMain
    call BiochemCropMain(noahmp)

    !---------------------------------------------------------------------
    ! Error check for energy and water balance
    !--------------------------------------------------------------------- 

    call BalanceWaterCheck(noahmp)
    call BalanceEnergyCheck(noahmp) 

    !---------------------------------------------------------------------
    ! End of all NoahMP column processes
    !--------------------------------------------------------------------- 

  end subroutine NoahmpMain

end module NoahmpMainMod
