module IrrigationMicroMod

!!! Estimate irrigation water depth (m) based on Micro irrigation method
!!! Reference: chapter 7 of NRCS, Part 623 National Engineering Handbook
!!! Irrigation water is applied under the canopy, within first layer
!!! (at ~5 cm depth) considering current soil moisture

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use IrrigationInfilPhilipMod, only : IrrigationInfilPhilip

  implicit none

contains

  subroutine IrrigationMicro(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: MICRO_IRRIGATION
! Original code: P. Valayamkunnath (NCAR) <prasanth@ucar.edu> (08/06/2020)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J              ! grid indices
    real(kind=kind_noahmp) :: InfilRateSfc(noahmp%config%domain%ITS:noahmp%config%domain%ITE,noahmp%config%domain%JTS:noahmp%config%domain%JTE)   ! surface infiltration rate [m/s]
    real(kind=kind_noahmp)           :: IrriRateTmp       ! temporary micro irrigation rate [m/timestep]

    !$acc data create(InfilRateSfc)

    ! estimate infiltration rate based on Philips Eq.
    call IrrigationInfilPhilip(noahmp, noahmp%config%domain%SoilTimeStep, InfilRateSfc)

   !$acc parallel loop collapse(2) gang vector present(noahmp, InfilRateSfc) private(IrriRateTmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

      if ( .not.((noahmp%config%domain%FlagCropland(I,J) .eqv. .true.) .and. (noahmp%water%state%IrrigationAmtFlood(I,J) > 0.0)) ) cycle
! --------------------------------------------------------------------
    associate(                                                               &
              SoilTimeStep        => noahmp%config%domain%SoilTimeStep      ,& ! in,    noahmp soil time step [s]
              DepthSoilLayer      => noahmp%config%domain%DepthSoilLayer    ,& ! in,    depth [m] of layer-bottom from soil surface
              IrrigationFracMicro => noahmp%water%state%IrrigationFracMicro(I,J),& ! in,    fraction of grid under micro irrigation (0 to 1)
              IrriMicroRate       => noahmp%water%param%IrriMicroRate(I,J)  ,& ! in,    micro irrigation rate [mm/hr]
              SoilLiqWater        => noahmp%water%state%SoilLiqWater        ,& ! inout, soil water content [m3/m3]
              IrrigationAmtMicro  => noahmp%water%state%IrrigationAmtMicro(I,J),& ! inout, micro irrigation water amount [m]
              IrrigationRateMicro => noahmp%water%flux%IrrigationRateMicro(I,J) & ! inout, micro irrigation water rate [m/timestep]
             )
! ----------------------------------------------------------------------
    
    ! initialize local variables
    IrriRateTmp = 0.0


    ! irrigation rate of micro irrigation
    IrriRateTmp         = IrriMicroRate * (1.0/1000.0) * SoilTimeStep/ 3600.0                   ! NRCS rate/time step - calibratable
    IrrigationRateMicro = min(0.5*InfilRateSfc(I,J)*SoilTimeStep, IrrigationAmtMicro, IrriRateTmp)   ! Limit irrigation rate to minimum of 0.5*infiltration rate
                                                                                                ! and to the NRCS recommended rate, (m)
    IrrigationRateMicro = IrrigationRateMicro * IrrigationFracMicro

    if ( IrrigationRateMicro >= IrrigationAmtMicro ) then
       IrrigationRateMicro = IrrigationAmtMicro
       IrrigationAmtMicro  = 0.0
    else
       IrrigationAmtMicro  = IrrigationAmtMicro - IrrigationRateMicro
    endif

    ! update soil moisture
    ! we implement drip in first layer of the Noah-MP. Change layer 1 moisture wrt to irrigation rate
    SoilLiqWater(I,1,J) = SoilLiqWater(I,1,J) + (IrrigationRateMicro / (-1.0*DepthSoilLayer(I,1,J)))

    end associate

      end do
    end do
   !$acc end parallel loop
   !$acc end data
  end subroutine IrrigationMicro

end module IrrigationMicroMod
