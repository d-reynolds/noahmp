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

    associate(                                                               &
              SoilTimeStep        => noahmp%config%domain%SoilTimeStep      ,& ! in,    noahmp soil time step [s]
              DepthSoilLayer      => noahmp%config%domain%DepthSoilLayer    ,& ! in,    depth [m] of layer-bottom from soil surface
              IrrigationFracMicro => noahmp%water%state%IrrigationFracMicro,& ! in,    fraction of grid under micro irrigation (0 to 1)
              IrriMicroRate       => noahmp%water%param%IrriMicroRate  ,& ! in,    micro irrigation rate [mm/hr]
              SoilLiqWater        => noahmp%water%state%SoilLiqWater        ,& ! inout, soil water content [m3/m3]
              IrrigationAmtMicro  => noahmp%water%state%IrrigationAmtMicro,& ! inout, micro irrigation water amount [m]
              IrrigationRateMicro => noahmp%water%flux%IrrigationRateMicro & ! inout, micro irrigation water rate [m/timestep]
             )

   !$acc parallel loop collapse(2) gang vector default(present) private(IrriRateTmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

      if ( .not.((noahmp%config%domain%FlagCropland(I,J) .eqv. .true.) .and. (noahmp%water%state%IrrigationAmtFlood(I,J) > 0.0)) ) cycle
    
    ! initialize local variables
    IrriRateTmp = 0.0


    ! irrigation rate of micro irrigation
    IrriRateTmp         = IrriMicroRate(I,J) * (1.0/1000.0) * SoilTimeStep/ 3600.0                   ! NRCS rate/time step - calibratable
    IrrigationRateMicro(I,J) = min(0.5*InfilRateSfc(I,J)*SoilTimeStep, IrrigationAmtMicro(I,J), IrriRateTmp)   ! Limit irrigation rate to minimum of 0.5*infiltration rate
                                                                                                ! and to the NRCS recommended rate, (m)
    IrrigationRateMicro(I,J) = IrrigationRateMicro(I,J) * IrrigationFracMicro(I,J)

    if ( IrrigationRateMicro(I,J) >= IrrigationAmtMicro(I,J) ) then
       IrrigationRateMicro(I,J) = IrrigationAmtMicro(I,J)
       IrrigationAmtMicro(I,J)  = 0.0
    else
       IrrigationAmtMicro(I,J)  = IrrigationAmtMicro(I,J) - IrrigationRateMicro(I,J)
    endif

    ! update soil moisture
    ! we implement drip in first layer of the Noah-MP. Change layer 1 moisture wrt to irrigation rate
    SoilLiqWater(I,1,J) = SoilLiqWater(I,1,J) + (IrrigationRateMicro(I,J) / (-1.0*DepthSoilLayer(I,1,J)))


      end do
    end do
   !$acc end parallel loop
   !$acc end data

    end associate

  end subroutine IrrigationMicro

end module IrrigationMicroMod
