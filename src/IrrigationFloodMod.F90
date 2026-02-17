module IrrigationFloodMod

!!! Estimate irrigation water depth (m) based on surface flooding irrigation method
!!! Reference: chapter 4 of NRCS, Part 623 National Engineering Handbook
!!! Irrigation water is applied on the surface based on present soil moisture and
!!! infiltration rate of the soil. Flooding or overland flow is based on infiltration excess

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use IrrigationInfilPhilipMod, only : IrrigationInfilPhilip

  implicit none

contains

  subroutine IrrigationFlood(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: FLOOD_IRRIGATION
! Original code: P. Valayamkunnath (NCAR) <prasanth@ucar.edu> (08/06/2020)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J           ! grid indices
    real(kind=kind_noahmp) :: InfilRateSfc(noahmp%config%domain%ITS:noahmp%config%domain%ITE,noahmp%config%domain%JTS:noahmp%config%domain%JTE)   ! surface infiltration rate [m/s]

    !$acc data create(InfilRateSfc)
    ! estimate infiltration rate based on Philips Eq.
    call IrrigationInfilPhilip(noahmp, noahmp%config%domain%SoilTimeStep, InfilRateSfc)

    associate(                                                               &
              SoilTimeStep        => noahmp%config%domain%SoilTimeStep ,& ! in,    noahmp soil time step [s]
              NumSoilTimeStep     => noahmp%config%domain%NumSoilTimeStep,& ! in,    number of time step for calculating soil processes
              IrriFloodRateFac    => noahmp%water%param%IrriFloodRateFac,& ! in,    flood application rate factor
              IrrigationFracFlood => noahmp%water%state%IrrigationFracFlood,& ! in,    fraction of grid under flood irrigation (0 to 1)
              IrrigationAmtFlood  => noahmp%water%state%IrrigationAmtFlood,& ! inout, flood irrigation water amount [m]
              SoilSfcInflowAcc    => noahmp%water%flux%SoilSfcInflowAcc,& ! inout, accumulated water flux into soil during soil timestep [m/s * dt_soil/dt_main]
              IrrigationRateFlood => noahmp%water%flux%IrrigationRateFlood & ! inout, flood irrigation water rate [m/timestep]
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if ( .not.((noahmp%config%domain%FlagCropland(I,J) .eqv. .true.) .and. (noahmp%water%state%IrrigationAmtFlood(I,J) > 0.0)) ) cycle



    ! irrigation rate of flood irrigation. It should be
    ! greater than infiltration rate to get infiltration
    ! excess runoff at the time of application
    IrrigationRateFlood(I,J) = InfilRateSfc(I,J) * SoilTimeStep * IrriFloodRateFac(I,J)   ! Limit irrigation rate to fac*infiltration rate 
    IrrigationRateFlood(I,J) = IrrigationRateFlood(I,J) * IrrigationFracFlood(I,J)

    if ( IrrigationRateFlood(I,J) >= IrrigationAmtFlood(I,J) ) then
       IrrigationRateFlood(I,J) = IrrigationAmtFlood(I,J)
       IrrigationAmtFlood(I,J)  = 0.0
    else
       IrrigationAmtFlood(I,J)  = IrrigationAmtFlood(I,J) - IrrigationRateFlood(I,J)
    endif

    ! update water flux going to surface soil
    SoilSfcInflowAcc(I,J) = SoilSfcInflowAcc(I,J) + (IrrigationRateFlood(I,J) / SoilTimeStep * NumSoilTimeStep)  ! [m/s * dt_soil/dt_main]


      end do
    end do
   !$acc end parallel loop
   !$acc end data


    end associate

  end subroutine IrrigationFlood

end module IrrigationFloodMod
