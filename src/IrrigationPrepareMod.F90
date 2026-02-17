module IrrigationPrepareMod

!!! Prepare dynamic irrigation variables and trigger irrigation based on conditions

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use IrrigationTriggerMod, only : IrrigationTrigger

  implicit none

contains

  subroutine IrrigationPrepare(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: None (embedded in NOAHMP_SFLX
! Original code: P. Valayamkunnath (NCAR) <prasanth@ucar.edu> (08/06/2020)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J     ! grid indices

    associate(                                                                           &
              FlagCropland            => noahmp%config%domain%FlagCropland         ,& ! in,    flag to identify croplands
              FlagSoilProcess         => noahmp%config%domain%FlagSoilProcess           ,& ! in,    flag to calculate soil processes
              OptIrrigationMethod     => noahmp%config%nmlist%OptIrrigationMethod       ,& ! in,    irrigation method option
              IrriFracThreshold       => noahmp%water%param%IrriFracThreshold      ,& ! in,    irrigation fraction threshold
              IrriStopPrecipThr       => noahmp%water%param%IrriStopPrecipThr      ,& ! in,    maximum precipitation to stop irrigation trigger
              IrrigationFracGrid      => noahmp%water%state%IrrigationFracGrid     ,& ! in,    total input irrigation fraction of a grid
              IrrigationAmtSprinkler  => noahmp%water%state%IrrigationAmtSprinkler ,& ! inout, irrigation water amount [m] to be applied, Sprinkler
              IrrigationAmtFlood      => noahmp%water%state%IrrigationAmtFlood     ,& ! inout, flood irrigation water amount [m]
              IrrigationAmtMicro      => noahmp%water%state%IrrigationAmtMicro     ,& ! inout, micro irrigation water amount [m]
              RainfallRefHeight       => noahmp%water%flux%RainfallRefHeight       ,& ! inout, rainfall [mm/s] at reference height
              IrrigationFracSprinkler => noahmp%water%state%IrrigationFracSprinkler,& ! out,   sprinkler irrigation fraction (0 to 1)
              IrrigationFracMicro     => noahmp%water%state%IrrigationFracMicro    ,& ! out,   fraction of grid under micro irrigation (0 to 1)
              IrrigationFracFlood     => noahmp%water%state%IrrigationFracFlood     & ! out,   fraction of grid under flood irrigation (0 to 1)
             )

   !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! if OptIrrigationMethod = 0 and if methods are unknown for certain area, then use sprinkler irrigation method
    if ( (OptIrrigationMethod == 0) .and. (IrrigationFracSprinkler(I,J) == 0.0) .and. (IrrigationFracMicro(I,J) == 0.0) &
         .and. (IrrigationFracFlood(I,J) == 0.0) .and. (IrrigationFracGrid(I,J) >= IrriFracThreshold(I,J)) ) then
       IrrigationFracSprinkler(I,J) = 1.0
    endif

    ! choose method based on user namelist choice
    if ( OptIrrigationMethod == 1 ) then
       IrrigationFracSprinkler(I,J) = 1.0
       IrrigationFracMicro(I,J)     = 0.0
       IrrigationFracFlood(I,J)     = 0.0
    elseif ( OptIrrigationMethod == 2 ) then
       IrrigationFracSprinkler(I,J) = 0.0
       IrrigationFracMicro(I,J)     = 1.0
       IrrigationFracFlood(I,J)     = 0.0
    elseif ( OptIrrigationMethod == 3 ) then
       IrrigationFracSprinkler(I,J) = 0.0
       IrrigationFracMicro(I,J)     = 0.0
       IrrigationFracFlood(I,J)     = 1.0
    endif

    ! trigger irrigation only at soil water timestep to be consistent for solving soil water
    if ( FlagSoilProcess .eqv. .true. ) then
       if ( (FlagCropland(I,J) .eqv. .true.) .and. (IrrigationFracGrid(I,J) >= IrriFracThreshold(I,J)) .and. &
            (RainfallRefHeight(I,J) < (IrriStopPrecipThr(I,J)/3600.0)) .and. &
            ((IrrigationAmtSprinkler(I,J)+IrrigationAmtMicro(I,J)+IrrigationAmtFlood(I,J)) == 0.0) ) then
          call IrrigationTrigger(noahmp, I, J)
       endif

       ! set irrigation off if larger than IrriStopPrecipThr mm/h for this time step and irr triggered last time step
       if ( (RainfallRefHeight(I,J) >= (IrriStopPrecipThr(I,J)/3600.0)) .or. (IrrigationFracGrid(I,J) < IrriFracThreshold(I,J)) ) then
          IrrigationAmtSprinkler(I,J) = 0.0
          IrrigationAmtMicro(I,J)     = 0.0
          IrrigationAmtFlood(I,J)     = 0.0
       endif
    endif


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine IrrigationPrepare

end module IrrigationPrepareMod
