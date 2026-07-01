module NoahmpDriverMainMod

  use Machine
  use NoahmpVarType
  use NoahmpIOVarType, only : NoahmpIO_type
  use ConfigVarInitMod
  use EnergyVarInitMod
  use ForcingVarInitMod
  use WaterVarInitMod
  use BiochemVarInitMod
  use ConfigVarInTransferMod
  use EnergyVarInTransferMod
  use ForcingVarInTransferMod
  use WaterVarInTransferMod
  use BiochemVarInTransferMod
  use ConfigVarOutTransferMod
  use ForcingVarOutTransferMod
  use EnergyVarOutTransferMod
  use WaterVarOutTransferMod
  use BiochemVarOutTransferMod
  use NoahmpMainMod
  use NoahmpMainGlacierMod
!   use module_ra_gfdleta,  only: cal_mon_day

  implicit none

  type(noahmp_type), allocatable :: noahmp
  logical :: noahmp_initialized = .false.

contains

  subroutine NoahmpDriverMain(NoahmpIO)
  
! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: noahmplsm
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! ------------------------------------------------------------------------- 
 
    implicit none 
    
    type(NoahmpIO_type), intent(inout)  :: NoahmpIO
    
    ! local variables
    integer                             :: I
    integer                             :: J
    integer                             :: K
    integer                             :: JMONTH, JDAY
    real(kind=kind_noahmp)              :: SOLAR_TIME 
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%nsoil ) :: SAND
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%nsoil ) :: CLAY
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%nsoil ) :: ORGM
    integer                          :: LoopInd    ! loop index for array section expansion
! ------------------------------------------------------------------------- 

    !---------------------------------------------------------------------
    !  Treatment of Noah-MP soil timestep
    !---------------------------------------------------------------------
    NoahmpIO%CALCULATE_SOIL    = .false.
    NoahmpIO%SOIL_UPDATE_STEPS = nint(NoahmpIO%SOILTSTEP / NoahmpIO%DTBL)
    NoahmpIO%SOIL_UPDATE_STEPS = max(NoahmpIO%SOIL_UPDATE_STEPS,1)

    if ( NoahmpIO%SOIL_UPDATE_STEPS == 1 ) then
       !$acc parallel loop collapse(2) gang vector default(present) private(I, J) private(LoopInd)
       do J = NoahmpIO%JTS, NoahmpIO%JTE
       do I = NoahmpIO%ITS, NoahmpIO%ITE
          NoahmpIO%ACC_SSOILXY(I,J)    = 0.0
          NoahmpIO%ACC_QINSURXY(I,J)   = 0.0
          NoahmpIO%ACC_QSEVAXY(I,J)    = 0.0
          !$acc loop seq
          do LoopInd = 1, NoahmpIO%NSOIL
             NoahmpIO%ACC_ETRANIXY(I,LoopInd,J) = 0.0
          enddo
          NoahmpIO%ACC_DWATERXY(I,J)   = 0.0
          NoahmpIO%ACC_PRCPXY(I,J)     = 0.0
          NoahmpIO%ACC_ECANXY(I,J)     = 0.0
          NoahmpIO%ACC_ETRANXY(I,J)    = 0.0
          NoahmpIO%ACC_EDIRXY(I,J)     = 0.0
          NoahmpIO%ACC_GLAFLWXY(I,J)   = 0.0
       enddo
       enddo
    endif

    if ( NoahmpIO%SOIL_UPDATE_STEPS > 1 ) then
       if ( mod(NoahmpIO%ITIMESTEP, NoahmpIO%SOIL_UPDATE_STEPS) == 1 ) then
          !$acc parallel loop collapse(2) gang vector default(present) private(I, J) firstprivate(LoopInd)
          do J = NoahmpIO%JTS, NoahmpIO%JTE
          do I = NoahmpIO%ITS, NoahmpIO%ITE
             NoahmpIO%ACC_SSOILXY(I,J)    = 0.0
             NoahmpIO%ACC_QINSURXY(I,J)   = 0.0
             NoahmpIO%ACC_QSEVAXY(I,J)    = 0.0
             !$acc loop seq
             do LoopInd = 1, NoahmpIO%NSOIL
                NoahmpIO%ACC_ETRANIXY(I,LoopInd,J) = 0.0
             enddo
             NoahmpIO%ACC_DWATERXY(I,J)   = 0.0
             NoahmpIO%ACC_PRCPXY(I,J)     = 0.0
             NoahmpIO%ACC_ECANXY(I,J)     = 0.0
             NoahmpIO%ACC_ETRANXY(I,J)    = 0.0
             NoahmpIO%ACC_EDIRXY(I,J)     = 0.0
             NoahmpIO%ACC_GLAFLWXY(I,J)   = 0.0
          enddo
          enddo
       endif
    endif

    !if ( mod(NoahmpIO%ITIMESTEP, NoahmpIO%SOIL_UPDATE_STEPS) == 0 ) NoahmpIO%CALCULATE_SOIL = .true.
    ! Prevent stale values of calculate_soil from leaking across cpu threads in if-statement above
    NoahmpIO%CALCULATE_SOIL = mod(NoahmpIO%ITIMESTEP, NoahmpIO%SOIL_UPDATE_STEPS) == 0

    !---------------------------------------------------------------------
    !  Prepare Noah-MP driver
    !---------------------------------------------------------------------
    
    ! find length of year for phenology (also S Hemisphere)
    NoahmpIO%YEARLEN = 365
    if (mod(NoahmpIO%YR,4) == 0)then
       NoahmpIO%YEARLEN = 366
       if (mod(NoahmpIO%YR,100) == 0)then
          NoahmpIO%YEARLEN = 365
          if (mod(NoahmpIO%YR,400) == 0)then
             NoahmpIO%YEARLEN = 366
          endif
       endif
    endif

    noahmp%config%domain%DayJulianInYear = NoahmpIO%JULIAN

    ! depth to soil interfaces (<0) [m]
    NoahmpIO%ZSOIL(1) = -NoahmpIO%DZS(1)
    do K = 2, NoahmpIO%NSOIL
       NoahmpIO%ZSOIL(K) = -NoahmpIO%DZS(K) + NoahmpIO%ZSOIL(K-1)
    enddo
    
    !$acc update device(NoahmpIO%YEARLEN, NoahmpIO%CALCULATE_SOIL, NoahmpIO%ZSOIL, noahmp%config%domain%DayJulianInYear)

    if ( NoahmpIO%ITIMESTEP == 1 ) then
       !$acc parallel loop collapse(2) gang vector default(present) private(I, J, K)
       JLOOP : do J = NoahmpIO%JTS, NoahmpIO%JTE
          ILOOP : do I = NoahmpIO%ITS, NoahmpIO%ITE
               if ( (NoahmpIO%XLAND(I,J)-1.5) >= 0.0 .and. &                                                                                                                                                                             
                    NoahmpIO%IVGTYP(I,J) == NoahmpIO%ISWATER_TABLE ) then  ! Open water point (real water, not SCHNAPS snow-cell flag)                                                                                                   
                if ( NoahmpIO%XICE(I,J) == 1.0 ) print*,' sea-ice at water point, I=',I,'J=',J
                NoahmpIO%SMSTAV(I,J) = 1.0
                NoahmpIO%SMSTOT(I,J) = 1.0
                do K = 1, NoahmpIO%NSOIL
                   NoahmpIO%SMOIS(I,K,J) = 1.0
                   NoahmpIO%TSLB(I,K,J)  = 273.16
                enddo
             else
                if ( NoahmpIO%XICE(I,J) == 1.0 ) then      ! Sea-ice case
                   NoahmpIO%SMSTAV(I,J) = 1.0
                   NoahmpIO%SMSTOT(I,J) = 1.0
                   do K = 1, NoahmpIO%NSOIL
                      NoahmpIO%SMOIS(I,K,J) = 1.0
                   enddo
                endif
             endif
            enddo ILOOP  ! I loop
         enddo  JLOOP    ! J loop
       endif  ! end of initialization over ocean

   !$acc parallel loop collapse(2) gang vector default(present) private(I, J, K)
   JLOOP2 : do J = NoahmpIO%JTS, NoahmpIO%JTE
      ILOOP2 : do I = NoahmpIO%ITS, NoahmpIO%ITE

          if ( NoahmpIO%XICE(I,J) >= NoahmpIO%XICE_THRESHOLD ) then  ! Sea-ice point
             NoahmpIO%ICE                        = 1
             NoahmpIO%SH2O(I,1:NoahmpIO%NSOIL,J) = 1.0
             NoahmpIO%LAI (I,J)                  = 0.01
         !     cycle ILOOP2                                             ! Skip any sea-ice points
         !  else
         !     if ( (NoahmpIO%XLAND(I,J)-1.5) >= 0.0 ) cycle ILOOP2     ! Skip any open water points
          endif
       enddo ILOOP2  ! I loop
    enddo  JLOOP2    ! J loop


             !------------------------------------------------------------------------------------
             !  Wait for async device creates from NoahmpDriverInit (first timestep only)
             !------------------------------------------------------------------------------------
            !  if (.not. noahmp_initialized) then
            !     !$acc wait(NOAHMP_ACC_QUEUE)
            !  endif

             !------------------------------------------------------------------------------------
             !  Transfer all the inputs from 2-D NoahmpIO to noahmp column variables
             !------------------------------------------------------------------------------------
             call ConfigVarInTransfer   (noahmp, NoahmpIO)
             call ForcingVarInTransfer  (noahmp, NoahmpIO)
             call EnergyVarInTransfer   (noahmp, NoahmpIO)
             call WaterVarInTransfer    (noahmp, NoahmpIO)
             call BiochemVarInTransfer  (noahmp, NoahmpIO)

    !$acc parallel loop collapse(2) gang vector default(present) private(I, J)
    do J = NoahmpIO%JTS, NoahmpIO%JTE
      do I = NoahmpIO%ITS, NoahmpIO%ITE
          if (.not.( NoahmpIO%XICE(I,J) >= NoahmpIO%XICE_THRESHOLD )) then  ! Sea-ice point
            !  if ( (NoahmpIO%XLAND(I,J)-1.5) >= 0.0 ) cycle ILOOP2     ! Skip any open water points

             ! glacier ice
             if (noahmp%config%domain%VegType(I,J) == noahmp%config%domain%IndexIcePoint ) then
                 noahmp%config%domain%IndicatorIceSfc(I,J) = -1  ! Land-ice point      
                 noahmp%forcing%TemperatureSoilBottom(I,J) = min(noahmp%forcing%TemperatureSoilBottom(I,J),263.15) ! set deep glaicer temp to >= -10C
             ! non-glacier land
             else
                 noahmp%config%domain%IndicatorIceSfc(I,J) = 0   ! land soil point.
             endif ! glacial split ends

          endif
      enddo
    enddo

             !---------------------------------------------------------------------
             !  hydrological processes for vegetation in urban model
             !  irrigate vegetaion only in urban area, MAY-SEP, 9-11pm
             ! need to be separated from Noah-MP into outside urban specific module 
             !---------------------------------------------------------------------
             !
             ! DR Jan 2026 -- SCHNAPS currently does not support Urban parameterization schemes, so commenting out
             !
            !  if ( (NoahmpIO%IVGTYP(I,J) == NoahmpIO%ISURBAN_TABLE) .or. &
            !       (NoahmpIO%IVGTYP(I,J) > NoahmpIO%URBTYPE_beg) ) then
            !     if ( (NoahmpIO%SF_URBAN_PHYSICS > 0) .and. (NoahmpIO%IRI_URBAN == 1) ) then
            !        SOLAR_TIME = (NoahmpIO%JULIAN - int(NoahmpIO%JULIAN))*24 + NoahmpIO%XLONG(I,J)/15.0
            !        if ( SOLAR_TIME < 0.0 ) SOLAR_TIME = SOLAR_TIME + 24.0
            !        call CAL_MON_DAY(int(NoahmpIO%JULIAN), NoahmpIO%YR, JMONTH, JDAY)
            !        if ( (SOLAR_TIME >= 21.0) .and. (SOLAR_TIME <= 23.0) .and. &
            !             (JMONTH >= 5) .and. (JMONTH <= 9) ) then
            !            noahmp%water%state%SoilMoisture(1) = &
            !                   max(noahmp%water%state%SoilMoisture(1),noahmp%water%param%SoilMoistureFieldCap(1))
            !            noahmp%water%state%SoilMoisture(2) = &
            !                   max(noahmp%water%state%SoilMoisture(2),noahmp%water%param%SoilMoistureFieldCap(2))
            !        endif
            !     endif
            !  endif

             !------------------------------------------------------------------------
             !  Call 1D Noah-MP LSM  
             !------------------------------------------------------------------------
         
             ! glacier ice
            !  if (noahmp%config%domain%VegType(I,J) == noahmp%config%domain%IndexIcePoint ) then
            !      noahmp%config%domain%IndicatorIceSfc(I,J) = -1  ! Land-ice point      
            !      noahmp%forcing%TemperatureSoilBottom(I,J) = min(noahmp%forcing%TemperatureSoilBottom(I,J),263.15) ! set deep glaicer temp to >= -10C
               !   call NoahmpMainGlacier(noahmp)
             ! non-glacier land
            !  else
               !   noahmp%config%domain%IndicatorIceSfc(I,J) = 0   ! land soil point.
                 call NoahmpMain(noahmp)
            !  endif ! glacial split ends

             !---------------------------------------------------------------------
             !  Transfer 1-D Noah-MP column variables to 2-D output variables
             !---------------------------------------------------------------------
             call ConfigVarOutTransfer (noahmp, NoahmpIO)
             call ForcingVarOutTransfer(noahmp, NoahmpIO)
             call EnergyVarOutTransfer (noahmp, NoahmpIO)
             call WaterVarOutTransfer  (noahmp, NoahmpIO)
             call BiochemVarOutTransfer(noahmp, NoahmpIO) 

   !        endif     ! land-sea split ends

   !     enddo ILOOP  ! I loop
   !  enddo  JLOOP    ! J loop
              
  end subroutine NoahmpDriverMain

  subroutine NoahmpDriverInit(NoahmpIO)

    implicit none

    type(NoahmpIO_type), intent(inout) :: NoahmpIO

    ! If re-initializing, clean up existing device data first
    if (noahmp_initialized) then
       call NoahmpDriverCleanup()
    endif

    allocate(noahmp)
    ! Copy the top-level noahmp structure to device (sync, needed for ConfigVarInTransfer)
    !$acc enter data copyin(noahmp)

    ! Config must be first — sets domain dimensions needed by other VarInit modules
    call ConfigVarInitDefault(noahmp)
    call ConfigVarInTransfer(noahmp, NoahmpIO)

    ! Remaining VarInit calls: allocate + batched async device create
    call ForcingVarInitDefault(noahmp)
    call EnergyVarInitDefault(noahmp)
    call WaterVarInitDefault(noahmp)
    call BiochemVarInitDefault(noahmp)
    
    noahmp_initialized = .true.

  end subroutine NoahmpDriverInit

  subroutine NoahmpDriverCleanup()
    ! Delete all device data created by VarInit modules and the noahmp structure.
    ! Called before re-initialization to prevent GPU memory leaks.

    implicit none

    call ConfigVarExitDevice(noahmp)
    call ForcingVarExitDevice(noahmp)
    call EnergyVarExitDevice(noahmp)
    call WaterVarExitDevice(noahmp)
    call BiochemVarExitDevice(noahmp)

    ! Delete the noahmp structure itself (was copyin'd in NoahmpDriverInit)
    !$acc exit data delete(noahmp)

    deallocate(noahmp)

    noahmp_initialized = .false.

  end subroutine NoahmpDriverCleanup

end module NoahmpDriverMainMod
