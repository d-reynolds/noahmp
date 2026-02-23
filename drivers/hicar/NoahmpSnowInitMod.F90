module NoahmpSnowInitMod

!  Module to initialize Noah-MP Snow variables

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type
  
  implicit none
  
contains

  subroutine NoahmpSnowInitMain(NoahmpIO)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SNOW_INIT
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! ------------------------------------------------------------------------- 

    implicit none 
    
    type(NoahmpIO_type), intent(inout) :: NoahmpIO
    
! local variables
    integer                                                               :: I,J,IZ,itf,jtf
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: DZSNO
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: DZSNSO

!------------------------------------------------------------------------------------------    
!   Initialize snow arrays for Noah-MP LSM, based in input SNOWDEP, NSNOW
!   ISNOWXY is an index array, indicating the index of the top snow layer.  Valid indices
!           for snow layers range from 0 (no snow) and -1 (shallow snow) to (-NSNOW)+1 (deep snow).
!   TSNOXY holds the temperature of the snow layer.  Snow layers are initialized with 
!          temperature = ground temperature [?].  Snow-free levels in the array have value 0.0
!   SNICEXY is the frozen content of a snow layer.  Initial estimate based on SNOWH and SNOW
!   SNLIQXY is the liquid content of a snow layer.  Initialized to 0.0
!   ZNSNOXY is the layer depth from the surface.  
!   SNRDSXY is the snow layer effective grain radius [microns, m-6]
!   SNFRXY  is the snow layer rate of snow freezing [mm/s]
!------------------------------------------------------------------------------------------

    itf = min0(NoahmpIO%ite, NoahmpIO%ide-1)
    jtf = min0(NoahmpIO%jte, NoahmpIO%jde-1)

    allocate(DZSNO(NoahmpIO%its:itf, -NoahmpIO%NSNOW+1:0, NoahmpIO%jts:jtf))
    allocate(DZSNSO(NoahmpIO%its:itf, -NoahmpIO%NSNOW+1:NoahmpIO%NSOIL, NoahmpIO%jts:jtf))
    !$acc data create(DZSNO, DZSNSO)

    !$acc parallel loop collapse(2) gang vector default(present) private(I,J,IZ)
    do J = NoahmpIO%jts, jtf
       do I = NoahmpIO%its, itf

          ! initialize snow layers and thickness
          ! no explicit snow layer
          if ( NoahmpIO%SNOWH(I,J) < 0.025 ) then
             NoahmpIO%ISNOWXY(I,J) = 0
             DZSNO(I,-NoahmpIO%NSNOW+1:0,J) = 0.0
          else
             ! 1 layer snow
             if ( (NoahmpIO%SNOWH(I,J) >= 0.025) .and. (NoahmpIO%SNOWH(I,J) <= 0.05) ) then
                NoahmpIO%ISNOWXY(I,J) = -1
                DZSNO(I,0,J)  = NoahmpIO%SNOWH(I,J)
             ! 2 layer snow
             elseif ( (NoahmpIO%SNOWH(I,J) > 0.05) .and. (NoahmpIO%SNOWH(I,J) <= 0.10) ) then
                NoahmpIO%ISNOWXY(I,J) = -2
                DZSNO(I,-1,J) = NoahmpIO%SNOWH(I,J) / 2.0
                DZSNO(I, 0,J) = NoahmpIO%SNOWH(I,J) / 2.0
             ! 2 layer thick snow
             elseif ( (NoahmpIO%SNOWH(I,J) > 0.10) .and. (NoahmpIO%SNOWH(I,J) <= 0.25) ) then
                NoahmpIO%ISNOWXY(I,J) = -2
                DZSNO(I,-1,J) = 0.05
                DZSNO(I, 0,J) = NoahmpIO%SNOWH(I,J) - DZSNO(I,-1,J)
             ! 3 layer snow
             elseif ( (NoahmpIO%SNOWH(I,J) > 0.25) .and. (NoahmpIO%SNOWH(I,J) <= 0.45) ) then
                NoahmpIO%ISNOWXY(I,J) = -3
                DZSNO(I,-2,J) = 0.05
                DZSNO(I,-1,J) = 0.5 * (NoahmpIO%SNOWH(I,J)-DZSNO(I,-2,J))
                DZSNO(I, 0,J) = 0.5 * (NoahmpIO%SNOWH(I,J)-DZSNO(I,-2,J))
             ! 3 layer thick snow
             elseif ( NoahmpIO%SNOWH(I,J) > 0.45 ) then
                NoahmpIO%ISNOWXY(I,J) = -3
                DZSNO(I,-2,J) = 0.05
                DZSNO(I,-1,J) = 0.20
                DZSNO(I, 0,J) = NoahmpIO%SNOWH(I,J) - DZSNO(I,-1,J) - DZSNO(I,-2,J)
             else
                print*, "Problem with the logic assigning snow layers."
                stop
             endif
          endif

          ! initialize snow temperatuer and ice/liquid content
          !$acc loop seq
          do IZ = -NoahmpIO%NSNOW+1, 0
            NoahmpIO%TSNOXY (I,IZ,J) = 0.0
            NoahmpIO%SNICEXY(I,IZ,J) = 0.0
            NoahmpIO%SNLIQXY(I,IZ,J) = 0.0
          enddo
          !$acc loop seq
          do IZ = NoahmpIO%ISNOWXY(I,J)+1, 0
             NoahmpIO%TSNOXY(I,IZ,J)  = NoahmpIO%TGXY(I,J)
             NoahmpIO%SNLIQXY(I,IZ,J) = 0.0
             NoahmpIO%SNICEXY(I,IZ,J) = 1.0 * DZSNO(I,IZ,J) * (NoahmpIO%SNOW(I,J)/NoahmpIO%SNOWH(I,J))
          enddo

          ! Assign local variable DZSNSO, the soil/snow layer thicknesses, for snow layers
          !$acc loop seq
          do IZ = NoahmpIO%ISNOWXY(I,J)+1, 0
             DZSNSO(I,IZ,J) = -DZSNO(I,IZ,J)
          enddo

          ! Assign local variable DZSNSO, the soil/snow layer thicknesses, for soil layers
          DZSNSO(I,1,J) = NoahmpIO%ZSOIL(1)
          !$acc loop seq
          do IZ = 2, NoahmpIO%NSOIL
             DZSNSO(I,IZ,J) = NoahmpIO%ZSOIL(IZ) - NoahmpIO%ZSOIL(IZ-1)
          enddo

          ! Assign ZSNSOXY, the layer depths, for soil and snow layers
          NoahmpIO%ZSNSOXY(I,NoahmpIO%ISNOWXY(I,J)+1,J) = DZSNSO(I,NoahmpIO%ISNOWXY(I,J)+1,J)
          !$acc loop seq
          do IZ = NoahmpIO%ISNOWXY(I,J)+2, NoahmpIO%NSOIL
             NoahmpIO%ZSNSOXY(I,IZ,J) = NoahmpIO%ZSNSOXY(I,IZ-1,J) + DZSNSO(I,IZ,J)
          enddo

          ! SNICAR
          if ( NoahmpIO%IOPT_ALB == 3 )then
            !$acc loop seq
            do IZ = -NoahmpIO%NSNOW+1, 0
               NoahmpIO%SNRDSXY(I,IZ,J)  = 0.0
               NoahmpIO%SNFRXY (I,IZ,J)  = 0.0
               NoahmpIO%BCPHIXY(I,IZ,J)  = 0.0
               NoahmpIO%BCPHOXY(I,IZ,J)  = 0.0
               NoahmpIO%OCPHIXY(I,IZ,J)  = 0.0
               NoahmpIO%OCPHOXY(I,IZ,J)  = 0.0
               NoahmpIO%DUST1XY(I,IZ,J)  = 0.0
               NoahmpIO%DUST2XY(I,IZ,J)  = 0.0
               NoahmpIO%DUST3XY(I,IZ,J)  = 0.0
               NoahmpIO%DUST4XY(I,IZ,J)  = 0.0
               NoahmpIO%DUST5XY(I,IZ,J)  = 0.0
            enddo
          endif

       enddo ! I
    enddo    ! J
    !$acc end data
    deallocate(DZSNO)
    deallocate(DZSNSO)

  end subroutine NoahmpSnowInitMain

end module NoahmpSnowInitMod

