module SoilWaterMainMod

!!! Main soil water module including all soil water processes & update soil moisture
!!! surface runoff, infiltration, soil water diffusion, subsurface runoff, tile drainage

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use RunoffSurfaceTopModelGrdMod,       only : RunoffSurfaceTopModelGrd
  use RunoffSurfaceTopModelEquiMod,      only : RunoffSurfaceTopModelEqui
  use RunoffSurfaceFreeDrainMod,         only : RunoffSurfaceFreeDrain
  use RunoffSurfaceBatsMod,              only : RunoffSurfaceBATS
  use RunoffSurfaceTopModelMmfMod,       only : RunoffSurfaceTopModelMMF
  use RunoffSurfaceVicMod,               only : RunoffSurfaceVIC
  use RunoffSurfaceXinAnJiangMod,        only : RunoffSurfaceXinAnJiang
  use RunoffSurfaceDynamicVicMod,        only : RunoffSurfaceDynamicVic
  use RunoffSurfaceWetlandMod,           only : RunoffSurfaceWetland
  use RunoffSubSurfaceEquiWaterTableMod, only : RunoffSubSurfaceEquiWaterTable
  use RunoffSubSurfaceGroundWaterMod,    only : RunoffSubSurfaceGroundWater
  use RunoffSubSurfaceDrainageMod,       only : RunoffSubSurfaceDrainage
  use RunoffSubSurfaceShallowMmfMod,     only : RunoffSubSurfaceShallowWaterMMF
  use SoilWaterDiffusionRichardsMod,     only : SoilWaterDiffusionRichards
  use SoilMoistureSolverMod,             only : SoilMoistureSolver
  use TileDrainageSimpleMod,             only : TileDrainageSimple
  use TileDrainageHooghoudtMod,          only : TileDrainageHooghoudt

  implicit none

contains

  subroutine SoilWaterMain(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SOILWATER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout)  :: noahmp

! local variables
    integer                           :: I, J                         ! grid indices
    integer                           :: LoopInd1, LoopInd2           ! loop index
    integer                           :: IndIter                      ! iteration index
    integer                           :: NumIterSoilWat               ! iteration times soil moisture
    real(kind=kind_noahmp)            :: TimeStepFine                 ! fine time step [s]
    real(kind=kind_noahmp)            :: SoilSatExcAcc                ! accumulation of soil saturation excess [m]
    real(kind=kind_noahmp)            :: SoilWatConductAcc            ! sum of SoilWatConductivity*ThicknessSnowSoilLayer
    real(kind=kind_noahmp)            :: WaterRemove                  ! water mass removal [mm]
    real(kind=kind_noahmp)            :: SoilWatRem                   ! temporary remaining soil water [mm]
    real(kind=kind_noahmp)            :: SoilWaterMin                 ! minimum soil water [mm]
    real(kind=kind_noahmp)            :: DrainSoilBotAcc              ! accumulated drainage water [mm] at fine time step
    real(kind=kind_noahmp)            :: RunoffSurfaceAcc             ! accumulated surface runoff [mm] at fine time step
    real(kind=kind_noahmp)            :: InfilSfcAcc                  ! accumulated infiltration rate [m/s]
    real(kind=kind_noahmp), parameter :: SoilImpPara = 4.0            ! soil impervious fraction parameter
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatRight ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatLeft1 ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatLeft2 ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatLeft3 ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: SoilLiqTmp   ! temporary soil liquid water [mm]
    ! 2D accumulator arrays that persist across parallel regions
    real(kind=kind_noahmp), allocatable, dimension(:,:)   :: SoilSatExcAcc2D
    real(kind=kind_noahmp), allocatable, dimension(:,:)   :: DrainSoilBotAcc2D
    real(kind=kind_noahmp), allocatable, dimension(:,:)   :: RunoffSurfaceAcc2D

    associate(                                                                      &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer ,& ! in,    number of soil layers
              SoilTimeStep           => noahmp%config%domain%SoilTimeStep ,& ! in,    noahmp soil time step [s]
              OptRunoffSurface       => noahmp%config%nmlist%OptRunoffSurface ,& ! in,    options for surface runoff
              OptRunoffSubsurface    => noahmp%config%nmlist%OptRunoffSubsurface ,& ! in,    options for subsurface runoff
              OptTileDrainage        => noahmp%config%nmlist%OptTileDrainage ,& ! in,    options for tile drainage
              OptWetlandModel        => noahmp%config%nmlist%OptWetlandModel ,& ! in,    options for wetland model
                  ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
                  SoilMoistureSat        => noahmp%water%param%SoilMoistureSat ,& ! in,    saturated value of soil moisture [m3/m3]
                  SoilIce                => noahmp%water%state%SoilIce ,& ! in,    soil ice content [m3/m3]
                  SoilLiqWater           => noahmp%water%state%SoilLiqWater ,& ! inout, soil water content [m3/m3]
                  RunoffSurface          => noahmp%water%flux%RunoffSurface ,& ! out,   surface runoff [mm per soil timestep]
                  RunoffSubsurface       => noahmp%water%flux%RunoffSubsurface ,& ! out,   subsurface runoff [mm per soil timestep]
                  InfilRateSfc           => noahmp%water%flux%InfilRateSfc ,& ! out,   infiltration rate at surface [m/s]
                  TileDrain              => noahmp%water%flux%TileDrain ,& ! out,   tile drainage [mm per soil timestep]
                  SoilEffPorosity        => noahmp%water%state%SoilEffPorosity ,& ! out,   soil effective porosity [m3/m3]
                  SoilIceFrac            => noahmp%water%state%SoilIceFrac ,& ! out,   ice fraction in frozen soil
                  SoilImpervFrac         => noahmp%water%state%SoilImpervFrac ,& ! out,   impervious fraction due to frozen soil
                  SoilIceMax             => noahmp%water%state%SoilIceMax ,& ! out,   maximum soil ice content [m3/m3]
                  SoilImpervFracMax      => noahmp%water%state%SoilImpervFracMax ,& ! out,   maximum soil imperviousness fraction
                  SoilLiqWaterMin        => noahmp%water%state%SoilLiqWaterMin ,& ! out,   minimum soil liquid water content [m3/m3]
                  SoilMoisture           => noahmp%water%state%SoilMoisture ,& ! inout, total soil moisture [m3/m3]
                  RechargeGwDeepWT       => noahmp%water%state%RechargeGwDeepWT ,& ! inout, recharge to or from water table when deep [m]
                  SoilWatConductivity    => noahmp%water%state%SoilWatConductivity  & ! in,    soil hydraulic conductivity [m/s]
             )


    ! allocate 3D matrix arrays (I, NumSoilLayer, J) for passing to subroutines
    if (.not. allocated(MatRight)) &
       allocate(MatRight(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                         1:NumSoilLayer, &
                         noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    if (.not. allocated(MatLeft1)) &
       allocate(MatLeft1(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                         1:NumSoilLayer, &
                         noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    if (.not. allocated(MatLeft2)) &
       allocate(MatLeft2(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                         1:NumSoilLayer, &
                         noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    if (.not. allocated(MatLeft3)) &
       allocate(MatLeft3(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                         1:NumSoilLayer, &
                         noahmp%config%domain%JTS:noahmp%config%domain%JTE))

    ! allocate 2D accumulator arrays
    if (.not. allocated(SoilSatExcAcc2D)) &
       allocate(SoilSatExcAcc2D(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                                noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    if (.not. allocated(DrainSoilBotAcc2D)) &
       allocate(DrainSoilBotAcc2D(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                                  noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    if (.not. allocated(RunoffSurfaceAcc2D)) &
       allocate(RunoffSurfaceAcc2D(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                                   noahmp%config%domain%JTS:noahmp%config%domain%JTE))

    !$acc data create(MatRight, MatLeft1, MatLeft2, MatLeft3, &
    !$acc             SoilSatExcAcc2D, DrainSoilBotAcc2D, RunoffSurfaceAcc2D)

    ! ===== Region 1: Initialization and soil property computation =====
    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd1, SoilSatExcAcc) private(InfilSfcAcc)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points


        ! initialization
        !$acc loop seq
        do LoopInd1 = 1, NumSoilLayer
           MatRight(I,LoopInd1,J) = 0.0
           MatLeft1(I,LoopInd1,J) = 0.0
           MatLeft2(I,LoopInd1,J) = 0.0
           MatLeft3(I,LoopInd1,J) = 0.0
        enddo
         RunoffSurface(I,J)    = 0.0
         RunoffSubsurface(I,J) = 0.0
         InfilRateSfc(I,J)     = 0.0
         SoilSatExcAcc    = 0.0
         InfilSfcAcc      = 1.0e-06

        ! for the case when snowmelt water is too large
        !$acc loop seq
        do LoopInd1 = 1, NumSoilLayer
           SoilEffPorosity(I,LoopInd1,J) = max(1.0e-4, (SoilMoistureSat(I,LoopInd1,J) - SoilIce(I,LoopInd1,J)))
           SoilSatExcAcc = SoilSatExcAcc + max(0.0, SoilLiqWater(I,LoopInd1,J) - &
                                SoilEffPorosity(I,LoopInd1,J)) * ThicknessSnowSoilLayer(I,LoopInd1,J)
           SoilLiqWater(I,LoopInd1,J) = min(SoilEffPorosity(I,LoopInd1,J), SoilLiqWater(I,LoopInd1,J))
        enddo
        SoilSatExcAcc2D(I,J) = SoilSatExcAcc

        ! impervious fraction due to frozen soil
        !$acc loop seq
        do LoopInd1 = 1, NumSoilLayer
           SoilIceFrac(I,LoopInd1,J)    = min(1.0, SoilIce(I,LoopInd1,J) / SoilMoistureSat(I,LoopInd1,J))
           SoilImpervFrac(I,LoopInd1,J) = max(0.0, exp(-SoilImpPara*(1.0-SoilIceFrac(I,LoopInd1,J))) - &
                                           exp(-SoilImpPara)) / (1.0 - exp(-SoilImpPara))
        enddo

        ! maximum soil ice content and minimum liquid water of all layers
        SoilIceMax(I,J)        = 0.0
        SoilImpervFracMax(I,J) = 0.0
        SoilLiqWaterMin(I,J)   = SoilMoistureSat(I,1,J)
        !$acc loop seq
        do LoopInd1 = 1, NumSoilLayer
           if ( SoilIce(I,LoopInd1,J) > SoilIceMax(I,J) )               SoilIceMax(I,J)        = SoilIce(I,LoopInd1,J)
           if ( SoilImpervFrac(I,LoopInd1,J) > SoilImpervFracMax(I,J) ) SoilImpervFracMax(I,J) = SoilImpervFrac(I,LoopInd1,J)
           if ( SoilLiqWater(I,LoopInd1,J) < SoilLiqWaterMin(I,J) )     SoilLiqWaterMin(I,J)   = SoilLiqWater(I,LoopInd1,J)
        enddo


      end do
    end do
    !$acc end parallel loop

    ! subsurface runoff for runoff scheme option 2
    if ( OptRunoffSubsurface == 2 ) call RunoffSubSurfaceEquiWaterTable(noahmp)

    ! jref impermable surface at urban
      !$acc loop gang vector collapse(2)
      do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
        if ( noahmp%config%domain%FlagUrban(I,J) .eqv. .true. ) then
          noahmp%water%state%SoilImpervFrac(I,1,J) = 0.95
        endif
      enddo
      enddo
    ! surface runoff and infiltration rate using different schemes
    if ( OptRunoffSurface == 1 ) call RunoffSurfaceTopModelGrd(noahmp)
    if ( OptRunoffSurface == 2 ) call RunoffSurfaceTopModelEqui(noahmp)
    if ( OptRunoffSurface == 3 ) call RunoffSurfaceFreeDrain(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 4 ) call RunoffSurfaceBATS(noahmp)
    if ( OptRunoffSurface == 5 ) call RunoffSurfaceTopModelMMF(noahmp)
    if ( OptRunoffSurface == 6 ) call RunoffSurfaceVIC(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 7 ) call RunoffSurfaceXinAnJiang(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 8 ) call RunoffSurfaceDynamicVic(noahmp,SoilTimeStep,InfilSfcAcc)

    ! special treatment for wetland points (due to subgrid wetland treatment, currently no flag control)
    !if ( (FlagWetland .eqv. .true.) .and. (OptWetlandModel > 0) ) call RunoffSurfaceWetland(noahmp)
    if ( OptWetlandModel > 0 ) call RunoffSurfaceWetland(noahmp)

    ! determine iteration times to solve soil water diffusion and moisture
    ! Use maximum iteration count for GPU (all grid points use same count)
    NumIterSoilWat = 6
    TimeStepFine = SoilTimeStep / NumIterSoilWat

    ! solve soil moisture
    InfilSfcAcc      = 1.0e-06
    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
        DrainSoilBotAcc2D(I,J)  = 0.0
        RunoffSurfaceAcc2D(I,J) = 0.0
      end do
    end do
    !$acc end parallel loop

    do IndIter = 1, NumIterSoilWat

       ! surface runoff update within iteration (subroutines have own parallel regions)
       if ( OptRunoffSurface == 3 ) call RunoffSurfaceFreeDrain(noahmp,TimeStepFine)
       if ( OptRunoffSurface == 6 ) call RunoffSurfaceVIC(noahmp,TimeStepFine)
       if ( OptRunoffSurface == 7 ) call RunoffSurfaceXinAnJiang(noahmp,TimeStepFine)
       if ( OptRunoffSurface == 8 ) call RunoffSurfaceDynamicVic(noahmp,TimeStepFine,InfilSfcAcc)

       call SoilWaterDiffusionRichards(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)
       call SoilMoistureSolver(noahmp, TimeStepFine, MatLeft1, MatLeft2, MatLeft3, MatRight)

       !$acc parallel loop collapse(2) gang vector default(present)
       do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
           SoilSatExcAcc2D(I,J)    = SoilSatExcAcc2D(I,J)    + noahmp%water%state%SoilSaturationExcess(I,J)
           DrainSoilBotAcc2D(I,J)  = DrainSoilBotAcc2D(I,J)  + noahmp%water%flux%DrainSoilBot(I,J)
           RunoffSurfaceAcc2D(I,J) = RunoffSurfaceAcc2D(I,J) + noahmp%water%flux%RunoffSurface(I,J)
         end do
       end do
       !$acc end parallel loop

    enddo  ! IndIter

    !$acc parallel loop collapse(2) gang vector default(present) firstprivate(NumIterSoilWat)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points
        noahmp%water%flux%DrainSoilBot(I,J) = DrainSoilBotAcc2D(I,J) / NumIterSoilWat
        noahmp%water%flux%RunoffSurface(I,J) = RunoffSurfaceAcc2D(I,J) / NumIterSoilWat
        noahmp%water%flux%RunoffSurface(I,J) = noahmp%water%flux%RunoffSurface(I,J) * 1000.0 + &
                                                SoilSatExcAcc2D(I,J) * 1000.0 / SoilTimeStep  ! m/s -> mm/s
        noahmp%water%flux%DrainSoilBot(I,J)  = noahmp%water%flux%DrainSoilBot(I,J) * 1000.0   ! m/s -> mm/s
      end do
    end do
    !$acc end parallel loop

    if ( (OptTileDrainage == 1) .and. (OptRunoffSurface == 3) ) then
       call TileDrainageSimple(noahmp)  ! simple tile drainage
    endif
    if ( (OptTileDrainage == 2) .and. (OptRunoffSurface == 3) ) then
       call TileDrainageHooghoudt(noahmp)  ! Hooghoudt tile drain
    endif


    allocate(SoilLiqTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                        1:NumSoilLayer, &
                        noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    !$acc data create(SoilLiqTmp)

    !$acc parallel loop collapse(2) gang vector default(present) &
    !$acc private(LoopInd1, LoopInd2, SoilWatConductAcc, WaterRemove, SoilWatRem, SoilWaterMin)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points


        !$acc loop seq
        do LoopInd1 = 1, NumSoilLayer
           SoilLiqTmp(I,LoopInd1,J)    = 0.0
        enddo
        ! removal of soil water due to subsurface runoff (option 2)
        if ( OptRunoffSubsurface == 2 ) then
           SoilWatConductAcc = 0.0
           !$acc loop seq
           do LoopInd1 = 1, NumSoilLayer
              SoilWatConductAcc = SoilWatConductAcc + SoilWatConductivity(I,LoopInd1,J) * &
                                  ThicknessSnowSoilLayer(I,LoopInd1,J)
           enddo
           !$acc loop seq
           do LoopInd1 = 1, NumSoilLayer
              WaterRemove = RunoffSubsurface(I,J) * SoilTimeStep * &
                            (SoilWatConductivity(I,LoopInd1,J)*ThicknessSnowSoilLayer(I,LoopInd1,J)) / SoilWatConductAcc
              SoilLiqWater(I,LoopInd1,J) = SoilLiqWater(I,LoopInd1,J) - &
                                            WaterRemove / (ThicknessSnowSoilLayer(I,LoopInd1,J)*1000.0)
           enddo
        endif

        ! Limit SoilLiqTmp to be greater than or equal to watmin.
        ! Get water needed to bring SoilLiqTmp equal SoilWaterMin from lower layer.
        if ( OptRunoffSubsurface /= 1 ) then
           !$acc loop seq
           do LoopInd2 = 1, NumSoilLayer
              SoilLiqTmp(I,LoopInd2,J) = SoilLiqWater(I,LoopInd2,J) * ThicknessSnowSoilLayer(I,LoopInd2,J) * 1000.0
           enddo

           SoilWaterMin = 0.01   ! mm
           !$acc loop seq
           do LoopInd2 = 1, NumSoilLayer-1
              if ( SoilLiqTmp(I,LoopInd2,J) < 0.0 ) then
                 SoilWatRem = SoilWaterMin - SoilLiqTmp(I,LoopInd2,J)
              else
                 SoilWatRem = 0.0
              endif
              SoilLiqTmp(I,LoopInd2  ,J) = SoilLiqTmp(I,LoopInd2  ,J) + SoilWatRem
              SoilLiqTmp(I,LoopInd2+1,J) = SoilLiqTmp(I,LoopInd2+1,J) - SoilWatRem
           enddo
           LoopInd2 = NumSoilLayer
           if ( SoilLiqTmp(I,LoopInd2,J) < SoilWaterMin ) then
               SoilWatRem = SoilWaterMin - SoilLiqTmp(I,LoopInd2,J)
           else
               SoilWatRem = 0.0
           endif
           SoilLiqTmp(I,LoopInd2,J) = SoilLiqTmp(I,LoopInd2,J) + SoilWatRem
           RunoffSubsurface(I,J)     = RunoffSubsurface(I,J) - SoilWatRem/SoilTimeStep

           if ( OptRunoffSubsurface == 5 ) RechargeGwDeepWT(I,J) = RechargeGwDeepWT(I,J) - SoilWatRem * 1.0e-3

           !$acc loop seq
           do LoopInd2 = 1, NumSoilLayer
              SoilLiqWater(I,LoopInd2,J) = SoilLiqTmp(I,LoopInd2,J) / &
                                            (ThicknessSnowSoilLayer(I,LoopInd2,J)*1000.0)
           enddo
        endif ! OptRunoffSubsurface /= 1


      end do
    end do
    !$acc end parallel loop

    !$acc end data
    deallocate(SoilLiqTmp)

    ! compute groundwater and subsurface runoff
    if ( OptRunoffSubsurface == 1 ) call RunoffSubSurfaceGroundWater(noahmp)

    ! compute subsurface runoff based on drainage rate
    if ( (OptRunoffSubsurface == 3) .or. (OptRunoffSubsurface == 4) .or. (OptRunoffSubsurface == 6) .or. &
         (OptRunoffSubsurface == 7) .or. (OptRunoffSubsurface == 8) ) then
         call RunoffSubSurfaceDrainage(noahmp)
    endif

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd2)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points


        ! update soil moisture
        !$acc loop seq
        do LoopInd2 = 1, NumSoilLayer
            SoilMoisture(I,LoopInd2,J) = SoilLiqWater(I,LoopInd2,J) + SoilIce(I,LoopInd2,J)
        enddo


      end do
    end do
    !$acc end parallel loop

    ! compute subsurface runoff and shallow water table for MMF scheme
    if ( OptRunoffSubsurface == 5 ) call RunoffSubSurfaceShallowWaterMMF(noahmp)

    ! accumulated water flux over soil timestep [mm]
    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points
        noahmp%water%flux%RunoffSurface(I,J)    = noahmp%water%flux%RunoffSurface(I,J)    * SoilTimeStep
        noahmp%water%flux%RunoffSubsurface(I,J) = noahmp%water%flux%RunoffSubsurface(I,J) * SoilTimeStep
        noahmp%water%flux%TileDrain(I,J)        = noahmp%water%flux%TileDrain(I,J)        * SoilTimeStep
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    ! deallocate local arrays
    deallocate(MatRight)
    deallocate(MatLeft1)
    deallocate(MatLeft2)
    deallocate(MatLeft3)
    deallocate(SoilSatExcAcc2D)
    deallocate(DrainSoilBotAcc2D)
    deallocate(RunoffSurfaceAcc2D)



    end associate

  end subroutine SoilWaterMain

end module SoilWaterMainMod
