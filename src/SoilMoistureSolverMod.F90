module SoilMoistureSolverMod

!!! Compute soil moisture content using based on Richards diffusion & tri-diagonal matrix
!!! Dependent on the output from SoilWaterDiffusionRichards subroutine

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use MatrixSolverTriDiagonalMod, only : MatrixSolverTriDiagonal

  implicit none

contains

  subroutine SoilMoistureSolver(noahmp, TimeStep, MatLeft1, MatLeft2, MatLeft3, MatRight)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: SSTEP
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)     , intent(inout) :: noahmp
    real(kind=kind_noahmp), intent(in)    :: TimeStep                               ! timestep (may not be the same as model timestep)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatRight    ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft1    ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft2    ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft3    ! left-hand side term of the matrix

! local variable
    integer                                           :: LoopInd                    ! soil layer loop index 
    real(kind=kind_noahmp)                            :: WatDefiTmp                 ! temporary water deficiency
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatRightTmp                ! temporary MatRight matrix coefficient
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: MatLeft3Tmp                ! temporary MatLeft3 matrix coefficient
    integer                                           :: I, J                        ! grid indices
    associate(                                                                       &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,    depth [m] of layer-bottom from soil surface
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
              OptRunoffSubsurface    => noahmp%config%nmlist%OptRunoffSubsurface    ,& ! in,    options for drainage and subsurface runoff
              SoilMoistureSat        => noahmp%water%param%SoilMoistureSat          ,& ! in,    saturated value of soil moisture [m3/m3]
              WaterTableDepth        => noahmp%water%state%WaterTableDepth          ,& ! in,    water table depth [m]
              SoilIce                => noahmp%water%state%SoilIce                  ,& ! in,    soil ice content [m3/m3]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! inout, soil water content [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! inout, total soil moisture [m3/m3]
              SoilMoistureToWT       => noahmp%water%state%SoilMoistureToWT         ,& ! inout, soil moisture between bottom of soil & water table
              RechargeGwDeepWT       => noahmp%water%state%RechargeGwDeepWT         ,& ! inout, recharge to or from the water table when deep [m]
              DrainSoilBot           => noahmp%water%flux%DrainSoilBot              ,& ! inout, soil bottom drainage (m/s)
              SoilEffPorosity        => noahmp%water%state%SoilEffPorosity               ,& ! out,   soil effective porosity (m3/m3)
              SoilSaturationExcess   => noahmp%water%state%SoilSaturationExcess      & ! out,   saturation excess of the total soil [m]
             )

    allocate(MatRightTmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                         1:NumSoilLayer, &
                         noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(MatLeft3Tmp(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                         1:NumSoilLayer, &
                         noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    !$acc data create(MatRightTmp, MatLeft3Tmp)

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd, &
    !$acc WatDefiTmp) firstprivate(TimeStep)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points

    ! initialization
    SoilSaturationExcess(I,J) = 0.0

    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       MatRightTmp(I,LoopInd,J) = 0.0
       MatLeft3Tmp(I,LoopInd,J) = 0.0
       SoilEffPorosity(I,LoopInd,J) = 0.0
    enddo

    ! update tri-diagonal matrix elements
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       MatRight(I,LoopInd,J) =       MatRight(I,LoopInd,J) * TimeStep
       MatLeft1(I,LoopInd,J) =       MatLeft1(I,LoopInd,J) * TimeStep
       MatLeft2(I,LoopInd,J) = 1.0 + MatLeft2(I,LoopInd,J) * TimeStep
       MatLeft3(I,LoopInd,J) =       MatLeft3(I,LoopInd,J) * TimeStep

       ! copy values for input variables before calling rosr12
       MatRightTmp(I,LoopInd,J) = MatRight(I,LoopInd,J)
       MatLeft3Tmp(I,LoopInd,J) = MatLeft3(I,LoopInd,J)
    enddo

    ! call ROSR12 to solve the tri-diagonal matrix
    call MatrixSolverTriDiagonal(MatLeft3(I,:,J),MatLeft1(I,:,J),MatLeft2(I,:,J),MatLeft3Tmp(I,:,J),MatRightTmp(I,:,J),MatRight(I,:,J),1,NumSoilLayer,0)

    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
        SoilLiqWater(I,LoopInd,J) = SoilLiqWater(I,LoopInd,J) + MatLeft3(I,LoopInd,J)
    enddo

    !  excessive water above saturation in a layer is moved to
    !  its unsaturated layer like in a bucket

    ! for MMF scheme, there is soil moisture below NumSoilLayer, to the water table
    if ( OptRunoffSubsurface == 5 ) then
       ! update SoilMoistureToWT
       if ( WaterTableDepth(I,J) < (DepthSoilLayer(I,NumSoilLayer,J)-ThicknessSnowSoilLayer(I,NumSoilLayer,J)) ) then
          ! accumulate soil drainage to update deep water table and soil moisture later
          RechargeGwDeepWT(I,J)           = RechargeGwDeepWT(I,J) + TimeStep * DrainSoilBot(I,J)
       else
          SoilMoistureToWT(I,J)           = SoilMoistureToWT(I,J) + &
                                       TimeStep * DrainSoilBot(I,J) / ThicknessSnowSoilLayer(I,NumSoilLayer,J)
          SoilSaturationExcess(I,J)       = max((SoilMoistureToWT(I,J) - SoilMoistureSat(I,NumSoilLayer,J)), 0.0) * &
                                       ThicknessSnowSoilLayer(I,NumSoilLayer,J)
          WatDefiTmp                 = max((1.0e-4 - SoilMoistureToWT(I,J)), 0.0) * ThicknessSnowSoilLayer(I,NumSoilLayer,J)
          SoilMoistureToWT(I,J)           = max(min(SoilMoistureToWT(I,J), SoilMoistureSat(I,NumSoilLayer,J)), 1.0e-4)
          SoilLiqWater(I,NumSoilLayer,J) = SoilLiqWater(I,NumSoilLayer,J) + &
                                       SoilSaturationExcess(I,J) / ThicknessSnowSoilLayer(I,NumSoilLayer,J)
          ! reduce fluxes at the bottom boundaries accordingly
          DrainSoilBot(I,J)               = DrainSoilBot(I,J) - SoilSaturationExcess(I,J)/TimeStep
          RechargeGwDeepWT(I,J)           = RechargeGwDeepWT(I,J) - WatDefiTmp
       endif
    endif

    !$acc loop seq
    do LoopInd = NumSoilLayer, 2, -1
       SoilEffPorosity(I,LoopInd,J) = max(1.0e-4, (SoilMoistureSat(I,LoopInd,J) - SoilIce(I,LoopInd,J)))
       SoilSaturationExcess(I,J)     = max((SoilLiqWater(I,LoopInd,J)-SoilEffPorosity(I,LoopInd,J)), 0.0) * &
                                  ThicknessSnowSoilLayer(I,LoopInd,J)
       SoilLiqWater(I,LoopInd,J)    = min(SoilEffPorosity(I,LoopInd,J), SoilLiqWater(I,LoopInd,J) )
       SoilLiqWater(I,LoopInd-1,J)  = SoilLiqWater(I,LoopInd-1,J) + SoilSaturationExcess(I,J) / ThicknessSnowSoilLayer(I,LoopInd-1,J)
    enddo

    SoilEffPorosity(I,1,J)   = max(1.0e-4, (SoilMoistureSat(I,1,J)-SoilIce(I,1,J)))
    SoilSaturationExcess(I,J) = max((SoilLiqWater(I,1,J)-SoilEffPorosity(I,1,J)), 0.0) * ThicknessSnowSoilLayer(I,1,J)
    SoilLiqWater(I,1,J)      = min(SoilEffPorosity(I,1,J), SoilLiqWater(I,1,J))

    if ( SoilSaturationExcess(I,J) > 0.0 ) then
       SoilLiqWater(I,2,J) = SoilLiqWater(I,2,J) + SoilSaturationExcess(I,J) / ThicknessSnowSoilLayer(I,2,J)
       !$acc loop seq
       do LoopInd = 2, NumSoilLayer-1
          SoilEffPorosity(I,LoopInd,J) = max(1.0e-4, (SoilMoistureSat(I,LoopInd,J) - SoilIce(I,LoopInd,J)))
          SoilSaturationExcess(I,J)     = max((SoilLiqWater(I,LoopInd,J)-SoilEffPorosity(I,LoopInd,J)), 0.0) * &
                                     ThicknessSnowSoilLayer(I,LoopInd,J)
          SoilLiqWater(I,LoopInd,J)    = min(SoilEffPorosity(I,LoopInd,J), SoilLiqWater(I,LoopInd,J))
          SoilLiqWater(I,LoopInd+1,J)  = SoilLiqWater(I,LoopInd+1,J) + SoilSaturationExcess(I,J) / ThicknessSnowSoilLayer(I,LoopInd+1,J)
       enddo
       SoilEffPorosity(I,NumSoilLayer,J) = max(1.0e-4, (SoilMoistureSat(I,NumSoilLayer,J) - SoilIce(I,NumSoilLayer,J)))
       SoilSaturationExcess(I,J)          = max((SoilLiqWater(I,NumSoilLayer,J)-SoilEffPorosity(I,NumSoilLayer,J)), 0.0) * &
                                       ThicknessSnowSoilLayer(I,NumSoilLayer,J)
       SoilLiqWater(I,NumSoilLayer,J)    = min(SoilEffPorosity(I,NumSoilLayer,J), SoilLiqWater(I,NumSoilLayer,J))
    endif

    SoilMoisture = SoilLiqWater + SoilIce


      end do
    end do

    !$acc end data
    deallocate(MatRightTmp)
    deallocate(MatLeft3Tmp)

    end associate

  end subroutine SoilMoistureSolver

end module SoilMoistureSolverMod
