module SoilSnowThermalDiffusionMod

!!! Solve soil and snow layer thermal diffusion
!!! Calculate the right hand side of the time tendency term of the soil
!!! and snow thermal diffusion equation. Currently snow and soil layers
!!! are coupled in solving the equations. Also compute/prepare the matrix
!!! coefficients for the tri-diagonal matrix of the implicit time scheme.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SoilSnowThermalDiffusion(noahmp, SoilTimeStep, MatLeft1, MatLeft2, MatLeft3, MatRight)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: HRT
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)     , intent(inout) :: noahmp
    real(kind=kind_noahmp), intent(in)    :: SoilTimeStep                             ! SoilTimeStep (may not be the same as model timestep)
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatRight  ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft1  ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft2  ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft3  ! left-hand side term of the matrix
! local variable
    integer                               :: I, J              ! grid indices
    integer                               :: LoopInd           ! loop index
    real(kind=kind_noahmp)                :: DepthSnowSoilTmp  ! temporary snow/soil layer depth [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: DepthSnowSoilInv   ! inverse of snow/soil layer depth [1/m]
    real(kind=kind_noahmp)                :: HeatCapacPerArea   ! Heat capacity per area [J/m2/K]
    real(kind=kind_noahmp)                :: TempGradDepth      ! temperature gradient [K/m]
    real(kind=kind_noahmp)                :: TempGradDepthPrev     ! previous temperature gradient
    real(kind=kind_noahmp)                :: EnergyExcess       ! energy flux excess [W/m2]

    associate(                                                                            &
              NumSoilLayer             => noahmp%config%domain%NumSoilLayer              ,& ! in,  number of soil layers
              NumSnowLayerMax          => noahmp%config%domain%NumSnowLayerMax           ,& ! in,  maximum number of snow layers
              NumSnowLayerNeg          => noahmp%config%domain%NumSnowLayerNeg      ,& ! in,  actual number of snow layers (negative)
              DepthSnowSoilLayer       => noahmp%config%domain%DepthSnowSoilLayer        ,& ! in,  depth of snow/soil layer-bottom [m]
              OptSoilTemperatureBottom => noahmp%config%nmlist%OptSoilTemperatureBottom  ,& ! in,  options for lower boundary condition of soil temp.
              OptSnowSoilTempTime      => noahmp%config%nmlist%OptSnowSoilTempTime       ,& ! in,  options for snow/soil temperature time scheme
              TemperatureSoilBottom    => noahmp%forcing%TemperatureSoilBottom      ,& ! in,  bottom boundary soil temperature [K]
              DepthSoilTempBotToSno    => noahmp%energy%state%DepthSoilTempBotToSno ,& ! in,  depth of lower boundary condition [m] from snow surface
              TemperatureSoilSnow      => noahmp%energy%state%TemperatureSoilSnow        ,& ! in,  snow and soil layer temperature [K]
              ThermConductSoilSnow     => noahmp%energy%state%ThermConductSoilSnow       ,& ! in,  thermal conductivity [W/m/K] for all soil & snow
              HeatCapacSoilSnow        => noahmp%energy%state%HeatCapacSoilSnow          ,& ! in,  heat capacity [J/m3/K] for all soil & snow
              HeatGroundTotMean        => noahmp%energy%flux%HeatGroundTotMean      ,& ! in,  total ground heat flux [W/m2] averaged during soil timestep
              RadSwPenetrateGrd        => noahmp%energy%flux%RadSwPenetrateGrd           ,& ! in,  light penetrating through soil/snow water [W/m2]
              HeatFromSoilBot          => noahmp%energy%flux%HeatFromSoilBot         & ! out, energy influx from soil bottom [W/m2]
             )

    allocate(DepthSnowSoilInv(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
                             -NumSnowLayerMax+1:NumSoilLayer, &
                             noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    !$acc data create(DepthSnowSoilInv)

    !$acc parallel loop collapse(2) gang vector default(present) private(DepthSnowSoilTmp, &
    !$acc EnergyExcess, HeatCapacPerArea, LoopInd, TempGradDepth, TempGradDepthPrev) firstprivate(SoilTimeStep)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if (noahmp%config%domain%IndicatorIceSfc(I,J) == -1) cycle  ! glacier points handled by GlacierThermalDiffusion

    ! initialization
    !$acc loop seq
    do LoopInd = -NumSnowLayerMax+1, NumSoilLayer
       MatRight(I,LoopInd,J) = 0.0
       MatLeft1(I,LoopInd,J) = 0.0
       MatLeft2(I,LoopInd,J) = 0.0
       MatLeft3(I,LoopInd,J) = 0.0
    enddo
    HeatCapacPerArea = 0.0
    TempGradDepth = 0.0
    EnergyExcess = 0.0
    TempGradDepthPrev = 0.0

    !$acc loop seq
     do LoopInd = -NumSnowLayerMax+1, NumSoilLayer
         DepthSnowSoilInv(I,LoopInd,J) = 0.0
     enddo

    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       if ( LoopInd == (NumSnowLayerNeg(I,J)+1) ) then
          DepthSnowSoilTmp          = - DepthSnowSoilLayer(I,LoopInd+1,J)
       elseif ( LoopInd < NumSoilLayer ) then
          DepthSnowSoilTmp          = DepthSnowSoilLayer(I,LoopInd-1,J) - DepthSnowSoilLayer(I,LoopInd+1,J)
       endif
          DepthSnowSoilInv(I,LoopInd,J) = 2.0 / DepthSnowSoilTmp
    enddo

    ! compute gradient and flux of soil/snow thermal diffusion
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       if ( LoopInd == (NumSnowLayerNeg(I,J)+1) ) then
          HeatCapacPerArea = - DepthSnowSoilLayer(I,LoopInd,J) * HeatCapacSoilSnow(I,LoopInd,J)
          DepthSnowSoilTmp = - DepthSnowSoilLayer(I,LoopInd+1,J)
          TempGradDepth = 2.0 * (TemperatureSoilSnow(I,LoopInd,J) - TemperatureSoilSnow(I,LoopInd+1,J)) / DepthSnowSoilTmp
          EnergyExcess = ThermConductSoilSnow(I,LoopInd,J) * TempGradDepth - &
                            HeatGroundTotMean(I,J) - RadSwPenetrateGrd(I,LoopInd,J)
       elseif ( LoopInd < NumSoilLayer ) then
          HeatCapacPerArea = (DepthSnowSoilLayer(I,LoopInd-1,J) - DepthSnowSoilLayer(I,LoopInd,J)) * HeatCapacSoilSnow(I,LoopInd,J)
          DepthSnowSoilTmp = DepthSnowSoilLayer(I,LoopInd-1,J) - DepthSnowSoilLayer(I,LoopInd+1,J)
          TempGradDepth = 2.0 * (TemperatureSoilSnow(I,LoopInd,J) - TemperatureSoilSnow(I,LoopInd+1,J)) / DepthSnowSoilTmp
          EnergyExcess = (ThermConductSoilSnow(I,LoopInd,J)*TempGradDepth - &
                            ThermConductSoilSnow(I,LoopInd-1,J) * TempGradDepthPrev ) - RadSwPenetrateGrd(I,LoopInd,J)
       elseif ( LoopInd == NumSoilLayer ) then
          HeatCapacPerArea = (DepthSnowSoilLayer(I,LoopInd-1,J) - DepthSnowSoilLayer(I,LoopInd,J)) * HeatCapacSoilSnow(I,LoopInd,J)
          DepthSnowSoilTmp = DepthSnowSoilLayer(I,LoopInd-1,J) - DepthSnowSoilLayer(I,LoopInd,J)
          if ( OptSoilTemperatureBottom == 1 ) then
             HeatFromSoilBot(I,J) = 0.0
          endif
          if ( OptSoilTemperatureBottom == 2 ) then
             TempGradDepth = (TemperatureSoilSnow(I,LoopInd,J) - TemperatureSoilBottom(I,J)) / &
                                (0.5*(DepthSnowSoilLayer(I,LoopInd-1,J)+DepthSnowSoilLayer(I,LoopInd,J)) - DepthSoilTempBotToSno(I,J))
             HeatFromSoilBot(I,J) = -ThermConductSoilSnow(I,LoopInd,J) * TempGradDepth
          endif
          EnergyExcess = (-HeatFromSoilBot(I,J) - ThermConductSoilSnow(I,LoopInd-1,J) * TempGradDepthPrev) - &
                            RadSwPenetrateGrd(I,LoopInd,J)
       endif

       ! prepare the matrix coefficients for the tri-diagonal matrix
       if ( LoopInd == (NumSnowLayerNeg(I,J)+1) ) then
          MatLeft1(I,LoopInd,J) = 0.0
          MatLeft3(I,LoopInd,J) = - ThermConductSoilSnow(I,LoopInd,J) * DepthSnowSoilInv(I,LoopInd,J) / HeatCapacPerArea
          if ( (OptSnowSoilTempTime == 1) .or. (OptSnowSoilTempTime == 3) ) then
             MatLeft2(I,LoopInd,J) = - MatLeft3(I,LoopInd,J)
          endif
          if ( OptSnowSoilTempTime == 2 ) then
             MatLeft2(I,LoopInd,J) = - MatLeft3(I,LoopInd,J) + ThermConductSoilSnow(I,LoopInd,J) / &
                            (0.5*DepthSnowSoilLayer(I,LoopInd,J)*DepthSnowSoilLayer(I,LoopInd,J)*HeatCapacSoilSnow(I,LoopInd,J))
          endif
       elseif ( LoopInd < NumSoilLayer ) then
          MatLeft1(I,LoopInd,J) = - ThermConductSoilSnow(I,LoopInd-1,J) * DepthSnowSoilInv(I,LoopInd-1,J) / HeatCapacPerArea
          MatLeft3(I,LoopInd,J) = - ThermConductSoilSnow(I,LoopInd  ,J) * DepthSnowSoilInv(I,LoopInd,J) / HeatCapacPerArea
          MatLeft2(I,LoopInd,J) = - (MatLeft1(I,LoopInd,J) + MatLeft3(I,LoopInd,J))
       elseif ( LoopInd == NumSoilLayer ) then
          MatLeft1(I,LoopInd,J) = - ThermConductSoilSnow(I,LoopInd-1,J) * DepthSnowSoilInv(I,LoopInd-1,J) / HeatCapacPerArea
          MatLeft3(I,LoopInd,J) = 0.0
          MatLeft2(I,LoopInd,J) = - (MatLeft1(I,LoopInd,J) + MatLeft3(I,LoopInd,J))
       endif
       MatRight(I,LoopInd,J) = EnergyExcess / (-HeatCapacPerArea)

       ! save for next iteration
       TempGradDepthPrev = TempGradDepth
    enddo
    ! accumulate soil bottom flux for soil timestep
    HeatFromSoilBot(I,J) = HeatFromSoilBot(I,J) * SoilTimeStep

      end do
    end do
    !$acc end parallel loop

    !$acc end data
    deallocate(DepthSnowSoilInv)

    end associate

  end subroutine SoilSnowThermalDiffusion

end module SoilSnowThermalDiffusionMod
