module RunoffSurfaceFreeDrainMod

!!! Calculate inflitration rate at soil surface and surface runoff for free drainage scheme

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SoilHydraulicPropertyMod, only : SoilDiffusivityConductivityOpt2

  implicit none

contains

  subroutine RunoffSurfaceFreeDrain(noahmp, TimeStep)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: INFIL
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

! IN & OUT variabls
    type(noahmp_type)     , intent(inout) :: noahmp
    real(kind=kind_noahmp), intent(in)    :: TimeStep                      ! timestep (may not be the same as model timestep)

! local variable
    integer                :: I, J                                         ! grid indices
    integer                :: IndSoilFrz                                   ! number of interaction
    integer                :: LoopInd1, LoopInd2,  LoopInd3                ! do-loop index
    integer, parameter     :: FrzSoilFac = 3                               ! frozen soil pre-factor
    real(kind=kind_noahmp) :: FracVoidRem                                  ! remaining fraction
    real(kind=kind_noahmp) :: SoilWatHoldMaxRem                            ! remaining accumulated maximum holdable soil water [m]
    real(kind=kind_noahmp) :: WaterInSfc                                   ! surface in water [m]
    real(kind=kind_noahmp) :: TimeStepDay                                  ! time indices
    real(kind=kind_noahmp) :: SoilWatHoldMaxAcc                            ! accumulated maximum holdable soil water [m]
    real(kind=kind_noahmp) :: SoilIceWatTmp                                ! maximum soil ice water [m]
    real(kind=kind_noahmp) :: SoilImpervFrac                               ! impervious fraction due to frozen soil
    real(kind=kind_noahmp) :: IndAcc                                       ! accumulation index
    real(kind=kind_noahmp) :: SoilIceCoeff                                 ! soil ice coefficient
    real(kind=kind_noahmp) :: SoilWatDiffusivity                           ! soil water diffusivity [m2/s]
    real(kind=kind_noahmp) :: SoilWatConductivity                          ! soil water conductivity [m/s]
    real(kind=kind_noahmp) :: SoilWatHoldCap                               ! soil moisture holding capacity [m3/m3]
    real(kind=kind_noahmp) :: InfilRateMax                                 ! maximum infiltration rate [m/s]
    real(kind=kind_noahmp) :: SoilWatMaxHold(1:noahmp%config%domain%NumSoilLayer)  ! maximum soil water that can hold [m]

! --------------------------------------------------------------------

    associate(                                                                  &
              NumSoilLayer        => noahmp%config%domain%NumSoilLayer          ,& ! in,  number of soil layers
              DepthSoilLayer      => noahmp%config%domain%DepthSoilLayer        ,& ! in,  depth [m] of layer-bottom from soil surface
              FlagUrban           => noahmp%config%domain%FlagUrban        ,& ! in,  logical flag for urban grid
              SoilLiqWater        => noahmp%water%state%SoilLiqWater            ,& ! in,  soil water content [m3/m3]
              SoilIce             => noahmp%water%state%SoilIce                 ,& ! in,  soil ice content [m3/m3]
              SoilIceMax          => noahmp%water%state%SoilIceMax         ,& ! in,  maximum soil ice content [m3/m3]
              SoilSfcInflowMean   => noahmp%water%flux%SoilSfcInflowMean   ,& ! in,  water input on soil surface [m/s]
              SoilMoistureSat     => noahmp%water%param%SoilMoistureSat         ,& ! in,  saturated value of soil moisture [m3/m3]
              SoilMoistureWilt    => noahmp%water%param%SoilMoistureWilt        ,& ! in,  wilting point soil moisture [m3/m3]
              SoilInfilMaxCoeff   => noahmp%water%param%SoilInfilMaxCoeff  ,& ! in,  parameter to calculate maximum infiltration rate
              SoilImpervFracCoeff => noahmp%water%param%SoilImpervFracCoeff,& ! in,  parameter to calculate frozen soil impermeable fraction
              RunoffSurface       => noahmp%water%flux%RunoffSurface       ,& ! out, surface runoff [m/s]
              InfilRateSfc        => noahmp%water%flux%InfilRateSfc         & ! out, infiltration rate at surface [m/s]
             )

   !$acc parallel loop collapse(2) gang vector default(present) &
   !$acc private(IndSoilFrz,LoopInd1,LoopInd2,LoopInd3,FracVoidRem,SoilWatHoldMaxRem,WaterInSfc) &
   !$acc private(TimeStepDay,SoilWatHoldMaxAcc,SoilIceWatTmp,SoilImpervFrac,IndAcc,SoilIceCoeff) &
   !$acc private(SoilWatDiffusivity,SoilWatConductivity,SoilWatHoldCap,InfilRateMax,SoilWatMaxHold) &
   !$acc firstprivate(TimeStep)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points


    ! initialize
    SoilWatMaxHold(1:4) = 0.0

    ! start infiltration for free drainage scheme
    if ( SoilSfcInflowMean(I,J) > 0.0 ) then

       TimeStepDay    = TimeStep / 86400.0
       SoilWatHoldCap = SoilMoistureSat(I,1,J) - SoilMoistureWilt(I,1,J)

       ! compute maximum infiltration rate
       SoilWatMaxHold(1) = -DepthSoilLayer(I,1,J) * SoilWatHoldCap
       SoilIceWatTmp     = -DepthSoilLayer(I,1,J) * SoilIce(I,1,J)
       SoilWatMaxHold(1) =  SoilWatMaxHold(1) * (1.0-(SoilLiqWater(I,1,J)+SoilIce(I,1,J)-SoilMoistureWilt(I,1,J)) / SoilWatHoldCap)
       SoilWatHoldMaxAcc =  SoilWatMaxHold(1)
       !$acc loop seq
       do LoopInd3 = 2, NumSoilLayer
          SoilIceWatTmp            = SoilIceWatTmp + (DepthSoilLayer(I,LoopInd3-1,J) - DepthSoilLayer(I,LoopInd3,J))*SoilIce(I,LoopInd3,J)
          SoilWatMaxHold(LoopInd3) = (DepthSoilLayer(I,LoopInd3-1,J) - DepthSoilLayer(I,LoopInd3,J)) * SoilWatHoldCap
          SoilWatMaxHold(LoopInd3) = SoilWatMaxHold(LoopInd3) * (1.0 - (SoilLiqWater(I,LoopInd3,J) + SoilIce(I,LoopInd3,J) - &
                                                                 SoilMoistureWilt(I,LoopInd3,J)) / SoilWatHoldCap)
          SoilWatHoldMaxAcc        = SoilWatHoldMaxAcc + SoilWatMaxHold(LoopInd3)
       enddo
       FracVoidRem       = 1.0 - exp(-1.0 * SoilInfilMaxCoeff(I,J) * TimeStepDay)
       SoilWatHoldMaxRem = SoilWatHoldMaxAcc * FracVoidRem
       WaterInSfc        = max(0.0, SoilSfcInflowMean(I,J) * TimeStep)
       InfilRateMax      = (WaterInSfc * (SoilWatHoldMaxRem/(WaterInSfc + SoilWatHoldMaxRem))) / TimeStep

       ! impermeable fraction due to frozen soil
       SoilImpervFrac = 1.0
       if ( SoilIceWatTmp > 1.0e-2 ) then
          SoilIceCoeff = FrzSoilFac * SoilImpervFracCoeff(I,J) / SoilIceWatTmp
          IndAcc       = 1.0
          IndSoilFrz   = FrzSoilFac - 1
          !$acc loop seq
          do LoopInd1 = 1, IndSoilFrz
             LoopInd3  = 1
             !$acc loop seq
             do LoopInd2 = LoopInd1+1, IndSoilFrz
                LoopInd3 = LoopInd3 * LoopInd2
             enddo
             IndAcc = IndAcc + (SoilIceCoeff ** (FrzSoilFac-LoopInd1)) / float(LoopInd3)
          enddo
          SoilImpervFrac = 1.0 - exp(-SoilIceCoeff) * IndAcc
       endif

       ! correction of infiltration limitation
       InfilRateMax = InfilRateMax * SoilImpervFrac
       ! jref for urban areas
       ! if ( FlagUrban .eqv. .true. ) InfilRateMax == InfilRateMax * 0.05

       ! soil hydraulic conductivity and diffusivity
       call SoilDiffusivityConductivityOpt2(noahmp, SoilWatDiffusivity, SoilWatConductivity, SoilLiqWater(I,1,J), SoilIceMax(I,J), 1, I, J)

       InfilRateMax = max(InfilRateMax, SoilWatConductivity)
       InfilRateMax = min(InfilRateMax, WaterInSfc/TimeStep)

       ! compute surface runoff and infiltration rate
       RunoffSurface(I,J) = max(0.0, SoilSfcInflowMean(I,J)-InfilRateMax)
       InfilRateSfc(I,J)  = SoilSfcInflowMean(I,J) - RunoffSurface(I,J)

    endif ! SoilSfcInflowMean(I,J) > 0.0


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine RunoffSurfaceFreeDrain

end module RunoffSurfaceFreeDrainMod
