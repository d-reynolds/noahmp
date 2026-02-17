module RunoffSurfaceXinAnJiangMod

!!! Compute surface infiltration rate and surface runoff based on XinAnJiang runoff scheme
!!! Reference: Knoben, W. J., et al., (2019): Modular Assessment of Rainfall-Runoff Models 
!!! Toolbox (MARRMoT) v1.2 an open-source, extendable framework providing implementations 
!!! of 46 conceptual hydrologic models as continuous state-space formulations.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine RunoffSurfaceXinAnJiang(noahmp, TimeStep)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: COMPUTE_XAJ_SURFRUNOFF
! Original code: Prasanth Valayamkunnath <prasanth@ucar.edu>
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

! IN & OUT variables
    type(noahmp_type)     , intent(inout) :: noahmp
    real(kind=kind_noahmp), intent(in)    :: TimeStep            ! timestep (may not be the same as model timestep)

! local variable
    integer                               :: I, J                ! grid indices
    integer                               :: LoopInd             ! do-loop index
    real(kind=kind_noahmp)                :: SoilWaterTmp        ! temporary soil water [m]
    real(kind=kind_noahmp)                :: SoilWaterMax        ! maximum soil water [m]
    real(kind=kind_noahmp)                :: SoilWaterFree       ! free soil water [m]
    real(kind=kind_noahmp)                :: SoilWaterFreeMax    ! maximum free soil water [m]
    real(kind=kind_noahmp)                :: RunoffSfcImp        ! impervious surface runoff [m]
    real(kind=kind_noahmp)                :: RunoffSfcPerv       ! pervious surface runoff [m]

    associate(                                                                 &
              NumSoilLayer         => noahmp%config%domain%NumSoilLayer       ,& ! in,  number of soil layers
              DepthSoilLayer       => noahmp%config%domain%DepthSoilLayer     ,& ! in,  depth [m] of layer-bottom from soil surface
              SoilMoisture         => noahmp%water%state%SoilMoisture         ,& ! in,  total soil moisture [m3/m3]
              SoilImpervFrac       => noahmp%water%state%SoilImpervFrac       ,& ! in,  fraction of imperviousness due to frozen soil
              SoilSfcInflowMean    => noahmp%water%flux%SoilSfcInflowMean     ,& ! in,  mean water input on soil surface [m/s]
              SoilMoistureSat      => noahmp%water%param%SoilMoistureSat      ,& ! in,  saturated value of soil moisture [m3/m3]
              SoilMoistureFieldCap => noahmp%water%param%SoilMoistureFieldCap ,& ! in,  reference soil moisture (field capacity) [m3/m3]
              TensionWatDistrInfl  => noahmp%water%param%TensionWatDistrInfl  ,& ! in,  Tension water distribution inflection parameter
              TensionWatDistrShp   => noahmp%water%param%TensionWatDistrShp   ,& ! in,  Tension water distribution shape parameter
              FreeWatDistrShp      => noahmp%water%param%FreeWatDistrShp      ,& ! in,  Free water distribution shape parameter
              RunoffSurface        => noahmp%water%flux%RunoffSurface         ,& ! out, surface runoff [m/s]
              InfilRateSfc         => noahmp%water%flux%InfilRateSfc           & ! out, infiltration rate at surface [m/s]
             )

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd, RunoffSfcImp, RunoffSfcPerv, &
    !$acc SoilWaterFree, SoilWaterFreeMax, SoilWaterMax, SoilWaterTmp) firstprivate(TimeStep)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points

    ! initialization 
    SoilWaterTmp     = 0.0
    SoilWaterMax     = 0.0
    SoilWaterFree    = 0.0
    SoilWaterFreeMax = 0.0
    RunoffSfcImp     = 0.0
    RunoffSfcPerv    = 0.0
    RunoffSurface(I,J)    = 0.0
    InfilRateSfc(I,J)     = 0.0

    do LoopInd = 1, NumSoilLayer-2
       if ( (SoilMoisture(I,LoopInd,J)-SoilMoistureFieldCap(I,LoopInd,J)) > 0.0 ) then   ! soil moisture greater than field capacity
          SoilWaterFree = SoilWaterFree + &
                          (SoilMoisture(I,LoopInd,J)-SoilMoistureFieldCap(I,LoopInd,J)) * (-1.0) * DepthSoilLayer(I,LoopInd,J)
          SoilWaterTmp  = SoilWaterTmp + SoilMoistureFieldCap(I,LoopInd,J) * (-1.0) * DepthSoilLayer(I,LoopInd,J)
       else
          SoilWaterTmp  = SoilWaterTmp + SoilMoisture(I,LoopInd,J) * (-1.0) * DepthSoilLayer(I,LoopInd,J)
       endif
       SoilWaterMax     = SoilWaterMax + SoilMoistureFieldCap(I,LoopInd,J) * (-1.0) * DepthSoilLayer(I,LoopInd,J)
       SoilWaterFreeMax = SoilWaterFreeMax + &
                          (SoilMoistureSat(I,LoopInd,J)-SoilMoistureFieldCap(I,LoopInd,J)) * (-1.0) * DepthSoilLayer(I,LoopInd,J)
    enddo
    SoilWaterTmp  = min(SoilWaterTmp, SoilWaterMax)      ! tension water [m] 
    SoilWaterFree = min(SoilWaterFree, SoilWaterFreeMax) ! free water [m]

    ! impervious surface runoff R_IMP    
    RunoffSfcImp = SoilImpervFrac(I,1,J) * SoilSfcInflowMean(I,J) * TimeStep

    ! solve pervious surface runoff (m) based on Eq. (310)
    if ( (SoilWaterTmp/SoilWaterMax) <= (0.5-TensionWatDistrInfl(I,J)) ) then
       RunoffSfcPerv = (1.0-SoilImpervFrac(I,1,J)) * SoilSfcInflowMean(I,J) * TimeStep * &
                       ((0.5-TensionWatDistrInfl(I,J))**(1.0-TensionWatDistrShp(I,J))) * &
                       ((SoilWaterTmp/SoilWaterMax)**TensionWatDistrShp(I,J))
    else
       RunoffSfcPerv = (1.0-SoilImpervFrac(I,1,J)) * SoilSfcInflowMean(I,J) * TimeStep * &
                       (1.0-(((0.5+TensionWatDistrInfl(I,J))**(1.0-TensionWatDistrShp(I,J))) * &
                       ((1.0-(SoilWaterTmp/SoilWaterMax))**TensionWatDistrShp(I,J))))
    endif

    ! estimate surface runoff based on Eq. (313)
    if ( SoilSfcInflowMean(I,J) == 0.0 ) then
      RunoffSurface(I,J) = 0.0
    else
      RunoffSurface(I,J) = RunoffSfcPerv * (1.0-((1.0-(SoilWaterFree/SoilWaterFreeMax))**FreeWatDistrShp(I,J))) + RunoffSfcImp
    endif
    RunoffSurface(I,J) = RunoffSurface(I,J) / TimeStep
    RunoffSurface(I,J) = max(0.0,RunoffSurface(I,J))
    RunoffSurface(I,J) = min(SoilSfcInflowMean(I,J), RunoffSurface(I,J))
    InfilRateSfc(I,J)  = SoilSfcInflowMean(I,J) - RunoffSurface(I,J)

      end do
    end do


    end associate

  end subroutine RunoffSurfaceXinAnJiang

end module RunoffSurfaceXinAnJiangMod
