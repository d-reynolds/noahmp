module WaterTableEquilibriumMod

!!! Calculate equilibrium water table depth (Niu et al., 2005)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine WaterTableEquilibrium(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: ZWTEQ
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: IndSoil                           ! do-loop index
    integer, parameter               :: NumSoilFineLy = 100               ! no. of fine soil layers of 6m soil
    real(kind=kind_noahmp)           :: WatDeficitCoarse                  ! water deficit from coarse (4-L) soil moisture profile
    real(kind=kind_noahmp)           :: WatDeficitFine                    ! water deficit from fine (100-L) soil moisture profile
    real(kind=kind_noahmp)           :: ThickSoilFineLy                   ! layer thickness of the 100-L soil layers to 6.0 m
    real(kind=kind_noahmp)           :: TmpVar                            ! temporary variable
    real(kind=kind_noahmp), dimension(1:NumSoilFineLy) :: DepthSoilFineLy ! layer-bottom depth of the 100-L soil layers to 6.0 m
    integer                          :: I, J                              ! grid indices
    associate(                                                                       &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,  number of soil layers
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,  depth [m] of layer-bottom from soil surface
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,  thickness of snow/soil layers [m]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! in,  soil water content [m3/m3]
              SoilMoistureSat        => noahmp%water%param%SoilMoistureSat          ,& ! in,  saturated value of soil moisture [m3/m3]
              SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat      ,& ! in,  saturated soil matric potential [m]
              SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB            ,& ! in,  soil B parameter
              WaterTableDepth        => noahmp%water%state%WaterTableDepth      & ! out, water table depth [m]
             )

    !$acc parallel loop collapse(2) gang vector default(present) private(DepthSoilFineLy) private(IndSoil, &
    !$acc ThickSoilFineLy, TmpVar, WatDeficitCoarse, WatDeficitFine)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points
    !$acc loop seq
    do IndSoil = 1, NumSoilFineLy
      DepthSoilFineLy(IndSoil) = 0.0
    enddo
    WatDeficitCoarse                 = 0.0
    !$acc loop seq
    do IndSoil = 1, NumSoilLayer
       WatDeficitCoarse = WatDeficitCoarse + (SoilMoistureSat(I,1,J) - SoilLiqWater(I,IndSoil,J)) * &
                                             ThicknessSnowSoilLayer(I,IndSoil,J)   ! [m]
    enddo

    ThickSoilFineLy = 3.0 * (-DepthSoilLayer(I,NumSoilLayer,J)) / NumSoilFineLy
    !$acc loop seq
    do IndSoil = 1, NumSoilFineLy
       DepthSoilFineLy(IndSoil) = float(IndSoil) * ThickSoilFineLy
    enddo

    WaterTableDepth(I,J) = -3.0 * DepthSoilLayer(I,NumSoilLayer,J) - 0.001              ! initial value [m]

    WatDeficitFine = 0.0
    !$acc loop seq
    do IndSoil = 1, NumSoilFineLy
       TmpVar         = 1.0 + (WaterTableDepth(I,J) - DepthSoilFineLy(IndSoil)) / SoilMatPotentialSat(I,1,J)
       WatDeficitFine = WatDeficitFine + SoilMoistureSat(I,1,J) * &
                                         (1.0 - TmpVar**(-1.0/SoilExpCoeffB(I,1,J))) * ThickSoilFineLy
       if ( abs(WatDeficitFine-WatDeficitCoarse) <= 0.01 ) then
          WaterTableDepth(I,J) = DepthSoilFineLy(IndSoil)
          exit
       endif
    enddo


   enddo
enddo


    end associate

  end subroutine WaterTableEquilibrium

end module WaterTableEquilibriumMod
