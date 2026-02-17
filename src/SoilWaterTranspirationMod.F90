module SoilWaterTranspirationMod

!!! compute soil water transpiration factor that will be used for 
!!! stomata resistance and evapotranspiration calculations

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SoilWaterTranspiration(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in ENERGY subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: I, J          ! grid indices
    integer                          :: IndSoil       ! loop index
    real(kind=kind_noahmp)           :: SoilWetFac    ! temporary variable
    real(kind=kind_noahmp)           :: MinThr        ! minimum threshold to prevent divided by zero


    associate(                                                                             &
              SurfaceType               => noahmp%config%domain%SurfaceType               ,& ! in,  surface type 1-soil; 2-lake
              ThicknessSnowSoilLayer    => noahmp%config%domain%ThicknessSnowSoilLayer    ,& ! in,  thickness of snow/soil layers [m]
              DepthSoilLayer            => noahmp%config%domain%DepthSoilLayer            ,& ! in,  depth [m] of layer-bottom from soil surface
              OptSoilWaterTranspiration => noahmp%config%nmlist%OptSoilWaterTranspiration ,& ! in,  option for soil moisture factor for stomatal resistance & ET
              NumSoilLayerRoot          => noahmp%water%param%NumSoilLayerRoot            ,& ! in,  number of soil layers with root present
              SoilMoistureWilt          => noahmp%water%param%SoilMoistureWilt            ,& ! in,  wilting point soil moisture [m3/m3]
              SoilMoistureFieldCap      => noahmp%water%param%SoilMoistureFieldCap        ,& ! in,  reference soil moisture (field capacity) [m3/m3]
              SoilMatPotentialWilt      => noahmp%water%param%SoilMatPotentialWilt        ,& ! in,  soil metric potential for wilting point [m]
              SoilMatPotentialSat       => noahmp%water%param%SoilMatPotentialSat         ,& ! in,  saturated soil matric potential [m]
              SoilMoistureSat           => noahmp%water%param%SoilMoistureSat             ,& ! in,  saturated value of soil moisture [m3/m3]
              SoilExpCoeffB             => noahmp%water%param%SoilExpCoeffB               ,& ! in,  soil B parameter              
              SoilLiqWater              => noahmp%water%state%SoilLiqWater                ,& ! in,  soil water content [m3/m3]
              SoilTranspFac             => noahmp%water%state%SoilTranspFac               ,& ! out, soil water transpiration factor (0 to 1)
              SoilTranspFacAcc          => noahmp%water%state%SoilTranspFacAcc            ,& ! out, accumulated soil water transpiration factor (0 to 1)
              SoilMatPotential          => noahmp%water%state%SoilMatPotential             & ! out, soil matrix potential [m]
             )

    !$acc parallel loop collapse(2) gang vector default(present) &
    !$acc private(IndSoil,SoilWetFac,MinThr)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! soil moisture factor controlling stomatal resistance and evapotranspiration
    MinThr         = 1.0e-6
    SoilTranspFacAcc(I,J) = 0.0

    ! only for soil point
    if ( SurfaceType(I,J) ==1 ) then
       !$acc loop seq
       do IndSoil = 1, NumSoilLayerRoot(I,J)
          if ( OptSoilWaterTranspiration == 1 ) then  ! Noah
             SoilWetFac = (SoilLiqWater(I,IndSoil,J) - SoilMoistureWilt(I,IndSoil,J)) / &
                          (SoilMoistureFieldCap(I,IndSoil,J) - SoilMoistureWilt(I,IndSoil,J))
          endif
          if ( OptSoilWaterTranspiration == 2 ) then  ! CLM
             SoilMatPotential(I,IndSoil,J) = max(SoilMatPotentialWilt(I,J), -SoilMatPotentialSat(I,IndSoil,J) * &
                                                (max(0.01,SoilLiqWater(I,IndSoil,J))/SoilMoistureSat(I,IndSoil,J)) ** &
                                                (-SoilExpCoeffB(I,IndSoil,J)))
             SoilWetFac                = (1.0 - SoilMatPotential(I,IndSoil,J)/SoilMatPotentialWilt(I,J)) / &
                                         (1.0 + SoilMatPotentialSat(I,IndSoil,J)/SoilMatPotentialWilt(I,J))
          endif
          if ( OptSoilWaterTranspiration == 3 ) then  ! SSiB
             SoilMatPotential(I,IndSoil,J) = max(SoilMatPotentialWilt(I,J), -SoilMatPotentialSat(I,IndSoil,J) * &
                                                (max(0.01,SoilLiqWater(I,IndSoil,J))/SoilMoistureSat(I,IndSoil,J)) ** &
                                                (-SoilExpCoeffB(I,IndSoil,J)))
             SoilWetFac = 1.0 - exp(-5.8*(log(SoilMatPotentialWilt(I,J)/SoilMatPotential(I,IndSoil,J))))
          endif
          SoilWetFac    = min(1.0, max(0.0,SoilWetFac))

          SoilTranspFac(I,IndSoil,J) = max(MinThr, ThicknessSnowSoilLayer(I,IndSoil,J) / &
                                          (-DepthSoilLayer(I,NumSoilLayerRoot(I,J),J)) * SoilWetFac)
          SoilTranspFacAcc(I,J)           = SoilTranspFacAcc(I,J) + SoilTranspFac(I,IndSoil,J)
       enddo

       SoilTranspFacAcc(I,J) = max(MinThr, SoilTranspFacAcc(I,J))
       !$acc loop seq
       do IndSoil = 1, NumSoilLayerRoot(I,J)
          SoilTranspFac(I,IndSoil,J) = SoilTranspFac(I,IndSoil,J) / SoilTranspFacAcc(I,J)
       enddo
    endif
    

      end do
    end do
    !$acc end parallel loop



    end associate

  end subroutine SoilWaterTranspiration

end module SoilWaterTranspirationMod
