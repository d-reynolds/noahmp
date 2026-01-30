module SoilWaterDiffusionRichardsMod

!!! Solve Richards equation for soil water movement/diffusion
!!! Compute the right hand side of the time tendency term of the soil
!!! water diffusion equation.  also to compute (prepare) the matrix
!!! coefficients for the tri-diagonal matrix of the implicit time scheme.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SoilHydraulicPropertyMod

  implicit none

contains

  subroutine SoilWaterDiffusionRichards(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: SRT
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)     , intent(inout) :: noahmp
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatRight     ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft1     ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft2     ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:,:,:), intent(inout) :: MatLeft3     ! left-hand side term of the matrix

! local variable
    integer                                           :: LoopInd                     ! loop index
    real(kind=kind_noahmp)                            :: DepthSnowSoilTmp            ! temporary snow/soil layer depth [m]
    real(kind=kind_noahmp)                            :: SoilMoistTmpToWT            ! temporary soil moisture between bottom of the soil and water table
    real(kind=kind_noahmp)                            :: SoilMoistBotTmp             ! temporary soil moisture below bottom to calculate flux
    real(kind=kind_noahmp)                            :: DepthSnowSoilInv(1:noahmp%config%domain%NumSoilLayer)            ! inverse of snow/soil layer depth [1/m]
    real(kind=kind_noahmp)                            :: SoilThickTmp(1:noahmp%config%domain%NumSoilLayer)                ! temporary soil thickness
    real(kind=kind_noahmp)                            :: SoilWaterGrad(1:noahmp%config%domain%NumSoilLayer)               ! temporary soil moisture vertical gradient
    real(kind=kind_noahmp)                            :: WaterExcess(1:noahmp%config%domain%NumSoilLayer)                 ! temporary excess water flux
    real(kind=kind_noahmp)                            :: SoilMoistureTmp(1:noahmp%config%domain%NumSoilLayer)             ! temporary soil moisture
    integer                                           :: I, J                        ! grid indices
    !$acc parallel loop collapse(2) gang vector present(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight) private(DepthSnowSoilTmp, SoilMoistTmpToWT, SoilMoistBotTmp, DepthSnowSoilInv, SoilThickTmp, SoilWaterGrad, WaterExcess, SoilMoistureTmp, LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
! --------------------------------------------------------------------
    associate(                                                                             &
              NumSoilLayer              => noahmp%config%domain%NumSoilLayer              ,& ! in,  number of soil layers
              DepthSoilLayer            => noahmp%config%domain%DepthSoilLayer            ,& ! in,  depth [m] of layer-bottom from soil surface
              OptSoilPermeabilityFrozen => noahmp%config%nmlist%OptSoilPermeabilityFrozen ,& ! in,  options for frozen soil permeability
              OptRunoffSubsurface       => noahmp%config%nmlist%OptRunoffSubsurface       ,& ! in,  options for drainage and subsurface runoff
              SoilDrainSlope            => noahmp%water%param%SoilDrainSlope(I,J)              ,& ! in,  slope index for soil drainage
              InfilRateSfc              => noahmp%water%flux%InfilRateSfc(I,J)                 ,& ! in,  infiltration rate at surface [m/s]
              EvapSoilSfcLiqMean        => noahmp%water%flux%EvapSoilSfcLiqMean(I,J)           ,& ! in,  mean evaporation from soil surface [m/s]
              TranspWatLossSoilMean     => noahmp%water%flux%TranspWatLossSoilMean        ,& ! in,  mean transpiration water loss from soil layers [m/s]
              SoilLiqWater              => noahmp%water%state%SoilLiqWater                ,& ! in,  soil water content [m3/m3]
              SoilMoisture              => noahmp%water%state%SoilMoisture                ,& ! in,  total soil moisture [m3/m3]
              WaterTableDepth           => noahmp%water%state%WaterTableDepth(I,J)             ,& ! in,  water table depth [m]
              SoilImpervFrac            => noahmp%water%state%SoilImpervFrac              ,& ! in,  fraction of imperviousness due to frozen soil
              SoilImpervFracMax         => noahmp%water%state%SoilImpervFracMax(I,J)           ,& ! in,  maximum soil imperviousness fraction
              SoilIceMax                => noahmp%water%state%SoilIceMax(I,J)                  ,& ! in,  maximum soil ice content [m3/m3]
              SoilMoistureToWT          => noahmp%water%state%SoilMoistureToWT(I,J)            ,& ! in,  soil moisture between bottom of the soil and the water table
              SoilWatConductivity       => noahmp%water%state%SoilWatConductivity         ,& ! out, soil hydraulic conductivity [m/s]
              SoilWatDiffusivity        => noahmp%water%state%SoilWatDiffusivity          ,& ! out, soil water diffusivity [m2/s]
              DrainSoilBot              => noahmp%water%flux%DrainSoilBot(I,J)                  & ! out, soil bottom drainage [m/s]
             )
! ----------------------------------------------------------------------

    ! initialization
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
      MatRight(I,LoopInd,J)     = 0.0
      MatLeft1(I,LoopInd,J)     = 0.0
      MatLeft2(I,LoopInd,J)     = 0.0
      MatLeft3(I,LoopInd,J)     = 0.0
      DepthSnowSoilInv(LoopInd) = 0.0
      SoilThickTmp(LoopInd)     = 0.0
      SoilWaterGrad(LoopInd)    = 0.0
      WaterExcess(LoopInd)      = 0.0
      SoilMoistureTmp(LoopInd)  = 0.0
    enddo
    ! compute soil hydraulic conductivity and diffusivity
    if ( OptSoilPermeabilityFrozen == 1 ) then
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          call SoilDiffusivityConductivityOpt1(noahmp,SoilWatDiffusivity(I,LoopInd,J),SoilWatConductivity(I,LoopInd,J),&
                                               SoilMoisture(I,LoopInd,J),SoilImpervFrac(I,LoopInd,J),LoopInd,I,J) 
          SoilMoistureTmp(LoopInd) = SoilMoisture(I,LoopInd,J)
       enddo
       if ( OptRunoffSubsurface == 5 ) SoilMoistTmpToWT = SoilMoistureToWT
    endif

    if ( OptSoilPermeabilityFrozen == 2 ) then
       !$acc loop seq
       do LoopInd = 1, NumSoilLayer
          call SoilDiffusivityConductivityOpt2(noahmp,SoilWatDiffusivity(I,LoopInd,J),SoilWatConductivity(I,LoopInd,J),&
                                               SoilLiqWater(I,LoopInd,J),SoilIceMax,LoopInd, I, J)
          SoilMoistureTmp(LoopInd) = SoilLiqWater(I,LoopInd,J)
       enddo
       if ( OptRunoffSubsurface == 5 ) &
          SoilMoistTmpToWT = SoilMoistureToWT * SoilLiqWater(I,NumSoilLayer,J) / SoilMoisture(I,NumSoilLayer,J)  !same liquid fraction as in the bottom layer
    endif

    ! compute gradient and flux of soil water diffusion terms
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       if ( LoopInd == 1 ) then
          SoilThickTmp(LoopInd)     = - DepthSoilLayer(I,LoopInd,J)
          DepthSnowSoilTmp          = - DepthSoilLayer(I,LoopInd+1,J)
          DepthSnowSoilInv(LoopInd) = 2.0 / DepthSnowSoilTmp
          SoilWaterGrad(LoopInd)    = 2.0 * (SoilMoistureTmp(LoopInd)-SoilMoistureTmp(LoopInd+1)) / DepthSnowSoilTmp
          WaterExcess(LoopInd)      = SoilWatDiffusivity(I,LoopInd,J)*SoilWaterGrad(LoopInd) + SoilWatConductivity(I,LoopInd,J) - &
                                      InfilRateSfc + TranspWatLossSoilMean(I,LoopInd,J) + EvapSoilSfcLiqMean
       else if ( LoopInd < NumSoilLayer ) then
          SoilThickTmp(LoopInd)     = (DepthSoilLayer(I,LoopInd-1,J) - DepthSoilLayer(I,LoopInd,J))
          DepthSnowSoilTmp          = (DepthSoilLayer(I,LoopInd-1,J) - DepthSoilLayer(I,LoopInd+1,J))
          DepthSnowSoilInv(LoopInd) = 2.0 / DepthSnowSoilTmp
          SoilWaterGrad(LoopInd)    = 2.0 * (SoilMoistureTmp(LoopInd) - SoilMoistureTmp(LoopInd+1)) / DepthSnowSoilTmp
          WaterExcess(LoopInd)      = SoilWatDiffusivity(I,LoopInd,J)*SoilWaterGrad(LoopInd) + SoilWatConductivity(I,LoopInd,J) - &
                                      SoilWatDiffusivity(I,LoopInd-1,J)*SoilWaterGrad(LoopInd-1) - SoilWatConductivity(I,LoopInd-1,J) + &
                                      TranspWatLossSoilMean(I,LoopInd,J)
       else
          SoilThickTmp(LoopInd) = (DepthSoilLayer(I,LoopInd-1,J) - DepthSoilLayer(I,LoopInd,J))
          if ( (OptRunoffSubsurface == 1) .or. (OptRunoffSubsurface == 2) ) then
             DrainSoilBot = 0.0
          endif
          if ( (OptRunoffSubsurface == 3) .or. (OptRunoffSubsurface == 6) .or. &
               (OptRunoffSubsurface == 7) .or. (OptRunoffSubsurface == 8) ) then
             DrainSoilBot = SoilDrainSlope * SoilWatConductivity(I,LoopInd,J)
          endif
          if ( OptRunoffSubsurface == 4 ) then
             DrainSoilBot = (1.0 - SoilImpervFracMax) * SoilWatConductivity(I,LoopInd,J)
          endif
          if ( OptRunoffSubsurface == 5 ) then   ! gmm new m-m&f water table dynamics formulation
             DepthSnowSoilTmp  = 2.0 * SoilThickTmp(LoopInd)
             if ( WaterTableDepth < (DepthSoilLayer(I,NumSoilLayer,J)-SoilThickTmp(NumSoilLayer)) ) then
                ! gmm interpolate from below, midway to the water table, 
                ! to the middle of the auxiliary layer below the soil bottom
                SoilMoistBotTmp = SoilMoistureTmp(LoopInd) - (SoilMoistureTmp(LoopInd)-SoilMoistTmpToWT) * &
                                  SoilThickTmp(LoopInd)*2.0 / (SoilThickTmp(LoopInd)+DepthSoilLayer(I,LoopInd,J)-WaterTableDepth)
             else
                SoilMoistBotTmp = SoilMoistTmpToWT
             endif
             SoilWaterGrad(LoopInd) = 2.0 * (SoilMoistureTmp(LoopInd) - SoilMoistBotTmp) / DepthSnowSoilTmp
             DrainSoilBot           = SoilWatDiffusivity(I,LoopInd,J) * SoilWaterGrad(LoopInd) + SoilWatConductivity(I,LoopInd,J)
          endif
          WaterExcess(LoopInd) = -(SoilWatDiffusivity(I,LoopInd-1,J)*SoilWaterGrad(LoopInd-1)) - SoilWatConductivity(I,LoopInd-1,J) + &
                                 TranspWatLossSoilMean(I,LoopInd,J) + DrainSoilBot
       endif
    enddo

    ! prepare the matrix coefficients for the tri-diagonal matrix
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
       if ( LoopInd == 1 ) then
          MatLeft1(I,LoopInd,J) =   0.0
          MatLeft2(I,LoopInd,J) =   SoilWatDiffusivity(I,LoopInd  ,J) * DepthSnowSoilInv(LoopInd  ) / SoilThickTmp(LoopInd)
          MatLeft3(I,LoopInd,J) = - MatLeft2(I,LoopInd,J)
       else if ( LoopInd < NumSoilLayer ) then
          MatLeft1(I,LoopInd,J) = - SoilWatDiffusivity(I,LoopInd-1,J) * DepthSnowSoilInv(LoopInd-1) / SoilThickTmp(LoopInd)
          MatLeft3(I,LoopInd,J) = - SoilWatDiffusivity(I,LoopInd  ,J) * DepthSnowSoilInv(LoopInd  ) / SoilThickTmp(LoopInd)
          MatLeft2(I,LoopInd,J) = - (MatLeft1(I,LoopInd,J) + MatLeft3(I,LoopInd,J))
       else
          MatLeft1(I,LoopInd,J) = - SoilWatDiffusivity(I,LoopInd-1,J) * DepthSnowSoilInv(LoopInd-1) / SoilThickTmp(LoopInd)
          MatLeft3(I,LoopInd,J) =   0.0
          MatLeft2(I,LoopInd,J) = - (MatLeft1(I,LoopInd,J) + MatLeft3(I,LoopInd,J))
       endif
       MatRight(I,LoopInd,J) = WaterExcess(LoopInd) / (-SoilThickTmp(LoopInd))
    enddo

    end associate
   end do
   end do

  end subroutine SoilWaterDiffusionRichards

end module SoilWaterDiffusionRichardsMod
