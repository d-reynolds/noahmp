module GeneralInitMod

!!! General initialization for variables (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine GeneralInit(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in NOAHMP_SFLX)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J      ! grid indices
    integer                          :: LoopInd   ! loop index
! --------------------------------------------------------------------
        associate(                                                                       &
                  LandUseDataName        => noahmp%config%domain%LandUseDataName        ,& ! in,  landuse data name (USGS or MODIS_IGBP)
                  VegType                => noahmp%config%domain%VegType           ,& ! in,  vegetation type
                  NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,  number of soil layers
                  DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,  depth [m] of layer-bottom from soil surface
                  NumSoilLayerRoot       => noahmp%water%param%NumSoilLayerRoot    ,& ! in,  number of soil layers with root present
                  NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg   ,& ! in,  actual number of snow layers (negative)
                  FlagCropland           => noahmp%config%domain%FlagCropland      ,& ! out, flag to identify croplands
                  FlagWetland            => noahmp%config%domain%FlagWetland       ,& ! out, flag to identify wetlands
                  TemperatureRootZone    => noahmp%energy%state%TemperatureRootZone & ! out, root-zone averaged temperature [K]
                 )

    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

  if (noahmp%config%domain%IndicatorIceSfc(I,J) == 0) then

    ! initialize snow/soil layer thickness
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       if ( LoopInd == NumSnowLayerNeg(I,J)+1 ) then
          noahmp%config%domain%ThicknessSnowSoilLayer(I,LoopInd,J) = - noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd,J)
       else
          noahmp%config%domain%ThicknessSnowSoilLayer(I,LoopInd,J) = noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd-1,J) - &
                                                                       noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd,J)
       endif
    enddo

    ! initialize root-zone soil temperature
    TemperatureRootZone(I,J) = 0.0
    !$acc loop seq
    do LoopInd = 1, NumSoilLayerRoot(I,J)
       TemperatureRootZone(I,J) = TemperatureRootZone(I,J) + &
                             noahmp%energy%state%TemperatureSoilSnow(I,LoopInd,J) * &
                             noahmp%config%domain%ThicknessSnowSoilLayer(I,LoopInd,J) / (-DepthSoilLayer(I,NumSoilLayerRoot(I,J),J))
    enddo

    ! initialize special land type flags
    FlagCropland(I,J) = .false.
    FlagWetland(I,J)  = .false.
    if ( LandUseDataName == ConstLU_USGS) then
       if ( (VegType(I,J) >= 3 ) .and. (VegType(I,J) <= 6 ) ) FlagCropland(I,J) = .true.
       if ( (VegType(I,J) >= 17) .and. (VegType(I,J) <= 18) ) FlagWetland(I,J)  = .true.
    elseif ( LandUseDataName == ConstLU_IGBP_MODIS_NOAH) then
       if ( (VegType(I,J) == 12) .or. (VegType(I,J) == 14) )  FlagCropland(I,J) = .true.
       if ( (VegType(I,J) == 11) )                       FlagWetland(I,J)  = .true.
    endif

  else if (noahmp%config%domain%IndicatorIceSfc(I,J) == -1) then
    ! initialize snow/soil layer thickness for glaciers
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
       if ( LoopInd == (NumSnowLayerNeg(I,J)+1) ) then
          noahmp%config%domain%ThicknessSnowSoilLayer(I,LoopInd,J) = - noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd,J)
       else
          noahmp%config%domain%ThicknessSnowSoilLayer(I,LoopInd,J) = noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd-1,J) - noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd,J)
       endif
    enddo
  endif

      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine GeneralInit

end module GeneralInitMod
