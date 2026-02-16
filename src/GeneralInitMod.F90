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
    !$acc parallel loop collapse(2) gang vector present(noahmp) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        associate(                                                                       &
                  LandUseDataName        => noahmp%config%domain%LandUseDataName        ,& ! in,  landuse data name (USGS or MODIS_IGBP)
                  VegType                => noahmp%config%domain%VegType(I,J)           ,& ! in,  vegetation type
                  NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,  number of soil layers
                  DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,  depth [m] of layer-bottom from soil surface
                  NumSoilLayerRoot       => noahmp%water%param%NumSoilLayerRoot(I,J)    ,& ! in,  number of soil layers with root present
                  NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg(I,J)   ,& ! in,  actual number of snow layers (negative)
                  FlagCropland           => noahmp%config%domain%FlagCropland(I,J)      ,& ! out, flag to identify croplands
                  FlagWetland            => noahmp%config%domain%FlagWetland(I,J)       ,& ! out, flag to identify wetlands
                  TemperatureRootZone    => noahmp%energy%state%TemperatureRootZone(I,J) & ! out, root-zone averaged temperature [K]
                 )
! ----------------------------------------------------------------------
  if (noahmp%config%domain%IndicatorIceSfc(I,J) == 0) then

    ! initialize snow/soil layer thickness
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       if ( LoopInd == NumSnowLayerNeg+1 ) then
          noahmp%config%domain%ThicknessSnowSoilLayer(I,LoopInd,J) = - noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd,J)
       else
          noahmp%config%domain%ThicknessSnowSoilLayer(I,LoopInd,J) = noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd-1,J) - &
                                                                       noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd,J)
       endif
    enddo

    ! initialize root-zone soil temperature
    TemperatureRootZone = 0.0
    !$acc loop seq
    do LoopInd = 1, NumSoilLayerRoot
       TemperatureRootZone = TemperatureRootZone + &
                             noahmp%energy%state%TemperatureSoilSnow(I,LoopInd,J) * &
                             noahmp%config%domain%ThicknessSnowSoilLayer(I,LoopInd,J) / (-DepthSoilLayer(I,NumSoilLayerRoot,J))
    enddo

    ! initialize special land type flags
    FlagCropland = .false.
    FlagWetland  = .false.
    if ( LandUseDataName == ConstLU_USGS) then
       if ( (VegType >= 3 ) .and. (VegType <= 6 ) ) FlagCropland = .true.
       if ( (VegType >= 17) .and. (VegType <= 18) ) FlagWetland  = .true.
    elseif ( LandUseDataName == ConstLU_IGBP_MODIS_NOAH) then
       if ( (VegType == 12) .or. (VegType == 14) )  FlagCropland = .true.
       if ( (VegType == 11) )                       FlagWetland  = .true.
    endif

  else if (noahmp%config%domain%IndicatorIceSfc(I,J) == -1) then
    ! initialize snow/soil layer thickness for glaciers
    !$acc loop seq
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       if ( LoopInd == (NumSnowLayerNeg+1) ) then
          noahmp%config%domain%ThicknessSnowSoilLayer(I,LoopInd,J) = - noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd,J)
       else
          noahmp%config%domain%ThicknessSnowSoilLayer(I,LoopInd,J) = noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd-1,J) - noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd,J)
       endif
    enddo
  endif
        end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine GeneralInit

end module GeneralInitMod
