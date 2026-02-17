module GroundThermalPropertyMod

!!! Compute snow and soil thermal conductivity and heat capacity (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowThermalPropertyMod, only : SnowThermalProperty
  use SoilThermalPropertyMod, only : SoilThermalProperty

  implicit none

contains

  subroutine GroundThermalProperty(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: THERMOPROP
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd    ! loop index
    integer                          :: I, J       ! grid indices

! --------------------------------------------------------------------
    ! initialize - done in parallel
        associate(                                                                           &
                  NumSoilLayer           => noahmp%config%domain%NumSoilLayer               ,& ! in,  number of soil layers
                  SurfaceType            => noahmp%config%domain%SurfaceType                ,& ! in,  surface type 1-soil; 2-lake
                  MainTimeStep           => noahmp%config%domain%MainTimeStep               ,& ! in,  main noahmp timestep [s]
                  ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer     ,& ! in,  thickness of snow/soil layers [m] (3D)
                  NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg       ,& ! in,  actual number of snow layers (negative)
                  FlagUrban              => noahmp%config%domain%FlagUrban             ,& ! in,  logical flag for urban grid
                  SnowDepth              => noahmp%water%state%SnowDepth               ,& ! in,  snow depth [m]
                  TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow         ,& ! in,  snow and soil layer temperature [K] (3D)
                  ThermConductSoilSnow   => noahmp%energy%state%ThermConductSoilSnow        ,& ! out, thermal conductivity [W/m/K] for all soil & snow (3D)
                  HeatCapacSoilSnow      => noahmp%energy%state%HeatCapacSoilSnow           ,& ! out, heat capacity [J/m3/K] for all soil & snow (3D)
                  PhaseChgFacSoilSnow    => noahmp%energy%state%PhaseChgFacSoilSnow         ,& ! out, energy factor for soil & snow phase change (3D)
                  HeatCapacVolSnow       => noahmp%energy%state%HeatCapacVolSnow            ,& ! out, snow layer volumetric specific heat [J/m3/K] (3D)
                  ThermConductSnow       => noahmp%energy%state%ThermConductSnow            ,& ! out, snow layer thermal conductivity [W/m/K] (3D)
                  HeatCapacVolSoil       => noahmp%energy%state%HeatCapacVolSoil            ,& ! out, soil layer volumetric specific heat [J/m3/K] (3D)
                  ThermConductSoil       => noahmp%energy%state%ThermConductSoil             & ! out, soil layer thermal conductivity [W/m/K] (3D)
                 )

    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
        noahmp%energy%state%HeatCapacSoilSnow(I,:,J)    = 0.0
        noahmp%energy%state%ThermConductSoilSnow(I,:,J) = 0.0
      end do
    end do
    !$acc end parallel loop

    ! compute snow thermal conductivity and heat capacity
    call SnowThermalProperty(noahmp)

    ! compute soil thermal properties
    call SoilThermalProperty(noahmp)

    ! combine snow and soil thermal properties, compute phase change factor
    !$acc parallel loop collapse(2) gang vector default(present) private(LoopInd)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


        ! copy snow thermal properties
        !$acc loop seq
        do LoopInd = NumSnowLayerNeg(I,J)+1, 0
           ThermConductSoilSnow(I,LoopInd,J) = ThermConductSnow(I,LoopInd,J)
           HeatCapacSoilSnow(I,LoopInd,J)    = HeatCapacVolSnow(I,LoopInd,J)
        enddo

        ! copy soil thermal properties
        !$acc loop seq
        do LoopInd = 1, NumSoilLayer
           ThermConductSoilSnow(I,LoopInd,J) = ThermConductSoil(I,LoopInd,J)
           HeatCapacSoilSnow(I,LoopInd,J)    = HeatCapacVolSoil(I,LoopInd,J)
        enddo

        ! urban override
        if ( FlagUrban(I,J) .eqv. .true. ) then
           !$acc loop seq
           do LoopInd = 1, NumSoilLayer
              ThermConductSoilSnow(I,LoopInd,J) = 3.24
           enddo
        endif

        ! heat flux reduction effect from the overlying green canopy, adapted from
        ! section 2.1.2 of Peters-Lidard et al. (1997, JGR, VOL 102(D4)).
        ! not in use because of the separation of the canopy layer from the ground.
        ! but this may represent the effects of leaf litter (Niu comments)
        ! ThermConductSoilSnow(I,1,J) = ThermConductSoilSnow(I,1,J) * EXP (SBETA * VegFracGreen(I,J))

        ! compute lake thermal properties (no consideration of turbulent mixing for this version)
        if ( SurfaceType(I,J) == 2 ) then
           !$acc loop seq
           do LoopInd = 1, NumSoilLayer
              if ( TemperatureSoilSnow(I,LoopInd,J) > ConstFreezePoint) then
                 HeatCapacSoilSnow(I,LoopInd,J)    = ConstHeatCapacWater
                 ThermConductSoilSnow(I,LoopInd,J) = ConstThermConductWater  !+ KEDDY * ConstHeatCapacWater
              else
                 HeatCapacSoilSnow(I,LoopInd,J)    = ConstHeatCapacIce
                 ThermConductSoilSnow(I,LoopInd,J) = ConstThermConductIce
              endif
           enddo
        endif

        ! combine a temporary variable used for melting/freezing of snow and frozen soil
        !$acc loop seq
        do LoopInd = NumSnowLayerNeg(I,J)+1, NumSoilLayer
           PhaseChgFacSoilSnow(I,LoopInd,J) = MainTimeStep / (HeatCapacSoilSnow(I,LoopInd,J) * ThicknessSnowSoilLayer(I,LoopInd,J))
        enddo

        ! snow/soil interface
        if ( NumSnowLayerNeg(I,J) == 0 ) then
           ThermConductSoilSnow(I,1,J) = (ThermConductSoilSnow(I,1,J)*ThicknessSnowSoilLayer(I,1,J) + 0.35*SnowDepth(I,J)) / &
                                     (SnowDepth(I,J) + ThicknessSnowSoilLayer(I,1,J))
        else
           ThermConductSoilSnow(I,1,J) = (ThermConductSoilSnow(I,1,J)*ThicknessSnowSoilLayer(I,1,J) + &
                                      ThermConductSoilSnow(I,0,J)*ThicknessSnowSoilLayer(I,0,J)) / &
                                     (ThicknessSnowSoilLayer(I,0,J) + ThicknessSnowSoilLayer(I,1,J))
        endif


      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine GroundThermalProperty

end module GroundThermalPropertyMod
