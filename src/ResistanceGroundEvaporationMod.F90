module ResistanceGroundEvaporationMod

!!! Compute soil surface resistance to ground evaporation/sublimation (2D GPU-optimized)
!!! It represents the resistance imposed by the molecular diffusion in soil
!!! surface (as opposed to aerodynamic resistance computed elsewhere in the model)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine ResistanceGroundEvaporation(noahmp)

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
    integer                          :: I, J                   ! grid indices
    real(kind=kind_noahmp)           :: SoilEvapFac            ! soil water evaporation factor (0- 1)
    real(kind=kind_noahmp)           :: DrySoilThickness       ! Dry-layer thickness [m] for computing RSURF (Sakaguchi and Zeng, 2009)
    real(kind=kind_noahmp)           :: VapDiffuseRed          ! Reduced vapor diffusivity [m2/s] in soil for computing RSURF (SZ09)
    real(kind=kind_noahmp)           :: SoilMatPotentialSfc    ! surface layer soil matric potential [m]

! --------------------------------------------------------------------
        associate(                                                                         &
                  SurfaceType             => noahmp%config%domain%SurfaceType             ,& ! in,  surface type 1-soil; 2-lake
                  DepthSoilLayer          => noahmp%config%domain%DepthSoilLayer          ,& ! in,  depth [m] of layer-bottom from soil surface (2D)
                  FlagUrban               => noahmp%config%domain%FlagUrban          ,& ! in,  logical flag for urban grid
                  OptGroundResistanceEvap => noahmp%config%nmlist%OptGroundResistanceEvap ,& ! in,  options for ground resistance to evaporation/sublimation
                  ResistanceSoilExp       => noahmp%energy%param%ResistanceSoilExp        ,& ! in,  exponent in the shape parameter for soil resistance
                  ResistanceSnowSfc       => noahmp%energy%param%ResistanceSnowSfc        ,& ! in,  surface resistance for snow [s/m]
                  SoilMoistureSat         => noahmp%water%param%SoilMoistureSat           ,& ! in,  saturated value of soil moisture [m3/m3] (3D)
                  SoilMoistureWilt        => noahmp%water%param%SoilMoistureWilt          ,& ! in,  wilting point soil moisture [m3/m3] (3D)
                  SoilExpCoeffB           => noahmp%water%param%SoilExpCoeffB             ,& ! in,  soil B parameter (3D)
                  SoilMatPotentialSat     => noahmp%water%param%SoilMatPotentialSat       ,& ! in,  saturated soil matric potential [m] (3D)
                  SoilLiqWater            => noahmp%water%state%SoilLiqWater              ,& ! in,  soil water content [m3/m3] (3D)
                  SnowCoverFrac           => noahmp%water%state%SnowCoverFrac        ,& ! in,  snow cover fraction
                  SnowDepth               => noahmp%water%state%SnowDepth            ,& ! in,  snow depth [m]
                  TemperatureGrd          => noahmp%energy%state%TemperatureGrd      ,& ! in,  ground temperature [K]
                  ResistanceGrdEvap       => noahmp%energy%state%ResistanceGrdEvap   ,& ! out, ground surface resistance [s/m] to evaporation
                  RelHumidityGrd          => noahmp%energy%state%RelHumidityGrd       & ! out, raltive humidity in surface soil/snow air space
                 )

    !$acc parallel loop collapse(2) gang vector default(present) private(DrySoilThickness, SoilEvapFac, &
    !$acc SoilMatPotentialSfc, VapDiffuseRed)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


        ! initialization
        SoilEvapFac = max(0.0, SoilLiqWater(I,1,J)/SoilMoistureSat(I,1,J))

        if ( SurfaceType(I,J) == 2 ) then         ! lake point
           ResistanceGrdEvap(I,J) = 1.0           ! avoid being divided by 0
           RelHumidityGrd(I,J)    = 1.0
        else                                 ! soil point
           ! Sakaguchi and Zeng, 2009
           if ( (OptGroundResistanceEvap == 1) .or. (OptGroundResistanceEvap == 4) ) then
              DrySoilThickness  = (-DepthSoilLayer(I,1,J)) * (exp((1.0 - min(1.0,SoilLiqWater(I,1,J)/SoilMoistureSat(I,1,J))) ** &
                                                          ResistanceSoilExp(I,J)) - 1.0) / (2.71828-1.0)
              VapDiffuseRed     = 2.2e-5 * SoilMoistureSat(I,1,J) * SoilMoistureSat(I,1,J) * &
                                  (1.0 - SoilMoistureWilt(I,1,J)/SoilMoistureSat(I,1,J)) ** (2.0 + 3.0/SoilExpCoeffB(I,1,J))
              ResistanceGrdEvap(I,J) = DrySoilThickness / VapDiffuseRed

           ! Sellers (1992) original
           elseif ( OptGroundResistanceEvap == 2 ) then
              ResistanceGrdEvap(I,J) = SnowCoverFrac(I,J) * 1.0 + (1.0 - SnowCoverFrac(I,J)) * exp(8.25 - 4.225*SoilEvapFac)

           ! Sellers (1992) adjusted to decrease ResistanceGrdEvap for wet soil
           elseif ( OptGroundResistanceEvap == 3 ) then
              ResistanceGrdEvap(I,J) = SnowCoverFrac(I,J) * 1.0 + (1.0 - SnowCoverFrac(I,J)) * exp(8.25 - 6.0*SoilEvapFac)
           endif

           ! SnowCoverFrac weighted; snow ResistanceGrdEvap set in MPTABLE v3.8
           if ( OptGroundResistanceEvap == 4 ) then
              ResistanceGrdEvap(I,J) = 1.0 / (SnowCoverFrac(I,J) * (1.0/ResistanceSnowSfc(I,J)) + &
                                         (1.0-SnowCoverFrac(I,J)) * (1.0/max(ResistanceGrdEvap(I,J),0.001)))
           endif
           if ( (SoilLiqWater(I,1,J) < 0.01) .and. (SnowDepth(I,J) == 0.0) ) ResistanceGrdEvap(I,J) = 1.0e6

           SoilMatPotentialSfc = -SoilMatPotentialSat(I,1,J) * &
                                 (max(0.01,SoilLiqWater(I,1,J)) / SoilMoistureSat(I,1,J)) ** (-SoilExpCoeffB(I,1,J))
           RelHumidityGrd(I,J)      = SnowCoverFrac(I,J) + &
                                 (1.0-SnowCoverFrac(I,J)) * exp(SoilMatPotentialSfc*ConstGravityAcc/(ConstGasWaterVapor*TemperatureGrd(I,J)))
        endif

        ! urban
        if ( (FlagUrban(I,J) .eqv. .true.) .and. (SnowDepth(I,J) == 0.0) ) then
           ResistanceGrdEvap(I,J) = 1.0e6
        endif

        if (noahmp%config%domain%IndicatorIceSfc(I,J) == -1) then
            ResistanceGrdEvap(I,J) = 1.0
            RelHumidityGrd(I,J)    = 1.0
        endif


      end do
    end do
    !$acc end parallel loop


        end associate

  end subroutine ResistanceGroundEvaporation

end module ResistanceGroundEvaporationMod
