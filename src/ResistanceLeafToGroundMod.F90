module ResistanceLeafToGroundMod

!!! Compute under-canopy aerodynamic resistance and leaf boundary layer resistance (2D GPU-optimized)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine ResistanceLeafToGround(noahmp, IndIter, VegAreaIndEff, HeatSenGrdTmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: RAGRB
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    integer               , intent(in   ) :: IndIter         ! iteration index
    real(kind=kind_noahmp), allocatable, intent(in   ) :: HeatSenGrdTmp(:,:)   ! temporary ground sensible heat flux (w/m2) in each iteration
    real(kind=kind_noahmp), allocatable, intent(in   ) :: VegAreaIndEff(:,:)   ! temporary effective vegetation area index with constraint (<=6.0)
    type(noahmp_type)     , intent(inout) :: noahmp

! local variable
    integer                               :: I, J            ! grid indices
    real(kind=kind_noahmp)                :: MPE             ! prevents overflow for division by zero
    real(kind=kind_noahmp)                :: KH              ! turbulent transfer coefficient, sensible heat, (m2/s)
    real(kind=kind_noahmp)                :: TMP1            ! temporary calculation
    real(kind=kind_noahmp)                :: TMP2            ! temporary calculation
    real(kind=kind_noahmp)                :: TMPRAH2         ! temporary calculation for aerodynamic resistances
    real(kind=kind_noahmp)                :: TMPRB           ! temporary calculation for rb
    real(kind=kind_noahmp)                :: FHGNEW          ! temporary vars

! --------------------------------------------------------------------
    !$acc parallel loop collapse(2) gang vector present(noahmp, HeatSenGrdTmp, VegAreaIndEff) &
    !$acc private(MPE, KH, TMP1, TMP2, TMPRAH2, TMPRB, FHGNEW)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        if ( .not. ((noahmp%energy%state%VegAreaIndEff(I,J) > 0.0 ) .and. (noahmp%energy%state%VegFrac(I,J) > 0)) ) cycle ! skip non-vegetated surface

        associate(                                                                      &
                  LeafDimLength          => noahmp%energy%param%LeafDimLength(I,J)          ,& ! in,    characteristic leaf dimension [m]
                  CanopyWindExtFac       => noahmp%energy%param%CanopyWindExtFac(I,J)       ,& ! in,    canopy wind extinction parameter
                  DensityAirRefHeight    => noahmp%energy%state%DensityAirRefHeight(I,J) ,& ! in,    density air [kg/m3]
                  TemperatureCanopyAir   => noahmp%energy%state%TemperatureCanopyAir(I,J) ,& ! in,    canopy air temperature [K]
                  ZeroPlaneDispSfc       => noahmp%energy%state%ZeroPlaneDispSfc(I,J)  ,& ! in,    zero plane displacement [m]
                  RoughLenMomGrd         => noahmp%energy%state%RoughLenMomGrd(I,J)    ,& ! in,    roughness length [m], momentum, ground
                  CanopyHeight           => noahmp%energy%state%CanopyHeight(I,J)      ,& ! in,    canopy height [m]
                  WindSpdCanopyTop       => noahmp%energy%state%WindSpdCanopyTop(I,J)  ,& ! in,    wind speed at top of canopy [m/s]
                  RoughLenShCanopy       => noahmp%energy%state%RoughLenShCanopy(I,J)  ,& ! in,    roughness length [m], sensible heat, canopy
                  RoughLenShVegGrd       => noahmp%energy%state%RoughLenShVegGrd(I,J)  ,& ! in,    roughness length [m], sensible heat ground, below canopy
                  FrictionVelVeg         => noahmp%energy%state%FrictionVelVeg(I,J)    ,& ! in,    friction velocity [m/s], vegetated
                  MoStabCorrShUndCan     => noahmp%energy%state%MoStabCorrShUndCan(I,J) ,& ! inout, stability correction ground, below canopy
                  WindExtCoeffCanopy     => noahmp%energy%state%WindExtCoeffCanopy(I,J) ,& ! out,   canopy wind extinction coefficient
                  MoStabParaUndCan       => noahmp%energy%state%MoStabParaUndCan(I,J)  ,& ! out,   Monin-Obukhov stability parameter ground, below canopy
                  MoLengthUndCan         => noahmp%energy%state%MoLengthUndCan(I,J)    ,& ! out,   Monin-Obukhov length [m], ground, below canopy
                  ResistanceMomUndCan    => noahmp%energy%state%ResistanceMomUndCan(I,J) ,& ! out,   ground aerodynamic resistance for momentum [s/m]
                  ResistanceShUndCan     => noahmp%energy%state%ResistanceShUndCan(I,J) ,& ! out,   ground aerodynamic resistance for sensible heat [s/m]
                  ResistanceLhUndCan     => noahmp%energy%state%ResistanceLhUndCan(I,J) ,& ! out,   ground aerodynamic resistance for water vapor [s/m]
                  ResistanceLeafBoundary => noahmp%energy%state%ResistanceLeafBoundary(I,J) & ! out,   bulk leaf boundary layer resistance [s/m]
                 )
! ----------------------------------------------------------------------

        ! initialization
        MPE              = 1.0e-6
        MoStabParaUndCan = 0.0
        MoLengthUndCan   = 0.0

        ! stability correction to below canopy resistance
        if ( IndIter > 1 ) then
           TMP1 = ConstVonKarman * (ConstGravityAcc / TemperatureCanopyAir) * HeatSenGrdTmp(I,J) / &
                  (DensityAirRefHeight * ConstHeatCapacAir)
           if ( abs(TMP1) <= MPE ) TMP1 = MPE
           MoLengthUndCan   = -1.0 * FrictionVelVeg**3 / TMP1
           MoStabParaUndCan = min((ZeroPlaneDispSfc-RoughLenMomGrd)/MoLengthUndCan, 1.0)
        endif
        if ( MoStabParaUndCan < 0.0 ) then
           FHGNEW = (1.0 - 15.0 * MoStabParaUndCan)**(-0.25)
        else
           FHGNEW = 1.0 + 4.7 * MoStabParaUndCan
        endif
        if ( IndIter == 1 ) then
           MoStabCorrShUndCan = FHGNEW
        else
           MoStabCorrShUndCan = 0.5 * (MoStabCorrShUndCan + FHGNEW)
        endif

        ! wind attenuation within canopy
        WindExtCoeffCanopy = (CanopyWindExtFac * VegAreaIndEff(I,J) * CanopyHeight * MoStabCorrShUndCan)**0.5
        TMP1               = exp(-WindExtCoeffCanopy * RoughLenShVegGrd / CanopyHeight)
        TMP2               = exp(-WindExtCoeffCanopy * (RoughLenShCanopy + ZeroPlaneDispSfc) / CanopyHeight)
        TMPRAH2            = CanopyHeight * exp(WindExtCoeffCanopy) / WindExtCoeffCanopy * (TMP1-TMP2)

        ! aerodynamic resistances raw and rah between heights ZeroPlaneDisp+RoughLenShVegGrd and RoughLenShVegGrd.
        KH                  = max(ConstVonKarman*FrictionVelVeg*(CanopyHeight-ZeroPlaneDispSfc), MPE)
        ResistanceMomUndCan = 0.0
        ResistanceShUndCan  = TMPRAH2 / KH
        ResistanceLhUndCan  = ResistanceShUndCan

        ! leaf boundary layer resistance
        TMPRB                  = WindExtCoeffCanopy * 50.0 / (1.0 - exp(-WindExtCoeffCanopy/2.0))
        ResistanceLeafBoundary = TMPRB * sqrt(LeafDimLength / WindSpdCanopyTop)
        ResistanceLeafBoundary = min(max(ResistanceLeafBoundary, 5.0), 50.0)      ! limit ResistanceLeafBoundary to 5-50, typically <50

        end associate

      enddo
    enddo

  end subroutine ResistanceLeafToGround

end module ResistanceLeafToGroundMod
