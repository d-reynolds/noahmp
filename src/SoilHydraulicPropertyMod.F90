module SoilHydraulicPropertyMod

!!! Two methods for calculating soil water diffusivity and soil hydraulic conductivity (2D GPU-optimized)
!!! Option 1: linear effects (more permeable, Niu and Yang,2006); Option 2: nonlinear effects (less permeable)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SoilDiffusivityConductivityOpt1(noahmp, SoilWatDiffusivity, SoilWatConductivity, &
                                             SoilMoisture, SoilImpervFrac, IndLayer, I, J)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: WDFCND1
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

! IN and OUT variables
    type(noahmp_type)     , intent(inout) :: noahmp
    integer               , intent(in)    :: I, J                   ! grid indices
    integer               , intent(in)    :: IndLayer               ! soil layer index
    real(kind=kind_noahmp), intent(in)    :: SoilMoisture           ! soil moisture [m3/m3]
    real(kind=kind_noahmp), intent(in)    :: SoilImpervFrac         ! impervious fraction due to frozen soil
    real(kind=kind_noahmp), intent(out)   :: SoilWatConductivity    ! soil water conductivity [m/s]
    real(kind=kind_noahmp), intent(out)   :: SoilWatDiffusivity     ! soil water diffusivity [m2/s]

! local variable
    real(kind=kind_noahmp)                :: SoilExpTmp             ! exponential local factor
    real(kind=kind_noahmp)                :: SoilPreFac             ! pre-factor
    !$acc routine seq

    SoilPreFac = max(0.01, SoilMoisture/noahmp%water%param%SoilMoistureSat(I,IndLayer,J))

    ! soil water diffusivity
    SoilExpTmp         = noahmp%water%param%SoilExpCoeffB(I,IndLayer,J) + 2.0
    SoilWatDiffusivity = noahmp%water%param%SoilWatDiffusivitySat(I,IndLayer,J) * SoilPreFac ** SoilExpTmp
    SoilWatDiffusivity = SoilWatDiffusivity * (1.0 - SoilImpervFrac)

    ! soil hydraulic conductivity
    SoilExpTmp          = 2.0 * noahmp%water%param%SoilExpCoeffB(I,IndLayer,J) + 3.0
    SoilWatConductivity = noahmp%water%param%SoilWatConductivitySat(I,IndLayer,J) * SoilPreFac ** SoilExpTmp
    SoilWatConductivity = SoilWatConductivity * (1.0 - SoilImpervFrac)

  end subroutine SoilDiffusivityConductivityOpt1


  subroutine SoilDiffusivityConductivityOpt2(noahmp, SoilWatDiffusivity, SoilWatConductivity, &
                                             SoilMoisture, SoilIce, IndLayer, I, J)

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: WDFCND2
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

! IN and OUT variables
    type(noahmp_type)     , intent(inout) :: noahmp
    integer               , intent(in)    :: I, J                  ! grid indices
    integer               , intent(in)    :: IndLayer              ! soil layer index
    real(kind=kind_noahmp), intent(in)    :: SoilMoisture          ! soil moisture [m3/m3]
    real(kind=kind_noahmp), intent(in)    :: SoilIce               ! soil ice content [m3/m3]
    real(kind=kind_noahmp), intent(out)   :: SoilWatConductivity   ! soil water conductivity [m/s]
    real(kind=kind_noahmp), intent(out)   :: SoilWatDiffusivity    ! soil water diffusivity [m2/s]

! local variable
    real(kind=kind_noahmp)                :: SoilExpTmp            ! exponential local factor
    real(kind=kind_noahmp)                :: SoilPreFac1           ! pre-factor
    real(kind=kind_noahmp)                :: SoilPreFac2           ! pre-factor
    real(kind=kind_noahmp)                :: SoilIceWgt            ! weights
    !$acc routine seq

    SoilPreFac1 = 0.05 / noahmp%water%param%SoilMoistureSat(I,IndLayer,J)
    SoilPreFac2 = max(0.01, SoilMoisture/noahmp%water%param%SoilMoistureSat(I,IndLayer,J))
    SoilPreFac1 = min(SoilPreFac1, SoilPreFac2)

    ! soil water diffusivity
    SoilExpTmp         = noahmp%water%param%SoilExpCoeffB(I,IndLayer,J) + 2.0
    SoilWatDiffusivity = noahmp%water%param%SoilWatDiffusivitySat(I,IndLayer,J) * SoilPreFac2 ** SoilExpTmp
    if ( SoilIce > 0.0 ) then
       SoilIceWgt         = 1.0 / (1.0 + (500.0 * SoilIce)**3.0)
       SoilWatDiffusivity = SoilIceWgt * SoilWatDiffusivity + &
                            (1.0-SoilIceWgt) * noahmp%water%param%SoilWatDiffusivitySat(I,IndLayer,J) * SoilPreFac1**SoilExpTmp
    endif

    ! soil hydraulic conductivity
    SoilExpTmp          = 2.0 * noahmp%water%param%SoilExpCoeffB(I,IndLayer,J) + 3.0
    SoilWatConductivity = noahmp%water%param%SoilWatConductivitySat(I,IndLayer,J) * SoilPreFac2 ** SoilExpTmp

  end subroutine SoilDiffusivityConductivityOpt2

end module SoilHydraulicPropertyMod
