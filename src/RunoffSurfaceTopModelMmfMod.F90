module RunoffSurfaceTopModelMmfMod

!!! Calculate surface runoff based on TOPMODEL with MMF groundwater scheme

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine RunoffSurfaceTopModelMMF(noahmp)

! ------------------------ Code history --------------------------------------------------
! Originally embeded in SOILWATER subroutine instead of as a separate subroutine
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    ! local variable
    integer                          :: I, J      ! grid indices

    associate(                                                           &
              SoilSfcInflowMean => noahmp%water%flux%SoilSfcInflowMean  ,& ! in,  mean water input on soil surface [m/s]
              RunoffDecayFac    => noahmp%water%param%RunoffDecayFac    ,& ! in,  runoff decay factor [1/m]
              SoilSfcSatFracMax => noahmp%water%param%SoilSfcSatFracMax ,& ! in,  maximum surface saturated fraction (global mean)
              SoilImpervFrac    => noahmp%water%state%SoilImpervFrac         ,& ! in,  impervious fraction due to frozen soil
              WaterTableDepth   => noahmp%water%state%WaterTableDepth   ,& ! in,  water table depth [m]
              SoilSaturateFrac  => noahmp%water%state%SoilSaturateFrac  ,& ! out, fractional saturated area for soil moisture
              RunoffSurface     => noahmp%water%flux%RunoffSurface      ,& ! out, surface runoff [m/s]
              InfilRateSfc      => noahmp%water%flux%InfilRateSfc        & ! out, infiltration rate at surface [m/s]
             )

    !$acc parallel loop collapse(2) gang vector default(present)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE 
         if ( noahmp%config%domain%IndicatorIceSfc(I,J) == -1 ) cycle  ! skip soil process for ice surface points

    ! set up key parameter
    RunoffDecayFac(I,J) = 6.0

    ! compute saturated area fraction
    SoilSaturateFrac(I,J) = SoilSfcSatFracMax(I,J) * exp(-0.5 * RunoffDecayFac(I,J) * max(-2.0-WaterTableDepth(I,J),0.0))

    ! compute surface runoff and infiltration  m/s
    if ( SoilSfcInflowMean(I,J) > 0.0 ) then
       RunoffSurface(I,J) = SoilSfcInflowMean(I,J) * ((1.0-SoilImpervFrac(I,1,J)) * SoilSaturateFrac(I,J) + SoilImpervFrac(I,1,J))
       InfilRateSfc(I,J)  = SoilSfcInflowMean(I,J) - RunoffSurface(I,J) 
    endif


      end do
    end do


    end associate

  end subroutine RunoffSurfaceTopModelMMF

end module RunoffSurfaceTopModelMmfMod
