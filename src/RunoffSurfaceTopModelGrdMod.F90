module RunoffSurfaceTopModelGrdMod

!!! Calculate surface runoff based on TOPMODEL with groundwater scheme (Niu et al., 2007)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine RunoffSurfaceTopModelGrd(noahmp)

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


    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
! --------------------------------------------------------------------
    associate(                                                                 &
              SoilSfcInflowMean => noahmp%water%flux%SoilSfcInflowMean(I,J)        ,& ! in,  mean water input on soil surface [m/s]
              RunoffDecayFac    => noahmp%water%param%RunoffDecayFac(I,J)          ,& ! in,  runoff decay factor [1/m]
              SoilSfcSatFracMax => noahmp%water%param%SoilSfcSatFracMax(I,J)       ,& ! in,  maximum surface saturated fraction (global mean)
              SoilExpCoeffB     => noahmp%water%param%SoilExpCoeffB                ,& ! in,  soil B parameter
              SoilImpervFrac    => noahmp%water%state%SoilImpervFrac               ,& ! in,  impervious fraction due to frozen soil
              WaterTableDepth   => noahmp%water%state%WaterTableDepth(I,J)         ,& ! in,  water table depth [m]
              SoilSaturateFrac  => noahmp%water%state%SoilSaturateFrac(I,J)        ,& ! out, fractional saturated area for soil moisture
              RunoffSurface     => noahmp%water%flux%RunoffSurface(I,J)            ,& ! out, surface runoff [m/s]
              InfilRateSfc      => noahmp%water%flux%InfilRateSfc(I,J)              & ! out, infiltration rate at surface [m/s]
             )
! ----------------------------------------------------------------------

    ! set up key parameter
    !RunoffDecayFac = 6.0
    RunoffDecayFac = SoilExpCoeffB(I,1,J) / 3.0 ! calibratable, GY Niu's update 2022

    ! compute saturated area fraction
    !SoilSaturateFrac = SoilSfcSatFracMax * exp(-0.5 * RunoffDecayFac * (WaterTableDepth-2.0))
    SoilSaturateFrac = SoilSfcSatFracMax * exp(-0.5 * RunoffDecayFac * WaterTableDepth) ! GY Niu's update 2022

    ! compute surface runoff and infiltration  m/s
    if ( SoilSfcInflowMean > 0.0 ) then
       RunoffSurface = SoilSfcInflowMean * ((1.0-SoilImpervFrac(I,1,J)) * SoilSaturateFrac + SoilImpervFrac(I,1,J))
       InfilRateSfc  = SoilSfcInflowMean - RunoffSurface 
    endif

    end associate

      end do
    end do
    !$acc end parallel loop


  end subroutine RunoffSurfaceTopModelGrd

end module RunoffSurfaceTopModelGrdMod
