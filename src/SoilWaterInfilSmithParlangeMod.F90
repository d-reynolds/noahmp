module SoilWaterInfilSmithParlangeMod

!!! Compute soil surface infiltration rate based on Smith-Parlange equation
!!! Reference: Smith, R.E. (2002), Infiltration Theory for Hydrologic Applications

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SoilHydraulicPropertyMod, only : SoilDiffusivityConductivityOpt2

  implicit none

contains

  subroutine SoilWaterInfilSmithParlange(noahmp, IndInfilMax, InfilSfcAcc, InfilSfcTmp, I, J)
  !$acc routine seq
! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: SMITH_PARLANGE_INFIL
! Original code: Prasanth Valayamkunnath <prasanth@ucar.edu>
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variabls
    type(noahmp_type)     , intent(inout) :: noahmp
    integer               , intent(in)    :: IndInfilMax           ! check for maximum infiltration at SoilMoistureWilt 
    real(kind=kind_noahmp), intent(inout) :: InfilSfcAcc           ! accumulated infiltration rate [m/s]
    real(kind=kind_noahmp), intent(out)   :: InfilSfcTmp           ! surface infiltration rate [m/s]
    integer               , intent(in)    :: I, J                  ! grid indices
! local variables
    integer                               :: IndSoil               ! soil layer index
    real(kind=kind_noahmp)                :: SoilWatDiffusivity    ! soil water diffusivity [m2/s]
    real(kind=kind_noahmp)                :: SoilWatConductivity   ! soil water conductivity [m/s]
    real(kind=kind_noahmp)                :: InfilFacTmp           ! temporary infiltrability variable
    real(kind=kind_noahmp)                :: WeighFac              ! smith-parlang weighing parameter
    real(kind=kind_noahmp)                :: IniSoilIce            ! zero soil ice

! ----------------------------------------------------------------------

    ! smith-parlang weighing parameter, Gamma
    WeighFac   = 0.82
    IniSoilIce = 0.0
    IndSoil    = 1

    ! check whether we are estimating infiltration for current SoilMoisture or SoilMoistureWilt
    if ( IndInfilMax == 1 ) then ! not active for now as the maximum infiltration is estimated based on table values

       ! estimate initial soil hydraulic conductivty (Ki in the equation) (m/s)
       call SoilDiffusivityConductivityOpt2(noahmp, SoilWatDiffusivity, SoilWatConductivity, &
                                            noahmp%water%param%SoilMoistureWilt(I,IndSoil,J), IniSoilIce, IndSoil, I, J)

       ! Maximum infiltrability based on the Eq. 6.25. (m/s)
       InfilFacTmp = noahmp%water%param%InfilCapillaryDynVic(I,J) * (noahmp%water%param%SoilMoistureSat(I,IndSoil,J) - noahmp%water%param%SoilMoistureWilt(I,IndSoil,J)) * &
                     (-1.0) * noahmp%config%domain%DepthSoilLayer(I,IndSoil,J)
       InfilSfcTmp = noahmp%water%param%SoilWatConductivitySat(I,IndSoil,J) + (WeighFac*(noahmp%water%param%SoilWatConductivitySat(I,IndSoil,J)-SoilWatConductivity) / &
                                                       (exp(WeighFac*1.0e-05/InfilFacTmp) - 1.0))

       ! infiltration rate at surface
       if ( noahmp%water%param%SoilWatConductivitySat(I,IndSoil,J) < noahmp%water%flux%SoilSfcInflowMean(I,J) ) then
          InfilSfcTmp = min(noahmp%water%flux%SoilSfcInflowMean(I,J), InfilSfcTmp)
       else
          InfilSfcTmp = noahmp%water%flux%SoilSfcInflowMean(I,J)
       endif
       if ( InfilSfcTmp < 0.0 ) InfilSfcTmp = SoilWatConductivity

    else

       ! estimate initial soil hydraulic conductivty (Ki in the equation) (m/s)
       call SoilDiffusivityConductivityOpt2(noahmp, SoilWatDiffusivity, SoilWatConductivity, &
                                            noahmp%water%state%SoilMoisture(I,IndSoil,J), noahmp%water%state%SoilIce(I,IndSoil,J), IndSoil, I, J)

       ! Maximum infiltrability based on the Eq. 6.25. (m/s)
       InfilFacTmp = noahmp%water%param%InfilCapillaryDynVic(I,J) * max(0.0, (noahmp%water%param%SoilMoistureSat(I,IndSoil,J) - noahmp%water%state%SoilMoisture(I,IndSoil,J))) * &
                     (-1.0) * noahmp%config%domain%DepthSoilLayer(I,IndSoil,J)
       if ( InfilFacTmp == 0.0 ) then ! infiltration at surface == saturated hydraulic conductivity
          InfilSfcTmp = SoilWatConductivity
       else
          InfilSfcTmp = noahmp%water%param%SoilWatConductivitySat(I,IndSoil,J) + (WeighFac*(noahmp%water%param%SoilWatConductivitySat(I,IndSoil,J)-SoilWatConductivity) / &
                                                          (exp(WeighFac*InfilSfcAcc/InfilFacTmp) - 1.0))
       endif

       ! infiltration rate at surface  
       if ( noahmp%water%param%SoilWatConductivitySat(I,IndSoil,J) < noahmp%water%flux%SoilSfcInflowMean(I,J) ) then
          InfilSfcTmp = min(noahmp%water%flux%SoilSfcInflowMean(I,J), InfilSfcTmp)
       else
          InfilSfcTmp = noahmp%water%flux%SoilSfcInflowMean(I,J)
       endif

       ! accumulated infiltration function
       InfilSfcAcc = InfilSfcAcc + InfilSfcTmp

    endif

  end subroutine SoilWaterInfilSmithParlange

end module SoilWaterInfilSmithParlangeMod
