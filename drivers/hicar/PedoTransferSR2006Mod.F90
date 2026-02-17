module PedoTransferSR2006Mod

!!! Compute soil water infiltration based on different soil composition

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type
  use NoahmpVarType

  implicit none

contains

  subroutine PedoTransferSR2006(NoahmpIO, noahmp, Sand, Clay, Orgm, I, J)
  !$acc routine seq
! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: PEDOTRANSFER_SR2006
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(NoahmpIO_type), intent(inout) :: NoahmpIO
    type(noahmp_type)  , intent(inout) :: noahmp
    integer, intent(in)                :: I, J

    real(kind=kind_noahmp), dimension(1:NoahmpIO%NSOIL), intent(inout) :: Sand
    real(kind=kind_noahmp), dimension(1:NoahmpIO%NSOIL), intent(inout) :: Clay
    real(kind=kind_noahmp), dimension(1:NoahmpIO%NSOIL), intent(inout) :: Orgm

! local
    integer                                                 :: k
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: theta_1500t
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: theta_1500
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: theta_33t
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: theta_33
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: theta_s33t
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: theta_s33
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: psi_et
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: psi_e
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: smcmax
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: smcref
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: smcwlt
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: smcdry
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: bexp
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: psisat
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: dksat
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: dwsat
    real(kind=kind_noahmp), dimension( 1:NoahmpIO%NSOIL )   :: quartz

! ------------------------------------------------------------------------------

    ! initialize
    smcmax  = 0.0
    smcref  = 0.0
    smcwlt  = 0.0
    smcdry  = 0.0
    bexp    = 0.0
    psisat  = 0.0
    dksat   = 0.0
    dwsat   = 0.0
    quartz  = 0.0

    do k = 1,4
      if(Sand(k) <= 0 .or. Clay(k) <= 0) then
         Sand(k) = 0.41
         Clay(k) = 0.18
      end if
      if(Orgm(k) <= 0 ) Orgm(k) = 0.0
    end do

    ! compute soil properties
    theta_1500t =   NoahmpIO%sr2006_theta_1500t_a_TABLE*Sand       &
                  + NoahmpIO%sr2006_theta_1500t_b_TABLE*Clay       &
                  + NoahmpIO%sr2006_theta_1500t_c_TABLE*Orgm       &
                  + NoahmpIO%sr2006_theta_1500t_d_TABLE*Sand*Orgm  &
                  + NoahmpIO%sr2006_theta_1500t_e_TABLE*Clay*Orgm  &
                  + NoahmpIO%sr2006_theta_1500t_f_TABLE*Sand*Clay  &
                  + NoahmpIO%sr2006_theta_1500t_g_TABLE

    theta_1500  =   theta_1500t                                     &
                  + NoahmpIO%sr2006_theta_1500_a_TABLE*theta_1500t  &
                  + NoahmpIO%sr2006_theta_1500_b_TABLE

    theta_33t   =   NoahmpIO%sr2006_theta_33t_a_TABLE*Sand       &
                  + NoahmpIO%sr2006_theta_33t_b_TABLE*Clay       &
                  + NoahmpIO%sr2006_theta_33t_c_TABLE*Orgm       &
                  + NoahmpIO%sr2006_theta_33t_d_TABLE*Sand*Orgm  &
                  + NoahmpIO%sr2006_theta_33t_e_TABLE*Clay*Orgm  &
                  + NoahmpIO%sr2006_theta_33t_f_TABLE*Sand*Clay  &
                  + NoahmpIO%sr2006_theta_33t_g_TABLE

    theta_33    =   theta_33t                                                     &
                  + NoahmpIO%sr2006_theta_33_a_TABLE*theta_33t*theta_33t  &
                  + NoahmpIO%sr2006_theta_33_b_TABLE*theta_33t            &
                  + NoahmpIO%sr2006_theta_33_c_TABLE

    theta_s33t  =   NoahmpIO%sr2006_theta_s33t_a_TABLE*Sand      &
                  + NoahmpIO%sr2006_theta_s33t_b_TABLE*Clay      &
                  + NoahmpIO%sr2006_theta_s33t_c_TABLE*Orgm      &
                  + NoahmpIO%sr2006_theta_s33t_d_TABLE*Sand*Orgm &
                  + NoahmpIO%sr2006_theta_s33t_e_TABLE*Clay*Orgm &
                  + NoahmpIO%sr2006_theta_s33t_f_TABLE*Sand*Clay &
                  + NoahmpIO%sr2006_theta_s33t_g_TABLE

    theta_s33   = theta_s33t                                      &
                  + NoahmpIO%sr2006_theta_s33_a_TABLE*theta_s33t  &
                  + NoahmpIO%sr2006_theta_s33_b_TABLE

    psi_et      =   NoahmpIO%sr2006_psi_et_a_TABLE*Sand           &
                  + NoahmpIO%sr2006_psi_et_b_TABLE*Clay           &
                  + NoahmpIO%sr2006_psi_et_c_TABLE*theta_s33      &
                  + NoahmpIO%sr2006_psi_et_d_TABLE*Sand*theta_s33 &
                  + NoahmpIO%sr2006_psi_et_e_TABLE*Clay*theta_s33 &
                  + NoahmpIO%sr2006_psi_et_f_TABLE*Sand*Clay      &
                  + NoahmpIO%sr2006_psi_et_g_TABLE

    psi_e       =   psi_et                                       &
                  + NoahmpIO%sr2006_psi_e_a_TABLE*psi_et*psi_et  &
                  + NoahmpIO%sr2006_psi_e_b_TABLE*psi_et         &
                  + NoahmpIO%sr2006_psi_e_c_TABLE

    theta_33    = max(10.0**-3.0,theta_33)  ! For numerical stability
    theta_1500  = max(10.0**-5.0,theta_1500)  ! For numerical stability

    ! assign property values
    smcwlt = theta_1500
    smcref = theta_33
    smcmax = theta_33                              &
             + theta_s33                           &
             + NoahmpIO%sr2006_smcmax_a_TABLE*Sand &
             + NoahmpIO%sr2006_smcmax_b_TABLE

    bexp   = 3.816712826 / (log(theta_33) - log(theta_1500) )
    psisat = psi_e
    dksat  = 1930.0 * (smcmax - theta_33) ** (3.0 - 1.0/bexp)
    quartz = Sand

    ! Units conversion
    psisat = max(0.1, psisat)               ! arbitrarily impose a limit of 0.1kpa
    psisat = 0.101997 * psisat              ! convert kpa to m
    dksat  = dksat / 3600000.0              ! convert mm/h to m/s
    dwsat  = dksat * psisat * bexp / smcmax ! units should be m*m/s
    smcdry = smcwlt

    ! Introducing somewhat arbitrary limits (based on NoahmpTable soil) to prevent bad things
    smcmax = max(0.32 ,min(smcmax,  0.50 ))
    smcref = max(0.17 ,min(smcref, smcmax))
    smcwlt = max(0.01 ,min(smcwlt, smcref))
    smcdry = max(0.01 ,min(smcdry, smcref))
    bexp   = max(2.50 ,min(bexp,    12.0 ))
    psisat = max(0.03 ,min(psisat,  1.00 ))
    dksat  = max(5.e-7,min(dksat,   1.e-5))
    dwsat  = max(1.e-6,min(dwsat,   3.e-5))
    quartz = max(0.05 ,min(quartz,  0.95 ))

    !$acc loop seq
    do k = 1,NoahmpIO%NSOIL
      noahmp%water%param%SoilMoistureWilt(I,k,J)       = smcwlt(k)
      noahmp%water%param%SoilMoistureFieldCap(I,k,J)   = smcref(k)
      noahmp%water%param%SoilMoistureSat(I,k,J)        = smcmax(k)
      noahmp%water%param%SoilMoistureDry(I,k,J)        = smcdry(k)
      noahmp%water%param%SoilExpCoeffB(I,k,J)          = bexp(k)
      noahmp%water%param%SoilMatPotentialSat(I,k,J)    = psisat(k)
      noahmp%water%param%SoilWatConductivitySat(I,k,J) = dksat(k)
      noahmp%water%param%SoilWatDiffusivitySat(I,k,J)  = dwsat(k)
      noahmp%energy%param%SoilQuartzFrac(I,k,J)        = quartz(k)
    enddo

  end subroutine PedoTransferSR2006

end module PedoTransferSR2006Mod
