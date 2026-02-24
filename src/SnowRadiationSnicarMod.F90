module SnowRadiationSnicarMod

!!! Determine reflectance of, and vertically-resolved solar absorption in, 
!!! snow with impurities.
!!! Original references on physical models of snow reflectance include: 
!!! Wiscombe and Warren [1980] and Warren and Wiscombe [1980],Journal of Atmospheric Sciences, 37,
!!! The multi-layer solution for multiple-scattering used here is from:
!!! Toon et al. [1989], Rapid calculation of radiative heating rates 
!!! and photodissociation rates in inhomogeneous multiple scattering atmospheres, 
!!! J. Geophys. Res., 94, D13, 16287-16301
!!! The implementation of the SNICAR model in CLM/CSIM is described in:
!!! Flanner, M., C. Zender, J. Randerson, and P. Rasch [2007], 
!!! Present-day climate forcing and response from black carbon in snow,
!!! J. Geophys. Res., 112, D11202, doi: 10.1029/2006JD008003
!!! Updated radiative transfer solver:
!!! The multi-layer solution for multiple-scattering used here is from:
!!! Briegleb, P. and Light, B.: A Delta-Eddington mutiple scattering
!!! parameterization for solar radiation in the sea ice component of the
!!! community climate system model, 2007.
!!! The implementation of the SNICAR-AD model in CLM is described in:
!!! Dang et al.2019, Inter-comparison and improvement of 2-stream shortwave
!!! radiative transfer models for unified treatment of cryospheric surfaces
!!! in ESMs; and Flanner et al. 2021, SNICAR-ADv3: a community tool for modeling 
!!! spectral snow albedo

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use PiecewiseLinearInterp1dMod

  implicit none

contains

  subroutine SnowRadiationSnicar(noahmp,FlagSwRadType)

! ------------------------ Code history -----------------------------------
! Implementation: T.-S. Lin, C. He, et al. (2025, JHM)
! Adapted from SNICAR module SNICAR_RT from CTSM
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout)   :: noahmp
    integer,           intent(in)      :: FlagSwRadType           ! flag: =1 for direct-beam incident flux,=2 for diffuse incident flux

! local variables
    ! GPU port parameters and variables
    integer, parameter                 :: NSNOW_MAX = 5
    integer, parameter                 :: NAER = 9
    integer                            :: II, JJ, k
    real(kind=kind_noahmp)             :: flx_slrd_val
    real(kind=kind_noahmp)             :: flx_slri_val
    real(kind=kind_noahmp)             :: albsfc_val
    real(kind=kind_noahmp)             :: flx_sum_wgt
    real(kind=kind_noahmp)             :: flx_sum_wgt_vis
    ! Promoted arrays (band-sized, allocated on device before parallel loop)
    real(kind=kind_noahmp), allocatable :: flx_wgt_3d(:,:,:)
    real(kind=kind_noahmp), allocatable :: albout_lcl_3d(:,:,:)
    real(kind=kind_noahmp), allocatable :: flx_abs_lcl_4d(:,:,:,:)
    ! general local variables
    integer                            :: i,idb,igb
    integer                            :: j                       ! aerosol number index [idx]
    integer                            :: n                       ! tridiagonal matrix index [idx]
    integer                            :: ng                      ! gaussian integration index
    integer, parameter                 :: ngmax = 8
    integer                            :: trip                    ! flag: =1 to redo RT calculation if result is unrealistic
    ! NAER defined as parameter above
    integer                            :: SnowLayerTop            ! top snow layer index [idx]
    integer                            :: SnowLayerBottom         ! bottom snow layer index [idx]
    integer                            :: LoopInd                 ! do loop/array indices
    integer                            :: nir_bnd_bgn             ! first band index in near-IR spectrum [idx] cenlin
    integer                            :: nir_bnd_end             ! ending near-IR band index [idx] cenlin
    integer                            :: flg_nosnl               ! flag: =1 if there is snow, but zero snow layers, =0 if at least 1 snow layer [flg]
    integer                            :: snl_lcl                 ! negative number of snow layers [nbr]
    integer                            :: flg_dover               ! defines conditions for RT redo (explained below)
    integer                            :: err_idx                 ! counter for number of times through error loop [nbr]
    integer                            :: APRX_TYP                ! two-stream approximation type (1=Eddington, 2=Quadrature, 3=Hemispheric Mean) [nbr]
    integer                            :: rds_idx                 ! snow effective radius index for retrieving, Mie parameters from lookup table [idx]
    integer                            :: snl_btm_itf             ! index of bottom snow layer interfaces (1) [idx]
    integer,  parameter                :: snw_rds_min_tbl = 30    ! minimium effective radius defined in Mie lookup table [microns]
    integer,  parameter                :: snw_rds_max_tbl = 1500  ! maximum effective radius defined in Mie lookup table [microns]
    integer                            :: snw_rds_lcl(-NSNOW_MAX+1:0)             ! snow effective radius [m^-6]
    real(kind=kind_noahmp)             :: tau_sum                 ! cumulative (snow+aerosol) optical depth [unitless]
    real(kind=kind_noahmp)             :: omega_sum               ! temporary summation of single-scatter albedo of all aerosols [frc]
    real(kind=kind_noahmp)             :: g_sum                   ! temporary summation of asymmetry parameter of all aerosols [frc]
    real(kind=kind_noahmp)             :: F_direct_btm            ! direct-beam radiation at bottom of snowpack [W/m^2]
    real(kind=kind_noahmp)             :: F_sfc_pls               ! upward radiative flux at snowpack top [W/m^2]
    real(kind=kind_noahmp)             :: F_btm_net               ! net flux at bottom of snowpack [W/m^2]
    real(kind=kind_noahmp)             :: F_sfc_net               ! net flux at top of snowpack [W/m^2]
    real(kind=kind_noahmp)             :: energy_sum              ! sum of all energy terms; should be 0.0 [W/m^2]
    real(kind=kind_noahmp)             :: albedo                  ! temporary snow albedo [frc]
    real(kind=kind_noahmp)             :: F_abs_sum               ! total absorbed energy in column [W/m^2]
    real(kind=kind_noahmp)             :: flx_sum                 ! temporary summation variable for NIR weighting

    ! local constant and coefficients used for SZA parameterization
    real(kind=kind_noahmp), parameter  :: sza_a0 =  0.085730_kind_noahmp
    real(kind=kind_noahmp), parameter  :: sza_a1 = -0.630883_kind_noahmp
    real(kind=kind_noahmp), parameter  :: sza_a2 =  1.303723_kind_noahmp
    real(kind=kind_noahmp), parameter  :: sza_b0 =  1.467291_kind_noahmp
    real(kind=kind_noahmp), parameter  :: sza_b1 = -3.338043_kind_noahmp
    real(kind=kind_noahmp), parameter  :: sza_b2 =  6.807489_kind_noahmp
    real(kind=kind_noahmp), parameter  :: puny   =  1.0e-11_kind_noahmp
    real(kind=kind_noahmp), parameter  :: mu_75  =  0.2588_kind_noahmp   ! cosine of 75 degree
    real(kind=kind_noahmp)             :: sza_c1                  ! coefficient, SZA parameteirzation
    real(kind=kind_noahmp)             :: sza_c0                  ! coefficient, SZA parameterization
    real(kind=kind_noahmp)             :: sza_factor              ! factor used to adjust NIR direct albedo
    real(kind=kind_noahmp)             :: flx_sza_adjust          ! direct NIR flux adjustment from sza_factor
    real(kind=kind_noahmp)             :: mu0                     ! incident solar zenith angle

    ! local constants used in algorithm
    real(kind=kind_noahmp)             :: mu_not                  ! cosine of solar zenith angle (used locally) [frc]
    real(kind=kind_noahmp), parameter  :: c0     = 0.0_kind_noahmp
    real(kind=kind_noahmp), parameter  :: c1     = 1.0_kind_noahmp
    real(kind=kind_noahmp), parameter  :: c3     = 3.0_kind_noahmp
    real(kind=kind_noahmp), parameter  :: c4     = 4.0_kind_noahmp
    real(kind=kind_noahmp), parameter  :: c6     = 6.0_kind_noahmp
    real(kind=kind_noahmp), parameter  :: cp01   = 0.01_kind_noahmp
    real(kind=kind_noahmp), parameter  :: cp5    = 0.5_kind_noahmp
    real(kind=kind_noahmp), parameter  :: cp75   = 0.75_kind_noahmp
    real(kind=kind_noahmp), parameter  :: c1p5   = 1.5_kind_noahmp
    real(kind=kind_noahmp), parameter  :: trmin  = 0.001_kind_noahmp
    real(kind=kind_noahmp), parameter  :: argmax = 10.0_kind_noahmp  ! maximum argument of exponential
    real(kind=kind_noahmp)             :: wvl_ct5(1:5)            ! band center wavelength (um) for 5-band case
    real(kind=kind_noahmp)             :: SnowWaterEquivMin       ! minimum snow mass required for SNICAR RT calculation [kg m-2] !samlin, may need to change this
    real(kind=kind_noahmp)             :: diam_ice                ! effective snow grain diameter (SSA-equivalent) unit: microns
    real(kind=kind_noahmp)             :: fs_sphd                 ! shape factor for spheroid snow
    real(kind=kind_noahmp)             :: fs_hex                  ! shape factor for reference hexagonal snow
    real(kind=kind_noahmp)             :: fs_hex0                 ! shape factor for hexagonal plate
    real(kind=kind_noahmp)             :: fs_koch                 ! shape factor for Koch snowflake
    real(kind=kind_noahmp)             :: AR_tmp                  ! aspect ratio temporary
    real(kind=kind_noahmp)             :: g_ice_Cg_tmp(1:7)       ! temporary asymmetry factor correction coeff
    real(kind=kind_noahmp)             :: gg_ice_F07_tmp(1:7)     ! temporary asymmetry factor related to geometric reflection & refraction
    real(kind=kind_noahmp)             :: g_Cg_intp               ! interpolated asymmetry factor correction coeff to target bands
    real(kind=kind_noahmp)             :: gg_F07_intp             ! interpolated asymmetry factor related to geometric reflection & refraction
    real(kind=kind_noahmp)             :: g_ice_F07               ! asymmetry factor for Fu 2007 parameterization value

    ! local variables used for nonspherical snow grain treatment (He et al. 2017 J of Climate):
    ! Constants and parameters for aspherical ice particles    
    ! asymmetry factor parameterization coefficients (6 bands) from Table 3 & Eqs. 6-7 in He et al. (2017)
    real(kind=kind_noahmp)             :: g_wvl(1:8)              ! wavelength (um) division point
    real(kind=kind_noahmp)             :: g_wvl_ct(1:7)           ! center point for wavelength band (um)
    real(kind=kind_noahmp)             :: g_b0(1:7)
    real(kind=kind_noahmp)             :: g_b1(1:7)
    real(kind=kind_noahmp)             :: g_b2(1:7)
    ! Tables 1 & 2 and Eqs. 3.1-3.4 from Fu, 2007 JAS
    real(kind=kind_noahmp)             :: g_F07_c2(1:7)
    real(kind=kind_noahmp)             :: g_F07_c1(1:7)
    real(kind=kind_noahmp)             :: g_F07_c0(1:7)
    real(kind=kind_noahmp)             :: g_F07_p2(1:7)
    real(kind=kind_noahmp)             :: g_F07_p1(1:7)
    real(kind=kind_noahmp)             :: g_F07_p0(1:7)
    ! variables used for BC-snow internal mixing (He et al. 2017 J of Climate):
    real(kind=kind_noahmp)             :: enh_omg_bcint            ! BC-induced enhancement in snow single-scattering co-albedo (1-omega)
    real(kind=kind_noahmp)             :: enh_omg_bcint_tmp(1:16)  ! temporary BC-induced enhancement in snow 1-omega
    real(kind=kind_noahmp)             :: enh_omg_bcint_tmp2(1:16) ! temporary BC-induced enhancement in snow 1-omega
    real(kind=kind_noahmp)             :: bcint_wvl(1:17)          ! Parameterization band (0.2-1.2um) for BC-induced enhancement in snow 1-omega
    real(kind=kind_noahmp)             :: bcint_wvl_ct(1:16)       ! Parameterization band center wavelength (um)
    real(kind=kind_noahmp)             :: bcint_d0(1:16)           ! Parameterization coefficients at each band center wavelength
    real(kind=kind_noahmp)             :: bcint_d1(1:16)           ! Parameterization coefficients at each band center wavelength
    real(kind=kind_noahmp)             :: bcint_d2(1:16)           ! Parameterization coefficients at each band center wavelength
    real(kind=kind_noahmp), parameter  :: den_bc = 1.49_kind_noahmp   ! target BC particle density (g/cm3) used in BC MAC adjustment
    real(kind=kind_noahmp), parameter  :: Re_bc = 0.045_kind_noahmp  ! target BC effective radius (um) used in BC MAC adjustment
    real(kind=kind_noahmp)             :: bcint_m(1:3)             ! Parameterization coefficients for BC size adjustment in BC-snow int mix
    real(kind=kind_noahmp)             :: bcint_n(1:3)             ! Parameterization coefficients for BC size adjustment in BC-snow int mix
    real(kind=kind_noahmp)             :: bcint_dd                 ! intermediate parameter
    real(kind=kind_noahmp)             :: bcint_dd2                ! intermediate parameter
    real(kind=kind_noahmp)             :: bcint_f                  ! intermediate parameter
    real(kind=kind_noahmp)             :: enh_omg_bcint_intp       ! BC-induced enhancement in snow 1-omega (logscale) interpolated to CLM wavelength
    real(kind=kind_noahmp)             :: enh_omg_bcint_intp2      ! BC-induced enhancement in snow 1-omega interpolated to CLM wavelength
    real(kind=kind_noahmp)             :: wvl_doint                ! wavelength doing BC-snow int mixing (<=1.2um)
    integer                            :: ibb                      ! loop index

    ! local variables used for dust-snow internal mixing (He et al. 2019 JAMES):
    real(kind=kind_noahmp)             :: enh_omg_dstint           ! dust-induced enhancement in snow single-scattering co-albedo (1-omega)
    real(kind=kind_noahmp)             :: enh_omg_dstint_tmp(1:6)  ! temporary dust-induced enhancement in snow 1-omega
    real(kind=kind_noahmp)             :: enh_omg_dstint_tmp2(1:6) ! temporary dust-induced enhancement in snow 1-omega
    real(kind=kind_noahmp)             :: dstint_wvl(1:7)          ! Parameterization band (0.2-1.2um) for dust-induced enhancement in snow 1-omega
    real(kind=kind_noahmp)             :: dstint_wvl_ct(1:6)       ! Parameterization band center wavelength (um)
    real(kind=kind_noahmp)             :: dstint_a1(1:6)           ! Parameterization coefficients at each band center wavelength
    real(kind=kind_noahmp)             :: dstint_a2(1:6)           ! Parameterization coefficients at each band center wavelength
    real(kind=kind_noahmp)             :: dstint_a3(1:6)           ! Parameterization coefficients at each band center wavelength
    real(kind=kind_noahmp)             :: enh_omg_dstint_intp      ! dust-induced enhancement in snow 1-omega (logscale) interpolated to CLM wavelength
    real(kind=kind_noahmp)             :: enh_omg_dstint_intp2     ! dust-induced enhancement in snow 1-omega interpolated to CLM wavelength
    real(kind=kind_noahmp)             :: tot_dst_snw_conc         ! total dust content in snow across all size bins (ppm=ug/g)
    real(kind=kind_noahmp)             :: sno_shp(-NSNOW_MAX+1:0)   ! Snow shape type: 1=sphere; 2=spheroid; 3=hexagonal plate; 4=koch snowflake
                                                                   ! currently only assuming same shapes for all snow layers
    real(kind=kind_noahmp)             :: sno_fs(-NSNOW_MAX+1:0)    ! Snow shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
                                                                   ! only activated when OptSnicarSnowShape > 1 (i.e. nonspherical)
                                                                   ! 0=use recommended default value (He et al. 2017);
                                                                   ! others(0<sno_fs<1)= user-specified value
    real(kind=kind_noahmp)             :: sno_AR(-NSNOW_MAX+1:0)    ! Snow grain aspect ratio: ratio of grain width to length
                                                                   ! only activated when snicar_snw_shape > 1 (i.e. nonspherical)
                                                                   ! 0=use recommended default value (He et al. 2017);
                                                                   ! others(0.1<fs<20)= use user-specified value
    ! other local snow and energy flux variables
    real(kind=kind_noahmp)             :: h2osno_liq_lcl(-NSNOW_MAX+1:0)       ! liquid water mass [kg/m2]
    real(kind=kind_noahmp)             :: h2osno_ice_lcl(-NSNOW_MAX+1:0)       ! ice mass [kg/m2]
    real(kind=kind_noahmp)             :: ss_alb_snw_lcl(-NSNOW_MAX+1:0)       ! single-scatter albedo of ice grains (lyr) [frc]
    real(kind=kind_noahmp)             :: asm_prm_snw_lcl(-NSNOW_MAX+1:0)      ! asymmetry parameter of ice grains (lyr) [frc]
    real(kind=kind_noahmp)             :: ext_cff_mss_snw_lcl(-NSNOW_MAX+1:0)  ! mass extinction coefficient of ice grains (lyr) [m2/kg]
    real(kind=kind_noahmp)             :: ss_alb_aer_lcl(1:NAER)       ! single-scatter albedo of aerosol species (aer_nbr) [frc]
    real(kind=kind_noahmp)             :: asm_prm_aer_lcl(1:NAER)      ! asymmetry parameter of aerosol species (aer_nbr) [frc]
    real(kind=kind_noahmp)             :: ext_cff_mss_aer_lcl(1:NAER)  ! mass extinction coefficient of aerosol species (aer_nbr) [m2/kg]
    real(kind=kind_noahmp)             :: F_direct(-NSNOW_MAX+1:0)             ! direct-beam radiation at bottom of layer interface (lyr) [W/m^2]
    real(kind=kind_noahmp)             :: F_net(-NSNOW_MAX+1:0)                ! net radiative flux at bottom of layer interface (lyr) [W/m^2]
    real(kind=kind_noahmp)             :: F_abs(-NSNOW_MAX+1:0)                ! net absorbed radiative energy (lyr) [W/m^2]
    real(kind=kind_noahmp)             :: L_snw(-NSNOW_MAX+1:0)                ! h2o mass (liquid+solid) in snow layer (lyr) [kg/m2]
    real(kind=kind_noahmp)             :: tau_snw(-NSNOW_MAX+1:0)              ! snow optical depth (lyr) [unitless]
    real(kind=kind_noahmp)             :: tau(-NSNOW_MAX+1:0)                  ! weighted optical depth of snow+aerosol layer (lyr) [unitless]
    real(kind=kind_noahmp)             :: omega(-NSNOW_MAX+1:0)                ! weighted single-scatter albedo of snow+aerosol layer (lyr) [frc]
    real(kind=kind_noahmp)             :: g(-NSNOW_MAX+1:0)                    ! weighted asymmetry parameter of snow+aerosol layer (lyr) [frc]
    real(kind=kind_noahmp)             :: tau_star(-NSNOW_MAX+1:0)             ! transformed (i.e. Delta-Eddington) optical depth of snow+aerosol layer! (lyr) [unitless]
    real(kind=kind_noahmp)             :: omega_star(-NSNOW_MAX+1:0)           ! transformed (i.e. Delta-Eddington) SSA of snow+aerosol layer (lyr) [frc]
    real(kind=kind_noahmp)             :: g_star(-NSNOW_MAX+1:0)               ! transformed (i.e. Delta-Eddington) asymmetry paramater of snow+aerosol layer! (lyr) [frc]
    real(kind=kind_noahmp)             :: tau_clm(-NSNOW_MAX+1:0)              ! column optical depth from layer bottom to snowpack top (lyr) [unitless]
    real(kind=kind_noahmp)             :: L_aer(-NSNOW_MAX+1:0, 1:NAER)              ! aerosol mass in snow layer (lyr,nbr_aer) [kg/m2]
    real(kind=kind_noahmp)             :: tau_aer(-NSNOW_MAX+1:0, 1:NAER)            ! aerosol optical depth (lyr,nbr_aer) [unitless]
    real(kind=kind_noahmp)             :: mss_cnc_aer_lcl(-NSNOW_MAX+1:0, 1:NAER)    ! aerosol mass concentration [kg/kg]

    ! local variables used for Toon et al. 1989 2-stream solver (Flanner et al. 2007):
    ! intermediate variables for radiative transfer approximation:
    real(kind=kind_noahmp)             :: gamma1(-NSNOW_MAX+1:0)             ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(kind=kind_noahmp)             :: gamma2(-NSNOW_MAX+1:0)             ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(kind=kind_noahmp)             :: gamma3(-NSNOW_MAX+1:0)             ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(kind=kind_noahmp)             :: gamma4(-NSNOW_MAX+1:0)             ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(kind=kind_noahmp)             :: lambda(-NSNOW_MAX+1:0)             ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(kind=kind_noahmp)             :: GAMMA(-NSNOW_MAX+1:0)              ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(kind=kind_noahmp)                              :: mu_one             ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(kind=kind_noahmp)             :: e1(-NSNOW_MAX+1:0)                 ! tri-diag intermediate variable from Toon et al. (lyr)
    real(kind=kind_noahmp)             :: e2(-NSNOW_MAX+1:0)                 ! tri-diag intermediate variable from Toon et al. (lyr)
    real(kind=kind_noahmp)             :: e3(-NSNOW_MAX+1:0)                 ! tri-diag intermediate variable from Toon et al. (lyr)
    real(kind=kind_noahmp)             :: e4(-NSNOW_MAX+1:0)                 ! tri-diag intermediate variable from Toon et al. (lyr)
    real(kind=kind_noahmp)             :: C_pls_btm(-NSNOW_MAX+1:0)          ! intermediate variable: upward flux at bottom interface (lyr) [W/m2]
    real(kind=kind_noahmp)             :: C_mns_btm(-NSNOW_MAX+1:0)          ! intermediate variable: downward flux at bottom interface (lyr) [W/m2]
    real(kind=kind_noahmp)             :: C_pls_top(-NSNOW_MAX+1:0)          ! intermediate variable: upward flux at top interface (lyr) [W/m2]
    real(kind=kind_noahmp)             :: C_mns_top(-NSNOW_MAX+1:0)          ! intermediate variable: downward flux at top interface (lyr) [W/m2]
    real(kind=kind_noahmp)             :: A(-2*NSNOW_MAX+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(kind=kind_noahmp)             :: B(-2*NSNOW_MAX+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(kind=kind_noahmp)             :: D(-2*NSNOW_MAX+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(kind=kind_noahmp)             :: E(-2*NSNOW_MAX+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(kind=kind_noahmp)             :: AS(-2*NSNOW_MAX+1:0)                 ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(kind=kind_noahmp)             :: DS(-2*NSNOW_MAX+1:0)                 ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(kind=kind_noahmp)             :: X(-2*NSNOW_MAX+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(kind=kind_noahmp)             :: Y(-2*NSNOW_MAX+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)

    ! local variables used for Adding-doubling 2-stream solver based on SNICAR-ADv3 version
    ! (Dang et al. 2019; Flanner et al. 2021)
    real(kind=kind_noahmp)             :: trndir(-NSNOW_MAX+1:1)             ! solar beam down transmission from top
    real(kind=kind_noahmp)             :: trntdr(-NSNOW_MAX+1:1)             ! total transmission to direct beam for layers above
    real(kind=kind_noahmp)             :: trndif(-NSNOW_MAX+1:1)             ! diffuse transmission to diffuse beam for layers above
    real(kind=kind_noahmp)             :: rupdir(-NSNOW_MAX+1:1)             ! reflectivity to direct radiation for layers below
    real(kind=kind_noahmp)             :: rupdif(-NSNOW_MAX+1:1)             ! reflectivity to diffuse radiation for layers below
    real(kind=kind_noahmp)             :: rdndif(-NSNOW_MAX+1:1)             ! reflectivity to diffuse radiation for layers above
    real(kind=kind_noahmp)             :: dfdir(-NSNOW_MAX+1:1)              ! down-up flux at interface due to direct beam at top surface
    real(kind=kind_noahmp)             :: dfdif(-NSNOW_MAX+1:1)              ! down-up flux at interface due to diffuse beam at top surface
    real(kind=kind_noahmp)             :: dftmp(-NSNOW_MAX+1:1)              ! temporary variable for down-up flux at interface
    real(kind=kind_noahmp)             :: rdir(-NSNOW_MAX+1:0)               ! layer reflectivity to direct radiation
    real(kind=kind_noahmp)             :: rdif_a(-NSNOW_MAX+1:0)             ! layer reflectivity to diffuse radiation from above
    real(kind=kind_noahmp)             :: rdif_b(-NSNOW_MAX+1:0)             ! layer reflectivity to diffuse radiation from below
    real(kind=kind_noahmp)             :: tdir(-NSNOW_MAX+1:0)               ! layer transmission to direct radiation (solar beam + diffuse)
    real(kind=kind_noahmp)             :: tdif_a(-NSNOW_MAX+1:0)             ! layer transmission to diffuse radiation from above
    real(kind=kind_noahmp)             :: tdif_b(-NSNOW_MAX+1:0)             ! layer transmission to diffuse radiation from below
    real(kind=kind_noahmp)             :: trnlay(-NSNOW_MAX+1:0)             ! solar beam transm for layer (direct beam only)
    real(kind=kind_noahmp)                              :: ts                 ! layer delta-scaled extinction optical depth
    real(kind=kind_noahmp)                              :: ws                 ! layer delta-scaled single scattering albedo
    real(kind=kind_noahmp)                              :: gs                 ! layer delta-scaled asymmetry parameter
    real(kind=kind_noahmp)                              :: extins             ! extinction
    real(kind=kind_noahmp)                              :: alp                ! temporary for alpha
    real(kind=kind_noahmp)                              :: gam                ! temporary for agamm
    real(kind=kind_noahmp)                              :: amg                ! alp - gam
    real(kind=kind_noahmp)                              :: apg                ! alp + gam
    real(kind=kind_noahmp)                              :: ue                 ! temporary for u
    real(kind=kind_noahmp)                              :: refk               ! interface multiple scattering
    real(kind=kind_noahmp)                              :: refkp1             ! interface multiple scattering for k+1
    real(kind=kind_noahmp)                              :: refkm1             ! interface multiple scattering for k-1
    real(kind=kind_noahmp)                              :: tdrrdir            ! direct tran times layer direct ref
    real(kind=kind_noahmp)                              :: tdndif             ! total down diffuse = tot tran - direct tran
    real(kind=kind_noahmp)                              :: taus               ! scaled extinction optical depth
    real(kind=kind_noahmp)                              :: omgs               ! scaled single particle scattering albedo
    real(kind=kind_noahmp)                              :: asys               ! scaled asymmetry parameter
    real(kind=kind_noahmp)                              :: lm                 ! temporary for el
    real(kind=kind_noahmp)                              :: mu                 ! cosine solar zenith for either snow or water
    real(kind=kind_noahmp)                              :: ne                 ! temporary for n
    real(kind=kind_noahmp)                              :: R1                 ! perpendicular polarization reflection amplitude
    real(kind=kind_noahmp)                              :: R2                 ! parallel polarization reflection amplitude
    real(kind=kind_noahmp)                              :: T1                 ! perpendicular polarization transmission amplitude
    real(kind=kind_noahmp)                              :: T2                 ! parallel polarization transmission amplitude
    real(kind=kind_noahmp)                              :: Rf_dir_a           ! fresnel reflection to direct radiation
    real(kind=kind_noahmp)                              :: Tf_dir_a           ! fresnel transmission to direct radiation
    real(kind=kind_noahmp)                              :: Rf_dif_a           ! fresnel reflection to diff radiation from above
    real(kind=kind_noahmp)                              :: Rf_dif_b           ! fresnel reflection to diff radiation from below
    real(kind=kind_noahmp)                              :: Tf_dif_a           ! fresnel transmission to diff radiation from above
    real(kind=kind_noahmp)                              :: Tf_dif_b           ! fresnel transmission to diff radiation from below
    real(kind=kind_noahmp)                              :: gwt                ! gaussian weight
    real(kind=kind_noahmp)                              :: swt                ! sum of weights
    real(kind=kind_noahmp)                              :: trn                ! layer transmission
    real(kind=kind_noahmp)                              :: rdr                ! rdir for gaussian integration
    real(kind=kind_noahmp)                              :: tdr                ! tdir for gaussian integration
    real(kind=kind_noahmp)                              :: smr                ! accumulator for rdif gaussian integration
    real(kind=kind_noahmp)                              :: smt                ! accumulator for tdif gaussian integration
    real(kind=kind_noahmp)                              :: exp_min            ! minimum exponential value
    real(kind=kind_noahmp)                              :: difgauspt(1:8)     ! Gaussian integration angle
    real(kind=kind_noahmp)                              :: difgauswt(1:8)     ! Gaussian integration coefficients/weights

! --------------------------------------------------------------------
    associate(                                                                            &  
              IndicatorIceSfc          => noahmp%config%domain%IndicatorIceSfc           ,& ! in,  indicator for ice surface/point (1=sea ice, 0=non-ice, -1=land ice)
              OptSnicarSnowShape       => noahmp%config%nmlist%OptSnicarSnowShape        ,& ! in,  Snow shape: 1=sphere; 2=spheroid; 3=hexagonal plate; 4=koch snowflake  
              OptSnicarRTSolver        => noahmp%config%nmlist%OptSnicarRTSolver         ,& ! in,  SNICAR radiative transfer solver
              FlagSnicarSnowBCIntmix   => noahmp%config%nmlist%FlagSnicarSnowBCIntmix    ,& ! in,  flag to activate BC-snow internal mixing in SNICAR (He et al. 2017 JC)
              FlagSnicarSnowDustIntmix => noahmp%config%nmlist%FlagSnicarSnowDustIntmix  ,& ! in,  flag to activate dust-snow internal mixing in SNICAR (He et al. 2017 JC)
              FlagSnicarUseAerosol     => noahmp%config%nmlist%FlagSnicarUseAerosol      ,& ! in,  flag to turn on/off aerosol deposition flux effect in snow in SNICAR
              FlagSnicarUseOC          => noahmp%config%nmlist%FlagSnicarUseOC           ,& ! in,  flag to activate OC in snow in SNICAR
              NumSnicarRadBand         => noahmp%config%domain%NumSnicarRadBand          ,& ! in,  wavelength bands used in SNICAR snow albedo calculation
              NumSwRadBand             => noahmp%config%domain%NumSwRadBand              ,& ! in,  number of shortwave radiation bands 
              NumSnowLayerMax          => noahmp%config%domain%NumSnowLayerMax           ,& ! in,  maximum number of snow layers
              NumSnowLayerNeg          => noahmp%config%domain%NumSnowLayerNeg           ,& ! in,  actual number of snow layers (negative)
              CosSolarZenithAngle      => noahmp%config%domain%CosSolarZenithAngle       ,& ! in,  cosine solar zenith angle
              SnowIce                  => noahmp%water%state%SnowIce                     ,& ! in,  snow layer ice [mm]
              SnowLiqWater             => noahmp%water%state%SnowLiqWater                ,& ! in,  snow layer liquid water [mm]
              SnowWaterEquiv           => noahmp%water%state%SnowWaterEquiv              ,& ! in,  snow water equivalent [mm]
              SnowRadius               => noahmp%water%state%SnowRadius                  ,& ! in,  effective grain radius [microns, m-6]
              AlbedoSoilDif            => noahmp%energy%state%AlbedoSoilDif              ,& ! in,  soil albedo (diffuse)
              AlbedoSoilDir            => noahmp%energy%state%AlbedoSoilDir              ,& ! in,  soil albedo (direct)
              AlbedoLandIce            => noahmp%energy%param%AlbedoLandIce              ,& ! in,  albedo land ice: 1=vis, 2=nir
              RadSwWgtDir              => noahmp%energy%param%RadSwWgtDir                ,& ! in,  downward solar radiation spectral weights (direct)
              RadSwWgtDif              => noahmp%energy%param%RadSwWgtDif                ,& ! in,  downward solar radiation spectral weights (diffuse)
              SsAlbSnwRadDir           => noahmp%energy%param%SsAlbSnwRadDir             ,& ! in,  Mie single scatter albedos for direct-beam ice
              AsyPrmSnwRadDir          => noahmp%energy%param%AsyPrmSnwRadDir            ,& ! in,  asymmetry parameter of direct-beam ice  
              ExtCffMassSnwRadDir      => noahmp%energy%param%ExtCffMassSnwRadDir        ,& ! in,  mass extinction coefficient for direct-beam ice [m2/kg]
              SsAlbSnwRadDif           => noahmp%energy%param%SsAlbSnwRadDif             ,& ! in,  Mie single scatter albedos for diffuse ice
              AsyPrmSnwRadDif          => noahmp%energy%param%AsyPrmSnwRadDif            ,& ! in,  asymmetry parameter of diffuse ice  
              ExtCffMassSnwRadDif      => noahmp%energy%param%ExtCffMassSnwRadDif        ,& ! in,  mass extinction coefficient for diffuse ice [m2/kg]
              SsAlbBCphi               => noahmp%energy%param%SsAlbBCphi                 ,& ! in,  Mie single scatter albedos for hydrophillic BC
              AsyPrmBCphi              => noahmp%energy%param%AsyPrmBCphi                ,& ! in,  asymmetry parameter for hydrophillic BC
              ExtCffMassBCphi          => noahmp%energy%param%ExtCffMassBCphi            ,& ! in,  mass extinction coefficient for hydrophillic BC [m2/kg]
              SsAlbBCpho               => noahmp%energy%param%SsAlbBCpho                 ,& ! in,  Mie single scatter albedos for hydrophobic BC
              AsyPrmBCpho              => noahmp%energy%param%AsyPrmBCpho                ,& ! in,  asymmetry parameter for hydrophobic BC
              ExtCffMassBCpho          => noahmp%energy%param%ExtCffMassBCpho            ,& ! in,  mass extinction coefficient for hydrophobic BC [m2/kg]
              SsAlbOCphi               => noahmp%energy%param%SsAlbOCphi                 ,& ! in,  Mie single scatter albedos for hydrophillic OC
              AsyPrmOCphi              => noahmp%energy%param%AsyPrmOCphi                ,& ! in,  asymmetry parameter for hydrophillic OC
              ExtCffMassOCphi          => noahmp%energy%param%ExtCffMassOCphi            ,& ! in,  mass extinction coefficient for hydrophillic OC [m2/kg]
              SsAlbOCpho               => noahmp%energy%param%SsAlbOCpho                 ,& ! in,  Mie single scatter albedos for hydrophobic OC
              AsyPrmOCpho              => noahmp%energy%param%AsyPrmOCpho                ,& ! in,  asymmetry parameter for hydrophobic OC
              ExtCffMassOCpho          => noahmp%energy%param%ExtCffMassOCpho            ,& ! in,  mass extinction coefficient for hydrophobic OC [m2/kg]
              SsAlbDustB1              => noahmp%energy%param%SsAlbDustB1                ,& ! in,  Mie single scatter albedos for dust species 1
              AsyPrmDustB1             => noahmp%energy%param%AsyPrmDustB1               ,& ! in,  asymmetry parameter for dust species 1
              ExtCffMassDustB1         => noahmp%energy%param%ExtCffMassDustB1           ,& ! in,  mass extinction coefficient for dust species 1 [m2/kg]
              SsAlbDustB2              => noahmp%energy%param%SsAlbDustB2                ,& ! in,  Mie single scatter albedos for dust species 2
              AsyPrmDustB2             => noahmp%energy%param%AsyPrmDustB2               ,& ! in,  asymmetry parameter for dust species 2
              ExtCffMassDustB2         => noahmp%energy%param%ExtCffMassDustB2           ,& ! in,  mass extinction coefficient for dust species 2 [m2/kg]
              SsAlbDustB3              => noahmp%energy%param%SsAlbDustB3                ,& ! in,  Mie single scatter albedos for dust species 3
              AsyPrmDustB3             => noahmp%energy%param%AsyPrmDustB3               ,& ! in,  asymmetry parameter for dust species 3
              ExtCffMassDustB3         => noahmp%energy%param%ExtCffMassDustB3           ,& ! in,  mass extinction coefficient for dust species 3 [m2/kg]
              SsAlbDustB4              => noahmp%energy%param%SsAlbDustB4                ,& ! in,  Mie single scatter albedos for dust species 4
              AsyPrmDustB4             => noahmp%energy%param%AsyPrmDustB4               ,& ! in,  asymmetry parameter for dust species 4
              ExtCffMassDustB4         => noahmp%energy%param%ExtCffMassDustB4           ,& ! in,  mass extinction coefficient for dust species 4 [m2/kg]
              SsAlbDustB5              => noahmp%energy%param%SsAlbDustB5                ,& ! in,  Mie single scatter albedos for dust species 5
              AsyPrmDustB5             => noahmp%energy%param%AsyPrmDustB5               ,& ! in,  asymmetry parameter for dust species 5
              ExtCffMassDustB5         => noahmp%energy%param%ExtCffMassDustB5           ,& ! in,  mass extinction coefficient for dust species 5 [m2/kg]
              MassConcBChydropho       => noahmp%water%state%MassConcBChydropho          ,& ! in,  mass concentration of hydrophobic Black Carbon in snow [kg/kg]
              MassConcBChydrophi       => noahmp%water%state%MassConcBChydrophi          ,& ! in,  mass concentration of hydrophillic Black Carbon in snow [kg/kg]
              MassConcOChydropho       => noahmp%water%state%MassConcOChydropho          ,& ! in,  mass concentration of hydrophobic Organic Carbon in snow [kg/kg]
              MassConcOChydrophi       => noahmp%water%state%MassConcOChydrophi          ,& ! in,  mass concentration of hydrophillic Organic Carbon in snow [kg/kg]
              MassConcDust1            => noahmp%water%state%MassConcDust1               ,& ! in,  mass concentration of dust species 1 in snow [kg/kg]
              MassConcDust2            => noahmp%water%state%MassConcDust2               ,& ! in,  mass concentration of dust species 2 in snow [kg/kg]
              MassConcDust3            => noahmp%water%state%MassConcDust3               ,& ! in,  mass concentration of dust species 3 in snow [kg/kg]
              MassConcDust4            => noahmp%water%state%MassConcDust4               ,& ! in,  mass concentration of dust species 4 in snow [kg/kg]
              MassConcDust5            => noahmp%water%state%MassConcDust5               ,& ! in,  mass concentration of dust species 5 in snow [kg/kg]
              AlbedoSnowDir            => noahmp%energy%state%AlbedoSnowDir              ,& ! out, snow albedo for direct(1=vis, 2=nir)
              AlbedoSnowDif            => noahmp%energy%state%AlbedoSnowDif              ,& ! out, snow albedo for diffuse(1=vis, 2=nir)
              FracRadSwAbsSnowDir      => noahmp%energy%flux%FracRadSwAbsSnowDir         ,& ! out, direct solar flux factor absorbed by snow [frc]
              FracRadSwAbsSnowDif      => noahmp%energy%flux%FracRadSwAbsSnowDif          & ! out, diffuse solar flux factor absorbed by snow [frc]
             )
! ----------------------------------------------------------------------

    ! Allocate promoted arrays (band-sized, on device)
    allocate(flx_wgt_3d(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
             1:NumSnicarRadBand, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(albout_lcl_3d(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
             1:NumSnicarRadBand, noahmp%config%domain%JTS:noahmp%config%domain%JTE))
    allocate(flx_abs_lcl_4d(noahmp%config%domain%ITS:noahmp%config%domain%ITE, &
             -NumSnowLayerMax+1:1, 1:NumSnicarRadBand, noahmp%config%domain%JTS:noahmp%config%domain%JTE))

    !$acc data create(flx_wgt_3d, albout_lcl_3d, flx_abs_lcl_4d)

    !$acc parallel loop collapse(2) gang vector default(present) &
    !$acc private(i, j, n, ng, trip, SnowLayerTop, SnowLayerBottom, LoopInd, &
    !$acc   nir_bnd_bgn, nir_bnd_end, flg_nosnl, snl_lcl, flg_dover, err_idx, &
    !$acc   APRX_TYP, rds_idx, snl_btm_itf, ibb, idb, igb, &
    !$acc   snw_rds_lcl, sno_shp, sno_fs, sno_AR, &
    !$acc   h2osno_liq_lcl, h2osno_ice_lcl, F_direct, F_net, F_abs, L_snw, &
    !$acc   tau_snw, tau, omega, g, tau_star, omega_star, g_star, tau_clm, &
    !$acc   ss_alb_snw_lcl, asm_prm_snw_lcl, ext_cff_mss_snw_lcl, &
    !$acc   ss_alb_aer_lcl, asm_prm_aer_lcl, ext_cff_mss_aer_lcl, &
    !$acc   L_aer, tau_aer, mss_cnc_aer_lcl, &
    !$acc   gamma1, gamma2, gamma3, gamma4, lambda, GAMMA, &
    !$acc   e1, e2, e3, e4, C_pls_btm, C_mns_btm, C_pls_top, C_mns_top, &
    !$acc   A, B, D, E, AS, DS, X, Y, &
    !$acc   trndir, trntdr, trndif, rupdir, rupdif, rdndif, &
    !$acc   dfdir, dfdif, dftmp, rdir, rdif_a, rdif_b, tdir, tdif_a, tdif_b, trnlay, &
    !$acc   difgauspt, difgauswt, wvl_ct5, &
    !$acc   g_wvl, g_wvl_ct, g_b0, g_b1, g_b2, &
    !$acc   g_F07_c2, g_F07_c1, g_F07_c0, g_F07_p2, g_F07_p1, g_F07_p0, &
    !$acc   g_ice_Cg_tmp, gg_ice_F07_tmp, &
    !$acc   bcint_wvl, bcint_wvl_ct, bcint_d0, bcint_d1, bcint_d2, bcint_m, bcint_n, &
    !$acc   enh_omg_bcint_tmp, enh_omg_bcint_tmp2, &
    !$acc   dstint_wvl, dstint_wvl_ct, dstint_a1, dstint_a2, dstint_a3, &
    !$acc   enh_omg_dstint_tmp, enh_omg_dstint_tmp2, &
    !$acc   tau_sum, omega_sum, g_sum, F_direct_btm, F_sfc_pls, F_btm_net, F_sfc_net, &
    !$acc   energy_sum, albedo, F_abs_sum, flx_sum, flx_slrd_val, flx_slri_val, &
    !$acc   albsfc_val, flx_sum_wgt, flx_sum_wgt_vis, &
    !$acc   mu_not, mu_one, SnowWaterEquivMin, diam_ice, &
    !$acc   fs_sphd, fs_hex, fs_hex0, fs_koch, AR_tmp, &
    !$acc   g_Cg_intp, gg_F07_intp, g_ice_F07, &
    !$acc   enh_omg_bcint, bcint_dd, bcint_dd2, bcint_f, &
    !$acc   enh_omg_bcint_intp, enh_omg_bcint_intp2, wvl_doint, &
    !$acc   enh_omg_dstint, enh_omg_dstint_intp, enh_omg_dstint_intp2, tot_dst_snw_conc, &
    !$acc   ts, ws, gs, extins, alp, gam, amg, apg, ue, refk, refkp1, refkm1, &
    !$acc   tdrrdir, tdndif, taus, omgs, asys, lm, mu, ne, &
    !$acc   R1, R2, T1, T2, Rf_dir_a, Tf_dir_a, Rf_dif_a, Rf_dif_b, Tf_dif_a, Tf_dif_b, &
    !$acc   gwt, swt, trn, rdr, tdr, smr, smt, exp_min, &
    !$acc   sza_c1, sza_c0, sza_factor, flx_sza_adjust)
    do JJ = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do II = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! Skip dark points (no sunlight)
    if (CosSolarZenithAngle(II,JJ) <= 0.0) cycle

    ! determin band start and end index
    if (NumSnicarRadBand == 5)   nir_bnd_bgn = 2
    if (NumSnicarRadBand == 480) nir_bnd_bgn = 51
    nir_bnd_end = NumSnicarRadBand

    ! initialize for adding-doubling solver parameter
    difgauspt(1:8) = & ! gaussian angles (radians)
                      (/ 0.9894009,  0.9445750, &
                         0.8656312,  0.7554044, &
                         0.6178762,  0.4580168, &
                         0.2816036,  0.0950125/)
    difgauswt(1:8) = & ! gaussian weights
                      (/ 0.0271525,  0.0622535, &
                         0.0951585,  0.1246290, &
                         0.1495960,  0.1691565, &
                         0.1826034,  0.1894506/)

    ! initialize for nonspherical snow grains
    !$acc loop seq
    do k = -NSNOW_MAX+1, 0
       sno_shp(k) = OptSnicarSnowShape
       sno_fs(k)  = 0.0
       sno_AR(k)  = 0.0
    enddo

    ! Table 3 of He et al 2017 JC
    g_wvl(1:8)    = (/ 0.25, 0.70, 1.41, 1.90, &
                       2.50, 3.50, 4.00, 5.00 /)
    !$acc loop seq
    do k = 1, 7
       g_wvl_ct(k) = g_wvl(k+1) / 2.0_kind_noahmp + g_wvl(k) / 2.0_kind_noahmp
    enddo
    g_b0(1:7)     = (/  9.76029E-1,  9.67798E-1,  1.00111, 1.00224,        &
                        9.64295E-1,  9.97475E-1,  9.97475E-1 /)
    g_b1(1:7)     = (/  5.21042E-1,  4.96181E-1,  1.83711E-1,  1.37082E-1, &
                        5.50598E-2,  8.48743E-2,  8.48743E-2 /)
    g_b2(1:7)     = (/ -2.66792E-4,  1.14088E-3,  2.37011E-4, -2.35905E-4, &
                        8.40449E-4, -4.71484E-4, -4.71484E-4 /)

    ! Tables 1 & 2 and Eqs. 3.1-3.4 from Fu, 2007 JAS
    g_F07_c2(1:7) = (/  1.349959E-1,  1.115697E-1,  9.853958E-2,  5.557793E-2, &
                       -1.233493E-1,  0.0        ,  0.0         /)
    g_F07_c1(1:7) = (/ -3.987320E-1, -3.723287E-1, -3.924784E-1, -3.259404E-1, &
                        4.429054E-2, -1.726586E-1, -1.726586E-1 /)
    g_F07_c0(1:7) = (/  7.938904E-1,  8.030084E-1,  8.513932E-1,  8.692241E-1, &
                        7.085850E-1,  6.412701E-1,  6.412701E-1 /)
    g_F07_p2(1:7) = (/  3.165543E-3,  2.014810E-3,  1.780838E-3,  6.987734E-4, &
                       -1.882932E-2, -2.277872E-2, -2.277872E-2 /)
    g_F07_p1(1:7) = (/  1.140557E-1,  1.143152E-1,  1.143814E-1,  1.071238E-1, &
                        1.353873E-1,  1.914431E-1,  1.914431E-1 /)
    g_F07_p0(1:7) = (/  5.292852E-1,  5.425909E-1,  5.601598E-1,  6.023407E-1, &
                        6.473899E-1,  4.634944E-1,  4.634944E-1 /)

    ! initialize for BC-snow internal mixing
    ! Eq. 8b & Table 4 in He et al., 2017 J. Climate (wavelength>1.2um, no BC-snow int mixing effect)
    bcint_wvl(1:17) = (/ 0.20, 0.25, 0.30, 0.33, 0.36, 0.40, 0.44, 0.48,       &
                         0.52, 0.57, 0.64, 0.69, 0.75, 0.78, 0.87, 1.0, 1.2 /)
    !$acc loop seq
    do k = 1, 16
       bcint_wvl_ct(k) = bcint_wvl(k+1) / 2.0_kind_noahmp + bcint_wvl(k) / 2.0_kind_noahmp
    enddo
    bcint_d0(1:16)  = (/ 2.48045   , 4.70305   , 4.68619   , 4.67369   , 4.65040   , &
                         2.40364   , 7.95408E-1, 2.92745E-1, 8.63396E-2, 2.76299E-2, &
                         1.40864E-2, 8.65705E-3, 6.12971E-3, 4.45697E-3, 3.06648E-2, &
                         7.96544E-1 /)
    bcint_d1(1:16)  = (/ 9.77209E-1, 9.73317E-1, 9.79650E-1, 9.84579E-1, 9.93537E-1, &
                         9.95955E-1, 9.95218E-1, 9.74284E-1, 9.81193E-1, 9.81239E-1, &
                         9.55515E-1, 9.10491E-1, 8.74196E-1, 8.27238E-1, 4.82870E-1, &
                         4.36649E-2 /)
    bcint_d2(1:16)  = (/ 3.95960E-1, 2.04820E-1, 2.07410E-1, 2.09390E-1, 2.13030E-1, &
                         4.18570E-1, 1.29682   , 3.75514   , 1.27372E+1, 3.93293E+1, &
                         8.78918E+1, 1.86969E+2, 3.45600E+2, 7.08637E+2, 1.41067E+3, &
                         2.57288E+2 /)
    ! Eq. 1a,1b and Table S1 in He et al. 2018 GRL
    bcint_m(1:3)    = (/ -0.8724, -0.1866, -0.0046 /)
    bcint_n(1:3)    = (/ -0.0072, -0.1918, -0.5177 /)

    ! initialize for dust-snow internal mixing
    ! Eq. 1 and Table 1 in He et al. 2019 JAMES (wavelength>1.2um, no dust-snow int mixing effect)
    dstint_wvl(1:7) = (/ 0.2, 0.2632, 0.3448, 0.4415, 0.625, 0.7782, 1.2422/)
    !$acc loop seq
    do k = 1, 6
       dstint_wvl_ct(k) = dstint_wvl(k+1) / 2.0_kind_noahmp + dstint_wvl(k) / 2.0_kind_noahmp
    enddo
    dstint_a1(1:6) = (/ -2.1307E+1, -1.5815E+1, -9.2880   , 1.1115   , 1.0307   , 1.0185    /)
    dstint_a2(1:6) = (/  1.1746E+2,  9.3241E+1,  4.0605E+1, 3.7389E-1, 1.4800E-2, 2.8921E-4 /)
    dstint_a3(1:6) = (/  9.9701E-1,  9.9781E-1,  9.9848E-1, 1.0035   , 1.0024   , 1.0356    /)

    ! SNICAR snow band center wavelength (um)
    wvl_ct5(1:5)  = (/ 0.5, 0.85, 1.1, 1.35, 3.25 /)  ! 5-band
    ! wvl_ct480 replaced by inline computation: (0.205_kind_noahmp + 0.01_kind_noahmp * (LoopInd - 1))

    ! Zero absorbed radiative fluxes:
    !$acc loop seq
    do LoopInd=-NSNOW_MAX+1,1,1
       !$acc loop seq
       do k=1,NumSnicarRadBand
          flx_abs_lcl_4d(II,LoopInd,k,JJ) = 0.0
       enddo
    enddo

    ! set SWE (mm) threshold for precision
    if (NumSnicarRadBand == 480) then 
       SnowWaterEquivMin = 1.0e-1
    elseif (NumSnicarRadBand == 5) then
       SnowWaterEquivMin = 1.0e-2
    endif

    ! Qualifier for computing snow RT: 
    ! minimum amount of snow on ground. 
    ! Otherwise, set snow albedo to zero
    if (SnowWaterEquiv(II,JJ) >= SnowWaterEquivMin) then

       ! If there is snow, but zero snow layers, we must create a layer locally.
       ! This layer is presumed to have the fresh snow effective radius.
          if (NumSnowLayerNeg(II,JJ) > -1) then
             flg_nosnl         =  1
             snl_lcl           =  -1
             h2osno_ice_lcl(0) =  SnowWaterEquiv(II,JJ)
             h2osno_liq_lcl(0) =  0.0
             snw_rds_lcl(0)    =  nint(SnowRadius(II,0,JJ))
          else
             flg_nosnl         =  0
             snl_lcl           =  NumSnowLayerNeg(II,JJ)
             !$acc loop seq
             do k = -NSNOW_MAX+1, 0
                h2osno_liq_lcl(k) = SnowLiqWater(II,k,JJ)
                h2osno_ice_lcl(k) = SnowIce(II,k,JJ)
                snw_rds_lcl(k)    = nint(SnowRadius(II,k,JJ))
             enddo
          endif

          SnowLayerBottom = 0
          SnowLayerTop = snl_lcl + 1

       ! Set local aerosol array
       if (FlagSnicarUseAerosol .eqv. .true.) then
          !$acc loop seq
          do k = -NSNOW_MAX+1, 0
             mss_cnc_aer_lcl(k,1) = MassConcBChydrophi(II,k,JJ)
             mss_cnc_aer_lcl(k,2) = MassConcBChydropho(II,k,JJ)
          enddo
          if (FlagSnicarUseOC .eqv. .true.) then
             !$acc loop seq
             do k = -NSNOW_MAX+1, 0
                mss_cnc_aer_lcl(k,3) = MassConcOChydrophi(II,k,JJ)
                mss_cnc_aer_lcl(k,4) = MassConcOChydropho(II,k,JJ)
             enddo
          else
             !$acc loop seq
             do k = -NSNOW_MAX+1, 0
                mss_cnc_aer_lcl(k,3) = 0.0
                mss_cnc_aer_lcl(k,4) = 0.0
             enddo
          endif
          !$acc loop seq
          do k = -NSNOW_MAX+1, 0
             mss_cnc_aer_lcl(k,5) = MassConcDust1(II,k,JJ)
             mss_cnc_aer_lcl(k,6) = MassConcDust2(II,k,JJ)
             mss_cnc_aer_lcl(k,7) = MassConcDust3(II,k,JJ)
             mss_cnc_aer_lcl(k,8) = MassConcDust4(II,k,JJ)
             mss_cnc_aer_lcl(k,9) = MassConcDust5(II,k,JJ)
          enddo
       else
          !$acc loop seq
          do k = -NSNOW_MAX+1, 0
             !$acc loop seq
             do j = 1, NAER
                mss_cnc_aer_lcl(k,j) = 0.0
             enddo
          enddo
       endif

       ! albsfc_val computed inline per band inside band loop

       ! Clamp snow grain size to valid range (GPU-safe)
       !$acc loop seq
       do i=SnowLayerTop,SnowLayerBottom,1
          snw_rds_lcl(i) = max(snw_rds_min_tbl, min(snw_rds_max_tbl, snw_rds_lcl(i)))
       enddo


       ! Incident flux weighting parameters
       !  - sum of all VIS bands must equal 1
       !  - sum of all NIR bands must equal 1
       !
       ! Spectral bands (5-band case)
       !  Band 1: 0.3-0.7um (VIS)
       !  Band 2: 0.7-1.0um (NIR)
       !  Band 3: 1.0-1.2um (NIR)
       !  Band 4: 1.2-1.5um (NIR)
       !  Band 5: 1.5-5.0um (NIR)
       !
       ! Hyperspectral (10-nm) bands (480-band case)
       ! Bands 1~50  : 0.2-0.7um (VIS)
       ! Bands 51~480: 0.7~5.0um (NIR)
       !
       ! The following weights are appropriate for surface-incident flux in a mid-latitude winter atmosphere
       !
       ! 3-band weights
       if (NumSnicarRadBand == 3) then
          ! Direct:
          if (FlagSwRadType == 1) then
             flx_wgt_3d(II,1,JJ) = 1.0
             flx_wgt_3d(II,2,JJ) = 0.66628670195247
             flx_wgt_3d(II,3,JJ) = 0.33371329804753
          ! Diffuse:
          elseif (FlagSwRadType == 2) then
             flx_wgt_3d(II,1,JJ) = 1.0
             flx_wgt_3d(II,2,JJ) = 0.77887652162877
             flx_wgt_3d(II,3,JJ) = 0.22112347837123
          endif
       else   ! works for both 5-band & 480-band, flux weights directly read from input data, cenlin
          ! Direct:
          if (FlagSwRadType == 1) then
             !$acc loop seq
             do k = 1, NumSnicarRadBand
                flx_wgt_3d(II,k,JJ) = RadSwWgtDir(II,k,JJ)
             enddo  ! VIS or NIR band sum is already normalized to 1.0 in input data
          ! Diffuse:
          elseif (FlagSwRadType == 2) then
             !$acc loop seq
             do k = 1, NumSnicarRadBand
                flx_wgt_3d(II,k,JJ) = RadSwWgtDif(II,k,JJ)
             enddo  ! VIS or NIR band sum is already normalized to 1.0 in input data
          endif
       endif

       exp_min = exp(-argmax)

       ! Loop over snow spectral bands
       !$acc loop seq
       do LoopInd = 1,NumSnicarRadBand

          ! Compute surface albedo for this band (replaces albsfc_lcl array)
          if (IndicatorIceSfc(II,JJ) == 0) then
             if (FlagSwRadType == 1) then
                if (LoopInd < nir_bnd_bgn) then
                   albsfc_val = AlbedoSoilDir(II,1,JJ)
                else
                   albsfc_val = AlbedoSoilDir(II,2,JJ)
                endif
             else
                if (LoopInd < nir_bnd_bgn) then
                   albsfc_val = AlbedoSoilDif(II,1,JJ)
                else
                   albsfc_val = AlbedoSoilDif(II,2,JJ)
                endif
             endif
          elseif (IndicatorIceSfc(II,JJ) == -1) then
             if (LoopInd < nir_bnd_bgn) then
                albsfc_val = AlbedoLandIce(II,1,JJ)
             else
                albsfc_val = AlbedoLandIce(II,2,JJ)
             endif
          endif

          ! Toon et al 2-stream
          if (OptSnicarRTSolver == 1) then
             mu_not = CosSolarZenithAngle(II,JJ)    ! must set here, because of error handling

          ! Adding-doubling 2-stream
          elseif (OptSnicarRTSolver == 2) then
             ! flg_dover is not used since this algorithm is stable for mu_not > 0.01
             ! mu_not is cosine solar zenith angle above the fresnel level; make
             ! sure mu_not is large enough for stable and meaningful radiation
             ! solution: .01 is like sun just touching horizon with its lower edge
             ! equivalent to mu0 in sea-ice shortwave model ice_shortwave.F90
             mu_not = max(CosSolarZenithAngle(II,JJ), cp01)
          endif

          flg_dover = 1    ! default is to redo
          err_idx   = 0    ! number of times through loop

          !$acc loop seq
          do while (flg_dover > 0)

             ! for Toon et al 2-stream solver:
             if (OptSnicarRTSolver == 1) then
                ! DEFAULT APPROXIMATIONS:
                !  VIS:       Delta-Eddington
                !  NIR (all): Delta-Hemispheric Mean
                !  WARNING:   DO NOT USE DELTA-EDDINGTON FOR NIR DIFFUSE - this sometimes results in negative albedo
                !  
                ! ERROR CONDITIONS:
                !  Conditions which cause "trip", resulting in redo of RT approximation:
                !   1. negative absorbed flux
                !   2. total absorbed flux greater than incident flux
                !   3. negative albedo
                !   NOTE: These errors have only been encountered in spectral bands 4 and 5
                !
                ! ERROR HANDLING
                !  1st error (flg_dover=2): switch approximation (Edd->HM or HM->Edd)
                !  2nd error (flg_dover=3): change zenith angle by 0.02 (this happens about 1 in 10^6 cases)
                !  3rd error (flg_dover=4): switch approximation with new zenith
                !  Subsequent errors: repeatedly change zenith and approximations...

                if (LoopInd < nir_bnd_bgn) then !VIS

                   if (flg_dover == 2) then
                      APRX_TYP = 3
                   elseif (flg_dover == 3) then
                      APRX_TYP = 1
                      if (CosSolarZenithAngle(II,JJ) > 0.5) then
                         mu_not = mu_not - 0.02
                      else
                         mu_not = mu_not + 0.02
                      endif
                   elseif (flg_dover == 4) then
                      APRX_TYP = 3
                   else
                      APRX_TYP = 1
                   endif

                else  ! NIR

                   if (flg_dover == 2) then
                      APRX_TYP = 1
                   elseif (flg_dover == 3) then
                      APRX_TYP = 3
                      if (CosSolarZenithAngle(II,JJ) > 0.5) then
                         mu_not = mu_not - 0.02
                      else
                         mu_not = mu_not + 0.02
                      endif
                   elseif (flg_dover == 4) then
                      APRX_TYP = 1
                   else
                      APRX_TYP = 3
                   endif

                endif ! end if < nir_bnd_bgn

             endif ! end if OptSnicarRTSolver == 1

             ! Set direct or diffuse incident irradiance to 1
             ! (This has to be within the bnd loop because mu_not is adjusted in rare cases)
             if (FlagSwRadType == 1) then
                flx_slrd_val = 1.0/(mu_not*ConstPI) ! this corresponds to incident irradiance of 1.0
                flx_slri_val = 0.0
             else
                flx_slrd_val = 0.0
                flx_slri_val = 1.0
             endif

             ! Pre-emptive error handling: aerosols can reap havoc on these absorptive bands.
             ! Since extremely high soot concentrations have a negligible effect on these bands, zero them.
             if ( (NumSnicarRadBand == 5).and.((LoopInd == 5).or.(LoopInd == 4)) ) then
                !$acc loop seq
                do k = -NSNOW_MAX+1, 0
                   !$acc loop seq
                   do j = 1, NAER
                      mss_cnc_aer_lcl(k,j) = 0.0
                   enddo
                enddo
             endif

             if ( (NumSnicarRadBand == 3).and.(LoopInd == 3) ) then
                !$acc loop seq
                do k = -NSNOW_MAX+1, 0
                   !$acc loop seq
                   do j = 1, NAER
                      mss_cnc_aer_lcl(k,j) = 0.0
                   enddo
                enddo
             endif

             if ( (NumSnicarRadBand == 480).and.(LoopInd > 100) ) then ! >1.2um 
                !$acc loop seq
                do k = -NSNOW_MAX+1, 0
                   !$acc loop seq
                   do j = 1, NAER
                      mss_cnc_aer_lcl(k,j) = 0.0
                   enddo
                enddo
             endif

             !--------------------------- Start snow & aerosol optics --------------------------------
             ! Define local Mie parameters based on snow grain size and aerosol species retrieved from a lookup table.
             ! Spherical snow: single-scatter albedo, mass extinction coefficient, asymmetry factor
             if (FlagSwRadType == 1) then
                !$acc loop seq
                do i=SnowLayerTop,SnowLayerBottom,1
                   rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                   ! snow optical properties (direct radiation)
                   ss_alb_snw_lcl(i)      = SsAlbSnwRadDir(II,rds_idx,LoopInd,JJ)
                   ext_cff_mss_snw_lcl(i) = ExtCffMassSnwRadDir(II,rds_idx,LoopInd,JJ)
                   if (sno_shp(i) == 1) asm_prm_snw_lcl(i) = AsyPrmSnwRadDir(II,rds_idx,LoopInd,JJ)
                enddo
             elseif (FlagSwRadType == 2) then
                !$acc loop seq
                do i=SnowLayerTop,SnowLayerBottom,1
                   rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                   ! snow optical properties (diffuse radiation)
                   ss_alb_snw_lcl(i)      = SsAlbSnwRadDif(II,rds_idx,LoopInd,JJ)
                   ext_cff_mss_snw_lcl(i) = ExtCffMassSnwRadDif(II,rds_idx,LoopInd,JJ)
                   if (sno_shp(i) == 1) asm_prm_snw_lcl(i) = AsyPrmSnwRadDif(II,rds_idx,LoopInd,JJ)
                enddo
             endif


             ! Nonspherical snow: shape-dependent asymmetry factors
             !$acc loop seq
             do i=SnowLayerTop,SnowLayerBottom,1
             
                ! spheroid
                if (sno_shp(i) == 2) then
                   diam_ice = 2.0 * snw_rds_lcl(i)   ! unit: microns
                   if (sno_fs(i) == 0.0) then
                      fs_sphd = 0.929  ! default; He et al. (2017), Table 1
                   else
                      fs_sphd = sno_fs(i) ! user specified value
                   endif
                   fs_hex = 0.788      ! reference shape factor
                   if (sno_AR(i) == 0.0) then
                      AR_tmp = 0.5     ! default; He et al. (2017), Table 1
                   else
                      AR_tmp = sno_AR(i)  ! user specified value
                   endif
                   !$acc loop seq
                   do igb = 1,7
                      g_ice_Cg_tmp(igb) = g_b0(igb) * ((fs_sphd/fs_hex)**g_b1(igb)) * (diam_ice**g_b2(igb))   ! Eq.7, He et al. (2017)
                      gg_ice_F07_tmp(igb) = g_F07_c0(igb) + g_F07_c1(igb)*AR_tmp + g_F07_c2(igb)*(AR_tmp**2.0)  ! Eqn. 3.1 in Fu (2007)
                   enddo

                !hexagonal plate
                elseif (sno_shp(i) ==3) then
                   diam_ice = 2.0 * snw_rds_lcl(i)   ! unit: microns
                   if (sno_fs(i) == 0.0) then
                       fs_hex0 = 0.788  ! default; He et al. (2017), Table 1
                   else
                       fs_hex0 = sno_fs(i) ! user specified value
                   endif
                   fs_hex = 0.788      ! reference shape factor
                   if (sno_AR(i) == 0.0) then
                      AR_tmp = 2.5     ! default; He et al. (2017), Table 1
                   else
                      AR_tmp = sno_AR(i)  ! user specified value
                   endif
                   !$acc loop seq
                   do igb = 1,7
                      g_ice_Cg_tmp(igb) = g_b0(igb) * ((fs_hex0/fs_hex)**g_b1(igb)) * (diam_ice**g_b2(igb))   ! Eq.7, He et al. (2017)
                      gg_ice_F07_tmp(igb) = g_F07_p0(igb)+g_F07_p1(igb)*LOG(AR_tmp)+g_F07_p2(igb)*((LOG(AR_tmp))**2.0) ! Eqn. 3.3 in Fu (2007)
                   enddo

                 ! Koch snowflake
                elseif (sno_shp(i) == 4) then
                   diam_ice = 2.0 * snw_rds_lcl(i) / 0.544  ! unit: microns
                   if (sno_fs(i) == 0.0) then
                      fs_koch = 0.712  ! default; He et al. (2017), Table 1
                   else
                      fs_koch = sno_fs(i) ! user specified value
                   endif
                   fs_hex = 0.788      ! reference shape factor
                   if (sno_AR(i) == 0.0) then
                      AR_tmp = 2.5     ! default; He et al. (2017), Table 1
                   else
                      AR_tmp = sno_AR(i)  ! user specified value
                   endif
                   !$acc loop seq
                   do igb = 1,7
                      g_ice_Cg_tmp(igb) = g_b0(igb) * ((fs_koch/fs_hex)**g_b1(igb)) * (diam_ice**g_b2(igb))   ! Eq.7, He et al. (2017)
                      gg_ice_F07_tmp(igb) = g_F07_p0(igb)+g_F07_p1(igb)*LOG(AR_tmp)+g_F07_p2(igb)*((LOG(AR_tmp))**2.0) ! Eqn. 3.3 in Fu (2007)
                   enddo

                endif !snowshape

                ! compute nonspherical snow asymmetry factor
                if (sno_shp(i) > 1) then
                   ! 7 wavelength bands for g_ice to be interpolated into targeted SNICAR bands here
                   ! use the piecewise linear interpolation subroutine created at the end of this module
                   ! tests showed the piecewise linear interpolation has similar results as pchip interpolation
                   if (NumSnicarRadBand == 5) then
                      call PiecewiseLinearInterp1d(7,g_wvl_ct,g_ice_Cg_tmp,wvl_ct5(LoopInd),g_Cg_intp)
                      call PiecewiseLinearInterp1d(7,g_wvl_ct,gg_ice_F07_tmp,wvl_ct5(LoopInd),gg_F07_intp)
                   endif
                   if (NumSnicarRadBand == 480) then
                      call PiecewiseLinearInterp1d(7,g_wvl_ct,g_ice_Cg_tmp,(0.205_kind_noahmp + 0.01_kind_noahmp * (LoopInd - 1)),g_Cg_intp)
                      call PiecewiseLinearInterp1d(7,g_wvl_ct,gg_ice_F07_tmp,(0.205_kind_noahmp + 0.01_kind_noahmp * (LoopInd - 1)),gg_F07_intp)
                   endif
                   g_ice_F07 = gg_F07_intp + (1.0 - gg_F07_intp) / ss_alb_snw_lcl(i) / 2.0  ! Eq.2.2 in Fu (2007)
                   asm_prm_snw_lcl(i) = g_ice_F07 * g_Cg_intp     ! Eq.6, He et al. (2017)
                endif

                if (asm_prm_snw_lcl(i) > 0.99) asm_prm_snw_lcl(i) = 0.99 !avoid unreasonable values (rarely occur in large-size spheroid cases)

             enddo !snow layer

             ! aerosol species 2 optical properties, hydrophobic BC
             ss_alb_aer_lcl(2)        = SsAlbBCpho(II,LoopInd,JJ)
             asm_prm_aer_lcl(2)       = AsyPrmBCpho(II,LoopInd,JJ)
             ext_cff_mss_aer_lcl(2)   = ExtCffMassBCpho(II,LoopInd,JJ)

             ! aerosol species 3 optical properties, hydrophilic OC
             ss_alb_aer_lcl(3)        = SsAlbOCphi(II,LoopInd,JJ)
             asm_prm_aer_lcl(3)       = AsyPrmOCphi(II,LoopInd,JJ)
             ext_cff_mss_aer_lcl(3)   = ExtCffMassOCphi(II,LoopInd,JJ)

             ! aerosol species 4 optical properties, hydrophobic OC
             ss_alb_aer_lcl(4)        = SsAlbOCpho(II,LoopInd,JJ)
             asm_prm_aer_lcl(4)       = AsyPrmOCpho(II,LoopInd,JJ)
             ext_cff_mss_aer_lcl(4)   = ExtCffMassOCpho(II,LoopInd,JJ)

             ! 1. snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
             ! 2. optical Depths (tau_snw, tau_aer)
             ! 3. weighted Mie properties (tau, omega, g)

             ! Weighted Mie parameters of each layer
             !$acc loop seq
             do i=SnowLayerTop,SnowLayerBottom,1

                ! Optics for BC/dust-snow external mixing:
                ! aerosol species 1 optical properties, hydrophilic BC
                ss_alb_aer_lcl(1)        = SsAlbBCphi(II,LoopInd,JJ)
                asm_prm_aer_lcl(1)       = AsyPrmBCphi(II,LoopInd,JJ)
                ext_cff_mss_aer_lcl(1)   = ExtCffMassBCphi(II,LoopInd,JJ)
                ! aerosol species 5 optical properties, dust size1
                ss_alb_aer_lcl(5)      = SsAlbDustB1(II,LoopInd,JJ)
                asm_prm_aer_lcl(5)     = AsyPrmDustB1(II,LoopInd,JJ)
                ext_cff_mss_aer_lcl(5) = ExtCffMassDustB1(II,LoopInd,JJ)
                ! aerosol species 6 optical properties, dust size2
                ss_alb_aer_lcl(6)      = SsAlbDustB2(II,LoopInd,JJ)
                asm_prm_aer_lcl(6)     = AsyPrmDustB2(II,LoopInd,JJ)
                ext_cff_mss_aer_lcl(6) = ExtCffMassDustB2(II,LoopInd,JJ)
                ! aerosol species 7 optical properties, dust size3
                ss_alb_aer_lcl(7)      = SsAlbDustB3(II,LoopInd,JJ)
                asm_prm_aer_lcl(7)     = AsyPrmDustB3(II,LoopInd,JJ)
                ext_cff_mss_aer_lcl(7) = ExtCffMassDustB3(II,LoopInd,JJ)
                ! aerosol species 8 optical properties, dust size4
                ss_alb_aer_lcl(8)      = SsAlbDustB4(II,LoopInd,JJ)
                asm_prm_aer_lcl(8)     = AsyPrmDustB4(II,LoopInd,JJ)
                ext_cff_mss_aer_lcl(8) = ExtCffMassDustB4(II,LoopInd,JJ)
                ! aerosol species 9 optical properties, dust size5
                ss_alb_aer_lcl(9)      = SsAlbDustB5(II,LoopInd,JJ)
                asm_prm_aer_lcl(9)     = AsyPrmDustB5(II,LoopInd,JJ)
                ext_cff_mss_aer_lcl(9) = ExtCffMassDustB5(II,LoopInd,JJ)

                ! Start BC/dust-snow internal mixing for wavelength<=1.2um
                if (NumSnicarRadBand == 5)   wvl_doint = wvl_ct5(LoopInd)
                if (NumSnicarRadBand == 480) wvl_doint = (0.205_kind_noahmp + 0.01_kind_noahmp * (LoopInd - 1))
                
                if (wvl_doint <= 1.2) then
                   ! BC-snow internal mixing applied to hydrophilic BC if activated
                   ! BC-snow internal mixing primarily affect snow single-scattering albedo
                   if ( FlagSnicarSnowBCIntmix .and. (mss_cnc_aer_lcl(i,1) > 0.0) ) then
                      ! result from Eq.8b in He et al.(2017) is based on BC Re=0.1um &
                      ! MAC=6.81 m2/g (@550 nm) & BC density=1.7g/cm3.
                      ! To be consistent with Bond et al. 2006 recommeded value (BC MAC=7.5 m2/g @550nm)
                      ! we made adjustments on BC size & density as follows to get MAC=7.5m2/g:
                      ! (1) We use BC Re=0.045um [geometric mean diameter=0.06um (Dentener et al.2006, 
                      ! Yu and Luo,2009) & geometric std=1.5 (Flanner et al.2007;Aoki et al., 2011)].
                      ! (2) We tune BC density from 1.7 to 1.49 g/cm3 (Aoki et al., 2011).
                      ! These adjustments also lead to consistent results with Flanner et al. 2012 (ACP) lookup table
                      ! for BC-snow internal mixing enhancement in albedo reduction (He et al. 2018 ACP)

                      !$acc loop seq
                      do ibb=1,16

                         enh_omg_bcint_tmp(ibb) = bcint_d0(ibb) * &
                         ( (mss_cnc_aer_lcl(i,1)*1.0E9*1.7/den_bc + bcint_d2(ibb)) **bcint_d1(ibb) )
                         ! adjust enhancment factor for BC effective size from 0.1um to Re_bc (He et al. 2018 GRL Eqs.1a,1b)
                         if (ibb < 3) then ! near-UV
                            bcint_dd  = (Re_bc/0.05)**bcint_m(1)
                            bcint_dd2 = (0.1/0.05)**bcint_m(1)
                            bcint_f  = (Re_bc/0.1)**bcint_n(1)
                         endif
                         if ( (ibb >= 3) .and. (ibb <= 11) ) then ! visible
                            bcint_dd  = (Re_bc/0.05)**bcint_m(2)
                            bcint_dd2 = (0.1/0.05)**bcint_m(2)
                            bcint_f  = (Re_bc/0.1)**bcint_n(2)
                         endif
                         if ( ibb > 11 ) then ! NIR
                            bcint_dd  = (Re_bc/0.05)**bcint_m(3)
                            bcint_dd2 = (0.1/0.05)**bcint_m(3)
                            bcint_f  = (Re_bc/0.1)**bcint_n(3)
                         endif
                         enh_omg_bcint_tmp2(ibb)=LOG10(max(1.0,bcint_dd*((enh_omg_bcint_tmp(ibb)/bcint_dd2)**bcint_f)))
                      enddo
                      ! piecewise linear interpolate into targeted SNICAR bands in a logscale space
                      call PiecewiseLinearInterp1d(16,bcint_wvl_ct,enh_omg_bcint_tmp2,wvl_doint,enh_omg_bcint_intp)
                      ! update snow single-scattering albedo
                      enh_omg_bcint_intp2 = 10.0 ** enh_omg_bcint_intp
                      enh_omg_bcint_intp2 = min(1.0E5, max(enh_omg_bcint_intp2,1.0)) ! constrain enhancement to a reasonable range
                      ss_alb_snw_lcl(i)   = 1.0 - (1.0 - ss_alb_snw_lcl(i)) * enh_omg_bcint_intp2
                      ss_alb_snw_lcl(i)   = max(0.5, min(ss_alb_snw_lcl(i),1.0))
                      ! reset hydrophilic BC property to 0 since it is accounted by updated snow ss_alb above
                      ss_alb_aer_lcl(1)       = 0.0
                      asm_prm_aer_lcl(1)      = 0.0
                      ext_cff_mss_aer_lcl(1)  = 0.0

                   endif ! end if BC-snow mixing type

                   ! Dust-snow internal mixing applied to all size bins if activated
                   ! Dust-snow internal mixing primarily affect snow single-scattering albedo
                   ! default optics of externally mixed dust at 4 size bins based on effective
                   ! radius of 1.38um and sigma=2.0 with truncation to each size bin (Flanner et al. 2021 GMD)
                   ! parameterized dust-snow int mix results based on effective radius of 1.1um and sigma=2.0
                   ! from (He et al. 2019 JAMES). Thus, the parameterization can be approximately applied to
                   ! all dust size bins here.

                   tot_dst_snw_conc = (mss_cnc_aer_lcl(i,5) + mss_cnc_aer_lcl(i,6) + &
                                       mss_cnc_aer_lcl(i,7) + mss_cnc_aer_lcl(i,8) + &
                                       mss_cnc_aer_lcl(i,9)) * 1.0E6 !kg/kg->ppm

                   if ( FlagSnicarSnowDustIntmix .and. (tot_dst_snw_conc > 0.0) ) then
                      !$acc loop seq
                      do idb=1,6
                         enh_omg_dstint_tmp(idb) = dstint_a1(idb)+dstint_a2(idb)*(tot_dst_snw_conc**dstint_a3(idb))
                         enh_omg_dstint_tmp2(idb) = LOG10(max(enh_omg_dstint_tmp(idb),1.0))
                      enddo

                      ! piecewise linear interpolate into targeted SNICAR bands in a logscale space
                      call PiecewiseLinearInterp1d(6,dstint_wvl_ct,enh_omg_dstint_tmp2,wvl_doint,enh_omg_dstint_intp)
                      ! update snow single-scattering albedo
                      enh_omg_dstint_intp2 = 10.0 ** enh_omg_dstint_intp
                      enh_omg_dstint_intp2 = min(1.0E5, max(enh_omg_dstint_intp2,1.0)) ! constrain enhancement to a reasonable range
                      ss_alb_snw_lcl(i) = 1.0 - (1.0 - ss_alb_snw_lcl(i)) * enh_omg_dstint_intp2
                      ss_alb_snw_lcl(i) = max(0.5, min(ss_alb_snw_lcl(i),1.0))

                      ! reset all dust optics to zero  since it is accounted by updated snow ss_alb above
                      !$acc loop seq
                      do k = 5, NAER
                         ss_alb_aer_lcl(k)      = 0.0
                         asm_prm_aer_lcl(k)     = 0.0
                         ext_cff_mss_aer_lcl(k) = 0.0
                      enddo
                   endif ! end if dust-snow internal mixing

                endif ! end if BC/dust-snow internal mixing (bands<1.2um)

                L_snw(i)   = h2osno_ice_lcl(i)+h2osno_liq_lcl(i)
                tau_snw(i) = L_snw(i)*ext_cff_mss_snw_lcl(i)

                !$acc loop seq
                do j=1,NAER
                   L_aer(i,j)   = L_snw(i)*mss_cnc_aer_lcl(i,j)
                   tau_aer(i,j) = L_aer(i,j)*ext_cff_mss_aer_lcl(j)
                enddo

                tau_sum   = 0.0
                omega_sum = 0.0
                g_sum     = 0.0

                !$acc loop seq
                do j=1,NAER
                   tau_sum    = tau_sum + tau_aer(i,j)
                   omega_sum  = omega_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j))
                   g_sum      = g_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j)*asm_prm_aer_lcl(j))
                enddo

                tau(i)    = tau_sum + tau_snw(i)
                omega(i)  = (1/tau(i))*(omega_sum+(ss_alb_snw_lcl(i)*tau_snw(i)))
                g(i)      = (1/(tau(i)*omega(i)))*(g_sum+ (asm_prm_snw_lcl(i)*ss_alb_snw_lcl(i)*tau_snw(i)))

             enddo ! end do snow layers

             ! DELTA transformations, requested
             !$acc loop seq
             do i=SnowLayerTop,SnowLayerBottom,1
                g_star(i)     = g(i)/(1+g(i))
                omega_star(i) = ((1-(g(i)**2))*omega(i)) / (1-(omega(i)*(g(i)**2)))
                tau_star(i)   = (1-(omega(i)*(g(i)**2)))*tau(i)
             enddo
             !--------------------------- End of snow & aerosol optics --------------------------------


             !--------------------------- Start Toon et al. RT solver  --------------------------------
             if (OptSnicarRTSolver == 1) then
                ! Total column optical depth:
                ! tau_clm(i) = total optical depth above the bottom of layer i
                tau_clm(SnowLayerTop) = 0.0

                !$acc loop seq
                do i=SnowLayerTop+1,SnowLayerBottom,1
                   tau_clm(i) = tau_clm(i-1)+tau_star(i-1)
                enddo

                ! Direct radiation at bottom of snowpack:
                F_direct_btm = albsfc_val*mu_not * &
                               exp(-(tau_clm(SnowLayerBottom)+tau_star(SnowLayerBottom))/mu_not)*ConstPI*flx_slrd_val

                ! Intermediates
                ! Gamma values are approximation-specific.

                ! Eddington
                if (APRX_TYP==1) then
                   !$acc loop seq
                   do i=SnowLayerTop,SnowLayerBottom,1
                      gamma1(i) = (7.0-(omega_star(i)*(4.0+(3.0*g_star(i)))))/4.0
                      gamma2(i) = -(1.0-(omega_star(i)*(4.0-(3.0*g_star(i)))))/4.0
                      gamma3(i) = (2.0-(3.0*g_star(i)*mu_not))/4.0
                      gamma4(i) = 1.0-gamma3(i)
                      mu_one    = 0.5
                   enddo

                ! Quadrature
                elseif (APRX_TYP==2) then
                   !$acc loop seq
                   do i=SnowLayerTop,SnowLayerBottom,1
                      gamma1(i) = (3.0**0.5)*(2.0-(omega_star(i)*(1.0+g_star(i))))/2.0
                      gamma2(i) = omega_star(i)*(3.0**0.5)*(1.0-g_star(i))/2.0
                      gamma3(i) = (1.0-((3.0**0.5)*g_star(i)*mu_not))/2.0
                      gamma4(i) = 1.0-gamma3(i)
                      mu_one    = 1.0/(3.0**0.5)
                   enddo

                ! Hemispheric Mean
                elseif (APRX_TYP==3) then
                   !$acc loop seq
                   do i=SnowLayerTop,SnowLayerBottom,1
                      gamma1(i) = 2.0 - (omega_star(i)*(1.0+g_star(i)))
                      gamma2(i) = omega_star(i)*(1.0-g_star(i))
                      gamma3(i) = (1.0-((3.0**0.5)*g_star(i)*mu_not))/2.0
                      gamma4(i) = 1.0-gamma3(i)
                      mu_one    = 0.5
                   enddo
                endif

                ! Intermediates for tri-diagonal solution
                !$acc loop seq
                do i=SnowLayerTop,SnowLayerBottom,1
                   lambda(i) = sqrt(abs((gamma1(i)**2) - (gamma2(i)**2)))
                   GAMMA(i)  = gamma2(i)/(gamma1(i)+lambda(i))

                   e1(i)     = 1+(GAMMA(i)*exp(-lambda(i)*tau_star(i)))
                   e2(i)     = 1-(GAMMA(i)*exp(-lambda(i)*tau_star(i)))
                   e3(i)     = GAMMA(i) + exp(-lambda(i)*tau_star(i))
                   e4(i)     = GAMMA(i) - exp(-lambda(i)*tau_star(i))

                enddo !Snow layer

                !$acc loop seq
                do i=SnowLayerTop,SnowLayerBottom,1
                   if (FlagSwRadType == 1) then
                      C_pls_btm(i) = (omega_star(i)*ConstPI*flx_slrd_val* &
                              exp(-(tau_clm(i)+tau_star(i))/mu_not)*   &
                              (((gamma1(i)-(1/mu_not))*gamma3(i))+     &
                              (gamma4(i)*gamma2(i))))/((lambda(i)**2)-(1/(mu_not**2)))
                      C_mns_btm(i) = (omega_star(i)*ConstPI*flx_slrd_val* &
                              exp(-(tau_clm(i)+tau_star(i))/mu_not)*   &
                              (((gamma1(i)+(1/mu_not))*gamma4(i))+     &
                              (gamma2(i)*gamma3(i))))/((lambda(i)**2)-(1/(mu_not**2)))
                      C_pls_top(i) = (omega_star(i)*ConstPI*flx_slrd_val* &
                              exp(-tau_clm(i)/mu_not)*(((gamma1(i)-(1/mu_not))* &
                              gamma3(i))+(gamma4(i)*gamma2(i))))/((lambda(i)**2)-(1/(mu_not**2)))
                      C_mns_top(i) = (omega_star(i)*ConstPI*flx_slrd_val* &
                              exp(-tau_clm(i)/mu_not)*(((gamma1(i)+(1/mu_not))* &
                              gamma4(i))+(gamma2(i)*gamma3(i))))/((lambda(i)**2)-(1/(mu_not**2)))

                   else
                      C_pls_btm(i) = 0.0
                      C_mns_btm(i) = 0.0
                      C_pls_top(i) = 0.0
                      C_mns_top(i) = 0.0
                   endif
                enddo !Snow layer

                ! Coefficients for tridiaganol matrix solution
                !$acc loop seq
                do i=2*snl_lcl+1,0,1
                   !Boundary values for i=1 and i=2*snl_lcl, specifics for i=odd and i=even    
                   if (i==(2*snl_lcl+1)) then
                      A(i) = 0.0
                      B(i) = e1(SnowLayerTop)
                      D(i) = -e2(SnowLayerTop)
                      E(i) = flx_slri_val-C_mns_top(SnowLayerTop)
                   elseif(i==0) then
                      A(i) = e1(SnowLayerBottom)-(albsfc_val*e3(SnowLayerBottom))
                      B(i) = e2(SnowLayerBottom)-(albsfc_val*e4(SnowLayerBottom))
                      D(i) = 0.0
                      E(i) = F_direct_btm-C_pls_btm(SnowLayerBottom)+(albsfc_val*C_mns_btm(SnowLayerBottom))
                   elseif(mod(i,2)==-1) then   ! If odd and i>=3 (n=1 for i=3)
                      n=floor(i/2.0)
                      A(i) = (e2(n)*e3(n))-(e4(n)*e1(n))
                      B(i) = (e1(n)*e1(n+1))-(e3(n)*e3(n+1))
                      D(i) = (e3(n)*e4(n+1))-(e1(n)*e2(n+1))
                      E(i) = (e3(n)*(C_pls_top(n+1)-C_pls_btm(n)))+(e1(n)*(C_mns_btm(n)-C_mns_top(n+1)))
                   elseif(mod(i,2)==0) then    ! If even and i<=2*snl_lcl
                      n=(i/2)
                      A(i) = (e2(n+1)*e1(n))-(e3(n)*e4(n+1))
                      B(i) = (e2(n)*e2(n+1))-(e4(n)*e4(n+1))
                      D(i) = (e1(n+1)*e4(n+1))-(e2(n+1)*e3(n+1))
                      E(i) = (e2(n+1)*(C_pls_top(n+1)-C_pls_btm(n)))+(e4(n+1)*(C_mns_top(n+1)-C_mns_btm(n)))
                   endif
                enddo

                AS(0) = A(0)/B(0)
                DS(0) = E(0)/B(0)

                !$acc loop seq
                do i=-1,(2*snl_lcl+1),-1
                   X(i)  = 1/(B(i)-(D(i)*AS(i+1)))
                   AS(i) = A(i)*X(i)
                   DS(i) = (E(i)-(D(i)*DS(i+1)))*X(i)
                enddo

                Y(2*snl_lcl+1) = DS(2*snl_lcl+1)
                !$acc loop seq
                do i=(2*snl_lcl+2),0,1
                   Y(i) = DS(i)-(AS(i)*Y(i-1))
                enddo

                ! Downward direct-beam and net flux (F_net) at the base of each layer:
                !$acc loop seq
                do i=SnowLayerTop,SnowLayerBottom,1
                   F_direct(i) = mu_not*ConstPI*flx_slrd_val*exp(-(tau_clm(i)+tau_star(i))/mu_not)
                   F_net(i)    = (Y(2*i-1)*(e1(i)-e3(i))) + (Y(2*i)*(e2(i)-e4(i))) + &
                                     C_pls_btm(i) - C_mns_btm(i) - F_direct(i)
                enddo

                ! Upward flux at snowpack top:
                F_sfc_pls = (Y(2*snl_lcl+1)*(exp(-lambda(SnowLayerTop)*tau_star(SnowLayerTop))+ &
                         GAMMA(SnowLayerTop))) + (Y(2*snl_lcl+2)*(exp(-lambda(SnowLayerTop)* &
                         tau_star(SnowLayerTop))-GAMMA(SnowLayerTop))) + C_pls_top(SnowLayerTop)

                ! Net flux at bottom = absorbed radiation by underlying surface:
                F_btm_net = -F_net(SnowLayerBottom)

                ! Bulk column albedo and surface net flux
                albedo    = F_sfc_pls/((mu_not*ConstPI*flx_slrd_val)+flx_slri_val)
                F_sfc_net = F_sfc_pls - ((mu_not*ConstPI*flx_slrd_val)+flx_slri_val)

                trip = 0
                ! Absorbed flux in each layer
                !$acc loop seq
                do i=SnowLayerTop,SnowLayerBottom,1
                   if(i==SnowLayerTop) then
                      F_abs(i) = F_net(i)-F_sfc_net
                   else
                      F_abs(i) = F_net(i)-F_net(i-1)
                   endif
                   flx_abs_lcl_4d(II,i,LoopInd,JJ) = F_abs(i)

                   ! ERROR check: negative absorption
                   if (flx_abs_lcl_4d(II,i,LoopInd,JJ) < -0.00001) then
                      trip = 1
                   endif
                enddo

                flx_abs_lcl_4d(II,1,LoopInd,JJ) = F_btm_net

                if (flg_nosnl == 1) then
                   ! If there are no snow layers (but still snow), all absorbed energy must be in top soil layer
                   !flx_abs_lcl(:,bnd_idx) = 0._r8
                   !flx_abs_lcl(1,bnd_idx) = F_abs(0) + F_btm_net

                   ! changed on 20070408:
                   ! OK to put absorbed energy in the fictitous snow layer because routine SurfaceRadiation
                   ! handles the case of no snow layers. Then, if a snow layer is addded between now and
                   ! SurfaceRadiation (called in CanopyHydrology), absorbed energy will be properly distributed.
                   flx_abs_lcl_4d(II,0,LoopInd,JJ) = F_abs(0)
                   flx_abs_lcl_4d(II,1,LoopInd,JJ) = F_btm_net
                endif

                !Underflow check (we've already tripped the error condition above)
                !$acc loop seq
                do i=SnowLayerTop,1,1
                   if (flx_abs_lcl_4d(II,i,LoopInd,JJ) < 0.0) then
                      flx_abs_lcl_4d(II,i,LoopInd,JJ) = 0.0
                   endif
                enddo

                F_abs_sum = 0.0
                !$acc loop seq
                do i=SnowLayerTop,SnowLayerBottom,1
                   F_abs_sum = F_abs_sum + F_abs(i)
                enddo

                !ERROR check: absorption greater than incident flux
                ! (should make condition more generic than "1._r8")
                if (F_abs_sum > 1.0) then
                   trip = 1
                endif

                !ERROR check:
                if ((albedo < 0.0).and.(trip==0)) then
                   trip = 1
                endif

                ! Set conditions for redoing RT calculation
                if ((trip == 1).and.(flg_dover == 1)) then
                   flg_dover = 2
                elseif ((trip == 1).and.(flg_dover == 2)) then
                   flg_dover = 3
                elseif ((trip == 1).and.(flg_dover == 3)) then
                   flg_dover = 4
                elseif((trip == 1).and.(flg_dover == 4).and.(err_idx < 20)) then
                   flg_dover = 3
                   err_idx = err_idx + 1
                elseif((trip == 1).and.(flg_dover == 4).and.(err_idx >= 20)) then
                   flg_dover = 0
#ifndef _OPENACC
                   write(*,*) "SNICAR ERROR: FOUND A WORMHOLE. STUCK IN INFINITE LOOP!"
                   write(*,*) "SNICAR STATS: L_snw(0)= ", L_snw(0)
                   write(*,*) "SNICAR STATS: snw_rds_lcl(0)= ", snw_rds_lcl(0)
                   write(*,*) "SNICAR STATS: h2osno= ", SnowWaterEquiv(II,JJ), " snl= ", snl_lcl
                   write(*,*) "SNICAR STATS: BCphi(0)= ", mss_cnc_aer_lcl(0,1)
                   write(*,*) "SNICAR STATS: BCpho(0)= ", mss_cnc_aer_lcl(0,2)
                   write(*,*) "SNICAR STATS: dust1(0)= ", mss_cnc_aer_lcl(0,5)
                   write(*,*) "SNICAR STATS: dust2(0)= ", mss_cnc_aer_lcl(0,6)
                   write(*,*) "SNICAR STATS: dust3(0)= ", mss_cnc_aer_lcl(0,7)
                   write(*,*) "SNICAR STATS: dust4(0)= ", mss_cnc_aer_lcl(0,8)
                   write(*,*) "SNICAR STATS: dust5(0)= ", mss_cnc_aer_lcl(0,9)
#endif
                else
                   flg_dover = 0
                endif

             endif ! end if OptSnicarRTSolver == 1


             !--------------------------- Start Adding-doubling RT solver  --------------------------------
             if (OptSnicarRTSolver == 2) then
                ! Given input vertical profiles of optical properties, evaluate the
                ! monochromatic Delta-Eddington adding-doubling solution

                ! trndir, trntdr, trndif, rupdir, rupdif, rdndif are variables at the layer interface,
                ! for snow with layers from snl_top to snl_btm there are snl_top to snl_btm+1 layer interface
                snl_btm_itf = SnowLayerBottom + 1

                ! initialization for layer interface
                !$acc loop seq
                do i = SnowLayerTop,snl_btm_itf,1
                   trndir(i) = c0
                   trntdr(i) = c0
                   trndif(i) = c0
                   rupdir(i) = c0
                   rupdif(i) = c0
                   rdndif(i) = c0
                enddo

                ! initialize top interface of top layer
                trndir(SnowLayerTop) = c1
                trntdr(SnowLayerTop) = c1
                trndif(SnowLayerTop) = c1
                rdndif(SnowLayerTop) = c0

                ! begin main level loop for snow layer interfaces except for the very bottom
                !$acc loop seq
                do i = SnowLayerTop,SnowLayerBottom,1

                   ! initialize all layer apparent optical properties to 0
                   rdir  (i) = c0
                   rdif_a(i) = c0
                   rdif_b(i) = c0
                   tdir  (i) = c0
                   tdif_a(i) = c0
                   tdif_b(i) = c0
                   trnlay(i) = c0

                   ! compute next layer Delta-eddington solution only if total transmission
                   ! of radiation to the interface just above the layer exceeds trmin.
                   if (trntdr(i) > trmin ) then

                      ! delta-transformed single-scattering properties of this layer
                      ts = tau_star(i)
                      ws = omega_star(i)
                      gs = g_star(i)

                      ! Delta-Eddington solution expressions, Eq. 50: Briegleb and Light 2007
                      lm = sqrt(c3*(c1-ws)*(c1 - ws*gs))
                      ue = c1p5*(c1 - ws*gs)/lm
                      extins = max(exp_min, exp(-lm*ts))
                      ne = ((ue+c1)*(ue+c1)/extins) - ((ue-c1)*(ue-c1)*extins)

                      ! first calculation of rdif, tdif using Delta-Eddington formulas
                      ! Eq.: Briegleb 1992; alpha and gamma for direct radiation
                      rdif_a(i) = (ue**2-c1)*(c1/extins - extins)/ne
                      tdif_a(i) = c4*ue/ne

                      ! evaluate rdir,tdir for direct beam
                      trnlay(i) = max(exp_min, exp(-ts/mu_not))


                      ! Delta-Eddington solution expressions
                      ! Eq. 50: Briegleb and Light 2007; alpha and gamma for direct radiation
                      if (c1 - lm*lm*mu_not*mu_not /= 0.0) then
                         alp = cp75*ws*mu_not*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu_not*mu_not))
                         gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu_not*mu_not)/(c1-lm*lm*mu_not*mu_not))
                      else
                         alp = 0.0
                         gam = 0.0
                      endif
                      apg = alp + gam
                      amg = alp - gam
                      rdir(i) = apg*rdif_a(i) +  amg*(tdif_a(i)*trnlay(i) - c1)
                      tdir(i) = apg*tdif_a(i) + (amg* rdif_a(i)-apg+c1)*trnlay(i)


                      ! recalculate rdif,tdif using direct angular integration over rdir,tdir,
                      ! since Delta-Eddington rdif formula is not well-behaved (it is usually
                      ! biased low and can even be negative); use ngmax angles and gaussian
                      ! integration for most accuracy:
                      R1 = rdif_a(i) ! use R1 as temporary
                      T1 = tdif_a(i) ! use T1 as temporary
                      swt = c0
                      smr = c0
                      smt = c0
                      ! gaussian angles for the AD integral

                      !$acc loop seq
                      do ng=1,ngmax
                         mu  = difgauspt(ng)
                         gwt = difgauswt(ng)
                         swt = swt + mu*gwt
                         trn = max(exp_min, exp(-ts/mu))
                         alp = cp75*ws*mu*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu*mu))
                         gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu*mu)/(c1-lm*lm*mu*mu))
                         apg = alp + gam
                         amg = alp - gam
                         rdr = apg*R1 + amg*T1*trn - amg
                         tdr = apg*T1 + amg*R1*trn - apg*trn + trn
                         smr = smr + mu*rdr*gwt
                         smt = smt + mu*tdr*gwt
                      enddo      ! ng

                      rdif_a(i) = smr/swt
                      tdif_a(i) = smt/swt

                      ! homogeneous layer
                      rdif_b(i) = rdif_a(i)
                      tdif_b(i) = tdif_a(i)

                   endif  ! trntdr(i) > trmin 

                   ! Calculate the solar beam transmission, total transmission, and
                   ! reflectivity for diffuse radiation from below at interface i,
                   ! the top of the current layer k:
                   !
                   !              layers       interface
                   !
                   !       ---------------------  i-1
                   !                i-1
                   !       ---------------------  i
                   !                 i
                   !       ---------------------


                   trndir(i+1) = trndir(i)*trnlay(i)            ! solar beam transmission from top
                   refkm1      = c1/(c1 - rdndif(i)*rdif_a(i))  ! interface multiple scattering for i-1
                   tdrrdir     = trndir(i)*rdir(i)              ! direct tran times layer direct ref
                   tdndif      = trntdr(i) - trndir(i)          ! total down diffuse = tot tran - direct tran
                   trntdr(i+1) = trndir(i)*tdir(i) + &          ! total transmission to direct beam for layers above
                                 (tdndif + tdrrdir*rdndif(i))*refkm1*tdif_a(i)
                   ! Eq. B4; Briegleb and Light 2007
                   rdndif(i+1) = rdif_b(i) + &                  ! reflectivity to diffuse radiation for layers above
                                 (tdif_b(i)*rdndif(i)*refkm1*tdif_a(i))
                   trndif(i+1) = trndif(i)*refkm1*tdif_a(i)     ! diffuse transmission to diffuse beam for layers above

                enddo  !snow layer


                ! compute reflectivity to direct and diffuse radiation for layers
                ! below by adding succesive layers starting from the underlying
                ! ground and working upwards:
                !
                !              layers       interface
                !
                !       ---------------------  i
                !                 i
                !       ---------------------  i+1
                !                i+1
                !       ---------------------

                ! set the underlying ground albedo == albedo of near-IR
                ! unless bnd_idx < nir_bnd_bgn, for visible
                if (IndicatorIceSfc(II,JJ) == 0) then
                   if (FlagSwRadType == 1) then
                      rupdir(snl_btm_itf) = AlbedoSoilDir(II,2,JJ)
                      rupdif(snl_btm_itf) = AlbedoSoilDir(II,2,JJ)
                      if (LoopInd < nir_bnd_bgn) then
                         rupdir(snl_btm_itf) = AlbedoSoilDir(II,1,JJ)
                         rupdif(snl_btm_itf) = AlbedoSoilDir(II,1,JJ)
                      endif
                   elseif (FlagSwRadType == 2) then
                      rupdir(snl_btm_itf) = AlbedoSoilDif(II,2,JJ)
                      rupdif(snl_btm_itf) = AlbedoSoilDif(II,2,JJ)
                      if (LoopInd < nir_bnd_bgn) then
                         rupdir(snl_btm_itf) = AlbedoSoilDif(II,1,JJ)
                         rupdif(snl_btm_itf) = AlbedoSoilDif(II,1,JJ)
                      endif
                   endif
                elseif (IndicatorIceSfc(II,JJ) == -1) then !land ice
                   rupdir(snl_btm_itf) = AlbedoLandIce(II,2,JJ)
                   rupdif(snl_btm_itf) = AlbedoLandIce(II,2,JJ)
                   if (LoopInd < nir_bnd_bgn) then
                      rupdir(snl_btm_itf) = AlbedoLandIce(II,1,JJ)
                      rupdif(snl_btm_itf) = AlbedoLandIce(II,1,JJ)
                   endif
                endif

                !$acc loop seq
                do i=SnowLayerBottom,SnowLayerTop,-1
                   ! interface scattering Eq. B5; Briegleb and Light 2007
                   refkp1 = c1/( c1 - rdif_b(i)*rupdif(i+1))
                   ! dir from top layer plus exp tran ref from lower layer, interface
                   ! scattered and tran thru top layer from below, plus diff tran ref
                   ! from lower layer with interface scattering tran thru top from below
                   rupdir(i) = rdir(i) &
                               + (        trnlay(i)  *rupdir(i+1) &
                               +  (tdir(i)-trnlay(i))*rupdif(i+1) ) * refkp1 * tdif_b(i)
                   ! dif from top layer from above, plus dif tran upwards reflected and
                   ! interface scattered which tran top from below
                   rupdif(i) = rdif_a(i) + tdif_a(i)*rupdif(i+1)*refkp1*tdif_b(i)
                enddo       ! i


                ! net flux (down-up) at each layer interface from the
                ! snow top (i = snl_top) to bottom interface above land (i = snl_btm_itf)
                ! the interface reflectivities and transmissivities required
                ! to evaluate interface fluxes are returned from solution_dEdd;
                ! now compute up and down fluxes for each interface, using the
                ! combined layer properties at each interface:
                !
                !              layers       interface
                !
                !       ---------------------  i
                !                 i
                !       ---------------------

 
                !$acc loop seq
                do i = SnowLayerTop, snl_btm_itf
                   ! interface scattering, Eq. 52; Briegleb and Light 2007
                   refk = c1/(c1 - rdndif(i)*rupdif(i))
                   ! dir tran ref from below times interface scattering, plus diff
                   ! tran and ref from below times interface scattering
                   ! fdirup(i) = (trndir(i)*rupdir(i) + &
                   !                 (trntdr(i)-trndir(i))  &
                   !                 *rupdif(i))*refk
                   ! dir tran plus total diff trans times interface scattering plus
                   ! dir tran with up dir ref and down dif ref times interface scattering
                   ! fdirdn(i) = trndir(i) + (trntdr(i) &
                   !               - trndir(i) + trndir(i)  &
                   !               *rupdir(i)*rdndif(i))*refk
                   ! diffuse tran ref from below times interface scattering
                   ! fdifup(i) = trndif(i)*rupdif(i)*refk
                   ! diffuse tran times interface scattering
                   ! fdifdn(i) = trndif(i)*refk

                   ! netflux, down - up
                   ! dfdir = fdirdn - fdirup
                   dfdir(i) = trndir(i) &
                            + (trntdr(i)-trndir(i)) * (c1 - rupdif(i)) * refk &
                            -  trndir(i)*rupdir(i)  * (c1 - rdndif(i)) * refk
                   if (dfdir(i) < puny) dfdir(i) = c0
                   ! dfdif = fdifdn - fdifup
                   dfdif(i) = trndif(i) * (c1 - rupdif(i)) * refk
                   if (dfdif(i) < puny) dfdif(i) = c0
                enddo  ! i


                ! SNICAR_AD_RT is called twice for direct and diffuse incident fluxes
                ! direct incident
                if (FlagSwRadType == 1) then
                   albedo = rupdir(SnowLayerTop)
                   !$acc loop seq
                   do k = -NSNOW_MAX+1, 1
                      dftmp(k) = dfdir(k)
                   enddo
                   refk   = c1/(c1 - rdndif(SnowLayerTop)*rupdif(SnowLayerTop))
                   F_sfc_pls = (trndir(SnowLayerTop)*rupdir(SnowLayerTop) + &
                                   (trntdr(SnowLayerTop)-trndir(SnowLayerTop))  &
                                   *rupdif(SnowLayerTop))*refk
                !diffuse incident
                else
                   albedo = rupdif(SnowLayerTop)
                   !$acc loop seq
                   do k = -NSNOW_MAX+1, 1
                      dftmp(k) = dfdif(k)
                   enddo
                   refk   = c1/(c1 - rdndif(SnowLayerTop)*rupdif(SnowLayerTop))
                   F_sfc_pls = trndif(SnowLayerTop)*rupdif(SnowLayerTop)*refk
                endif

                ! Absorbed flux in each layer
                !$acc loop seq
                do i=SnowLayerTop,SnowLayerBottom,1
                   F_abs(i) = dftmp(i)-dftmp(i+1)
                   flx_abs_lcl_4d(II,i,LoopInd,JJ) = F_abs(i)

                   ! ERROR check: negative absorption
                   if (flx_abs_lcl_4d(II,i,LoopInd,JJ) < -0.0001) then !original -0.00001, but not work for Koch snowflake
#ifndef _OPENACC
                      write (*,"(a,e13.6,i0,i0,i0,i0)") "SNICAR ERROR: negative absoption : ", &
                            flx_abs_lcl_4d(II,i,LoopInd,JJ),i,LoopInd,SnowLayerTop,SnowLayerBottom
                      write(*,*) "SNICAR_AD STATS: L_snw(0)= ", L_snw(0)
                      write(*,*) "SNICAR_AD STATS: snw_rds_lcl(0)= ", snw_rds_lcl(0)
                      write(*,*) "SNICAR_AD STATS: coszen= ",  CosSolarZenithAngle(II,JJ)
                      write(*,*) 'SNICAR_AD STATS: wavelength=', (0.205_kind_noahmp + 0.01_kind_noahmp * (LoopInd - 1))
                      write(*,*) "SNICAR_AD STATS: h2osno= ", SnowWaterEquiv(II,JJ), " snl= ", snl_lcl
                      write(*,*) "SNICAR_AD STATS: BCphi(0)= ", mss_cnc_aer_lcl(0,1)
                      write(*,*) "SNICAR_AD STATS: BCpho(0)= ", mss_cnc_aer_lcl(0,2)
                      write(*,*) "SNICAR_AD STATS: OCphi(0)= ", mss_cnc_aer_lcl(0,3)
                      write(*,*) "SNICAR_AD STATS: OCpho(0)= ", mss_cnc_aer_lcl(0,4)
                      write(*,*) "SNICAR_AD STATS: dust1(0)= ", mss_cnc_aer_lcl(0,5)
                      write(*,*) "SNICAR_AD STATS: dust2(0)= ", mss_cnc_aer_lcl(0,6)
                      write(*,*) "SNICAR_AD STATS: dust3(0)= ", mss_cnc_aer_lcl(0,7)
                      write(*,*) "SNICAR_AD STATS: dust4(0)= ", mss_cnc_aer_lcl(0,8)
                      write(*,*) "SNICAR_AD STATS: dust5(0)= ", mss_cnc_aer_lcl(0,9)
                      stop "ERROR in SNICAR absorption"
#endif
                    endif
                enddo

                ! absobed flux by the underlying ground
                F_btm_net = dftmp(snl_btm_itf)

                ! note here, snl_btm_itf = 1 by snow column set up in CLM
                flx_abs_lcl_4d(II,1,LoopInd,JJ) = F_btm_net

                if (flg_nosnl == 1) then
                   ! If there are no snow layers (but still snow), all absorbed energy must be in top soil layer
                   !flx_abs_lcl(:,LoopInd) = 0.0
                   !flx_abs_lcl(1,LoopInd) = F_abs(0) + F_btm_net

                   ! changed on 20070408:
                   ! OK to put absorbed energy in the fictitous snow layer because routine SurfaceRadiation
                   ! handles the case of no snow layers. Then, if a snow layer is addded between now and
                   ! SurfaceRadiation (called in CanopyHydrology), absorbed energy will be properly distributed.
                   flx_abs_lcl_4d(II,0,LoopInd,JJ) = F_abs(0)
                   flx_abs_lcl_4d(II,1,LoopInd,JJ) = F_btm_net
                endif

                !Underflow check (we've already tripped the error condition above)
                !$acc loop seq
                do i=SnowLayerTop,1,1
                   if (flx_abs_lcl_4d(II,i,LoopInd,JJ) < 0.0) then
                      flx_abs_lcl_4d(II,i,LoopInd,JJ) = 0.0
                   endif
                enddo

                F_abs_sum = 0.0
                !$acc loop seq
                do i=SnowLayerTop,SnowLayerBottom,1
                   F_abs_sum = F_abs_sum + F_abs(i)
                enddo

                ! no need to repeat calculations for adding-doubling solver
                flg_dover = 0

             endif ! end if OptSnicarRTSolver == 2
             !--------------------------- End of Adding-doubling RT solver  --------------------------------

          enddo !enddo while (flg_dover > 0)

          ! Energy conservation check:
          ! Incident direct+diffuse radiation equals (absorbed+bulk_transmitted+bulk_reflected)
          energy_sum = (mu_not*ConstPI*flx_slrd_val) + flx_slri_val - (F_abs_sum + F_btm_net + F_sfc_pls)

          if (abs(energy_sum) > 0.00001) then
#ifndef _OPENACC
              write (*,*) "SNICAR ERROR: Energy conservation error of : ", energy_sum
              write (*,*) "Snow Top layer",SnowLayerTop
              write(*,*) "F_abs_sum: ",F_abs_sum
              write(*,*) "F_btm_net: ",F_btm_net
              write(*,*) "F_sfc_pls: ",F_sfc_pls
              write(*,*) "mu_not*pi*flx_slrd_val: ", mu_not*ConstPI*flx_slrd_val
              write(*,*) "flx_slri_val", flx_slri_val
              write(*,*) "bnd_idx", LoopInd
              write(*,*) "F_abs", F_abs
              write(*,*) "albedo", albedo
              write(*,*) "direct soil albedo",AlbedoSoilDir(II,1,JJ),AlbedoSoilDir(II,2,JJ)
              write(*,*) "diffuse soil albedo",AlbedoSoilDif(II,1,JJ),AlbedoSoilDif(II,2,JJ)
              stop "ERROR in SNICAR energy conservation"
#endif
          endif

          albout_lcl_3d(II,LoopInd,JJ) = albedo

          ! Check that albedo is less than 1
          if (albout_lcl_3d(II,LoopInd,JJ) > 1.0) then

#ifndef _OPENACC
              write (*,*) "SNICAR ERROR: Albedo > 1.0"
              write (*,*) "SNICAR STATS: bnd_idx= ",LoopInd
              write (*,*) "SNICAR STATS: albout_lcl_3d(II,bnd,JJ)= ",albout_lcl_3d(II,LoopInd,JJ), &
                       " albsfc_lcl(bnd_idx)= ",albsfc_val
              write (*,*) "SNICAR STATS: h2osno_total= ", SnowWaterEquiv(II,JJ), " snl= ", snl_lcl
              write (*,*) "SNICAR STATS: coszen= ", CosSolarZenithAngle(II,JJ), " flg_slr= ", FlagSwRadType
              write (*,*) "SNICAR STATS: BCphi(-2)= ", mss_cnc_aer_lcl(-2,1)
              write (*,*) "SNICAR STATS: BCphi(-1)= ", mss_cnc_aer_lcl(-1,1)
              write (*,*) "SNICAR STATS: BCphi(0)= ", mss_cnc_aer_lcl(0,1)

              write (*,*) "SNICAR STATS: L_snw(-2)= ", L_snw(-2)
              write (*,*) "SNICAR STATS: L_snw(-1)= ", L_snw(-1)
              write (*,*) "SNICAR STATS: L_snw(0)= ", L_snw(0)

              write (*,*) "SNICAR STATS: snw_rds(-2)= ", SnowRadius(II,-2,JJ)
              write (*,*) "SNICAR STATS: snw_rds(-1)= ", SnowRadius(II,-1,JJ)
              write (*,*) "SNICAR STATS: snw_rds(0)= ", SnowRadius(II,0,JJ)
              stop "ERROR in SNICAR too large albedo"
#endif

          endif

       enddo ! loop over all snow spectral bands

       ! Pre-compute spectral weight sums for post-processing
       flx_sum_wgt = 0.0_kind_noahmp
       !$acc loop seq
       do k = nir_bnd_bgn, nir_bnd_end
          flx_sum_wgt = flx_sum_wgt + flx_wgt_3d(II,k,JJ)
       enddo
       flx_sum_wgt_vis = 0.0_kind_noahmp
       !$acc loop seq
       do k = 1, nir_bnd_bgn - 1
          flx_sum_wgt_vis = flx_sum_wgt_vis + flx_wgt_3d(II,k,JJ)
       enddo

       ! Weight output NIR albedo appropriately
       ! for 5- and 3-band cases
       if (NumSnicarRadBand <= 5) then
          if (FlagSwRadType == 1) then
              AlbedoSnowDir(II,1,JJ) = albout_lcl_3d(II,1,JJ)
          elseif (FlagSwRadType == 2)then
              AlbedoSnowDif(II,1,JJ) = albout_lcl_3d(II,1,JJ)
          endif

          flx_sum         = 0.0
          !$acc loop seq
          do LoopInd= nir_bnd_bgn,nir_bnd_end
              flx_sum = flx_sum + flx_wgt_3d(II,LoopInd,JJ)*albout_lcl_3d(II,LoopInd,JJ)
          end do

          if (FlagSwRadType == 1) then
              AlbedoSnowDir(II,2,JJ) = flx_sum / flx_sum_wgt
          elseif (FlagSwRadType == 2)then
              AlbedoSnowDif(II,2,JJ) = flx_sum / flx_sum_wgt
          endif

       end if

       ! for 480-band case
       if (NumSnicarRadBand == 480) then
          ! average for VIS band
          flx_sum         = 0.0
          !$acc loop seq
          do LoopInd= 1, (nir_bnd_bgn-1)
             flx_sum = flx_sum + flx_wgt_3d(II,LoopInd,JJ)*albout_lcl_3d(II,LoopInd,JJ)
          end do

          if (FlagSwRadType == 1) then
             AlbedoSnowDir(II,1,JJ) = flx_sum / flx_sum_wgt_vis
          elseif (FlagSwRadType == 2)then
             AlbedoSnowDif(II,1,JJ) = flx_sum / flx_sum_wgt_vis
          endif

          ! average for NIR band
          flx_sum         = 0.0
          !$acc loop seq
          do LoopInd= nir_bnd_bgn,nir_bnd_end
             flx_sum = flx_sum + flx_wgt_3d(II,LoopInd,JJ)*albout_lcl_3d(II,LoopInd,JJ)
          end do

          if (FlagSwRadType == 1) then
              AlbedoSnowDir(II,2,JJ) = flx_sum / flx_sum_wgt
          elseif (FlagSwRadType == 2) then
              AlbedoSnowDif(II,2,JJ) = flx_sum / flx_sum_wgt
          endif

       end if

       if (NumSnicarRadBand <= 5) then
          if (FlagSwRadType == 1) then
             !$acc loop seq
             do k = -NSNOW_MAX+1, 1
                FracRadSwAbsSnowDir(II,k,1,JJ) = flx_abs_lcl_4d(II,k,1,JJ)
             enddo
          elseif (FlagSwRadType == 2) then
             !$acc loop seq
             do k = -NSNOW_MAX+1, 1
                FracRadSwAbsSnowDif(II,k,1,JJ) = flx_abs_lcl_4d(II,k,1,JJ)
             enddo
          endif

          !$acc loop seq
          do i=SnowLayerTop,1,1

             flx_sum = 0.0
             !$acc loop seq
             do LoopInd= nir_bnd_bgn,nir_bnd_end
                flx_sum = flx_sum + flx_wgt_3d(II,LoopInd,JJ)*flx_abs_lcl_4d(II,i,LoopInd,JJ)
             enddo

          if (FlagSwRadType == 1) then
             FracRadSwAbsSnowDir(II,i,2,JJ) = flx_sum / flx_sum_wgt
          elseif (FlagSwRadType == 2) then
             FracRadSwAbsSnowDif(II,i,2,JJ) = flx_sum / flx_sum_wgt
          endif
 
          end do

       endif

       ! for 480-band case
       if (NumSnicarRadBand == 480) then
          !$acc loop seq
          do i=SnowLayerTop,1,1

             ! average for VIS band
             flx_sum = 0.0
             !$acc loop seq
             do LoopInd= 1,(nir_bnd_bgn-1)
                flx_sum = flx_sum + flx_wgt_3d(II,LoopInd,JJ)*flx_abs_lcl_4d(II,i,LoopInd,JJ)
             enddo

             if (FlagSwRadType == 1) then
                FracRadSwAbsSnowDir(II,i,1,JJ)=flx_sum / flx_sum_wgt_vis
             elseif (FlagSwRadType == 2) then
                FracRadSwAbsSnowDif(II,i,1,JJ)=flx_sum / flx_sum_wgt_vis
             endif

             ! average for NIR band
             flx_sum = 0.0
             !$acc loop seq
             do LoopInd= nir_bnd_bgn,nir_bnd_end
                flx_sum = flx_sum + flx_wgt_3d(II,LoopInd,JJ)*flx_abs_lcl_4d(II,i,LoopInd,JJ)
             enddo

             if (FlagSwRadType == 1) then
                FracRadSwAbsSnowDir(II,i,2,JJ)=flx_sum / flx_sum_wgt
             elseif (FlagSwRadType == 2) then
                FracRadSwAbsSnowDif(II,i,2,JJ)=flx_sum / flx_sum_wgt
             endif
           end do

       end if

       ! high solar zenith angle adjustment for Adding-doubling solver results

       if ( OptSnicarRTSolver==2) then
          ! near-IR direct albedo/absorption adjustment for high solar zenith angles
          ! solar zenith angle parameterization
          ! calculate the scaling factor for NIR direct albedo if SZA>75 degree
          if ((mu_not < mu_75) .and. (FlagSwRadType == 1)) then
             sza_c1 = sza_a0 + sza_a1 * mu_not + sza_a2 * mu_not**2
             sza_c0 = sza_b0 + sza_b1 * mu_not + sza_b2 * mu_not**2
             sza_factor = sza_c1 * (log10(snw_rds_lcl(SnowLayerTop) * c1) - c6) + sza_c0
             flx_sza_adjust  = AlbedoSnowDir(II,2,JJ) * (sza_factor-c1) * flx_sum_wgt
             AlbedoSnowDir(II,2,JJ) = AlbedoSnowDir(II,2,JJ) * sza_factor
             FracRadSwAbsSnowDir(II,SnowLayerTop,2,JJ) = FracRadSwAbsSnowDir(II,SnowLayerTop,2,JJ) - flx_sza_adjust
          endif
       endif ! end of  OptSnicarRTSolver==2


    ! If snow < minimum_snow, but > 0, and there is sun, set albedo to underlying surface albedo
    elseif ((SnowWaterEquiv(II,JJ) < SnowWaterEquivMin) .and. (SnowWaterEquiv(II,JJ) > 0.0) ) then

       if (IndicatorIceSfc(II,JJ) == 0) then
          if (FlagSwRadType == 1) then
             AlbedoSnowDir(II,1,JJ) = AlbedoSoilDir(II,1,JJ) 
             AlbedoSnowDir(II,2,JJ) = AlbedoSoilDir(II,2,JJ) 
          elseif (FlagSwRadType == 2) then
             AlbedoSnowDif(II,1,JJ) = AlbedoSoilDif(II,1,JJ)
             AlbedoSnowDif(II,2,JJ) = AlbedoSoilDif(II,2,JJ)
          endif
       elseif (IndicatorIceSfc(II,JJ) == -1) then !land ice
          AlbedoSnowDif(II,1,JJ) = AlbedoLandIce(II,1,JJ)
          AlbedoSnowDif(II,2,JJ) = AlbedoLandIce(II,2,JJ)
       endif
    ! There is either zero snow, or no sun
    else

       if (FlagSwRadType == 1) then
          AlbedoSnowDir(II,1,JJ) = 0.0
          AlbedoSnowDir(II,2,JJ) = 0.0
       elseif (FlagSwRadType == 2) then
          AlbedoSnowDif(II,1,JJ) = 0.0
          AlbedoSnowDif(II,2,JJ) = 0.0
       endif

    endif ! if column has mim snow

    if (FlagSwRadType == 1) then
       if (AlbedoSnowDir(II,1,JJ)<0.0 .or. AlbedoSnowDir(II,2,JJ)<0.0 .or. AlbedoSnowDir(II,1,JJ)>1.0 .or. AlbedoSnowDir(II,2,JJ)>1.0)then
#ifndef _OPENACC
          print *,'Error in SNICAR direct snow albedo: ',FlagSwRadType,AlbedoSnowDir(II,1,JJ),AlbedoSnowDir(II,2,JJ),CosSolarZenithAngle(II,JJ)
          stop "Error in SNICAR direct snow albedo"
#endif
       endif
    endif

    if (FlagSwRadType == 2) then
       if (AlbedoSnowDif(II,1,JJ)<0.0 .or. AlbedoSnowDif(II,2,JJ)<0.0 .or. AlbedoSnowDif(II,1,JJ)>1.0 .or. AlbedoSnowDif(II,2,JJ)>1.0)then
#ifndef _OPENACC
          print *,'Error in SNICAR diffuse snow albedo',FlagSwRadType,AlbedoSnowDif(II,1,JJ),AlbedoSnowDif(II,2,JJ),CosSolarZenithAngle(II,JJ)
          stop "Error in SNICAR diffuse snow albedo"
#endif
       endif
    endif

      enddo ! II
    enddo ! JJ

    !$acc end data
    deallocate(flx_wgt_3d, albout_lcl_3d, flx_abs_lcl_4d)

    end associate

  end subroutine SnowRadiationSnicar

end module SnowRadiationSnicarMod
