module ResistanceBareGroundChen97Mod

!!! Compute bare ground resistance and exchange coefficient for momentum and heat
!!! based on Chen et al. (1997, BLM)
!!! This scheme can handle both over open water and over solid surface

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine ResistanceBareGroundChen97(noahmp, IndIter)
   !$acc routine gang
! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SFCDIF2 for bare ground portion
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

! in & out variables
    integer               , intent(in   ) :: IndIter   ! iteration index
    type(noahmp_type)     , intent(inout) :: noahmp

! local variables
    integer                               :: I, J      ! grid indices
    integer                               :: ILECH, ITR
    real(kind=kind_noahmp)                :: ZZ, PSLMU, PSLMS, PSLHU, PSLHS
    real(kind=kind_noahmp)                :: XX, PSPMU, YY, PSPMS, PSPHU, PSPHS
    real(kind=kind_noahmp)                :: ZILFC, ZU, ZT, RDZ, CXCH, DTHV, DU2
    real(kind=kind_noahmp)                :: BTGH, ZSLU, ZSLT, RLOGU, RLOGT, RLMA
    real(kind=kind_noahmp)                :: ZETALT, ZETALU, ZETAU, ZETAT, XLU4
    real(kind=kind_noahmp)                :: XLT4, XU4, XT4, XLU, XLT, XU, XT
    real(kind=kind_noahmp)                :: PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN
! local parameters
    integer               , parameter     :: ITRMX  = 5
    real(kind=kind_noahmp), parameter     :: WWST   = 1.2
    real(kind=kind_noahmp), parameter     :: WWST2  = WWST * WWST
    real(kind=kind_noahmp), parameter     :: VKRM   = 0.40
    real(kind=kind_noahmp), parameter     :: EXCM   = 0.001
    real(kind=kind_noahmp), parameter     :: BETA   = 1.0 / 270.0
    real(kind=kind_noahmp), parameter     :: BTG    = BETA * ConstGravityAcc
    real(kind=kind_noahmp), parameter     :: ELFC   = VKRM * BTG
    real(kind=kind_noahmp), parameter     :: WOLD   = 0.15
    real(kind=kind_noahmp), parameter     :: WNEW   = 1.0 - WOLD
    real(kind=kind_noahmp), parameter     :: PIHF   = 3.14159265 / 2.0
    real(kind=kind_noahmp), parameter     :: EPSU2  = 1.0e-4
    real(kind=kind_noahmp), parameter     :: EPSUST = 0.07
    real(kind=kind_noahmp), parameter     :: EPSIT  = 1.0e-4
    real(kind=kind_noahmp), parameter     :: EPSA   = 1.0e-8
    real(kind=kind_noahmp), parameter     :: ZTMIN  = -5.0
    real(kind=kind_noahmp), parameter     :: ZTMAX  = 1.0
    real(kind=kind_noahmp), parameter     :: HPBL   = 1000.0
    real(kind=kind_noahmp), parameter     :: SQVISC = 258.2
    real(kind=kind_noahmp), parameter     :: RIC    = 0.183
    real(kind=kind_noahmp), parameter     :: RRIC   = 1.0 / RIC
    real(kind=kind_noahmp), parameter     :: FHNEU  = 0.8
    real(kind=kind_noahmp), parameter     :: RFC    = 0.191
    real(kind=kind_noahmp), parameter     :: RFAC   = RIC / ( FHNEU * RFC * RFC )
! local statement functions
    ! LECH'S surface functions
    PSLMU(ZZ) = -0.96 * log(1.0 - 4.5 * ZZ)
    PSLMS(ZZ) = ZZ / RFC - 2.076 * (1.0 - 1.0/(ZZ + 1.0))
    PSLHU(ZZ) = -0.96 * log(1.0 - 4.5 * ZZ)
    PSLHS(ZZ) = ZZ * RFAC - 2.076 * (1.0 - exp(-1.2 * ZZ))

    ! PAULSON'S surface functions
    PSPMU(XX) = -2.0*log( (XX+1.0)*0.5 ) - log( (XX*XX+1.0)*0.5 ) + 2.0*atan(XX) - PIHF
    PSPMS(YY) = 5.0 * YY
    PSPHU(XX) = -2.0 * log( (XX*XX + 1.0)*0.5 )
    PSPHS(YY) = 5.0 * YY

! --------------------------------------------------------------------
    associate(                                                                           &
              SnowDepth               => noahmp%water%state%SnowDepth               ,& ! in,    snow depth [m]
              ZilitinkevichCoeff      => noahmp%energy%param%ZilitinkevichCoeff     ,& ! in,    Calculate roughness length of heat
              RefHeightAboveGrd       => noahmp%energy%state%RefHeightAboveGrd      ,& ! in,    reference height [m] above ground
              TemperaturePotRefHeight => noahmp%energy%state%TemperaturePotRefHeight,& ! in,    potential temp at reference height [K]
              WindSpdRefHeight        => noahmp%energy%state%WindSpdRefHeight       ,& ! in,    wind speed [m/s] at reference height
              RoughLenMomGrd          => noahmp%energy%state%RoughLenMomGrd         ,& ! in,    roughness length [m], momentum, ground
              TemperatureGrdBare      => noahmp%energy%state%TemperatureGrdBare     ,& ! in,    bare ground temperature [K]
              ExchCoeffMomBare        => noahmp%energy%state%ExchCoeffMomBare       ,& ! inout, exchange coeff [m/s] momentum, above ZeroPlaneDisp, bare ground
              ExchCoeffShBare         => noahmp%energy%state%ExchCoeffShBare        ,& ! inout, exchange coeff [m/s] for heat, above ZeroPlaneDisp, bare ground
              MoStabParaBare          => noahmp%energy%state%MoStabParaBare         ,& ! inout, Monin-Obukhov stability (z/L), above ZeroPlaneDisp, bare ground
              FrictionVelVertBare     => noahmp%energy%state%FrictionVelVertBare    ,& ! inout, friction velocity [m/s] in vertical direction, bare ground
              FrictionVelBare         => noahmp%energy%state%FrictionVelBare        ,& ! inout, friction velocity [m/s], bare ground
              ResistanceMomBareGrd    => noahmp%energy%state%ResistanceMomBareGrd   ,& ! out,   aerodynamic resistance for momentum [s/m], bare ground
              ResistanceShBareGrd     => noahmp%energy%state%ResistanceShBareGrd    ,& ! out,   aerodynamic resistance for sensible heat [s/m], bare ground
              ResistanceLhBareGrd     => noahmp%energy%state%ResistanceLhBareGrd     & ! out,   aerodynamic resistance for water vapor [s/m], bare ground
             )

   !$acc loop gang vector collapse(2) &
   !$acc private(ILECH,ZILFC,ZU,ZT,RDZ,CXCH,DTHV,DU2,BTGH,ZSLU,ZSLT,RLOGU,RLOGT,RLMA) &
   !$acc private(ZETALT,ZETALU,ZETAU,ZETAT,XLU4,XLT4,XU4,XT4,XLU,XLT,XU,XT) private(PSMZ,SIMM,PSHZ,SIMH,USTARK,RLMN)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1
    ! C......ZTFC=0.1
    ! ZilitinkevichCoeff: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
    ILECH = 0
    ZILFC = -ZilitinkevichCoeff(I,J) * VKRM * SQVISC
    ZU    = RoughLenMomGrd(I,J)
    RDZ   = 1.0 / RefHeightAboveGrd(I,J)
    CXCH  = EXCM * RDZ
    DTHV  = TemperaturePotRefHeight(I,J) - TemperatureGrdBare(I,J)

    ! BELJARS correction of friction velocity u*
    DU2   = max(WindSpdRefHeight(I,J)*WindSpdRefHeight(I,J), EPSU2)
    BTGH  = BTG * HPBL
    if ( IndIter == 1 ) then
       if ( (BTGH*ExchCoeffShBare(I,J)*DTHV) /= 0.0 ) then
          FrictionVelVertBare(I,J) = WWST2 * abs(BTGH*ExchCoeffShBare(I,J)*DTHV)**(2.0/3.0)
       else
          FrictionVelVertBare(I,J) = 0.0
       endif
       FrictionVelBare(I,J) = max(sqrt(ExchCoeffMomBare(I,J)*sqrt(DU2+FrictionVelVertBare(I,J))), EPSUST)
       MoStabParaBare(I,J)  = ELFC * ExchCoeffShBare(I,J) * DTHV / FrictionVelBare(I,J)**3
    endif

    ! ZILITINKEVITCH approach for ZT
    ZT    = max(1.0e-6, exp(ZILFC*sqrt(FrictionVelBare(I,J)*RoughLenMomGrd(I,J)))*RoughLenMomGrd(I,J))
    ZSLU  = RefHeightAboveGrd(I,J) + ZU
    ZSLT  = RefHeightAboveGrd(I,J) + ZT
    RLOGU = log(ZSLU / ZU)
    RLOGT = log(ZSLT / ZT)

    ! Monin-Obukhov length scale
    ZETALT         = max(ZSLT*MoStabParaBare(I,J), ZTMIN)
    MoStabParaBare(I,J) = ZETALT / ZSLT
    ZETALU         = ZSLU * MoStabParaBare(I,J)
    ZETAU          = ZU * MoStabParaBare(I,J)
    ZETAT          = ZT * MoStabParaBare(I,J)
    if ( ILECH == 0 ) then
       if ( MoStabParaBare(I,J) < 0.0 ) then
          XLU4 = 1.0 - 16.0 * ZETALU
          XLT4 = 1.0 - 16.0 * ZETALT
          XU4  = 1.0 - 16.0 * ZETAU
          XT4  = 1.0 - 16.0 * ZETAT
          XLU  = sqrt(sqrt(XLU4))
          XLT  = sqrt(sqrt(XLT4))
          XU   = sqrt(sqrt(XU4))
          XT   = sqrt(sqrt(XT4))
          PSMZ = PSPMU(XU)
          SIMM = PSPMU(XLU) - PSMZ + RLOGU
          PSHZ = PSPHU(XT)
          SIMH = PSPHU(XLT) - PSHZ + RLOGT
       else
          ZETALU = min(ZETALU, ZTMAX)
          ZETALT = min(ZETALT, ZTMAX)
          ZETAU  = min(ZETAU, ZTMAX/(ZSLU/ZU))   ! Barlage: add limit on ZETAU/ZETAT
          ZETAT  = min(ZETAT, ZTMAX/(ZSLT/ZT))   ! Barlage: prevent SIMM/SIMH < 0
          PSMZ   = PSPMS(ZETAU)
          SIMM   = PSPMS(ZETALU) - PSMZ + RLOGU
          PSHZ   = PSPHS(ZETAT)
          SIMH   = PSPHS(ZETALT) - PSHZ + RLOGT
       endif
    else ! LECH's functions
       if ( MoStabParaBare(I,J) < 0.0 ) then
          PSMZ = PSLMU(ZETAU)
          SIMM = PSLMU(ZETALU) - PSMZ + RLOGU
          PSHZ = PSLHU(ZETAT)
          SIMH = PSLHU(ZETALT) - PSHZ + RLOGT
       else
          ZETALU = min(ZETALU, ZTMAX)
          ZETALT = min(ZETALT, ZTMAX)
          PSMZ   = PSLMS(ZETAU)
          SIMM   = PSLMS(ZETALU) - PSMZ + RLOGU
          PSHZ   = PSLHS(ZETAT)
          SIMH   = PSLHS(ZETALT) - PSHZ + RLOGT
       endif
    endif

    ! BELJARS correction of friction velocity u*
    FrictionVelBare(I,J) = max(sqrt(ExchCoeffMomBare(I,J)*sqrt(DU2+FrictionVelVertBare(I,J))), EPSUST)

    ! ZILITINKEVITCH fix for ZT
    ZT     = max(1.0e-6, exp(ZILFC*sqrt(FrictionVelBare(I,J)*RoughLenMomGrd(I,J)))*RoughLenMomGrd(I,J))
    ZSLT   = RefHeightAboveGrd(I,J) + ZT
    RLOGT  = log(ZSLT/ZT)
    USTARK = FrictionVelBare(I,J) * VKRM

    ! avoid tangent linear problems near zero
    if ( SIMM < 1.0e-6 ) SIMM = 1.0e-6         ! Limit stability function
    ExchCoeffMomBare(I,J) = max(USTARK/SIMM, CXCH)
    if ( SIMH < 1.0e-6 ) SIMH = 1.0e-6         ! Limit stability function
    ExchCoeffShBare(I,J)  = max(USTARK/SIMH, CXCH)

    ! update vertical friction velocity w*
    if ( BTGH*ExchCoeffShBare(I,J)*DTHV /= 0.0 ) then
       FrictionVelVertBare(I,J) = WWST2 * abs(BTGH*ExchCoeffShBare(I,J)*DTHV)**(2.0/3.0)
    else
       FrictionVelVertBare(I,J) = 0.0
    endif

    ! update M-O stability parameter
    RLMN = ELFC * ExchCoeffShBare(I,J) * DTHV / FrictionVelBare(I,J)**3
    RLMA = MoStabParaBare(I,J) * WOLD + RLMN * WNEW
    MoStabParaBare(I,J) = RLMA

    ! Undo the multiplication by windspeed that applies to exchange coeff
    ExchCoeffShBare(I,J)  = ExchCoeffShBare(I,J) / WindSpdRefHeight(I,J)
    ExchCoeffMomBare(I,J) = ExchCoeffMomBare(I,J) / WindSpdRefHeight(I,J)
    if ( SnowDepth(I,J) > 0.0 ) then
       ExchCoeffMomBare(I,J) = min(0.01, ExchCoeffMomBare(I,J))    ! exch coeff is too large, causing
       ExchCoeffShBare(I,J)  = min(0.01, ExchCoeffShBare(I,J))     ! computational instability
    endif

    ! compute aerodynamic resistance
    ResistanceMomBareGrd(I,J) = max(1.0, 1.0/(ExchCoeffMomBare(I,J)*WindSpdRefHeight(I,J)))
    ResistanceShBareGrd(I,J)  = max(1.0, 1.0/(ExchCoeffShBare(I,J)*WindSpdRefHeight(I,J)))
    ResistanceLhBareGrd(I,J)  = ResistanceShBareGrd(I,J)


      end do
    end do


    end associate

  end subroutine ResistanceBareGroundChen97

end module ResistanceBareGroundChen97Mod
