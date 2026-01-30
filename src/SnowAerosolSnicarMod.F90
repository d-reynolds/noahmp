module SnowAerosolSnicarMod

!!! compute aerosol content in snow and its evolution prepared for SNICAR snow albedo calculation

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowAerosolSnicar(noahmp)

! --------------------------------- Code history -----------------------------------
! Implementation: T.-S. Lin, C. He, et al. (2025, JHM)
! Adapted from CTSM, AerosolFluxes, AerosolMasses, CalcAndApplyAerosolFluxes modules
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! ----------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd                       ! do loop/array indices
    real(kind=kind_noahmp)           :: SnowMass                      ! liquid+ice snow mass in a layer [kg/m2]
    real(kind=kind_noahmp)           :: FluxInBChydrophi              ! flux of hydrophilic BC into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutBChydrophi             ! flux of hydrophilic BC out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInBChydropho              ! flux of hydrophobic BC into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutBChydropho             ! flux of hydrophobic BC out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInOChydrophi              ! flux of hydrophilic OC into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutOChydrophi             ! flux of hydrophilic OC out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInOChydropho              ! flux of hydrophobic OC into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutOChydropho             ! flux of hydrophobic OC out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInDust1                   ! flux of dust species 1 into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutDust1                  ! flux of dust species 1 out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInDust2                   ! flux of dust species 2 into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutDust2                  ! flux of dust species 2 out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInDust3                   ! flux of dust species 3 into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutDust3                  ! flux of dust species 3 out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInDust4                   ! flux of dust species 4 into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutDust4                  ! flux of dust species 4 out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInDust5                   ! flux of dust species 5 into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutDust5                  ! flux of dust species 5 out of layer [kg/s]
    integer                          :: I, J                          ! grid indices

    !$acc parallel loop collapse(2) gang vector present(noahmp) &
    !$acc private(LoopInd, SnowMass, FluxInBChydrophi, FluxOutBChydrophi, FluxInBChydropho, FluxOutBChydropho) &
    !$acc private(FluxInOChydrophi, FluxOutOChydrophi, FluxInOChydropho, FluxOutOChydropho) &
    !$acc private(FluxInDust1, FluxOutDust1, FluxInDust2, FluxOutDust2, FluxInDust3, FluxOutDust3) &
    !$acc private(FluxInDust4, FluxOutDust4, FluxInDust5, FluxOutDust5)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE
! --------------------------------------------------------------------
    associate(                                                                          &
              MainTimeStep           => noahmp%config%domain%MainTimeStep              ,& ! in,    noahmp main time step [s]
              NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax           ,& ! in,    maximum number of snow layers
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg(I,J)           ,& ! in,    actual number of snow layers (negative)
              DepBChydropho          => noahmp%forcing%DepBChydropho(I,J)                   ,& ! in,    hydrophobic Black Carbon deposition [kg m-2 s-1] 
              DepBChydrophi          => noahmp%forcing%DepBChydrophi(I,J)                   ,& ! in,    hydrophillic Black Carbon deposition [kg m-2 s-1]
              DepOChydropho          => noahmp%forcing%DepOChydropho(I,J)                   ,& ! in,    hydrophobic Organic Carbon deposition [kg m-2 s-1]
              DepOChydrophi          => noahmp%forcing%DepOChydrophi(I,J)                   ,& ! in,    hydrophillic Organic Carbon deposition [kg m-2 s-1]
              DepDust1               => noahmp%forcing%DepDust1(I,J)                        ,& ! in,    dust species 1 deposition [kg m-2 s-1]
              DepDust2               => noahmp%forcing%DepDust2(I,J)                        ,& ! in,    dust species 2 deposition [kg m-2 s-1]
              DepDust3               => noahmp%forcing%DepDust3(I,J)                        ,& ! in,    dust species 3 deposition [kg m-2 s-1]
              DepDust4               => noahmp%forcing%DepDust4(I,J)                        ,& ! in,    dust species 4 deposition [kg m-2 s-1]
              DepDust5               => noahmp%forcing%DepDust5(I,J)                        ,& ! in,    dust species 5 deposition [kg m-2 s-1]
              ScavEffMeltScale       => noahmp%water%param%ScavEffMeltScale(I,J)            ,& ! in,    Scaling factor modifying scavenging factors for aerosol in meltwater (-)
              ScavEffMeltBCphi       => noahmp%water%param%ScavEffMeltBCphi(I,J)            ,& ! in,    scavenging factor for hydrophillic BC inclusion in meltwater [frc]
              ScavEffMeltBCpho       => noahmp%water%param%ScavEffMeltBCpho(I,J)            ,& ! in,    scavenging factor for hydrophobic BC inclusion in meltwater  [frc]
              ScavEffMeltOCphi       => noahmp%water%param%ScavEffMeltOCphi(I,J)            ,& ! in,    scavenging factor for hydrophillic OC inclusion in meltwater [frc]
              ScavEffMeltOCpho       => noahmp%water%param%ScavEffMeltOCpho(I,J)            ,& ! in,    scavenging factor for hydrophobic OC inclusion in meltwater  [frc]
              ScavEffMeltDust1       => noahmp%water%param%ScavEffMeltDust1(I,J)            ,& ! in,    scavenging factor for dust species 1 inclusion in meltwater  [frc]
              ScavEffMeltDust2       => noahmp%water%param%ScavEffMeltDust2(I,J)            ,& ! in,    scavenging factor for dust species 2 inclusion in meltwater  [frc]
              ScavEffMeltDust3       => noahmp%water%param%ScavEffMeltDust3(I,J)            ,& ! in,    scavenging factor for dust species 3 inclusion in meltwater  [frc]
              ScavEffMeltDust4       => noahmp%water%param%ScavEffMeltDust4(I,J)            ,& ! in,    scavenging factor for dust species 4 inclusion in meltwater  [frc]
              ScavEffMeltDust5       => noahmp%water%param%ScavEffMeltDust5(I,J)            ,& ! in,    scavenging factor for dust species 5 inclusion in meltwater  [frc]
              SnowIce                => noahmp%water%state%SnowIce                     ,& ! in,    snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater                ,& ! in,    snow layer liquid water [mm]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv(I,J)              ,& ! in,    snow water equivalent [mm]
              OutflowSnowLayer       => noahmp%water%flux%OutflowSnowLayer             ,& ! in,    water flow out of each snow layer [mm/s]
              MassBChydropho         => noahmp%water%state%MassBChydropho              ,& ! inout, mass of hydrophobic Black Carbon in snow [kg m-2]
              MassBChydrophi         => noahmp%water%state%MassBChydrophi              ,& ! inout, mass of hydrophillic Black Carbon in snow [kg m-2]
              MassOChydropho         => noahmp%water%state%MassOChydropho              ,& ! inout, mass of hydrophobic Organic Carbon in snow [kg m-2]
              MassOChydrophi         => noahmp%water%state%MassOChydrophi              ,& ! inout, mass of hydrophillic Organic Carbon in snow [kg m-2]
              MassDust1              => noahmp%water%state%MassDust1                   ,& ! inout, mass of dust species 1 in snow [kg m-2]
              MassDust2              => noahmp%water%state%MassDust2                   ,& ! inout, mass of dust species 2 in snow [kg m-2]
              MassDust3              => noahmp%water%state%MassDust3                   ,& ! inout, mass of dust species 3 in snow [kg m-2]
              MassDust4              => noahmp%water%state%MassDust4                   ,& ! inout, mass of dust species 4 in snow [kg m-2]
              MassDust5              => noahmp%water%state%MassDust5                   ,& ! inout, mass of dust species 5 in snow [kg m-2]
              MassConcBChydropho     => noahmp%water%state%MassConcBChydropho          ,& ! inout, mass concentration of hydrophobic Black Carbon in snow [kg/kg]
              MassConcBChydrophi     => noahmp%water%state%MassConcBChydrophi          ,& ! inout, mass concentration of hydrophillic Black Carbon in snow [kg/kg]
              MassConcOChydropho     => noahmp%water%state%MassConcOChydropho          ,& ! inout, mass concentration of hydrophobic Organic Carbon in snow [kg/kg]
              MassConcOChydrophi     => noahmp%water%state%MassConcOChydrophi          ,& ! inout, mass concentration of hydrophillic Organic Carbon in snow [kg/kg]
              MassConcDust1          => noahmp%water%state%MassConcDust1               ,& ! inout, mass concentration of dust species 1 in snow [kg/kg]
              MassConcDust2          => noahmp%water%state%MassConcDust2               ,& ! inout, mass concentration of dust species 2 in snow [kg/kg]
              MassConcDust3          => noahmp%water%state%MassConcDust3               ,& ! inout, mass concentration of dust species 3 in snow [kg/kg]
              MassConcDust4          => noahmp%water%state%MassConcDust4               ,& ! inout, mass concentration of dust species 4 in snow [kg/kg]
              MassConcDust5          => noahmp%water%state%MassConcDust5                & ! inout, mass concentration of dust species 5 in snow [kg/kg]
             )
! ----------------------------------------------------------------------

    ! initialize
    FluxInBChydropho  = 0.0
    FluxInBChydrophi  = 0.0
    FluxInOChydropho  = 0.0
    FluxInOChydrophi  = 0.0
    FluxInDust1       = 0.0
    FluxInDust2       = 0.0
    FluxInDust3       = 0.0
    FluxInDust4       = 0.0
    FluxInDust5       = 0.0
    FluxOutBChydropho = 0.0
    FluxOutBChydrophi = 0.0
    FluxOutOChydropho = 0.0
    FluxOutOChydrophi = 0.0
    FluxOutDust1      = 0.0
    FluxOutDust2      = 0.0
    FluxOutDust3      = 0.0
    FluxOutDust4      = 0.0
    FluxOutDust5      = 0.0

    ! compute aerosol mass in snow for each layer from interlayer flux
    !$acc loop seq
    do LoopInd = -NumSnowLayerMax+1, 0

       SnowMass = SnowLiqWater(I,LoopInd,J) + SnowIce(I,LoopInd,J)

       if (LoopInd >= NumSnowLayerNeg+1 .and. SnowMass > 0.0) then

          MassBChydropho(I,LoopInd,J) =  MassBChydropho (I,LoopInd,J) + FluxInBChydropho * MainTimeStep
          MassBChydrophi(I,LoopInd,J) =  MassBChydrophi (I,LoopInd,J) + FluxInBChydrophi * MainTimeStep
          MassOChydropho(I,LoopInd,J) =  MassOChydropho (I,LoopInd,J) + FluxInOChydropho * MainTimeStep
          MassOChydrophi(I,LoopInd,J) =  MassOChydrophi (I,LoopInd,J) + FluxInOChydrophi * MainTimeStep
          MassDust1(I,LoopInd,J)      =  MassDust1 (I,LoopInd,J)      + FluxInDust1 * MainTimeStep
          MassDust2(I,LoopInd,J)      =  MassDust2 (I,LoopInd,J)      + FluxInDust2 * MainTimeStep
          MassDust3(I,LoopInd,J)      =  MassDust3 (I,LoopInd,J)      + FluxInDust3 * MainTimeStep
          MassDust4(I,LoopInd,J)      =  MassDust4 (I,LoopInd,J)      + FluxInDust4 * MainTimeStep
          MassDust5(I,LoopInd,J)      =  MassDust5 (I,LoopInd,J)      + FluxInDust5 * MainTimeStep

          !BCPHO
          FluxOutBChydropho = OutflowSnowLayer(I,LoopInd,J) * ScavEffMeltScale * &
                              ScavEffMeltBCpho * (MassBChydropho(I,LoopInd,J) / SnowMass)
          if (FluxOutBChydropho * MainTimeStep > MassBChydropho(I,LoopInd,J)) then
             FluxOutBChydropho = MassBChydropho(I,LoopInd,J) / MainTimeStep
             MassBChydropho(I,LoopInd,J) = 0.0
          else
             MassBChydropho(I,LoopInd,J) = MassBChydropho(I,LoopInd,J) - FluxOutBChydropho * MainTimeStep
          end if
          FluxInBChydropho = FluxOutBChydropho

          !BCPHI
          FluxOutBChydrophi = OutflowSnowLayer(I,LoopInd,J) * ScavEffMeltScale * &
                              ScavEffMeltBCphi * (MassBChydrophi(I,LoopInd,J) / SnowMass)
          if (FluxOutBChydrophi * MainTimeStep > MassBChydrophi(I,LoopInd,J)) then
             FluxOutBChydrophi = MassBChydrophi(I,LoopInd,J) / MainTimeStep
             MassBChydrophi(I,LoopInd,J) = 0.0
          else
             MassBChydrophi(I,LoopInd,J) = MassBChydrophi(I,LoopInd,J) - FluxOutBChydrophi * MainTimeStep
          end if
          FluxInBChydrophi = FluxOutBChydrophi

          !OCPHO
          FluxOutOChydropho = OutflowSnowLayer(I,LoopInd,J) * ScavEffMeltScale * &
                              ScavEffMeltOCpho * (MassOChydropho(I,LoopInd,J) / SnowMass)
          if (FluxOutOChydropho * MainTimeStep > MassOChydropho(I,LoopInd,J)) then
             FluxOutOChydropho = MassOChydropho(I,LoopInd,J) / MainTimeStep
             MassOChydropho(I,LoopInd,J) = 0.0
          else
             MassOChydropho(I,LoopInd,J) = MassOChydropho(I,LoopInd,J) - FluxOutOChydropho * MainTimeStep
          end if
          FluxInOChydropho = FluxOutOChydropho

          !OCPHI
          FluxOutOChydrophi = OutflowSnowLayer(I,LoopInd,J) * ScavEffMeltScale * &
                              ScavEffMeltOCphi * (MassOChydrophi(I,LoopInd,J) / SnowMass)
          if (FluxOutOChydrophi * MainTimeStep > MassOChydrophi(I,LoopInd,J)) then
             FluxOutOChydrophi = MassOChydrophi(I,LoopInd,J) / MainTimeStep
             MassOChydrophi(I,LoopInd,J) = 0.0
          else
             MassOChydrophi(I,LoopInd,J) = MassOChydrophi(I,LoopInd,J) - FluxOutOChydrophi * MainTimeStep
          end if
          FluxInOChydrophi = FluxOutOChydrophi

          !Dust 1
          FluxOutDust1 = OutflowSnowLayer(I,LoopInd,J) * ScavEffMeltScale * &
                         ScavEffMeltDust1 * (MassDust1(I,LoopInd,J) / SnowMass)
          if (FluxOutDust1 * MainTimeStep > MassDust1(I,LoopInd,J)) then
             FluxOutDust1 = MassDust1(I,LoopInd,J) / MainTimeStep
             MassDust1(I,LoopInd,J) = 0.0
          else
             MassDust1(I,LoopInd,J) = MassDust1(I,LoopInd,J) - FluxOutDust1 * MainTimeStep
          end if
          FluxInDust1 = FluxOutDust1

          !Dust 2
          FluxOutDust2 = OutflowSnowLayer(I,LoopInd,J) * ScavEffMeltScale * &
                         ScavEffMeltDust2 * (MassDust2(I,LoopInd,J) / SnowMass)
          if (FluxOutDust2 * MainTimeStep > MassDust2(I,LoopInd,J)) then
             FluxOutDust2 = MassDust2(I,LoopInd,J) / MainTimeStep
             MassDust2(I,LoopInd,J) = 0.0
          else
             MassDust2(I,LoopInd,J) = MassDust2(I,LoopInd,J) - FluxOutDust2 * MainTimeStep
          end if
          FluxInDust2 = FluxOutDust2

          !Dust 3
          FluxOutDust3 = OutflowSnowLayer(I,LoopInd,J) * ScavEffMeltScale * &
                         ScavEffMeltDust3 * (MassDust3(I,LoopInd,J) / SnowMass)
          if (FluxOutDust3 * MainTimeStep > MassDust3(I,LoopInd,J)) then
             FluxOutDust3 = MassDust3(I,LoopInd,J) / MainTimeStep
             MassDust3(I,LoopInd,J) = 0.0
          else
             MassDust3(I,LoopInd,J) = MassDust3(I,LoopInd,J) - FluxOutDust3 * MainTimeStep
          end if
          FluxInDust3 = FluxOutDust3

          !Dust 4
          FluxOutDust4 = OutflowSnowLayer(I,LoopInd,J) * ScavEffMeltScale * &
                         ScavEffMeltDust4 * (MassDust4(I,LoopInd,J) / SnowMass)
          if (FluxOutDust4 * MainTimeStep > MassDust4(I,LoopInd,J)) then
             FluxOutDust4 = MassDust4(I,LoopInd,J) / MainTimeStep
             MassDust4(I,LoopInd,J) = 0.0
          else
             MassDust4(I,LoopInd,J) = MassDust4(I,LoopInd,J) - FluxOutDust4 * MainTimeStep
          end if
          FluxInDust4 = FluxOutDust4

          !Dust 5
          FluxOutDust5 = OutflowSnowLayer(I,LoopInd,J) * ScavEffMeltScale * &
                         ScavEffMeltDust5 * (MassDust5(I,LoopInd,J) / SnowMass)
          if (FluxOutDust5 * MainTimeStep > MassDust5(I,LoopInd,J)) then
             FluxOutDust5 = MassDust5(I,LoopInd,J) / MainTimeStep
             MassDust5(I,LoopInd,J) = 0.0
          else
             MassDust5(I,LoopInd,J) = MassDust5(I,LoopInd,J) - FluxOutDust5 * MainTimeStep
          end if
          FluxInDust5 = FluxOutDust5

       else ! SnowMass <=0 or non-existence snow layer

          MassBChydropho(I,LoopInd,J) = 0.0
          MassBChydrophi(I,LoopInd,J) = 0.0
          MassOChydropho(I,LoopInd,J) = 0.0
          MassOChydrophi(I,LoopInd,J) = 0.0
          MassDust1(I,LoopInd,J)      = 0.0
          MassDust2(I,LoopInd,J)      = 0.0
          MassDust3(I,LoopInd,J)      = 0.0
          MassDust4(I,LoopInd,J)      = 0.0
          MassDust5(I,LoopInd,J)      = 0.0
          FluxInBChydropho        = 0.0
          FluxInBChydrophi        = 0.0
          FluxInOChydropho        = 0.0
          FluxInOChydrophi        = 0.0
          FluxInDust1             = 0.0
          FluxInDust2             = 0.0
          FluxInDust3             = 0.0
          FluxInDust4             = 0.0
          FluxInDust5             = 0.0
          FluxOutBChydropho       = 0.0
          FluxOutBChydrophi       = 0.0
          FluxOutOChydropho       = 0.0
          FluxOutOChydrophi       = 0.0
          FluxOutDust1            = 0.0
          FluxOutDust2            = 0.0
          FluxOutDust3            = 0.0
          FluxOutDust4            = 0.0
          FluxOutDust5            = 0.0

       endif
    enddo

    ! update aerosol mass for the top snow layer from atmos deposition flux
    if (NumSnowLayerNeg < 0 ) then
       MassBChydropho(I,NumSnowLayerNeg+1,J) =  MassBChydropho(I,NumSnowLayerNeg+1,J) +  DepBChydropho * MainTimeStep
       MassBChydrophi(I,NumSnowLayerNeg+1,J) =  MassBChydrophi(I,NumSnowLayerNeg+1,J) +  DepBChydrophi * MainTimeStep
       MassOChydropho(I,NumSnowLayerNeg+1,J) =  MassOChydropho(I,NumSnowLayerNeg+1,J) +  DepOChydropho * MainTimeStep
       MassOChydrophi(I,NumSnowLayerNeg+1,J) =  MassOChydrophi(I,NumSnowLayerNeg+1,J) +  DepOChydrophi * MainTimeStep
       MassDust1(I,NumSnowLayerNeg+1,J)      =  MassDust1(I,NumSnowLayerNeg+1,J)      +  DepDust1 * MainTimeStep
       MassDust2(I,NumSnowLayerNeg+1,J)      =  MassDust2(I,NumSnowLayerNeg+1,J)      +  DepDust2 * MainTimeStep
       MassDust3(I,NumSnowLayerNeg+1,J)      =  MassDust3(I,NumSnowLayerNeg+1,J)      +  DepDust3 * MainTimeStep
       MassDust4(I,NumSnowLayerNeg+1,J)      =  MassDust4(I,NumSnowLayerNeg+1,J)      +  DepDust4 * MainTimeStep
       MassDust5(I,NumSnowLayerNeg+1,J)      =  MassDust5(I,NumSnowLayerNeg+1,J)      +  DepDust5 * MainTimeStep
    endif

    ! update aerosol mass concentration in snow for each layer
    !$acc loop seq
    do LoopInd = -NumSnowLayerMax+1, 0

       SnowMass = SnowLiqWater(I,LoopInd,J) + SnowIce(I,LoopInd,J) 

       if (LoopInd >= NumSnowLayerNeg+1 .and. SnowMass > 0.0) then
          MassConcBChydropho(I,LoopInd,J) =  MassBChydropho(I,LoopInd,J) / SnowMass
          MassConcBChydrophi(I,LoopInd,J) =  MassBChydrophi(I,LoopInd,J) / SnowMass
          MassConcOChydropho(I,LoopInd,J) =  MassOChydropho(I,LoopInd,J) / SnowMass
          MassConcOChydrophi(I,LoopInd,J) =  MassOChydrophi(I,LoopInd,J) / SnowMass
          MassConcDust1(I,LoopInd,J)      =  MassDust1(I,LoopInd,J) / SnowMass
          MassConcDust2(I,LoopInd,J)      =  MassDust2(I,LoopInd,J) / SnowMass
          MassConcDust3(I,LoopInd,J)      =  MassDust3(I,LoopInd,J) / SnowMass
          MassConcDust4(I,LoopInd,J)      =  MassDust4(I,LoopInd,J) / SnowMass
          MassConcDust5(I,LoopInd,J)      =  MassDust5(I,LoopInd,J) / SnowMass
       else
          MassConcBChydropho(I,LoopInd,J) =  0.0
          MassConcBChydrophi(I,LoopInd,J) =  0.0
          MassConcOChydropho(I,LoopInd,J) =  0.0
          MassConcOChydrophi(I,LoopInd,J) =  0.0
          MassConcDust1(I,LoopInd,J)      =  0.0
          MassConcDust2(I,LoopInd,J)      =  0.0
          MassConcDust3(I,LoopInd,J)      =  0.0
          MassConcDust4(I,LoopInd,J)      =  0.0
          MassConcDust5(I,LoopInd,J)      =  0.0
          MassBChydropho(I,LoopInd,J)     =  0.0
          MassBChydrophi(I,LoopInd,J)     =  0.0
          MassOChydropho(I,LoopInd,J)     =  0.0
          MassOChydrophi(I,LoopInd,J)     =  0.0
          MassDust1(I,LoopInd,J)          =  0.0
          MassDust2(I,LoopInd,J)          =  0.0
          MassDust3(I,LoopInd,J)          =  0.0
          MassDust4(I,LoopInd,J)          =  0.0
          MassDust5(I,LoopInd,J)          =  0.0
       endif

    enddo

    ! special treatment for very shallow snowpack (NumSnowLayerNeg = 0 and SnowMass > 0.0)
    if ( NumSnowLayerNeg == 0 ) then

       SnowMass = SnowWaterEquiv

       if ( SnowMass > 0.1 ) then ! set minimum threshold (0.1 mm SWE) for computing aerosol-snow albedo
          MassBChydropho(I,0,J)     =  MassBChydropho (I,0,J) +  DepBChydropho * MainTimeStep
          MassBChydrophi(I,0,J)     =  MassBChydrophi (I,0,J) +  DepBChydrophi * MainTimeStep
          MassOChydropho(I,0,J)     =  MassOChydropho (I,0,J) +  DepOChydropho * MainTimeStep
          MassOChydrophi(I,0,J)     =  MassOChydrophi (I,0,J) +  DepOChydrophi * MainTimeStep
          MassDust1(I,0,J)          =  MassDust1 (I,0,J)      +  DepDust1 * MainTimeStep
          MassDust2(I,0,J)          =  MassDust2 (I,0,J)      +  DepDust2 * MainTimeStep
          MassDust3(I,0,J)          =  MassDust3 (I,0,J)      +  DepDust3 * MainTimeStep
          MassDust4(I,0,J)          =  MassDust4 (I,0,J)      +  DepDust4 * MainTimeStep
          MassDust5(I,0,J)          =  MassDust5 (I,0,J)      +  DepDust5 * MainTimeStep
          MassConcBChydropho(I,0,J) =  MassBChydropho(I,0,J) / SnowMass
          MassConcBChydrophi(I,0,J) =  MassBChydrophi(I,0,J) / SnowMass
          MassConcOChydropho(I,0,J) =  MassOChydropho(I,0,J) / SnowMass
          MassConcOChydrophi(I,0,J) =  MassOChydrophi(I,0,J) / SnowMass
          MassConcDust1(I,0,J)      =  MassDust1(I,0,J) / SnowMass
          MassConcDust2(I,0,J)      =  MassDust2(I,0,J) / SnowMass
          MassConcDust3(I,0,J)      =  MassDust3(I,0,J) / SnowMass
          MassConcDust4(I,0,J)      =  MassDust4(I,0,J) / SnowMass
          MassConcDust5(I,0,J)      =  MassDust5(I,0,J) / SnowMass
       else
          MassBChydropho(I,0,J)     =  0.0
          MassBChydrophi(I,0,J)     =  0.0
          MassOChydropho(I,0,J)     =  0.0
          MassOChydrophi(I,0,J)     =  0.0
          MassDust1(I,0,J)          =  0.0
          MassDust2(I,0,J)          =  0.0
          MassDust3(I,0,J)          =  0.0
          MassDust4(I,0,J)          =  0.0
          MassDust5(I,0,J)          =  0.0
          MassConcBChydropho(I,0,J) =  0.0
          MassConcBChydrophi(I,0,J) =  0.0
          MassConcOChydropho(I,0,J) =  0.0
          MassConcOChydrophi(I,0,J) =  0.0
          MassConcDust1(I,0,J)      =  0.0
          MassConcDust2(I,0,J)      =  0.0
          MassConcDust3(I,0,J)      =  0.0
          MassConcDust4(I,0,J)      =  0.0
          MassConcDust5(I,0,J)      =  0.0
       endif

    endif

    end associate

      end do
    end do
    !$acc end parallel loop

  end subroutine SnowAerosolSnicar

end module SnowAerosolSnicarMod
