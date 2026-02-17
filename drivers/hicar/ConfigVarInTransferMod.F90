module ConfigVarInTransferMod

!!! Transfer input 2-D NoahmpIO Configuration variables to 1-D column variable
!!! 1-D variables should be first defined in /src/ConfigVarType.F90
!!! 2-D variables should be first defined in NoahmpIOVarType.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type
  use NoahmpVarType

  implicit none

contains

!=== initialize with input/restart data or table values

  subroutine ConfigVarInTransfer(noahmp, NoahmpIO)

    implicit none

    type(NoahmpIO_type) , intent(inout) :: NoahmpIO
    type(noahmp_type),    intent(inout) :: noahmp

    integer :: I, J, LoopInd

    ! config namelist variable
    noahmp%config%nmlist%OptDynamicVeg               = NoahmpIO%IOPT_DVEG
    noahmp%config%nmlist%OptRainSnowPartition        = NoahmpIO%IOPT_SNF
    noahmp%config%nmlist%OptSoilWaterTranspiration   = NoahmpIO%IOPT_BTR
    noahmp%config%nmlist%OptGroundResistanceEvap     = NoahmpIO%IOPT_RSF
    noahmp%config%nmlist%OptSurfaceDrag              = NoahmpIO%IOPT_SFC
    noahmp%config%nmlist%OptStomataResistance        = NoahmpIO%IOPT_CRS
    noahmp%config%nmlist%OptSnowAlbedo               = NoahmpIO%IOPT_ALB
    noahmp%config%nmlist%OptCanopyRadiationTransfer  = NoahmpIO%IOPT_RAD
    noahmp%config%nmlist%OptSnowSoilTempTime         = NoahmpIO%IOPT_STC
    noahmp%config%nmlist%OptSnowThermConduct         = NoahmpIO%IOPT_TKSNO
    noahmp%config%nmlist%OptSoilTemperatureBottom    = NoahmpIO%IOPT_TBOT
    noahmp%config%nmlist%OptSoilSupercoolWater       = NoahmpIO%IOPT_FRZ
    noahmp%config%nmlist%OptSoilPermeabilityFrozen   = NoahmpIO%IOPT_INF
    noahmp%config%nmlist%OptDynVicInfiltration       = NoahmpIO%IOPT_INFDV
    noahmp%config%nmlist%OptTileDrainage             = NoahmpIO%IOPT_TDRN
    noahmp%config%nmlist%OptIrrigation               = NoahmpIO%IOPT_IRR
    noahmp%config%nmlist%OptIrrigationMethod         = NoahmpIO%IOPT_IRRM
    noahmp%config%nmlist%OptCropModel                = NoahmpIO%IOPT_CROP
    noahmp%config%nmlist%OptSoilProperty             = NoahmpIO%IOPT_SOIL
    noahmp%config%nmlist%OptPedotransfer             = NoahmpIO%IOPT_PEDO
    noahmp%config%nmlist%OptRunoffSurface            = NoahmpIO%IOPT_RUNSRF
    noahmp%config%nmlist%OptRunoffSubsurface         = NoahmpIO%IOPT_RUNSUB
    noahmp%config%nmlist%OptGlacierTreatment         = NoahmpIO%IOPT_GLA
    noahmp%config%nmlist%OptSnowCompaction           = NoahmpIO%IOPT_COMPACT
    noahmp%config%nmlist%OptWetlandModel             = NoahmpIO%IOPT_WETLAND
    noahmp%config%nmlist%OptSnowCoverGround          = NoahmpIO%IOPT_SCF

    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then ! SNICAR namelist
       noahmp%config%nmlist%OptSnicarSnowShape          = NoahmpIO%SNICAR_SNOWSHAPE_OPT
       noahmp%config%nmlist%OptSnicarRTSolver           = NoahmpIO%SNICAR_RTSOLVER_OPT
       noahmp%config%nmlist%OptSnicarBandNum            = NoahmpIO%SNICAR_BANDNUMBER_OPT 
       noahmp%config%nmlist%OptSnicarSolarSpec          = NoahmpIO%SNICAR_SOLARSPEC_OPT
       noahmp%config%nmlist%OptSnicarSnwOptic           = NoahmpIO%SNICAR_SNOWOPTICS_OPT
       noahmp%config%nmlist%OptSnicarDustOptic          = NoahmpIO%SNICAR_DUSTOPTICS_OPT
       noahmp%config%nmlist%FlagSnicarSnowBCIntmix      = NoahmpIO%SNICAR_SNOWBC_INTMIX
       noahmp%config%nmlist%FlagSnicarSnowDustIntmix    = NoahmpIO%SNICAR_SNOWDUST_INTMIX
       noahmp%config%nmlist%FlagSnicarUseAerosol        = NoahmpIO%SNICAR_USE_AEROSOL
       noahmp%config%nmlist%FlagSnicarUseOC             = NoahmpIO%SNICAR_USE_OC
       noahmp%config%nmlist%FlagSnicarAerosolReadTable  = NoahmpIO%SNICAR_AEROSOL_READTABLE
    endif

    ! config domain variable
    noahmp%config%domain%NumSwRadBand                = NoahmpIO%NUMRAD
    noahmp%config%domain%NumCropGrowStage            = 8
    noahmp%config%domain%FlagSoilProcess             = NoahmpIO%calculate_soil
    noahmp%config%domain%NumSoilTimeStep             = NoahmpIO%soil_update_steps
    noahmp%config%domain%NumSnowLayerMax             = NoahmpIO%NSNOW
    noahmp%config%domain%NumSoilLayer                = NoahmpIO%NSOIL
    noahmp%config%domain%MainTimeStep                = NoahmpIO%DTBL
    noahmp%config%domain%SoilTimeStep                = NoahmpIO%DTBL * NoahmpIO%soil_update_steps
    noahmp%config%domain%GridSize                    = sqrt(max(10.0,NoahmpIO%DX) * max(10.0,NoahmpIO%DY))
    noahmp%config%domain%LandUseDataName             = NoahmpIO%LLANDUSE
    noahmp%config%domain%DayJulianInYear             = NoahmpIO%JULIAN
    noahmp%config%domain%NumDayInYear                = NoahmpIO%YEARLEN
    noahmp%config%domain%IndexWaterPoint             = NoahmpIO%ISWATER_TABLE
    noahmp%config%domain%IndexBarrenPoint            = NoahmpIO%ISBARREN_TABLE
    noahmp%config%domain%IndexIcePoint               = NoahmpIO%ISICE_TABLE
    noahmp%config%domain%IndexCropPoint              = NoahmpIO%ISCROP_TABLE
    noahmp%config%domain%IndexEBLForest              = NoahmpIO%EBLFOREST_TABLE
    noahmp%config%domain%RunoffSlopeType             = NoahmpIO%SLOPETYP

    noahmp%config%domain%ITS = NoahmpIO%ITS
    noahmp%config%domain%ITE = NoahmpIO%ITE
    noahmp%config%domain%JTS = NoahmpIO%JTS
    noahmp%config%domain%JTE = NoahmpIO%JTE

    ! config domain variable
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then ! SNICAR variables
       noahmp%config%domain%NumTempSnwAgeSnicar      = NoahmpIO%idx_T_max
       noahmp%config%domain%NumTempGradSnwAgeSnicar  = NoahmpIO%idx_Tgrd_max
       noahmp%config%domain%NumDensitySnwAgeSnicar   = NoahmpIO%idx_rhos_max
       noahmp%config%domain%NumSnicarRadBand         = NoahmpIO%snicar_numrad_snw
       noahmp%config%domain%NumRadiusSnwMieSnicar    = NoahmpIO%idx_Mie_snw_mx
    endif

    associate(                                      &
              NumSnowLayerMax => NoahmpIO%NSNOW    ,&
              NumSoilLayer    => NoahmpIO%NSOIL    ,&
              ITS             => noahmp%config%domain%ITS              ,&
              ITE             => noahmp%config%domain%ITE              ,&
              JTS             => noahmp%config%domain%JTS              ,&
              JTE             => noahmp%config%domain%JTE              &
             )
    ! the following initialization cannot be done in ConfigVarInitMod
    ! because the NumSoilLayer and NumSnowLayerMax are initialized with input values in this module
    if ( .not. allocated(noahmp%config%domain%DepthSoilLayer) ) then
       allocate( noahmp%config%domain%DepthSoilLayer(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%DepthSoilLayer)
    endif
    if ( .not. allocated(noahmp%config%domain%ThicknessSoilLayer) ) then
       allocate( noahmp%config%domain%ThicknessSoilLayer(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%ThicknessSoilLayer)
    endif
    if ( .not. allocated(noahmp%config%domain%SoilType) ) then
       allocate( noahmp%config%domain%SoilType(ITS:ITE,1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%SoilType)
    endif
    if ( .not. allocated(noahmp%config%domain%ThicknessSnowSoilLayer) ) then
       allocate( noahmp%config%domain%ThicknessSnowSoilLayer(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%ThicknessSnowSoilLayer)
    endif
    if ( .not. allocated(noahmp%config%domain%DepthSnowSoilLayer) ) then
       allocate( noahmp%config%domain%DepthSnowSoilLayer(ITS:ITE,-NumSnowLayerMax+1:NumSoilLayer,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%DepthSnowSoilLayer)
    endif
    if ( .not. allocated(noahmp%config%domain%VegType) ) then
       allocate( noahmp%config%domain%VegType(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%VegType)
    endif
    if ( .not. allocated(noahmp%config%domain%CropType) ) then
       allocate( noahmp%config%domain%CropType(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%CropType)
    endif
    if ( .not. allocated(noahmp%config%domain%Latitude) ) then
       allocate( noahmp%config%domain%Latitude(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%Latitude)
    endif
    if ( .not. allocated(noahmp%config%domain%SurfaceType) ) then
       allocate( noahmp%config%domain%SurfaceType(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%SurfaceType)
    endif
    if ( .not. allocated(noahmp%config%domain%FlagUrban) ) then
       allocate( noahmp%config%domain%FlagUrban(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%FlagUrban)
    endif
    if ( .not. allocated(noahmp%config%domain%SoilColor) ) then
       allocate( noahmp%config%domain%SoilColor(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%SoilColor)
    endif
    if ( .not. allocated(noahmp%config%domain%FlagCropland) ) then
       allocate( noahmp%config%domain%FlagCropland(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%FlagCropland)
    endif
    if ( .not. allocated(noahmp%config%domain%FlagWetland) ) then
       allocate( noahmp%config%domain%FlagWetland(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%FlagWetland)
    endif
    if ( .not. allocated(noahmp%config%domain%FlagDynamicCrop) ) then
       allocate( noahmp%config%domain%FlagDynamicCrop(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%FlagDynamicCrop)
    endif
    if ( .not. allocated(noahmp%config%domain%FlagDynamicVeg) ) then
       allocate( noahmp%config%domain%FlagDynamicVeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%FlagDynamicVeg)
    endif
    if ( .not. allocated(noahmp%config%domain%NumSnowLayerNeg) ) then
       allocate( noahmp%config%domain%NumSnowLayerNeg(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%NumSnowLayerNeg)
    endif
    if ( .not. allocated(noahmp%config%domain%IndicatorIceSfc) ) then
       allocate( noahmp%config%domain%IndicatorIceSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%IndicatorIceSfc)
    endif
    if ( .not. allocated(noahmp%config%domain%DepthSoilTempBottom) ) then
       allocate( noahmp%config%domain%DepthSoilTempBottom(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%DepthSoilTempBottom)
    endif
    if ( .not. allocated(noahmp%config%domain%RefHeightAboveSfc) ) then
       allocate( noahmp%config%domain%RefHeightAboveSfc(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%RefHeightAboveSfc)
    endif
    if ( .not. allocated(noahmp%config%domain%ThicknessAtmosBotLayer) ) then
       allocate( noahmp%config%domain%ThicknessAtmosBotLayer(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%ThicknessAtmosBotLayer)
    endif
    if ( .not. allocated(noahmp%config%domain%CosSolarZenithAngle) ) then
       allocate( noahmp%config%domain%CosSolarZenithAngle(ITS:ITE,JTS:JTE) )
       !$acc enter data create(noahmp%config%domain%CosSolarZenithAngle)
    endif
   end associate
   
    !$acc update device(noahmp%config)

! --------------------------------------------------------------------- 
    associate(                                      &
              NumSnowLayerMax => NoahmpIO%NSNOW    ,&
              NumSoilLayer    => NoahmpIO%NSOIL    ,&
              ITS             => noahmp%config%domain%ITS              ,&
              ITE             => noahmp%config%domain%ITE              ,&
              JTS             => noahmp%config%domain%JTS              ,&
              JTE             => noahmp%config%domain%JTE              &
             )
! ---------------------------------------------------------------------

    !$acc parallel loop collapse(2) gang vector default(present) private(I, J) private(LoopInd)
      do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! initial setting default values of 2D config arrays -- moved from ConfigVarInitMod,
    ! since these arrays are only allocated above and not earlier

      noahmp%config%domain%FlagUrban(I,J)              = .false.
      noahmp%config%domain%FlagCropland(I,J)           = .false.
      noahmp%config%domain%FlagWetland(I,J)            = .false.
      noahmp%config%domain%FlagDynamicCrop(I,J)        = .false.
      noahmp%config%domain%FlagDynamicVeg(I,J)         = .false.
      noahmp%config%domain%NumSnowLayerNeg(I,J)        = undefined_int
      noahmp%config%domain%VegType(I,J)                = undefined_int
      noahmp%config%domain%CropType(I,J)               = undefined_int
      noahmp%config%domain%SurfaceType(I,J)            = undefined_int
      noahmp%config%domain%SoilColor(I,J)              = undefined_int
      noahmp%config%domain%IndicatorIceSfc(I,J)        = undefined_int
      noahmp%config%domain%CosSolarZenithAngle(I,J)    = undefined_real
      noahmp%config%domain%RefHeightAboveSfc(I,J)      = undefined_real
      noahmp%config%domain%ThicknessAtmosBotLayer(I,J) = undefined_real
      noahmp%config%domain%Latitude(I,J)               = undefined_real
      noahmp%config%domain%DepthSoilTempBottom(I,J)    = NoahmpIO%ZBOT_TABLE

    noahmp%config%domain%SoilColor(I,J)                   = 4
    noahmp%config%domain%NumSnowLayerNeg(I,J)             = NoahmpIO%ISNOWXY(I,J)
    noahmp%config%domain%VegType(I,J)                     = NoahmpIO%IVGTYP(I,J)
    noahmp%config%domain%CropType(I,J)                    = NoahmpIO%CROPCAT(I,J)
    noahmp%config%domain%IndicatorIceSfc(I,J)             = NoahmpIO%ICE
    noahmp%config%domain%Latitude(I,J)                    = NoahmpIO%XLAT(I,J)
    noahmp%config%domain%RefHeightAboveSfc(I,J)           = NoahmpIO%DZ8W(I,1,J)*0.5
    noahmp%config%domain%ThicknessAtmosBotLayer(I,J)      = NoahmpIO%DZ8W(I,1,J)
    noahmp%config%domain%CosSolarZenithAngle(I,J)         = NoahmpIO%COSZEN(I,J) 
    noahmp%config%domain%SurfaceType(I,J)                 = 1
    if ( (NoahmpIO%IVGTYP(I,J) == NoahmpIO%ISWATER_TABLE)) then
       noahmp%config%domain%SurfaceType(I,J) = 2
    endif
    
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
      noahmp%config%domain%SoilType              (I,LoopInd,J)   = undefined_int
      noahmp%config%domain%DepthSoilLayer        (I,LoopInd,J)   = undefined_real
      noahmp%config%domain%ThicknessSoilLayer    (I,LoopInd,J)   = undefined_real
    enddo
    !$acc loop seq
    do LoopInd = -NumSnowLayerMax+1, NumSoilLayer
      noahmp%config%domain%ThicknessSnowSoilLayer(I,LoopInd,J)   = undefined_real
      noahmp%config%domain%DepthSnowSoilLayer    (I,LoopInd,J)   = undefined_real
    enddo
    if ( noahmp%config%nmlist%OptSoilProperty == 1 .or. noahmp%config%nmlist%OptSoilProperty == 3) then
      !$acc loop seq
      do LoopInd = 1, NumSoilLayer
         noahmp%config%domain%SoilType              (I,LoopInd,J)   = NoahmpIO%ISLTYP(I,J)  ! soil type same in all layers
      enddo
    elseif ( noahmp%config%nmlist%OptSoilProperty == 2 ) then
       noahmp%config%domain%SoilType(I,1,J) = nint(NoahmpIO%SOILCL1(I,J))        ! soil type in layer1
       noahmp%config%domain%SoilType(I,2,J) = nint(NoahmpIO%SOILCL2(I,J))        ! soil type in layer2
       noahmp%config%domain%SoilType(I,3,J) = nint(NoahmpIO%SOILCL3(I,J))        ! soil type in layer3
       noahmp%config%domain%SoilType(I,4,J) = nint(NoahmpIO%SOILCL4(I,J))        ! soil type in layer4
    endif 
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
      noahmp%config%domain%DepthSoilLayer              (I,LoopInd,J)   = NoahmpIO%ZSOIL(LoopInd)  ! soil type same in all layers
    enddo
    !$acc loop seq
    do LoopInd = -NumSnowLayerMax+1, NumSoilLayer
      noahmp%config%domain%DepthSnowSoilLayer(I,LoopInd,J)   = NoahmpIO%ZSNSOXY(I,LoopInd,J)
    enddo

    ! treatment for urban point
    if ( (NoahmpIO%IVGTYP(I,J) == NoahmpIO%ISURBAN_TABLE) .or. (NoahmpIO%IVGTYP(I,J) > NoahmpIO%URBTYPE_beg) ) then
       if ( NoahmpIO%SF_URBAN_PHYSICS == 0 ) then
           noahmp%config%domain%VegType(I,J) = NoahmpIO%ISURBAN_TABLE  ! treat as bulk urban point
           noahmp%config%domain%FlagUrban(I,J) = .true.
       else
           noahmp%config%domain%VegType(I,J) = NoahmpIO%NATURAL_TABLE  ! set rural vegetation type based on table natural
                                                                  ! urban is handled by explicit urban scheme outside Noah-MP
           NoahmpIO%GVFMAX(I,J)         = 0.96 * 100.0            ! unit: %
       endif
    endif

    ! treatment for crop point
    noahmp%config%domain%CropType(I,J) = 0
    if ( (NoahmpIO%IOPT_CROP > 0) .and. (NoahmpIO%IVGTYP(I,J) == NoahmpIO%ISCROP_TABLE) ) &
       noahmp%config%domain%CropType(I,J) = NoahmpIO%DEFAULT_CROP_TABLE   
       
    if ( (NoahmpIO%IOPT_CROP > 0) .and. (NoahmpIO%CROPCAT(I,J) > 0) ) then
       noahmp%config%domain%CropType(I,J) = NoahmpIO%CROPCAT(I,J)
       noahmp%config%domain%VegType(I,J)  = NoahmpIO%ISCROP_TABLE
       NoahmpIO%VEGFRA(I,J)          = 0.95 * 100.0              ! unit: %
       NoahmpIO%GVFMAX(I,J)          = 0.95 * 100.0              ! unit: %
    endif

    ! correct inconsistent soil type
    !$acc loop seq
    do LoopInd = 1, NumSoilLayer
      if ( (noahmp%config%domain%SoilType(I,LoopInd,J) == 14) .and. (NoahmpIO%XICE(I,J) == 0.0) ) then
#ifndef _OPENACC
       write(*,*) "SOIL TYPE FOUND TO BE WATER AT A LAND-POINT"
       write(*,*) "RESET SOIL type to SANDY CLAY LOAM at grid = ", I, J
#endif
         noahmp%config%domain%SoilType(I,LoopInd,J) = 7
      endif
     enddo

    ! set warning message for inconsistent surface and subsurface runoff option
    ! for now, only the same options for surface and subsurface runoff have been tested
    if ( noahmp%config%nmlist%OptRunoffSurface /= noahmp%config%nmlist%OptRunoffSubsurface ) then
#ifndef _OPENACC
       write(*,*) "Warning: Surface and subsurface runoff options are inconsistent! They may be incompatible!"
       write(*,*) "Warning: Currently only the same options for surface and subsurface runoff are tested."
#endif
    endif

         enddo
      enddo 
    end associate

  end subroutine ConfigVarInTransfer

end module ConfigVarInTransferMod
