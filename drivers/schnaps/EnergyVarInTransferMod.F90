module EnergyVarInTransferMod

!!! Transfer input 2-D NoahmpIO Energy variables to 1-D column variable
!!! 1-D variables should be first defined in /src/EnergyVarType.F90
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

!=== initialize with input data or table values

  subroutine EnergyVarInTransfer(noahmp, NoahmpIO)

    implicit none

    type(NoahmpIO_type), intent(inout) :: NoahmpIO
    type(noahmp_type),   intent(inout) :: noahmp

    ! local loop index
    integer :: I, J
    integer :: SoilLayerIndex, LoopInd

    associate(                                                                    &
              VegType               => noahmp%config%domain%VegType              ,&
              SoilType              => noahmp%config%domain%SoilType                  ,&
              CropType              => noahmp%config%domain%CropType             ,&
              SoilColor             => noahmp%config%domain%SoilColor            ,&
              FlagUrban             => noahmp%config%domain%FlagUrban            ,&
              NumSnowLayerMax       => noahmp%config%domain%NumSnowLayerMax      ,&
              NumSoilLayer          => noahmp%config%domain%NumSoilLayer         ,&
              NumSwRadBand          => noahmp%config%domain%NumSwRadBand         ,&
              NumSnicarRadBand      => noahmp%config%domain%NumSnicarRadBand     ,&
              NumRadiusSnwMieSnicar => noahmp%config%domain%NumRadiusSnwMieSnicar &
             )

    !$acc parallel loop collapse(2) default(present) private(LoopInd, SoilLayerIndex)
      do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
         do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    ! energy state variables
    noahmp%energy%state%LeafAreaIndex(I,J)                             = NoahmpIO%LAI     (I,J)
    noahmp%energy%state%StemAreaIndex(I,J)                             = NoahmpIO%XSAIXY  (I,J)
    noahmp%energy%state%SpecHumiditySfcMean(I,J)                       = NoahmpIO%QSFC    (I,J)
    noahmp%energy%state%TemperatureGrd(I,J)                            = NoahmpIO%TGXY    (I,J)
    noahmp%energy%state%TemperatureCanopy(I,J)                         = NoahmpIO%TVXY    (I,J)
    noahmp%energy%state%SnowAgeNondim(I,J)                             = NoahmpIO%TAUSSXY (I,J)
    noahmp%energy%state%AlbedoSnowPrev(I,J)                            = NoahmpIO%ALBOLDXY(I,J)
    noahmp%energy%state%PressureVaporCanAir(I,J)                       = NoahmpIO%EAHXY   (I,J)
    noahmp%energy%state%TemperatureCanopyAir(I,J)                      = NoahmpIO%TAHXY   (I,J)
    noahmp%energy%state%ExchCoeffShSfc(I,J)                            = NoahmpIO%CHXY    (I,J) 
    noahmp%energy%state%ExchCoeffMomSfc(I,J)                           = NoahmpIO%CMXY    (I,J)
    !$acc loop seq
    do LoopInd = -NumSnowLayerMax+1, 0
       noahmp%energy%state%TemperatureSoilSnow(I,LoopInd,J) = NoahmpIO%TSNOXY  (I,LoopInd,J)
    enddo
    !$acc loop seq 
    do LoopInd = 1, NumSoilLayer
      noahmp%energy%state%TemperatureSoilSnow(I,LoopInd,J) = NoahmpIO%TSLB    (I,LoopInd,J)
    enddo

    noahmp%energy%state%PressureAtmosCO2(I,J)                          = NoahmpIO%CO2_TABLE * noahmp%forcing%PressureAirRefHeight(I,J)
    noahmp%energy%state%PressureAtmosO2(I,J)                           = NoahmpIO%O2_TABLE  * noahmp%forcing%PressureAirRefHeight(I,J)
    !$acc loop seq
    do LoopInd = 1, NumSwRadBand
       noahmp%energy%state%AlbedoSoilDir(I,LoopInd,J) = NoahmpIO%ALBSOILDIRXY(I,LoopInd,J)
       noahmp%energy%state%AlbedoSoilDif(I,LoopInd,J) = NoahmpIO%ALBSOILDIFXY(I,LoopInd,J)
    enddo

    ! vegetation treatment for USGS land types (playa, lava, sand to bare)
    if ( (VegType(I,J) == 25) .or. (VegType(I,J) == 26) .or. (VegType(I,J) == 27) ) then
       noahmp%energy%state%VegFrac(I,J)       = 0.0
       noahmp%energy%state%LeafAreaIndex(I,J) = 0.0
    endif

    ! energy flux variables
    noahmp%energy%flux%HeatGroundTotAcc(I,J)                           = NoahmpIO%ACC_SSOILXY(I,J)

    ! energy parameter variables
    noahmp%energy%param%SoilHeatCapacity(I,J)                          = NoahmpIO%CSOIL_TABLE
    noahmp%energy%param%SnowAgeFacBats(I,J)                            = NoahmpIO%TAU0_TABLE
    noahmp%energy%param%SnowGrowVapFacBats(I,J)                        = NoahmpIO%GRAIN_GROWTH_TABLE
    noahmp%energy%param%SnowSootFacBats(I,J)                           = NoahmpIO%DIRT_SOOT_TABLE
    noahmp%energy%param%SnowGrowFrzFacBats(I,J)                        = NoahmpIO%EXTRA_GROWTH_TABLE
    noahmp%energy%param%SolarZenithAdjBats(I,J)                        = NoahmpIO%BATS_COSZ_TABLE
    noahmp%energy%param%FreshSnoAlbVisBats(I,J)                        = NoahmpIO%BATS_VIS_NEW_TABLE
    noahmp%energy%param%FreshSnoAlbNirBats(I,J)                        = NoahmpIO%BATS_NIR_NEW_TABLE
    noahmp%energy%param%SnoAgeFacDifVisBats(I,J)                       = NoahmpIO%BATS_VIS_AGE_TABLE
    noahmp%energy%param%SnoAgeFacDifNirBats(I,J)                       = NoahmpIO%BATS_NIR_AGE_TABLE
    noahmp%energy%param%SzaFacDirVisBats(I,J)                          = NoahmpIO%BATS_VIS_DIR_TABLE
    noahmp%energy%param%SzaFacDirNirBats(I,J)                          = NoahmpIO%BATS_NIR_DIR_TABLE
    noahmp%energy%param%SnowAlbRefClass(I,J)                           = NoahmpIO%CLASS_ALB_REF_TABLE
    noahmp%energy%param%SnowAgeFacClass(I,J)                           = NoahmpIO%CLASS_SNO_AGE_TABLE
    noahmp%energy%param%SnowAlbFreshClass(I,J)                         = NoahmpIO%CLASS_ALB_NEW_TABLE
    noahmp%energy%param%UpscatterCoeffSnowDir(I,J)                     = NoahmpIO%BETADS_TABLE
    noahmp%energy%param%UpscatterCoeffSnowDif(I,J)                     = NoahmpIO%BETAIS_TABLE
    noahmp%energy%param%ZilitinkevichCoeff(I,J)                        = NoahmpIO%CZIL_TABLE
    noahmp%energy%param%EmissivitySnow(I,J)                            = NoahmpIO%SNOW_EMIS_TABLE
    !$acc loop seq
    do LoopInd = 1, 2
      noahmp%energy%param%EmissivitySoilLake(I,LoopInd,J)              = NoahmpIO%EG_TABLE(LoopInd)
    enddo
    !$acc loop seq
    do LoopInd = 1, NumSwRadBand
      noahmp%energy%param%AlbedoLandIce(I,LoopInd,J)                   = NoahmpIO%ALBICE_TABLE(LoopInd)
    enddo
    noahmp%energy%param%RoughLenMomSnow(I,J)                           = NoahmpIO%Z0SNO_TABLE
    noahmp%energy%param%RoughLenMomSoil(I,J)                           = NoahmpIO%Z0SOIL_TABLE
    noahmp%energy%param%RoughLenMomLake(I,J)                           = NoahmpIO%Z0LAKE_TABLE
    noahmp%energy%param%EmissivityIceSfc(I,J)                          = NoahmpIO%EICE_TABLE
    noahmp%energy%param%ResistanceSoilExp(I,J)                         = NoahmpIO%RSURF_EXP_TABLE
    noahmp%energy%param%ResistanceSnowSfc(I,J)                         = NoahmpIO%RSURF_SNOW_TABLE
    noahmp%energy%param%VegFracAnnMax(I,J)                             = NoahmpIO%GVFMAX(I,J) / 100.0
    noahmp%energy%param%VegFracGreen(I,J)                              = NoahmpIO%VEGFRA(I,J) / 100.0
    noahmp%energy%param%TreeCrownRadius(I,J)                           = NoahmpIO%RC_TABLE    (VegType(I,J))
    noahmp%energy%param%HeightCanopyTop(I,J)                           = NoahmpIO%HVT_TABLE   (VegType(I,J))
    noahmp%energy%param%HeightCanopyBot(I,J)                           = NoahmpIO%HVB_TABLE   (VegType(I,J))
    noahmp%energy%param%RoughLenMomVeg(I,J)                            = NoahmpIO%Z0MVT_TABLE (VegType(I,J))
    noahmp%energy%param%CanopyWindExtFac(I,J)                          = NoahmpIO%CWPVT_TABLE (VegType(I,J))
    noahmp%energy%param%TreeDensity(I,J)                               = NoahmpIO%DEN_TABLE   (VegType(I,J))
    noahmp%energy%param%CanopyOrientIndex(I,J)                         = NoahmpIO%XL_TABLE    (VegType(I,J))
    noahmp%energy%param%ConductanceLeafMin(I,J)                        = NoahmpIO%BP_TABLE    (VegType(I,J))
    noahmp%energy%param%Co2MmConst25C(I,J)                             = NoahmpIO%KC25_TABLE  (VegType(I,J))
    noahmp%energy%param%O2MmConst25C(I,J)                              = NoahmpIO%KO25_TABLE  (VegType(I,J))
    noahmp%energy%param%Co2MmConstQ10(I,J)                             = NoahmpIO%AKC_TABLE   (VegType(I,J))
    noahmp%energy%param%O2MmConstQ10(I,J)                              = NoahmpIO%AKO_TABLE   (VegType(I,J))
    noahmp%energy%param%RadiationStressFac(I,J)                        = NoahmpIO%RGL_TABLE   (VegType(I,J))
    noahmp%energy%param%ResistanceStomataMin(I,J)                      = NoahmpIO%RS_TABLE    (VegType(I,J))
    noahmp%energy%param%ResistanceStomataMax(I,J)                      = NoahmpIO%RSMAX_TABLE (VegType(I,J))
    noahmp%energy%param%AirTempOptimTransp(I,J)                        = NoahmpIO%TOPT_TABLE  (VegType(I,J))
    noahmp%energy%param%VaporPresDeficitFac(I,J)                       = NoahmpIO%HS_TABLE    (VegType(I,J))
    noahmp%energy%param%LeafDimLength(I,J)                             = NoahmpIO%DLEAF_TABLE (VegType(I,J))
    noahmp%energy%param%HeatCapacCanFac(I,J)                           = NoahmpIO%CBIOM_TABLE (VegType(I,J))
    !$acc loop seq
    do LoopInd = 1, 12
      noahmp%energy%param%LeafAreaIndexMon (I,LoopInd,J)                   = NoahmpIO%LAIM_TABLE  (VegType(I,J),LoopInd)
      noahmp%energy%param%StemAreaIndexMon (I,LoopInd,J)                   = NoahmpIO%SAIM_TABLE  (VegType(I,J),LoopInd)
    end do
    !$acc loop seq
    do LoopInd = 1, NumSwRadBand
      noahmp%energy%param%ReflectanceLeaf  (I,LoopInd,J)         = NoahmpIO%RHOL_TABLE  (VegType(I,J),LoopInd)
      noahmp%energy%param%ReflectanceStem  (I,LoopInd,J)         = NoahmpIO%RHOS_TABLE  (VegType(I,J),LoopInd)
      noahmp%energy%param%TransmittanceLeaf(I,LoopInd,J)         = NoahmpIO%TAUL_TABLE  (VegType(I,J),LoopInd)
      noahmp%energy%param%TransmittanceStem(I,LoopInd,J)         = NoahmpIO%TAUS_TABLE  (VegType(I,J),LoopInd)
      noahmp%energy%param%AlbedoSoilSat    (I,LoopInd,J)         = NoahmpIO%ALBSAT_TABLE(SoilColor(I,J),LoopInd)
      noahmp%energy%param%AlbedoSoilDry    (I,LoopInd,J)         = NoahmpIO%ALBDRY_TABLE(SoilColor(I,J),LoopInd)
      noahmp%energy%param%AlbedoLakeFrz    (I,LoopInd,J)         = NoahmpIO%ALBLAK_TABLE(LoopInd)
      noahmp%energy%param%ScatterCoeffSnow (I,LoopInd,J)         = NoahmpIO%OMEGAS_TABLE(LoopInd)
    end do
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then ! SNICAR variables
       !$acc loop seq
       do LoopInd = 1, NumSnicarRadBand
         noahmp%energy%param%RadSwWgtDif        (I,LoopInd,J) = NoahmpIO%flx_wgt_dif(LoopInd)
         noahmp%energy%param%RadSwWgtDir        (I,LoopInd,J) = NoahmpIO%flx_wgt_dir(LoopInd)
         noahmp%energy%param%SsAlbBCphi         (I,LoopInd,J) = NoahmpIO%ss_alb_bc1       (LoopInd) 
         noahmp%energy%param%AsyPrmBCphi        (I,LoopInd,J) = NoahmpIO%asm_prm_bc1      (LoopInd)
         noahmp%energy%param%ExtCffMassBCphi    (I,LoopInd,J) = NoahmpIO%ext_cff_mss_bc1  (LoopInd)
         noahmp%energy%param%SsAlbBCpho         (I,LoopInd,J) = NoahmpIO%ss_alb_bc2       (LoopInd)
         noahmp%energy%param%AsyPrmBCpho        (I,LoopInd,J) = NoahmpIO%asm_prm_bc2      (LoopInd)
         noahmp%energy%param%ExtCffMassBCpho    (I,LoopInd,J) = NoahmpIO%ext_cff_mss_bc2  (LoopInd)
         noahmp%energy%param%SsAlbOCphi         (I,LoopInd,J) = NoahmpIO%ss_alb_oc1       (LoopInd)
         noahmp%energy%param%AsyPrmOCphi        (I,LoopInd,J) = NoahmpIO%asm_prm_oc1      (LoopInd)
         noahmp%energy%param%ExtCffMassOCphi    (I,LoopInd,J) = NoahmpIO%ext_cff_mss_oc1  (LoopInd)
         noahmp%energy%param%SsAlbOCpho         (I,LoopInd,J) = NoahmpIO%ss_alb_oc2       (LoopInd)
         noahmp%energy%param%AsyPrmOCpho        (I,LoopInd,J) = NoahmpIO%asm_prm_oc2      (LoopInd)
         noahmp%energy%param%ExtCffMassOCpho    (I,LoopInd,J) = NoahmpIO%ext_cff_mss_oc2  (LoopInd)
         noahmp%energy%param%SsAlbDustB1        (I,LoopInd,J) = NoahmpIO%ss_alb_dst1      (LoopInd)
         noahmp%energy%param%AsyPrmDustB1       (I,LoopInd,J) = NoahmpIO%asm_prm_dst1     (LoopInd)
         noahmp%energy%param%ExtCffMassDustB1   (I,LoopInd,J) = NoahmpIO%ext_cff_mss_dst1 (LoopInd)
         noahmp%energy%param%SsAlbDustB2        (I,LoopInd,J) = NoahmpIO%ss_alb_dst2      (LoopInd)
         noahmp%energy%param%AsyPrmDustB2       (I,LoopInd,J) = NoahmpIO%asm_prm_dst2     (LoopInd)
         noahmp%energy%param%ExtCffMassDustB2   (I,LoopInd,J) = NoahmpIO%ext_cff_mss_dst2 (LoopInd)
         noahmp%energy%param%SsAlbDustB3        (I,LoopInd,J) = NoahmpIO%ss_alb_dst3      (LoopInd)
         noahmp%energy%param%AsyPrmDustB3       (I,LoopInd,J) = NoahmpIO%asm_prm_dst3     (LoopInd)
         noahmp%energy%param%ExtCffMassDustB3   (I,LoopInd,J) = NoahmpIO%ext_cff_mss_dst3 (LoopInd)
         noahmp%energy%param%SsAlbDustB4        (I,LoopInd,J) = NoahmpIO%ss_alb_dst4      (LoopInd)
         noahmp%energy%param%AsyPrmDustB4       (I,LoopInd,J) = NoahmpIO%asm_prm_dst4     (LoopInd)
         noahmp%energy%param%ExtCffMassDustB4   (I,LoopInd,J) = NoahmpIO%ext_cff_mss_dst4 (LoopInd)
         noahmp%energy%param%SsAlbDustB5        (I,LoopInd,J) = NoahmpIO%ss_alb_dst5      (LoopInd)
         noahmp%energy%param%AsyPrmDustB5       (I,LoopInd,J) = NoahmpIO%asm_prm_dst5     (LoopInd)
         noahmp%energy%param%ExtCffMassDustB5   (I,LoopInd,J) = NoahmpIO%ext_cff_mss_dst5 (LoopInd)
         do SoilLayerIndex = 1, NumRadiusSnwMieSnicar
            noahmp%energy%param%SsAlbSnwRadDir     (I,SoilLayerIndex,LoopInd,J) = &
                        NoahmpIO%ss_alb_snw_drc     (SoilLayerIndex,LoopInd)
            noahmp%energy%param%AsyPrmSnwRadDir    (I,SoilLayerIndex,LoopInd,J) = &
                        NoahmpIO%asm_prm_snw_drc    (SoilLayerIndex,LoopInd)
            noahmp%energy%param%ExtCffMassSnwRadDir(I,SoilLayerIndex,LoopInd,J) = &
                        NoahmpIO%ext_cff_mss_snw_drc(SoilLayerIndex,LoopInd)
            noahmp%energy%param%SsAlbSnwRadDif     (I,SoilLayerIndex,LoopInd,J) = &
                        NoahmpIO%ss_alb_snw_dfs     (SoilLayerIndex,LoopInd)
            noahmp%energy%param%AsyPrmSnwRadDif    (I,SoilLayerIndex,LoopInd,J) = &
                        NoahmpIO%asm_prm_snw_dfs    (SoilLayerIndex,LoopInd)
            noahmp%energy%param%ExtCffMassSnwRadDif(I,SoilLayerIndex,LoopInd,J) = &
                        NoahmpIO%ext_cff_mss_snw_dfs(SoilLayerIndex,LoopInd)
       end do 
       end do
    endif

   !$acc loop seq
    do SoilLayerIndex = 1, NumSoilLayer
       noahmp%energy%param%SoilQuartzFrac(I,SoilLayerIndex,J) = NoahmpIO%QUARTZ_TABLE(SoilType(I,SoilLayerIndex,J))
    enddo

    ! spatial varying soil input
    if ( noahmp%config%nmlist%OptSoilProperty == 4 ) then
       !$acc loop seq
       do SoilLayerIndex = 1, NumSoilLayer
          noahmp%energy%param%SoilQuartzFrac(I,SoilLayerIndex,J) = NoahmpIO%QUARTZ_3D(I,SoilLayerIndex,J)
       enddo
    endif

    if ( FlagUrban(I,J) .eqv. .true. ) noahmp%energy%param%SoilHeatCapacity(I,J) = 3.0e6

    if ( CropType(I,J) > 0 ) then
       noahmp%energy%param%ConductanceLeafMin(I,J)             = NoahmpIO%BPI_TABLE  (CropType(I,J))
       noahmp%energy%param%Co2MmConst25C(I,J)                  = NoahmpIO%KC25I_TABLE(CropType(I,J))
       noahmp%energy%param%O2MmConst25C(I,J)                   = NoahmpIO%KO25I_TABLE(CropType(I,J))
       noahmp%energy%param%Co2MmConstQ10(I,J)                  = NoahmpIO%AKCI_TABLE (CropType(I,J))
       noahmp%energy%param%O2MmConstQ10(I,J)                   = NoahmpIO%AKOI_TABLE (CropType(I,J))
    endif


         end do
      end do
    !$acc end parallel loop


    end associate

  end subroutine EnergyVarInTransfer

end module EnergyVarInTransferMod
