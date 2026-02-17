module PhenologyMainMod

!!! Main Phenology module to estimate vegetation phenology
!!! considering vegeation canopy being buries by snow and evolution in time

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine PhenologyMain (noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: PHENOLOGY
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: I, J                    ! grid indices
    integer                          :: IntpMonth1,IntpMonth2   ! interpolation months
    real(kind=kind_noahmp)           :: ThicknessCanBury        ! thickness of canopy buried by snow [m]
    real(kind=kind_noahmp)           :: SnowDepthVegBury        ! critical snow depth at which short vege is fully covered by snow
    real(kind=kind_noahmp)           :: DayCurrent              ! current day of year (0<=DayCurrent<NumDayInYear)
    real(kind=kind_noahmp)           :: IntpWgt1,IntpWgt2       ! interpolation weights
    real(kind=kind_noahmp)           :: MonthCurrent            ! current month (1.00, ..., 12.00)

!---------------------------------------------------------------------
    associate(                                                                       &
              OptDynamicVeg          => noahmp%config%nmlist%OptDynamicVeg          ,& ! in,    dynamic vegetation option
              OptCropModel           => noahmp%config%nmlist%OptCropModel           ,& ! in,    crop model option
              VegType                => noahmp%config%domain%VegType                ,& ! in,    vegetation type
              CropType               => noahmp%config%domain%CropType               ,& ! in,    crop type
              IndexIcePoint          => noahmp%config%domain%IndexIcePoint          ,& ! in,    land ice flag
              IndexBarrenPoint       => noahmp%config%domain%IndexBarrenPoint       ,& ! in,    bare soil flag
              IndexWaterPoint        => noahmp%config%domain%IndexWaterPoint        ,& ! in,    water point flag
              FlagUrban              => noahmp%config%domain%FlagUrban              ,& ! in,    urban point flag
              FlagDynamicVeg         => noahmp%config%domain%FlagDynamicVeg         ,& ! in,    flag to activate dynamic vegetation model
              FlagDynamicCrop        => noahmp%config%domain%FlagDynamicCrop        ,& ! in,    flag to activate dynamic crop model
              Latitude               => noahmp%config%domain%Latitude               ,& ! in,    latitude [deg]
              NumDayInYear           => noahmp%config%domain%NumDayInYear                ,& ! in,    Number of days in the particular year
              DayJulianInYear        => noahmp%config%domain%DayJulianInYear             ,& ! in,    Julian day of year
              HeightCanopyTop        => noahmp%energy%param%HeightCanopyTop         ,& ! in,    top of canopy [m]
              HeightCanopyBot        => noahmp%energy%param%HeightCanopyBot         ,& ! in,    bottom of canopy [m]
              LeafAreaIndexMon       => noahmp%energy%param%LeafAreaIndexMon        ,& ! in,    monthly leaf area index, one-sided
              StemAreaIndexMon       => noahmp%energy%param%StemAreaIndexMon        ,& ! in,    monthly stem area index, one-sided
              VegFracAnnMax          => noahmp%energy%param%VegFracAnnMax           ,& ! in,    annual maximum vegetation fraction
              VegFracGreen           => noahmp%energy%param%VegFracGreen            ,& ! in,    green vegetation fraction
              TemperatureMinPhotosyn => noahmp%biochem%param%TemperatureMinPhotosyn ,& ! in,    minimum temperature for photosynthesis [K]
              PlantGrowStage         => noahmp%biochem%state%PlantGrowStage         ,& ! in,    plant growing stage
              SnowDepth              => noahmp%water%state%SnowDepth                ,& ! in,    snow depth [m]
              TemperatureCanopy      => noahmp%energy%state%TemperatureCanopy       ,& ! in,    vegetation temperature [K]
              LeafAreaIndex          => noahmp%energy%state%LeafAreaIndex           ,& ! inout, LeafAreaIndex, unadjusted for burying by snow
              StemAreaIndex          => noahmp%energy%state%StemAreaIndex           ,& ! inout, StemAreaIndex, unadjusted for burying by snow
              LeafAreaIndEff         => noahmp%energy%state%LeafAreaIndEff          ,& ! out,   leaf area index, after burying by snow
              StemAreaIndEff         => noahmp%energy%state%StemAreaIndEff          ,& ! out,   stem area index, after burying by snow
              VegFrac                => noahmp%energy%state%VegFrac                 ,& ! out,   green vegetation fraction
              CanopyFracSnowBury     => noahmp%energy%state%CanopyFracSnowBury      ,& ! out,   fraction of canopy buried by snow
              IndexGrowSeason        => noahmp%biochem%state%IndexGrowSeason         & ! out,   growing season index (0=off, 1=on)
             )

   !$acc parallel loop collapse(2) gang vector default(present) &
   !$acc private(IntpMonth1, IntpMonth2, ThicknessCanBury, SnowDepthVegBury, DayCurrent, IntpWgt1, IntpWgt2, MonthCurrent)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE


    ! compute LeafAreaIndex based on dynamic vegetation option
    if ( CropType(I,J) == 0 ) then

       ! no dynamic vegetation, use table LeafAreaIndex
       if ( (OptDynamicVeg == 1) .or. (OptDynamicVeg == 3) .or. (OptDynamicVeg == 4) ) then
          if ( Latitude(I,J) >= 0.0 ) then
            ! Northern Hemisphere
            DayCurrent = DayJulianInYear
          else
            ! Southern Hemisphere.  DayCurrent is shifted by 1/2 year.
            DayCurrent = mod(DayJulianInYear+(0.5*NumDayInYear), real(NumDayInYear))
          endif
          ! interpolate from monthly data to target time point
          MonthCurrent = 12.0 * DayCurrent / real(NumDayInYear)
          IntpMonth1   = MonthCurrent + 0.5
          IntpMonth2   = IntpMonth1 + 1
          IntpWgt1     = (IntpMonth1 + 0.5) - MonthCurrent
          IntpWgt2     = 1.0 - IntpWgt1
          if ( IntpMonth1 <  1 ) IntpMonth1 = 12
          if ( IntpMonth2 > 12 ) IntpMonth2 = 1
          LeafAreaIndex(I,J) = IntpWgt1 * LeafAreaIndexMon(I,IntpMonth1,J) + IntpWgt2 * LeafAreaIndexMon(I,IntpMonth2,J)
          StemAreaIndex(I,J) = IntpWgt1 * StemAreaIndexMon(I,IntpMonth1,J) + IntpWgt2 * StemAreaIndexMon(I,IntpMonth2,J)
       endif

       ! no dynamic vegetation, use input LeafAreaIndex time series
       if ( (OptDynamicVeg == 7) .or. (OptDynamicVeg == 8) .or. (OptDynamicVeg == 9) ) then
          StemAreaIndex(I,J) = max(0.05, 0.1*LeafAreaIndex(I,J))                                 ! set StemAreaIndex(I,J) to 10% LeafAreaIndex(I,J), but not below 0.05 MB: v3.8
          if ( LeafAreaIndex(I,J) < 0.05 ) StemAreaIndex(I,J) = 0.0                              ! if LeafAreaIndex(I,J) below minimum, make sure StemAreaIndex(I,J) = 0
       endif
       if ( StemAreaIndex(I,J) < 0.05 ) StemAreaIndex(I,J) = 0.0                                 ! MB: StemAreaIndex(I,J) CHECK, change to 0.05 v3.6
       if ( (LeafAreaIndex(I,J) < 0.05) .or. (StemAreaIndex(I,J) == 0.0) ) LeafAreaIndex(I,J) = 0.0   ! MB: LeafAreaIndex(I,J) CHECK

       ! for non-vegetation point
       if ( (VegType(I,J) == IndexWaterPoint) .or. (VegType(I,J) == IndexBarrenPoint) .or. &
            (VegType(I,J) == IndexIcePoint  ) .or. (FlagUrban(I,J) .eqv. .true.) ) then
          LeafAreaIndex(I,J) = 0.0
          StemAreaIndex(I,J) = 0.0
       endif

    endif   ! CropType(I,J) == 0

    ! vegetation fraction buried by snow
    ThicknessCanBury      = min(max(SnowDepth(I,J)-HeightCanopyBot(I,J),0.0), (HeightCanopyTop(I,J)-HeightCanopyBot(I,J)))
    CanopyFracSnowBury(I,J)    = ThicknessCanBury / max(1.0e-06, (HeightCanopyTop(I,J)-HeightCanopyBot(I,J)))           ! snow buried fraction
    if ( (HeightCanopyTop(I,J) > 0.0) .and. (HeightCanopyTop(I,J) <= 1.0) ) then                                   ! MB: change to 1.0 & 0.2 to reflect changes to HeightCanopyTop(I,J) in MPTABLE
       SnowDepthVegBury   = HeightCanopyTop(I,J) * exp(-min(SnowDepth(I,J),10.0) / 0.2)     ! CH: add limit to avoid numerical floating issue
       CanopyFracSnowBury(I,J) = min(SnowDepth(I,J), SnowDepthVegBury) / SnowDepthVegBury
    endif

    ! adjust LeafAreaIndex and StemAreaIndex bused on snow bury
    LeafAreaIndEff(I,J) = LeafAreaIndex(I,J) * (1.0 - CanopyFracSnowBury(I,J))
    StemAreaIndEff(I,J) = StemAreaIndex(I,J) * (1.0 - CanopyFracSnowBury(I,J))
    if ( (StemAreaIndEff(I,J) < 0.05) .and. (CropType(I,J) == 0) ) StemAreaIndEff(I,J) = 0.0                  ! MB: StemAreaIndEff(I,J) CHECK, change to 0.05 v3.6
    if ( ((LeafAreaIndEff(I,J) < 0.05) .or. (StemAreaIndEff(I,J) == 0.0)) .and. (CropType(I,J) == 0) ) &
       LeafAreaIndEff(I,J) = 0.0                                                                    ! MB: LeafAreaIndex(I,J) CHECK

    ! set growing season flag
    if ( ((TemperatureCanopy(I,J) > TemperatureMinPhotosyn(I,J)) .and. (CropType(I,J) == 0)) .or. &
         ((PlantGrowStage(I,J) > 2) .and. (PlantGrowStage(I,J) < 7) .and. (CropType(I,J) > 0))) then
       IndexGrowSeason(I,J) = 1.0
    else
       IndexGrowSeason(I,J) = 0.0
    endif 

    ! compute vegetation fraction
    ! input green vegetation fraction should be consistent with LeafAreaIndex
    if ( (OptDynamicVeg == 1) .or. (OptDynamicVeg == 6) .or. (OptDynamicVeg == 7) ) then      ! use VegFrac(I,J) = VegFracGreen(I,J) from input
       VegFrac(I,J) = VegFracGreen(I,J)
    elseif ( (OptDynamicVeg == 2) .or. (OptDynamicVeg == 3) .or. (OptDynamicVeg == 8) ) then  ! computed VegFrac(I,J) from LeafAreaIndex(I,J) & StemAreaIndex(I,J)
       VegFrac(I,J) = 1.0 - exp(-0.52 * (LeafAreaIndex(I,J) + StemAreaIndex(I,J)))
    elseif ( (OptDynamicVeg == 4) .or. (OptDynamicVeg == 5) .or. (OptDynamicVeg == 9) ) then  ! use yearly maximum vegetation fraction
       VegFrac(I,J) = VegFracAnnMax(I,J)
    else                                                                                      ! outside existing vegetation options
#ifdef _OPENACC
       ! GPU: Cannot write/stop in parallel region
#else
       write(*,*) "Un-recognized dynamic vegetation option (OptDynamicVeg)... "
       stop "Error: Namelist parameter OptDynamicVeg unknown"
#endif
    endif
    ! use maximum vegetation fraction for crop run
    if ( (OptCropModel > 0) .and. (CropType(I,J) > 0) ) then
       VegFrac(I,J) = VegFracAnnMax(I,J)
    endif

    ! adjust unreasonable vegetation fraction
    if ( VegFrac(I,J) <= 0.05 ) VegFrac(I,J) = 0.05
    if ( (FlagUrban(I,J) .eqv. .true.) .or. (VegType(I,J) == IndexBarrenPoint) ) VegFrac(I,J) = 0.0
    if ( (LeafAreaIndEff(I,J)+StemAreaIndEff(I,J)) == 0.0 ) VegFrac(I,J) = 0.0

    ! determine if activate dynamic vegetation or crop run
    FlagDynamicCrop(I,J) = .false.
    FlagDynamicVeg(I,J)  = .false.
    if ( (OptDynamicVeg == 2) .or. (OptDynamicVeg == 5) .or. (OptDynamicVeg == 6) ) &
       FlagDynamicVeg(I,J) = .true.
    if ( (OptCropModel > 0) .and. (CropType(I,J) > 0) ) then
       FlagDynamicCrop(I,J) = .true.
       FlagDynamicVeg(I,J)  = .false.
    endif


      end do
    end do
   !$acc end parallel loop


    end associate

  end subroutine PhenologyMain

end module PhenologyMainMod
