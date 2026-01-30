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
   !$acc parallel loop collapse(2) gang vector present(noahmp) &
   !$acc private(IntpMonth1, IntpMonth2, ThicknessCanBury, SnowDepthVegBury, DayCurrent, IntpWgt1, IntpWgt2, MonthCurrent)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

    associate(                                                                       &
              OptDynamicVeg          => noahmp%config%nmlist%OptDynamicVeg          ,& ! in,    dynamic vegetation option
              OptCropModel           => noahmp%config%nmlist%OptCropModel           ,& ! in,    crop model option
              VegType                => noahmp%config%domain%VegType(I,J)                ,& ! in,    vegetation type
              CropType               => noahmp%config%domain%CropType(I,J)               ,& ! in,    crop type
              IndexIcePoint          => noahmp%config%domain%IndexIcePoint          ,& ! in,    land ice flag
              IndexBarrenPoint       => noahmp%config%domain%IndexBarrenPoint       ,& ! in,    bare soil flag
              IndexWaterPoint        => noahmp%config%domain%IndexWaterPoint        ,& ! in,    water point flag
              FlagUrban              => noahmp%config%domain%FlagUrban(I,J)              ,& ! in,    urban point flag
              FlagDynamicVeg         => noahmp%config%domain%FlagDynamicVeg(I,J)         ,& ! in,    flag to activate dynamic vegetation model
              FlagDynamicCrop        => noahmp%config%domain%FlagDynamicCrop(I,J)        ,& ! in,    flag to activate dynamic crop model
              Latitude               => noahmp%config%domain%Latitude(I,J)               ,& ! in,    latitude [deg]
              NumDayInYear           => noahmp%config%domain%NumDayInYear                ,& ! in,    Number of days in the particular year
              DayJulianInYear        => noahmp%config%domain%DayJulianInYear             ,& ! in,    Julian day of year
              HeightCanopyTop        => noahmp%energy%param%HeightCanopyTop(I,J)         ,& ! in,    top of canopy [m]
              HeightCanopyBot        => noahmp%energy%param%HeightCanopyBot(I,J)         ,& ! in,    bottom of canopy [m]
              LeafAreaIndexMon       => noahmp%energy%param%LeafAreaIndexMon        ,& ! in,    monthly leaf area index, one-sided
              StemAreaIndexMon       => noahmp%energy%param%StemAreaIndexMon        ,& ! in,    monthly stem area index, one-sided
              VegFracAnnMax          => noahmp%energy%param%VegFracAnnMax(I,J)           ,& ! in,    annual maximum vegetation fraction
              VegFracGreen           => noahmp%energy%param%VegFracGreen(I,J)            ,& ! in,    green vegetation fraction
              TemperatureMinPhotosyn => noahmp%biochem%param%TemperatureMinPhotosyn(I,J) ,& ! in,    minimum temperature for photosynthesis [K]
              PlantGrowStage         => noahmp%biochem%state%PlantGrowStage(I,J)         ,& ! in,    plant growing stage
              SnowDepth              => noahmp%water%state%SnowDepth(I,J)                ,& ! in,    snow depth [m]
              TemperatureCanopy      => noahmp%energy%state%TemperatureCanopy(I,J)       ,& ! in,    vegetation temperature [K]
              LeafAreaIndex          => noahmp%energy%state%LeafAreaIndex(I,J)           ,& ! inout, LeafAreaIndex, unadjusted for burying by snow
              StemAreaIndex          => noahmp%energy%state%StemAreaIndex(I,J)           ,& ! inout, StemAreaIndex, unadjusted for burying by snow
              LeafAreaIndEff         => noahmp%energy%state%LeafAreaIndEff(I,J)          ,& ! out,   leaf area index, after burying by snow
              StemAreaIndEff         => noahmp%energy%state%StemAreaIndEff(I,J)          ,& ! out,   stem area index, after burying by snow
              VegFrac                => noahmp%energy%state%VegFrac(I,J)                 ,& ! out,   green vegetation fraction
              CanopyFracSnowBury     => noahmp%energy%state%CanopyFracSnowBury(I,J)      ,& ! out,   fraction of canopy buried by snow
              IndexGrowSeason        => noahmp%biochem%state%IndexGrowSeason(I,J)         & ! out,   growing season index (0=off, 1=on)
             )
!----------------------------------------------------------------------

    ! compute LeafAreaIndex based on dynamic vegetation option
    if ( CropType == 0 ) then

       ! no dynamic vegetation, use table LeafAreaIndex
       if ( (OptDynamicVeg == 1) .or. (OptDynamicVeg == 3) .or. (OptDynamicVeg == 4) ) then
          if ( Latitude >= 0.0 ) then
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
          LeafAreaIndex = IntpWgt1 * LeafAreaIndexMon(I,IntpMonth1,J) + IntpWgt2 * LeafAreaIndexMon(I,IntpMonth2,J)
          StemAreaIndex = IntpWgt1 * StemAreaIndexMon(I,IntpMonth1,J) + IntpWgt2 * StemAreaIndexMon(I,IntpMonth2,J)
       endif

       ! no dynamic vegetation, use input LeafAreaIndex time series
       if ( (OptDynamicVeg == 7) .or. (OptDynamicVeg == 8) .or. (OptDynamicVeg == 9) ) then
          StemAreaIndex = max(0.05, 0.1*LeafAreaIndex)                                 ! set StemAreaIndex to 10% LeafAreaIndex, but not below 0.05 MB: v3.8
          if ( LeafAreaIndex < 0.05 ) StemAreaIndex = 0.0                              ! if LeafAreaIndex below minimum, make sure StemAreaIndex = 0
       endif
       if ( StemAreaIndex < 0.05 ) StemAreaIndex = 0.0                                 ! MB: StemAreaIndex CHECK, change to 0.05 v3.6
       if ( (LeafAreaIndex < 0.05) .or. (StemAreaIndex == 0.0) ) LeafAreaIndex = 0.0   ! MB: LeafAreaIndex CHECK

       ! for non-vegetation point
       if ( (VegType == IndexWaterPoint) .or. (VegType == IndexBarrenPoint) .or. &
            (VegType == IndexIcePoint  ) .or. (FlagUrban .eqv. .true.) ) then
          LeafAreaIndex = 0.0
          StemAreaIndex = 0.0
       endif

    endif   ! CropType == 0

    ! vegetation fraction buried by snow
    ThicknessCanBury      = min(max(SnowDepth-HeightCanopyBot,0.0), (HeightCanopyTop-HeightCanopyBot))
    CanopyFracSnowBury    = ThicknessCanBury / max(1.0e-06, (HeightCanopyTop-HeightCanopyBot))           ! snow buried fraction
    if ( (HeightCanopyTop > 0.0) .and. (HeightCanopyTop <= 1.0) ) then                                   ! MB: change to 1.0 & 0.2 to reflect changes to HeightCanopyTop in MPTABLE
       SnowDepthVegBury   = HeightCanopyTop * exp(-min(SnowDepth,10.0) / 0.2)     ! CH: add limit to avoid numerical floating issue
       CanopyFracSnowBury = min(SnowDepth, SnowDepthVegBury) / SnowDepthVegBury
    endif

    ! adjust LeafAreaIndex and StemAreaIndex bused on snow bury
    LeafAreaIndEff = LeafAreaIndex * (1.0 - CanopyFracSnowBury)
    StemAreaIndEff = StemAreaIndex * (1.0 - CanopyFracSnowBury)
    if ( (StemAreaIndEff < 0.05) .and. (CropType == 0) ) StemAreaIndEff = 0.0                  ! MB: StemAreaIndEff CHECK, change to 0.05 v3.6
    if ( ((LeafAreaIndEff < 0.05) .or. (StemAreaIndEff == 0.0)) .and. (CropType == 0) ) &
       LeafAreaIndEff = 0.0                                                                    ! MB: LeafAreaIndex CHECK

    ! set growing season flag
    if ( ((TemperatureCanopy > TemperatureMinPhotosyn) .and. (CropType == 0)) .or. &
         ((PlantGrowStage > 2) .and. (PlantGrowStage < 7) .and. (CropType > 0))) then
       IndexGrowSeason = 1.0
    else
       IndexGrowSeason = 0.0
    endif 

    ! compute vegetation fraction
    ! input green vegetation fraction should be consistent with LeafAreaIndex
    if ( (OptDynamicVeg == 1) .or. (OptDynamicVeg == 6) .or. (OptDynamicVeg == 7) ) then      ! use VegFrac = VegFracGreen from input
       VegFrac = VegFracGreen
    elseif ( (OptDynamicVeg == 2) .or. (OptDynamicVeg == 3) .or. (OptDynamicVeg == 8) ) then  ! computed VegFrac from LeafAreaIndex & StemAreaIndex
       VegFrac = 1.0 - exp(-0.52 * (LeafAreaIndex + StemAreaIndex))
    elseif ( (OptDynamicVeg == 4) .or. (OptDynamicVeg == 5) .or. (OptDynamicVeg == 9) ) then  ! use yearly maximum vegetation fraction
       VegFrac = VegFracAnnMax
    else                                                                                      ! outside existing vegetation options
#ifdef _OPENACC
       ! GPU: Cannot write/stop in parallel region
#else
       write(*,*) "Un-recognized dynamic vegetation option (OptDynamicVeg)... "
       stop "Error: Namelist parameter OptDynamicVeg unknown"
#endif
    endif
    ! use maximum vegetation fraction for crop run
    if ( (OptCropModel > 0) .and. (CropType > 0) ) then
       VegFrac = VegFracAnnMax
    endif

    ! adjust unreasonable vegetation fraction
    if ( VegFrac <= 0.05 ) VegFrac = 0.05
    if ( (FlagUrban .eqv. .true.) .or. (VegType == IndexBarrenPoint) ) VegFrac = 0.0
    if ( (LeafAreaIndEff+StemAreaIndEff) == 0.0 ) VegFrac = 0.0

    ! determine if activate dynamic vegetation or crop run
    FlagDynamicCrop = .false.
    FlagDynamicVeg  = .false.
    if ( (OptDynamicVeg == 2) .or. (OptDynamicVeg == 5) .or. (OptDynamicVeg == 6) ) &
       FlagDynamicVeg = .true.
    if ( (OptCropModel > 0) .and. (CropType > 0) ) then
       FlagDynamicCrop = .true.
       FlagDynamicVeg  = .false.
    endif

    end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine PhenologyMain

end module PhenologyMainMod
