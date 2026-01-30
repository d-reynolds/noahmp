module SnowLayerCombineMod

!!! Snowpack layer combination process (2D GPU-optimized)
!!! Update snow ice, snow water, snow thickness, snow temperature

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowLayerWaterComboMod, only: SnowLayerWaterCombo

  implicit none

contains

  subroutine SnowLayerCombine(noahmp, II, JJ)
!$acc routine seq
! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: COMBINE
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp
    integer, intent(in)          :: II, JJ            ! grid indices
! local variable
    integer                          :: I, J, K, L        ! node indices
    integer                          :: NumSnowLayerOld   ! number of snow layer
    integer                          :: IndLayer          ! node index
    integer                          :: IndNeighbor       ! adjacent node selected for combination
    real(kind=kind_noahmp)           :: SnowIceTmp        ! total ice mass in snow
    real(kind=kind_noahmp)           :: SnowLiqTmp        ! total liquid water in snow
    real(kind=kind_noahmp)           :: SnowThickMin(3)   ! minimum thickness of each snow layer
    data SnowThickMin /0.025, 0.025, 0.1/                 ! MB: change limit
    !data SnowThickMin /0.045, 0.05, 0.2/

      if ( noahmp%config%domain%NumSnowLayerNeg(II,JJ) >= 0 ) return  ! no snow layers

! --------------------------------------------------------------------
    associate(                                                                       &
              OptSnowAlbedo          => noahmp%config%nmlist%OptSnowAlbedo          ,& ! in,    options for ground snow surface albedo
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg(II,JJ)        ,& ! inout, actual number of snow layers (negative)
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m]
              TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow     ,& ! inout, snow and soil layer temperature [K]
              SnowDepth              => noahmp%water%state%SnowDepth(II,JJ)                ,& ! inout, snow depth [m]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv(II,JJ)           ,& ! inout, snow water equivalent [mm]
              SnowIce                => noahmp%water%state%SnowIce                  ,& ! inout, snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater             ,& ! inout, snow layer liquid water [mm]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! inout, soil liquid moisture [m3/m3]
              SoilIce                => noahmp%water%state%SoilIce                  ,& ! inout, soil ice moisture [m3/m3]
              MassBChydropho         => noahmp%water%state%MassBChydropho           ,& ! inout, mass of hydrophobic Black Carbon in snow [kg m-2]
              MassBChydrophi         => noahmp%water%state%MassBChydrophi           ,& ! inout, mass of hydrophillic Black Carbon in snow [kg m-2]
              MassOChydropho         => noahmp%water%state%MassOChydropho           ,& ! inout, mass of hydrophobic Organic Carbon in snow [kg m-2]
              MassOChydrophi         => noahmp%water%state%MassOChydrophi           ,& ! inout, mass of hydrophillic Organic Carbon in snow [kg m-2]
              MassDust1              => noahmp%water%state%MassDust1                ,& ! inout, mass of dust species 1 in snow [kg m-2]
              MassDust2              => noahmp%water%state%MassDust2                ,& ! inout, mass of dust species 2 in snow [kg m-2]
              MassDust3              => noahmp%water%state%MassDust3                ,& ! inout, mass of dust species 3 in snow [kg m-2]
              MassDust4              => noahmp%water%state%MassDust4                ,& ! inout, mass of dust species 4 in snow [kg m-2]
              MassDust5              => noahmp%water%state%MassDust5                ,& ! inout, mass of dust species 5 in snow [kg m-2]
              SnowRadius             => noahmp%water%state%SnowRadius               ,& ! inout, effective grain radius [microns, m-6]
              PondSfcThinSnwComb     => noahmp%water%state%PondSfcThinSnwComb(II,JJ)       ,& ! out,   surface ponding [mm] from liquid in thin snow layer combination
              PondSfcThinSnwTrans    => noahmp%water%state%PondSfcThinSnwTrans(II,JJ)       & ! out,   surface ponding [mm] from thin snow when changing from multilayer to no layer
             )
! ----------------------------------------------------------------------

    ! check and combine small ice content layer
    NumSnowLayerOld = NumSnowLayerNeg

   !$acc loop seq
    do J = NumSnowLayerOld+1,0
       if ( SnowIce(II,J,JJ) <= 0.1 ) then
          if ( J /= 0 ) then
             SnowLiqWater(II,J+1,JJ)           = SnowLiqWater(II,J+1,JJ) + SnowLiqWater(II,J,JJ)
             SnowIce(II,J+1,JJ)                = SnowIce(II,J+1,JJ) + SnowIce(II,J,JJ)
             ThicknessSnowSoilLayer(II,J+1,JJ) = ThicknessSnowSoilLayer(II,J+1,JJ) + ThicknessSnowSoilLayer(II,J,JJ)

             if ( OptSnowAlbedo == 3 ) then
                MassBChydropho(II,J+1,JJ)      = MassBChydropho(II,J+1,JJ) +  MassBChydropho(II,J,JJ)
                MassBChydrophi(II,J+1,JJ)      = MassBChydrophi(II,J+1,JJ) +  MassBChydrophi(II,J,JJ)
                MassOChydropho(II,J+1,JJ)      = MassOChydropho(II,J+1,JJ) +  MassOChydropho(II,J,JJ)
                MassOChydrophi(II,J+1,JJ)      = MassOChydrophi(II,J+1,JJ) +  MassOChydrophi(II,J,JJ)
                MassDust1(II,J+1,JJ)           = MassDust1(II,J+1,JJ) +  MassDust1(II,J,JJ)
                MassDust2(II,J+1,JJ)           = MassDust2(II,J+1,JJ) +  MassDust2(II,J,JJ)
                MassDust3(II,J+1,JJ)           = MassDust3(II,J+1,JJ) +  MassDust3(II,J,JJ)
                MassDust4(II,J+1,JJ)           = MassDust4(II,J+1,JJ) +  MassDust4(II,J,JJ)
                MassDust5(II,J+1,JJ)           = MassDust5(II,J+1,JJ) +  MassDust5(II,J,JJ)
             endif

          else
             if ( NumSnowLayerNeg < -1 ) then    ! MB/KM: change to NumSnowLayerNeg
                SnowLiqWater(II,J-1,JJ)           = SnowLiqWater(II,J-1,JJ) + SnowLiqWater(II,J,JJ)
                SnowIce(II,J-1,JJ)                = SnowIce(II,J-1,JJ) + SnowIce(II,J,JJ)
                ThicknessSnowSoilLayer(II,J-1,JJ) = ThicknessSnowSoilLayer(II,J-1,JJ) + ThicknessSnowSoilLayer(II,J,JJ)

                if ( OptSnowAlbedo == 3 ) then
                   MassBChydropho(II,J-1,JJ)      = MassBChydropho(II,J-1,JJ) +  MassBChydropho(II,J,JJ)
                   MassBChydrophi(II,J-1,JJ)      = MassBChydrophi(II,J-1,JJ) +  MassBChydrophi(II,J,JJ)
                   MassOChydropho(II,J-1,JJ)      = MassOChydropho(II,J-1,JJ) +  MassOChydropho(II,J,JJ)
                   MassOChydrophi(II,J-1,JJ)      = MassOChydrophi(II,J-1,JJ) +  MassOChydrophi(II,J,JJ)
                   MassDust1(II,J-1,JJ)           = MassDust1(II,J-1,JJ) +  MassDust1(II,J,JJ)
                   MassDust2(II,J-1,JJ)           = MassDust2(II,J-1,JJ) +  MassDust2(II,J,JJ)
                   MassDust3(II,J-1,JJ)           = MassDust3(II,J-1,JJ) +  MassDust3(II,J,JJ)
                   MassDust4(II,J-1,JJ)           = MassDust4(II,J-1,JJ) +  MassDust4(II,J,JJ)
                   MassDust5(II,J-1,JJ)           = MassDust5(II,J-1,JJ) +  MassDust5(II,J,JJ)
                endif

             else
                if ( SnowIce(II,J,JJ) >= 0.0 ) then
                   PondSfcThinSnwComb = SnowLiqWater(II,J,JJ)                ! NumSnowLayerNeg WILL GET SET TO ZERO BELOW; PondSfcThinSnwComb WILL GET
                   SnowWaterEquiv     = SnowIce(II,J,JJ)                     ! ADDED TO PONDING FROM PHASECHANGE PONDING SHOULD BE
                   SnowDepth          = ThicknessSnowSoilLayer(II,J,JJ)      ! ZERO HERE BECAUSE IT WAS CALCULATED FOR THIN SNOW
                else  ! SnowIce OVER-SUBLIMATED EARLIER
                   PondSfcThinSnwComb = SnowLiqWater(II,J,JJ) + SnowIce(II,J,JJ)
                   if ( PondSfcThinSnwComb < 0.0 ) then                ! IF SnowIce AND SnowLiqWater SUBLIMATES REMOVE FROM SOIL
                      SoilIce(II,1,JJ) = SoilIce(II,1,JJ) + PondSfcThinSnwComb/(ThicknessSnowSoilLayer(II,1,JJ)*1000.0) ! negative SoilIce from oversublimation is adjusted below
                      PondSfcThinSnwComb = 0.0
                   endif
                   SnowWaterEquiv = 0.0
                   SnowDepth      = 0.0
                endif ! if(SnowIce(II,J,JJ) >= 0.0)
                SnowLiqWater(II,J,JJ)   = 0.0
                SnowIce(II,J,JJ)        = 0.0
                ThicknessSnowSoilLayer(II,J,JJ) = 0.0

                ! SNICAR, aerosol flux may infiltrate into top soil like PondSfcThinSnwComb, it
                ! would be more thorough to do so later
                if ( OptSnowAlbedo == 3 ) then
                   MassBChydropho(II,J,JJ) = 0.0
                   MassBChydrophi(II,J,JJ) = 0.0
                   MassOChydropho(II,J,JJ) = 0.0
                   MassOChydrophi(II,J,JJ) = 0.0
                   MassDust1(II,J,JJ)      = 0.0
                   MassDust2(II,J,JJ)      = 0.0
                   MassDust3(II,J,JJ)      = 0.0
                   MassDust4(II,J,JJ)      = 0.0
                   MassDust5(II,J,JJ)      = 0.0
                endif

             endif ! if(NumSnowLayerNeg < -1)
          endif ! if(J /= 0)

          ! shift all elements above this down by one.
          if ( (J > NumSnowLayerNeg+1) .and. (NumSnowLayerNeg < -1) ) then
            !$acc loop seq 
            do I = J, NumSnowLayerNeg+2, -1
                TemperatureSoilSnow(II,I,JJ)    = TemperatureSoilSnow(I-II,1,JJ)
                SnowLiqWater(II,I,JJ)           = SnowLiqWater(I-II,1,JJ)
                SnowIce(II,I,JJ)                = SnowIce(I-II,1,JJ)
                ThicknessSnowSoilLayer(II,I,JJ) = ThicknessSnowSoilLayer(I-II,1,JJ)

                if ( OptSnowAlbedo == 3 ) then
                   MassBChydropho(II,I,JJ)      = MassBChydropho(I-II,1,JJ)
                   MassBChydrophi(II,I,JJ)      = MassBChydrophi(I-II,1,JJ)
                   MassOChydropho(II,I,JJ)      = MassOChydropho(I-II,1,JJ)
                   MassOChydrophi(II,I,JJ)      = MassOChydrophi(I-II,1,JJ)
                   MassDust1(II,I,JJ)           = MassDust1(I-II,1,JJ)
                   MassDust2(II,I,JJ)           = MassDust2(I-II,1,JJ)
                   MassDust3(II,I,JJ)           = MassDust3(I-II,1,JJ)
                   MassDust4(II,I,JJ)           = MassDust4(I-II,1,JJ)
                   MassDust5(II,I,JJ)           = MassDust5(I-II,1,JJ)
                   SnowRadius(II,I,JJ)          = SnowRadius(I-II,1,JJ)
                endif

             enddo
          endif
          NumSnowLayerNeg = NumSnowLayerNeg + 1

       endif ! if(SnowIce(J) <= 0.1)
    enddo ! do J

    ! to conserve water in case of too large surface sublimation
    if ( SoilIce(II,1,JJ) < 0.0) then
       SoilLiqWater(II,1,JJ) = SoilLiqWater(II,1,JJ) + SoilIce(II,1,JJ)
       SoilIce(II,1,JJ)      = 0.0
    endif

    if ( NumSnowLayerNeg ==0 ) return   ! MB: get out if no longer multi-layer

    SnowWaterEquiv = 0.0
    SnowDepth      = 0.0
    SnowIceTmp     = 0.0
    SnowLiqTmp     = 0.0

    !$acc loop seq
    do J = NumSnowLayerNeg+1, 0
       SnowWaterEquiv = SnowWaterEquiv + SnowIce(II,J,JJ) + SnowLiqWater(II,J,JJ)
       SnowDepth      = SnowDepth + ThicknessSnowSoilLayer(II,J,JJ)
       SnowIceTmp     = SnowIceTmp + SnowIce(II,J,JJ)
       SnowLiqTmp     = SnowLiqTmp + SnowLiqWater(II,J,JJ)
    enddo

    ! check the snow depth - all snow gone, the liquid water assumes ponding on soil surface.
    ! if ( (SnowDepth < 0.05) .and. (NumSnowLayerNeg < 0) ) then
    if ( (SnowDepth < 0.025) .and. (NumSnowLayerNeg < 0) ) then ! MB: change limit
       NumSnowLayerNeg     = 0
       SnowWaterEquiv      = SnowIceTmp
       PondSfcThinSnwTrans = SnowLiqTmp                ! LIMIT OF NumSnowLayerNeg < 0 MEANS INPUT PONDING
       if ( SnowWaterEquiv <= 0.0 ) SnowDepth = 0.0    ! SHOULD BE ZERO; SEE ABOVE
    endif

    ! check the snow depth - snow layers combined
    if ( NumSnowLayerNeg < -1 ) then
       NumSnowLayerOld = NumSnowLayerNeg
       IndLayer        = 1
       !$acc loop seq
       do I = NumSnowLayerOld+1, 0
          if ( ThicknessSnowSoilLayer(II,I,JJ) < SnowThickMin(IndLayer) ) then
             if ( I == NumSnowLayerNeg+1 ) then
                IndNeighbor = I + 1
             else if ( I == 0 ) then
                IndNeighbor = I - 1
             else
                IndNeighbor = I + 1
                if ( (ThicknessSnowSoilLayer(I-II,1,JJ)+ThicknessSnowSoilLayer(II,I,JJ)) < &
                     (ThicknessSnowSoilLayer(I+II,1,JJ)+ThicknessSnowSoilLayer(II,I,JJ)) ) IndNeighbor = I-1
             endif
             ! Node l and j are combined and stored as node j.
             if ( IndNeighbor > I ) then
                J = IndNeighbor
                L = I
             else
                J = I
                L = IndNeighbor
             endif

             if ( OptSnowAlbedo == 3 ) then
                MassBChydropho(II,J,JJ) = MassBChydropho(II,J,JJ) +  MassBChydropho(II,L,JJ)
                MassBChydrophi(II,J,JJ) = MassBChydrophi(II,J,JJ) +  MassBChydrophi(II,L,JJ)
                MassOChydropho(II,J,JJ) = MassOChydropho(II,J,JJ) +  MassOChydropho(II,L,JJ)
                MassOChydrophi(II,J,JJ) = MassOChydrophi(II,J,JJ) +  MassOChydrophi(II,L,JJ)
                MassDust1(II,J,JJ)      = MassDust1(II,J,JJ) +  MassDust1(II,L,JJ)
                MassDust2(II,J,JJ)      = MassDust2(II,J,JJ) +  MassDust2(II,L,JJ)
                MassDust3(II,J,JJ)      = MassDust3(II,J,JJ) +  MassDust3(II,L,JJ)
                MassDust4(II,J,JJ)      = MassDust4(II,J,JJ) +  MassDust4(II,L,JJ)
                MassDust5(II,J,JJ)      = MassDust5(II,J,JJ) +  MassDust5(II,L,JJ)
                SnowRadius(II,J,JJ)     = (SnowRadius(II,J,JJ)*(SnowLiqWater(II,J,JJ)+SnowIce(II,J,JJ)) + SnowRadius(II,L,JJ)*(SnowLiqWater(II,L,JJ)+SnowIce(II,L,JJ))) / &
                                    (SnowLiqWater(II,J,JJ) + SnowIce(II,J,JJ) + SnowLiqWater(II,L,JJ) + SnowIce(II,L,JJ))
             endif

             ! update combined snow water & temperature
             call SnowLayerWaterCombo(ThicknessSnowSoilLayer(II,J,JJ), SnowLiqWater(II,J,JJ), SnowIce(II,J,JJ), TemperatureSoilSnow(II,J,JJ), &
                                      ThicknessSnowSoilLayer(II,L,JJ), SnowLiqWater(II,L,JJ), SnowIce(II,L,JJ), TemperatureSoilSnow(II,L,JJ) )


             ! Now shift all elements above this down one.
             if ( (J-1) > (NumSnowLayerNeg+1) ) then
                !$acc loop seq
                do K = J-1, NumSnowLayerNeg+2, -1
                   TemperatureSoilSnow(II,K,JJ)    = TemperatureSoilSnow(II,K-1,JJ)
                   SnowIce(II,K,JJ)                = SnowIce(II,K-1,JJ)
                   SnowLiqWater(II,K,JJ)           = SnowLiqWater(II,K-1,JJ)
                   ThicknessSnowSoilLayer(II,K,JJ) = ThicknessSnowSoilLayer(II,K-1,JJ)

                   if ( OptSnowAlbedo == 3 ) then
                      MassBChydropho(II,K,JJ)      = MassBChydropho(II,K-1,JJ)
                      MassBChydrophi(II,K,JJ)      = MassBChydrophi(II,K-1,JJ)
                      MassOChydropho(II,K,JJ)      = MassOChydropho(II,K-1,JJ)
                      MassOChydrophi(II,K,JJ)      = MassOChydrophi(II,K-1,JJ)
                      MassDust1(II,K,JJ)           = MassDust1(II,K-1,JJ)
                      MassDust2(II,K,JJ)           = MassDust2(II,K-1,JJ)
                      MassDust3(II,K,JJ)           = MassDust3(II,K-1,JJ)
                      MassDust4(II,K,JJ)           = MassDust4(II,K-1,JJ)
                      MassDust5(II,K,JJ)           = MassDust5(II,K-1,JJ)
                      SnowRadius(II,K,JJ)          = SnowRadius(II,K-1,JJ)
                   endif

                enddo
             endif
             ! Decrease the number of snow layers
             NumSnowLayerNeg = NumSnowLayerNeg + 1
             if ( NumSnowLayerNeg >= -1 ) Exit
          else 
             ! The layer thickness is greater than the prescribed minimum value
             IndLayer = IndLayer + 1
          endif
       enddo
    endif

    end associate

  end subroutine SnowLayerCombine

end module SnowLayerCombineMod
