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

    ! check and combine small ice content layer
    NumSnowLayerOld = noahmp%config%domain%NumSnowLayerNeg(II,JJ)

   !$acc loop seq
    do J = NumSnowLayerOld+1,0
       if ( noahmp%water%state%SnowIce(II,J,JJ) <= 0.1 ) then
          if ( J /= 0 ) then
             noahmp%water%state%SnowLiqWater(II,J+1,JJ)           = noahmp%water%state%SnowLiqWater(II,J+1,JJ) + noahmp%water%state%SnowLiqWater(II,J,JJ)
             noahmp%water%state%SnowIce(II,J+1,JJ)                = noahmp%water%state%SnowIce(II,J+1,JJ) + noahmp%water%state%SnowIce(II,J,JJ)
             noahmp%config%domain%ThicknessSnowSoilLayer(II,J+1,JJ) = noahmp%config%domain%ThicknessSnowSoilLayer(II,J+1,JJ) + noahmp%config%domain%ThicknessSnowSoilLayer(II,J,JJ)

             if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
                noahmp%water%state%MassBChydropho(II,J+1,JJ)      = noahmp%water%state%MassBChydropho(II,J+1,JJ) +  noahmp%water%state%MassBChydropho(II,J,JJ)
                noahmp%water%state%MassBChydrophi(II,J+1,JJ)      = noahmp%water%state%MassBChydrophi(II,J+1,JJ) +  noahmp%water%state%MassBChydrophi(II,J,JJ)
                noahmp%water%state%MassOChydropho(II,J+1,JJ)      = noahmp%water%state%MassOChydropho(II,J+1,JJ) +  noahmp%water%state%MassOChydropho(II,J,JJ)
                noahmp%water%state%MassOChydrophi(II,J+1,JJ)      = noahmp%water%state%MassOChydrophi(II,J+1,JJ) +  noahmp%water%state%MassOChydrophi(II,J,JJ)
                noahmp%water%state%MassDust1(II,J+1,JJ)           = noahmp%water%state%MassDust1(II,J+1,JJ) +  noahmp%water%state%MassDust1(II,J,JJ)
                noahmp%water%state%MassDust2(II,J+1,JJ)           = noahmp%water%state%MassDust2(II,J+1,JJ) +  noahmp%water%state%MassDust2(II,J,JJ)
                noahmp%water%state%MassDust3(II,J+1,JJ)           = noahmp%water%state%MassDust3(II,J+1,JJ) +  noahmp%water%state%MassDust3(II,J,JJ)
                noahmp%water%state%MassDust4(II,J+1,JJ)           = noahmp%water%state%MassDust4(II,J+1,JJ) +  noahmp%water%state%MassDust4(II,J,JJ)
                noahmp%water%state%MassDust5(II,J+1,JJ)           = noahmp%water%state%MassDust5(II,J+1,JJ) +  noahmp%water%state%MassDust5(II,J,JJ)
             endif

          else
             if ( noahmp%config%domain%NumSnowLayerNeg(II,JJ) < -1 ) then    ! MB/KM: change to NumSnowLayerNeg
                noahmp%water%state%SnowLiqWater(II,J-1,JJ)           = noahmp%water%state%SnowLiqWater(II,J-1,JJ) + noahmp%water%state%SnowLiqWater(II,J,JJ)
                noahmp%water%state%SnowIce(II,J-1,JJ)                = noahmp%water%state%SnowIce(II,J-1,JJ) + noahmp%water%state%SnowIce(II,J,JJ)
                noahmp%config%domain%ThicknessSnowSoilLayer(II,J-1,JJ) = noahmp%config%domain%ThicknessSnowSoilLayer(II,J-1,JJ) + noahmp%config%domain%ThicknessSnowSoilLayer(II,J,JJ)

                if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
                   noahmp%water%state%MassBChydropho(II,J-1,JJ)      = noahmp%water%state%MassBChydropho(II,J-1,JJ) +  noahmp%water%state%MassBChydropho(II,J,JJ)
                   noahmp%water%state%MassBChydrophi(II,J-1,JJ)      = noahmp%water%state%MassBChydrophi(II,J-1,JJ) +  noahmp%water%state%MassBChydrophi(II,J,JJ)
                   noahmp%water%state%MassOChydropho(II,J-1,JJ)      = noahmp%water%state%MassOChydropho(II,J-1,JJ) +  noahmp%water%state%MassOChydropho(II,J,JJ)
                   noahmp%water%state%MassOChydrophi(II,J-1,JJ)      = noahmp%water%state%MassOChydrophi(II,J-1,JJ) +  noahmp%water%state%MassOChydrophi(II,J,JJ)
                   noahmp%water%state%MassDust1(II,J-1,JJ)           = noahmp%water%state%MassDust1(II,J-1,JJ) +  noahmp%water%state%MassDust1(II,J,JJ)
                   noahmp%water%state%MassDust2(II,J-1,JJ)           = noahmp%water%state%MassDust2(II,J-1,JJ) +  noahmp%water%state%MassDust2(II,J,JJ)
                   noahmp%water%state%MassDust3(II,J-1,JJ)           = noahmp%water%state%MassDust3(II,J-1,JJ) +  noahmp%water%state%MassDust3(II,J,JJ)
                   noahmp%water%state%MassDust4(II,J-1,JJ)           = noahmp%water%state%MassDust4(II,J-1,JJ) +  noahmp%water%state%MassDust4(II,J,JJ)
                   noahmp%water%state%MassDust5(II,J-1,JJ)           = noahmp%water%state%MassDust5(II,J-1,JJ) +  noahmp%water%state%MassDust5(II,J,JJ)
                endif

             else
                if ( noahmp%water%state%SnowIce(II,J,JJ) >= 0.0 ) then
                   noahmp%water%state%PondSfcThinSnwComb(II,JJ) = noahmp%water%state%SnowLiqWater(II,J,JJ)                ! NumSnowLayerNeg WILL GET SET TO ZERO BELOW; PondSfcThinSnwComb WILL GET
                   noahmp%water%state%SnowWaterEquiv(II,JJ)     = noahmp%water%state%SnowIce(II,J,JJ)                     ! ADDED TO PONDING FROM PHASECHANGE PONDING SHOULD BE
                   noahmp%water%state%SnowDepth(II,JJ)          = noahmp%config%domain%ThicknessSnowSoilLayer(II,J,JJ)      ! ZERO HERE BECAUSE IT WAS CALCULATED FOR THIN SNOW
                else  ! SnowIce OVER-SUBLIMATED EARLIER
                   noahmp%water%state%PondSfcThinSnwComb(II,JJ) = noahmp%water%state%SnowLiqWater(II,J,JJ) + noahmp%water%state%SnowIce(II,J,JJ)
                   if ( noahmp%water%state%PondSfcThinSnwComb(II,JJ) < 0.0 ) then                ! IF SnowIce AND SnowLiqWater SUBLIMATES REMOVE FROM SOIL
                      noahmp%water%state%SoilIce(II,1,JJ) = noahmp%water%state%SoilIce(II,1,JJ) + noahmp%water%state%PondSfcThinSnwComb(II,JJ)/(noahmp%config%domain%ThicknessSnowSoilLayer(II,1,JJ)*1000.0) ! negative SoilIce from oversublimation is adjusted below
                      noahmp%water%state%PondSfcThinSnwComb(II,JJ) = 0.0
                   endif
                   noahmp%water%state%SnowWaterEquiv(II,JJ) = 0.0
                   noahmp%water%state%SnowDepth(II,JJ)      = 0.0
                endif ! if(SnowIce(II,J,JJ) >= 0.0)
                noahmp%water%state%SnowLiqWater(II,J,JJ)   = 0.0
                noahmp%water%state%SnowIce(II,J,JJ)        = 0.0
                noahmp%config%domain%ThicknessSnowSoilLayer(II,J,JJ) = 0.0

                ! SNICAR, aerosol flux may infiltrate into top soil like PondSfcThinSnwComb, it
                ! would be more thorough to do so later
                if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
                   noahmp%water%state%MassBChydropho(II,J,JJ) = 0.0
                   noahmp%water%state%MassBChydrophi(II,J,JJ) = 0.0
                   noahmp%water%state%MassOChydropho(II,J,JJ) = 0.0
                   noahmp%water%state%MassOChydrophi(II,J,JJ) = 0.0
                   noahmp%water%state%MassDust1(II,J,JJ)      = 0.0
                   noahmp%water%state%MassDust2(II,J,JJ)      = 0.0
                   noahmp%water%state%MassDust3(II,J,JJ)      = 0.0
                   noahmp%water%state%MassDust4(II,J,JJ)      = 0.0
                   noahmp%water%state%MassDust5(II,J,JJ)      = 0.0
                endif

             endif ! if(NumSnowLayerNeg < -1)
          endif ! if(J /= 0)

          ! shift all elements above this down by one.
          if ( (J > noahmp%config%domain%NumSnowLayerNeg(II,JJ)+1) .and. (noahmp%config%domain%NumSnowLayerNeg(II,JJ) < -1) ) then
            !$acc loop seq
            do I = J, noahmp%config%domain%NumSnowLayerNeg(II,JJ)+2, -1
                noahmp%energy%state%TemperatureSoilSnow(II,I,JJ)    = noahmp%energy%state%TemperatureSoilSnow(II,I-1,JJ)
                noahmp%water%state%SnowLiqWater(II,I,JJ)           = noahmp%water%state%SnowLiqWater(II,I-1,JJ)
                noahmp%water%state%SnowIce(II,I,JJ)                = noahmp%water%state%SnowIce(II,I-1,JJ)
                noahmp%config%domain%ThicknessSnowSoilLayer(II,I,JJ) = noahmp%config%domain%ThicknessSnowSoilLayer(II,I-1,JJ)

                if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
                   noahmp%water%state%MassBChydropho(II,I,JJ)      = noahmp%water%state%MassBChydropho(II,I-1,JJ)
                   noahmp%water%state%MassBChydrophi(II,I,JJ)      = noahmp%water%state%MassBChydrophi(II,I-1,JJ)
                   noahmp%water%state%MassOChydropho(II,I,JJ)      = noahmp%water%state%MassOChydropho(II,I-1,JJ)
                   noahmp%water%state%MassOChydrophi(II,I,JJ)      = noahmp%water%state%MassOChydrophi(II,I-1,JJ)
                   noahmp%water%state%MassDust1(II,I,JJ)           = noahmp%water%state%MassDust1(II,I-1,JJ)
                   noahmp%water%state%MassDust2(II,I,JJ)           = noahmp%water%state%MassDust2(II,I-1,JJ)
                   noahmp%water%state%MassDust3(II,I,JJ)           = noahmp%water%state%MassDust3(II,I-1,JJ)
                   noahmp%water%state%MassDust4(II,I,JJ)           = noahmp%water%state%MassDust4(II,I-1,JJ)
                   noahmp%water%state%MassDust5(II,I,JJ)           = noahmp%water%state%MassDust5(II,I-1,JJ)
                   noahmp%water%state%SnowRadius(II,I,JJ)          = noahmp%water%state%SnowRadius(II,I-1,JJ)
                endif

             enddo
          endif
          noahmp%config%domain%NumSnowLayerNeg(II,JJ) = noahmp%config%domain%NumSnowLayerNeg(II,JJ) + 1

       endif ! if(SnowIce(J) <= 0.1)
    enddo ! do J

    ! to conserve water in case of too large surface sublimation
    if ( noahmp%water%state%SoilIce(II,1,JJ) < 0.0) then
       noahmp%water%state%SoilLiqWater(II,1,JJ) = noahmp%water%state%SoilLiqWater(II,1,JJ) + noahmp%water%state%SoilIce(II,1,JJ)
       noahmp%water%state%SoilIce(II,1,JJ)      = 0.0
    endif

    if ( noahmp%config%domain%NumSnowLayerNeg(II,JJ) ==0 ) return   ! MB: get out if no longer multi-layer

    noahmp%water%state%SnowWaterEquiv(II,JJ) = 0.0
    noahmp%water%state%SnowDepth(II,JJ)      = 0.0
    SnowIceTmp     = 0.0
    SnowLiqTmp     = 0.0

    !$acc loop seq
    do J = noahmp%config%domain%NumSnowLayerNeg(II,JJ)+1, 0
       noahmp%water%state%SnowWaterEquiv(II,JJ) = noahmp%water%state%SnowWaterEquiv(II,JJ) + noahmp%water%state%SnowIce(II,J,JJ) + noahmp%water%state%SnowLiqWater(II,J,JJ)
       noahmp%water%state%SnowDepth(II,JJ)      = noahmp%water%state%SnowDepth(II,JJ) + noahmp%config%domain%ThicknessSnowSoilLayer(II,J,JJ)
       SnowIceTmp     = SnowIceTmp + noahmp%water%state%SnowIce(II,J,JJ)
       SnowLiqTmp     = SnowLiqTmp + noahmp%water%state%SnowLiqWater(II,J,JJ)
    enddo

    ! check the snow depth - all snow gone, the liquid water assumes ponding on soil surface.
    ! if ( (SnowDepth < 0.05) .and. (NumSnowLayerNeg < 0) ) then
    if ( (noahmp%water%state%SnowDepth(II,JJ) < 0.025) .and. (noahmp%config%domain%NumSnowLayerNeg(II,JJ) < 0) ) then ! MB: change limit
       noahmp%config%domain%NumSnowLayerNeg(II,JJ)     = 0
       noahmp%water%state%SnowWaterEquiv(II,JJ)      = SnowIceTmp
       noahmp%water%state%PondSfcThinSnwTrans(II,JJ) = SnowLiqTmp                ! LIMIT OF NumSnowLayerNeg < 0 MEANS INPUT PONDING
       if ( noahmp%water%state%SnowWaterEquiv(II,JJ) <= 0.0 ) noahmp%water%state%SnowDepth(II,JJ) = 0.0    ! SHOULD BE ZERO; SEE ABOVE
    endif

    ! check the snow depth - snow layers combined
    if ( noahmp%config%domain%NumSnowLayerNeg(II,JJ) < -1 ) then
       NumSnowLayerOld = noahmp%config%domain%NumSnowLayerNeg(II,JJ)
       IndLayer        = 1
       !$acc loop seq
       do I = NumSnowLayerOld+1, 0
          if ( noahmp%config%domain%ThicknessSnowSoilLayer(II,I,JJ) < SnowThickMin(IndLayer) ) then
             if ( I == noahmp%config%domain%NumSnowLayerNeg(II,JJ)+1 ) then
                IndNeighbor = I + 1
             else if ( I == 0 ) then
                IndNeighbor = I - 1
             else
                IndNeighbor = I + 1
                if ( (noahmp%config%domain%ThicknessSnowSoilLayer(II,I-1,JJ)+noahmp%config%domain%ThicknessSnowSoilLayer(II,I,JJ)) < &
                     (noahmp%config%domain%ThicknessSnowSoilLayer(II,I+1,JJ)+noahmp%config%domain%ThicknessSnowSoilLayer(II,I,JJ)) ) IndNeighbor = I-1
             endif
             ! Node l and j are combined and stored as node j.
             if ( IndNeighbor > I ) then
                J = IndNeighbor
                L = I
             else
                J = I
                L = IndNeighbor
             endif

             if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
                noahmp%water%state%MassBChydropho(II,J,JJ) = noahmp%water%state%MassBChydropho(II,J,JJ) +  noahmp%water%state%MassBChydropho(II,L,JJ)
                noahmp%water%state%MassBChydrophi(II,J,JJ) = noahmp%water%state%MassBChydrophi(II,J,JJ) +  noahmp%water%state%MassBChydrophi(II,L,JJ)
                noahmp%water%state%MassOChydropho(II,J,JJ) = noahmp%water%state%MassOChydropho(II,J,JJ) +  noahmp%water%state%MassOChydropho(II,L,JJ)
                noahmp%water%state%MassOChydrophi(II,J,JJ) = noahmp%water%state%MassOChydrophi(II,J,JJ) +  noahmp%water%state%MassOChydrophi(II,L,JJ)
                noahmp%water%state%MassDust1(II,J,JJ)      = noahmp%water%state%MassDust1(II,J,JJ) +  noahmp%water%state%MassDust1(II,L,JJ)
                noahmp%water%state%MassDust2(II,J,JJ)      = noahmp%water%state%MassDust2(II,J,JJ) +  noahmp%water%state%MassDust2(II,L,JJ)
                noahmp%water%state%MassDust3(II,J,JJ)      = noahmp%water%state%MassDust3(II,J,JJ) +  noahmp%water%state%MassDust3(II,L,JJ)
                noahmp%water%state%MassDust4(II,J,JJ)      = noahmp%water%state%MassDust4(II,J,JJ) +  noahmp%water%state%MassDust4(II,L,JJ)
                noahmp%water%state%MassDust5(II,J,JJ)      = noahmp%water%state%MassDust5(II,J,JJ) +  noahmp%water%state%MassDust5(II,L,JJ)
                noahmp%water%state%SnowRadius(II,J,JJ)     = (noahmp%water%state%SnowRadius(II,J,JJ)*(noahmp%water%state%SnowLiqWater(II,J,JJ)+noahmp%water%state%SnowIce(II,J,JJ)) + noahmp%water%state%SnowRadius(II,L,JJ)*(noahmp%water%state%SnowLiqWater(II,L,JJ)+noahmp%water%state%SnowIce(II,L,JJ))) / &
                                    (noahmp%water%state%SnowLiqWater(II,J,JJ) + noahmp%water%state%SnowIce(II,J,JJ) + noahmp%water%state%SnowLiqWater(II,L,JJ) + noahmp%water%state%SnowIce(II,L,JJ))
             endif

             ! update combined snow water & temperature
             call SnowLayerWaterCombo(noahmp%config%domain%ThicknessSnowSoilLayer(II,J,JJ), noahmp%water%state%SnowLiqWater(II,J,JJ), noahmp%water%state%SnowIce(II,J,JJ), noahmp%energy%state%TemperatureSoilSnow(II,J,JJ), &
                                      noahmp%config%domain%ThicknessSnowSoilLayer(II,L,JJ), noahmp%water%state%SnowLiqWater(II,L,JJ), noahmp%water%state%SnowIce(II,L,JJ), noahmp%energy%state%TemperatureSoilSnow(II,L,JJ) )


             ! Now shift all elements above this down one.
             if ( (J-1) > (noahmp%config%domain%NumSnowLayerNeg(II,JJ)+1) ) then
                !$acc loop seq
                do K = J-1, noahmp%config%domain%NumSnowLayerNeg(II,JJ)+2, -1
                   noahmp%energy%state%TemperatureSoilSnow(II,K,JJ)    = noahmp%energy%state%TemperatureSoilSnow(II,K-1,JJ)
                   noahmp%water%state%SnowIce(II,K,JJ)                = noahmp%water%state%SnowIce(II,K-1,JJ)
                   noahmp%water%state%SnowLiqWater(II,K,JJ)           = noahmp%water%state%SnowLiqWater(II,K-1,JJ)
                   noahmp%config%domain%ThicknessSnowSoilLayer(II,K,JJ) = noahmp%config%domain%ThicknessSnowSoilLayer(II,K-1,JJ)

                   if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
                      noahmp%water%state%MassBChydropho(II,K,JJ)      = noahmp%water%state%MassBChydropho(II,K-1,JJ)
                      noahmp%water%state%MassBChydrophi(II,K,JJ)      = noahmp%water%state%MassBChydrophi(II,K-1,JJ)
                      noahmp%water%state%MassOChydropho(II,K,JJ)      = noahmp%water%state%MassOChydropho(II,K-1,JJ)
                      noahmp%water%state%MassOChydrophi(II,K,JJ)      = noahmp%water%state%MassOChydrophi(II,K-1,JJ)
                      noahmp%water%state%MassDust1(II,K,JJ)           = noahmp%water%state%MassDust1(II,K-1,JJ)
                      noahmp%water%state%MassDust2(II,K,JJ)           = noahmp%water%state%MassDust2(II,K-1,JJ)
                      noahmp%water%state%MassDust3(II,K,JJ)           = noahmp%water%state%MassDust3(II,K-1,JJ)
                      noahmp%water%state%MassDust4(II,K,JJ)           = noahmp%water%state%MassDust4(II,K-1,JJ)
                      noahmp%water%state%MassDust5(II,K,JJ)           = noahmp%water%state%MassDust5(II,K-1,JJ)
                      noahmp%water%state%SnowRadius(II,K,JJ)          = noahmp%water%state%SnowRadius(II,K-1,JJ)
                   endif

                enddo
             endif
             ! Decrease the number of snow layers
             noahmp%config%domain%NumSnowLayerNeg(II,JJ) = noahmp%config%domain%NumSnowLayerNeg(II,JJ) + 1
             if ( noahmp%config%domain%NumSnowLayerNeg(II,JJ) >= -1 ) Exit
          else
             ! The layer thickness is greater than the prescribed minimum value
             IndLayer = IndLayer + 1
          endif
       enddo
    endif

  end subroutine SnowLayerCombine

end module SnowLayerCombineMod
