module GroundAlbedoMod

!!! Compute ground albedo based on soil and snow albedo

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine GroundAlbedo(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: GROUNDALB
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I, J                 ! grid indices
    integer                          :: IndSwBnd             ! solar radiation band index
    real(kind=kind_noahmp)           :: AlbedoSoilAdjWet     ! soil water correction factor for soil albedo

   !$acc parallel loop collapse(2) gang vector present(noahmp) private(IndSwBnd,AlbedoSoilAdjWet)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

! --------------------------------------------------------------------
    associate(                                                                      &
              NumSwRadBand        => noahmp%config%domain%NumSwRadBand             ,& ! in,  number of solar radiation wave bands
              SurfaceType         => noahmp%config%domain%SurfaceType(I,J)         ,& ! in,  surface type 1-soil; 2-lake
              CosSolarZenithAngle => noahmp%config%domain%CosSolarZenithAngle(I,J) ,& ! in,  cosine solar zenith angle
              SnowCoverFrac       => noahmp%water%state%SnowCoverFrac(I,J)         ,& ! in,  snow cover fraction
              SoilMoisture        => noahmp%water%state%SoilMoisture          ,& ! in,  total soil moisture [m3/m3]
              AlbedoSoilSat       => noahmp%energy%param%AlbedoSoilSat        ,& ! in,  saturated soil albedos: 1=vis, 2=nir
              AlbedoSoilDry       => noahmp%energy%param%AlbedoSoilDry        ,& ! in,  dry soil albedos: 1=vis, 2=nir
              AlbedoLakeFrz       => noahmp%energy%param%AlbedoLakeFrz        ,& ! in,  albedo frozen lakes: 1=vis, 2=nir
              TemperatureGrd      => noahmp%energy%state%TemperatureGrd(I,J)       ,& ! in,  ground temperature [K]
              AlbedoSnowDir       => noahmp%energy%state%AlbedoSnowDir        ,& ! in,  snow albedo for direct(1=vis, 2=nir)
              AlbedoSnowDif       => noahmp%energy%state%AlbedoSnowDif        ,& ! in,  snow albedo for diffuse(1=vis, 2=nir)
              AlbedoSoilDir       => noahmp%energy%state%AlbedoSoilDir        ,& ! out, soil albedo (direct)
              AlbedoSoilDif       => noahmp%energy%state%AlbedoSoilDif        ,& ! out, soil albedo (diffuse)
              AlbedoGrdDir        => noahmp%energy%state%AlbedoGrdDir         ,& ! out, ground albedo (direct beam: vis, nir)
              AlbedoGrdDif        => noahmp%energy%state%AlbedoGrdDif          & ! out, ground albedo (diffuse: vis, nir)
             )
! ----------------------------------------------------------------------

    ! solar radiation process is only done if there is light
    if ( CosSolarZenithAngle < 0 ) cycle

    !$acc loop seq
    do IndSwBnd = 1, NumSwRadBand

       AlbedoSoilAdjWet = max(0.11-0.40*SoilMoisture(I,1,J), 0.0)

       if ( SurfaceType == 1 )  then                      ! soil
          AlbedoSoilDir(I,IndSwBnd,J) = min(AlbedoSoilSat(I,IndSwBnd,J)+AlbedoSoilAdjWet, AlbedoSoilDry(I,IndSwBnd,J))
          AlbedoSoilDif(I,IndSwBnd,J) = AlbedoSoilDir(I,IndSwBnd,J)
       elseif ( TemperatureGrd > ConstFreezePoint ) then  ! unfrozen lake, wetland
          AlbedoSoilDir(I,IndSwBnd,J) = 0.06 / (max(0.01, CosSolarZenithAngle)**1.7+0.15)
          AlbedoSoilDif(I,IndSwBnd,J) = 0.06
       else                                               ! frozen lake, wetland
          AlbedoSoilDir(I,IndSwBnd,J) = AlbedoLakeFrz(I,IndSwBnd,J)
          AlbedoSoilDif(I,IndSwBnd,J) = AlbedoSoilDir(I,IndSwBnd,J)
       endif

       AlbedoGrdDir(I,IndSwBnd,J) = AlbedoSoilDir(I,IndSwBnd,J)*(1.0-SnowCoverFrac) + AlbedoSnowDir(I,IndSwBnd,J)*SnowCoverFrac
       AlbedoGrdDif(I,IndSwBnd,J) = AlbedoSoilDif(I,IndSwBnd,J)*(1.0-SnowCoverFrac) + AlbedoSnowDif(I,IndSwBnd,J)*SnowCoverFrac

    enddo

        end associate

      end do
    end do
   !$acc end parallel loop

  end subroutine GroundAlbedo

end module GroundAlbedoMod
