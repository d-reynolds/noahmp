module SnowAlbedoSnicarMod

!!! Compute snow albedo based on SNICAR scheme (Flanner et al. (2021) GMD)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  ! use SnowRadiationSnicarMod, only : SnowRadiationSnicar

  implicit none

contains

  subroutine SnowAlbedoSnicar(noahmp)

! ------------------------ Code history -----------------------------------
! Implementation: T.-S. Lin, C. He, et al. (2025, JHM)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: FlagSwRadType  ! flag: 1 for direct-beam incident flux, 2 for diffuse incident flux
    integer                          :: I, J           ! grid indices
    integer                          :: LoopInd        ! loop index for solar radiation bands

    !$acc parallel loop collapse(2) gang vector present(noahmp)
    do J = noahmp%config%domain%JTS, noahmp%config%domain%JTE
      do I = noahmp%config%domain%ITS, noahmp%config%domain%ITE

        ! solar radiation process is only done if there is light
        if ( noahmp%config%domain%CosSolarZenithAngle(I,J) <= 0 ) cycle

! --------------------------------------------------------------------
    associate(                                                          &
              NumSwRadBand  => noahmp%config%domain%NumSwRadBand ,& ! in,  number of solar radiation wave bands
              AlbedoSnowDir => noahmp%energy%state%AlbedoSnowDir ,& ! out, snow albedo for direct (1=vis, 2=nir)
              AlbedoSnowDif => noahmp%energy%state%AlbedoSnowDif  & ! out, snow albedo for diffuse (1=vis, 2=nir)
             )
! ----------------------------------------------------------------------

    ! initialization
    !$acc loop seq
    do LoopInd = 1, NumSwRadBand
      AlbedoSnowDir(I,LoopInd,J) = 0.0
      AlbedoSnowDif(I,LoopInd,J) = 0.0
    enddo

    end associate
      end do
    end do

    ! FlagSwRadType = 1 ! Direct
    ! call SnowRadiationSnicar(noahmp,FlagSwRadType) 

    ! FlagSwRadType = 2 ! Diffuse
    ! call SnowRadiationSnicar(noahmp,FlagSwRadType)
    

  end subroutine SnowAlbedoSnicar

end module SnowAlbedoSnicarMod
