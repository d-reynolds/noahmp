module NoahmpVarType

!!! Define 2D Noah-MP model variable data types (GPU-optimized)

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! GPU port (2D arrays): Full SoA transformation for OpenACC (2026)
! -------------------------------------------------------------------------

  use ForcingVarType
  use ConfigVarType
  use EnergyVarType
  use WaterVarType
  use BiochemVarType

  implicit none
  save
  private

  type, public :: noahmp_type

    ! Tile bounds for 2D domain
    integer :: ITS, ITE  ! Start/end indices in I direction
    integer :: JTS, JTE  ! Start/end indices in J direction

    ! Domain layer configuration
    integer :: NumSoilLayer     ! Number of soil layers (e.g., 4)
    integer :: NumSnowLayerMax  ! Maximum snow layers (e.g., 3)

    ! define specific variable types for Noah-MP (now with 2D arrays)
    type(forcing_type)  :: forcing
    type(config_type)   :: config
    type(energy_type)   :: energy
    type(water_type)    :: water
    type(biochem_type)  :: biochem

  end type noahmp_type

end module NoahmpVarType
