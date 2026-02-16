module LanduseConvertMod

  use ConstantDefineMod

  implicit none


contains

  function LanduseConvert(LU_string) result(LU_code)
! ------------------------ Code history -----------------------------------
! OpenACC Port: D. Reynolds (2026)
! -------------------------------------------------------------------------

    character(len=*), intent(in)  :: LU_string
    integer                       :: LU_code

    if (trim(LU_string) == 'USGS') then
      LU_code = ConstLU_USGS
    else if (trim(LU_string) == 'MODIFIED_IGBP_MODIS_NOAH') then
      LU_code = ConstLU_IGBP_MODIS_NOAH
    else
      LU_code = -1  ! Invalid landuse type
    end if

  end function LanduseConvert

  function LanduseConvertToString(LU_code) result(LU_string)
! ------------------------ Code history -----------------------------------
! OpenACC Port: D. Reynolds (2026)
! -------------------------------------------------------------------------

    integer,          intent(in)  :: LU_code
    character(len=32) :: LU_string

    select case (LU_code)
      case (ConstLU_USGS)
        LU_string = 'USGS'
      case (ConstLU_IGBP_MODIS_NOAH)
        LU_string = 'MODIFIED_IGBP_MODIS_NOAH'
      case default
        LU_string = 'UNKNOWN'
    end select

  end function LanduseConvertToString

end module LanduseConvertMod
