cmake_minimum_required(VERSION 3.20)
project(NoahMP LANGUAGES Fortran)


set(NOAHMP_MASTER_PROJECT OFF)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(NOAHMP_MASTER_PROJECT ON)
endif()

if(NOT DEFINED NETCDF_INCLUDES AND NOT DEFINED NETCDF_LIBRARIES)
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)
    set(CMAKE_PREFIX_PATH ${NETCDF_DIR} ${CMAKE_PREFIX_PATH})
    set(NETCDF_F90 "YES")
    find_package(NetCDF REQUIRED)
endif()

# Collect all Fortran source files
file(GLOB_RECURSE NOAHMP_SOURCES 
    "${CMAKE_CURRENT_SOURCE_DIR}/../../src/*.F90" 
    "${CMAKE_CURRENT_SOURCE_DIR}/../../utility/*.F90" 
    "${CMAKE_CURRENT_SOURCE_DIR}/*.F90"
)

# Create a static library
add_library(noahmp STATIC ${NOAHMP_SOURCES})

# Set preprocessor to handle .F90 files (uppercase = needs preprocessing)
set_source_files_properties(${NOAHMP_SOURCES} PROPERTIES Fortran_PREPROCESS ON)

# Compiler flags are inherited from parent project via CMAKE_Fortran_FLAGS
# which is set in the parent CMakeLists.txt before add_subdirectory()
# No need to explicitly set target_compile_options here

# Set module directory
set_target_properties(noahmp PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

# Make modules available to parent project
target_include_directories(noahmp INTERFACE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

# Set include directories matching the Makefile
target_include_directories(noahmp PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../../utility
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add NetCDF include directories (equivalent to NETCDFMOD)
if(DEFINED NETCDF_INCLUDES)
    target_include_directories(noahmp PUBLIC ${NETCDF_INCLUDES})
endif()


# Link NetCDF libraries if available
if(DEFINED NETCDF_LIBRARIES)
    target_link_libraries(noahmp PUBLIC ${NETCDF_LIBRARIES})
endif()

# Install the library (optional)
install(TARGETS noahmp DESTINATION ${CMAKE_INSTALL_LIBDIR})

# install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/ 
#         DESTINATION ${CMAKE_INSTALL_EXEDIR}
#         FILES_MATCHING PATTERN "*.mod")