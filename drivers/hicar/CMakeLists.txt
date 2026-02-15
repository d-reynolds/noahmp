cmake_minimum_required(VERSION 3.20)
project(NoahMP LANGUAGES Fortran)


set(NOAHMP_MASTER_PROJECT OFF)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(NOAHMP_MASTER_PROJECT ON)
endif()

if(NOT DEFINED NETCDF_INCLUDES AND NOT DEFINED NETCDF_LIBRARIES)
    if (NOT DEFINED NETCDF_DIR OR NETCDF_DIR STREQUAL "")
        set(NETCDF_DIR $ENV{NETCDF_DIR} CACHE PATH "Path to NETCDF installation")
    endif()
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

#Remove SnowInputSnicarMod.F90 and SnowRadiationSnicarMod.F90 from the build
list(FILTER NOAHMP_SOURCES EXCLUDE REGEX "SnowInputSnicarMod\\.F90$")
list(FILTER NOAHMP_SOURCES EXCLUDE REGEX "SnowRadiationSnicarMod\\.F90$")

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

if (NETCDF_FOUND)
    # Check for parallel NetCDF support
    execute_process(COMMAND nc-config --has-parallel
                    OUTPUT_VARIABLE NETCDF_HAS_PARALLEL
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    
    if(NETCDF_HAS_PARALLEL STREQUAL "yes")
        message(STATUS "NetCDF built with parallel support")
        add_definitions(-DNETCDF_PARALLEL)
    else()
        message(STATUS "NetCDF built WITHOUT parallel support")
    endif()
    
	find_program(NF_CONFIG nf-config)
	if(NF_CONFIG)
		execute_process(COMMAND ${NF_CONFIG} --flibs 
						OUTPUT_VARIABLE NETCDF_LIBRARIES 
						OUTPUT_STRIP_TRAILING_WHITESPACE
						ERROR_QUIET)
		separate_arguments(NETCDF_LIBRARIES)

		# Clean up any duplicate libraries and references to -lnetcdf
		list(REMOVE_DUPLICATES NETCDF_LIBRARIES)
		# clean references to base netcdf-c library
		# list(REMOVE_ITEM NETCDF_LIBRARIES -lnetcdf)

		message(STATUS "NetCDF libraries (with dependencies): ${NETCDF_LIBRARIES}")
	else()
		message(WARNING "nf-config not found, falling back to manual dependency resolution")
		find_package(HDF5 COMPONENTS C HL)
		if(HDF5_FOUND)
			target_link_libraries(noahmp PUBLIC ${HDF5_LIBRARIES})
		endif()
	endif()
    
endif()

#set debug flags
set(noahmp_compiler_flags -O3 -fopenmp -ffast-math -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops -march=native -fno-protect-parens)
target_compile_options(noahmp PUBLIC ${noahmp_compiler_flags})


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