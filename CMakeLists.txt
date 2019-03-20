
# CMake build script written by: Matthew V. Bilskie, PhD
# 03-20-2019

######################################################################
# GENERAL OPTIONS
######################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT("GRD2DEM")
######################################################################

######################################################################
#...Langauge Specifications
######################################################################
ENABLE_LANGUAGE(Fortran)
######################################################################

#... Make sure the defalt build is a RELEASE
IF(NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are: None, Debug & Release."
        FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#... FFLAGS depend on the compiler
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF(Fortran_Compiler_NAME MATCHES "gfortran.*")
    # gfortran
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
ELSEIF(Fortran_COMPILER_NAME MATCHES "ifort.*")
    # ifort
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -assume buffered_io")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
ENDIF()

# Build Executables
ADD_EXECUTABLE(GRD2DEM src/GRD2DEM.F90)

# Install Executables
INSTALL(TARGETS GRD2DEM RUNTIME DESTINATION bin)

