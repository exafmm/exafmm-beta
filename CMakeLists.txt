PROJECT(exafmm)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
       message(FATAL_ERROR "Rejecting in-place build. Please use a build directory")
ENDIF()

## Setup configuration parameters ##

option(USE_STOKES "Use Stokes kernels" OFF)
IF(USE_STOKES)
    SET(USE_CUDA ON)
    SET(USE_SPHERICAL ON)
    ADD_DEFINITIONS(-DSTOKES)
ENDIF()

# Use periodic images
OPTION(USE_IMAGES "Use periodic images" OFF)
IF(USE_IMAGES)
    ADD_DEFINITIONS(-DIMAGES)
ENDIF()

# Use topdown tree construction
OPTION(USE_TOPDOWN "Use topdown tree construction" OFF)
IF(USE_TOPDOWN)
    ADD_DEFINITIONS(-DTOPDOWN)
ENDIF()

# Optimization for spark processors
OPTION(USE_SPARC_SIMD "Use Spark optimizations" OFF)
IF(USE_TOPDOWN)
    ADD_DEFINITIONS(-DSPARC_SIMD)
ENDIF()

# Use spherical harmonics
OPTION(USE_SPHERICAL "Use Spherical Harmonics Expansions" OFF)
IF(USE_SPHERICAL)
    ADD_DEFINITIONS(-DSPHERICAL)
ENDIF()

OPTION(USE_TREECODE "Use treecode methods" OFF)
OPTION(USE_HYBRID "Use a hybrid algorithm (FMM+Treecode)" OFF)

IF(USE_TREECODE)
    add_definitions(-DTREECODE)
ENDIF()

IF(USE_HYBRID)
    add_definitions(-DHYBRID)
ENDIF()
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)
ADD_DEFINITIONS(-D_FILE_OFFSET_BITS=64)

#FindVTK
OPTION(USE_VTK "Use VTK" OFF)
IF(USE_VTK)
    FIND_PACKAGE(VTK REQUIRED)
    ADD_DEFINITIONS(-DVTK)
    INCLUDE(${VTK_USE_FILE})
    SET(VTK_LIBS vtkIO QVTK vtkRendering vtkGraphics vtkCommon)
    MESSAGE(STATUS "Enabling VTK")
ENDIF()

#FindMPI
OPTION(USE_MPI "Use MPI" OFF)
IF(USE_MPI)
    FIND_PACKAGE(MPI REQUIRED)
    SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
    MESSAGE(STATUS "Enabling MPI")
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
ENDIF()

OPTION(USE_GPU "Build GPU kernels" OFF)
IF(USE_GPU)
    SET(USE_CUDA ON)
    ADD_DEFINITIONS(-DGPU)
ENDIF()

#FindCUDA
IF(USE_CUDA)
    FIND_PACKAGE(CUDA REQUIRED)
    SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode arch=compute_20,code=sm_20 --ptxas-options=-v -use_fast_math)
    IF(NOT CMAKE_COMPILER_IS_GNUCXX)
        SET(CUDA_PROPAGATE_HOST_FLAGS OFF)
    ENDIF()
    MESSAGE(STATUS "Enabling CUDA ${CUDA_VERSION}")
    IF(${CUDA_VERSION} MATCHES "4.1")
        ADD_DEFINITIONS(-DCUDA_ENABLE_MEMORY_MAPPING)
    ENDIF()
    INCLUDE_DIRECTORIES(${CUDA_SDK_ROOT_DIR}/C/common/inc)
ENDIF()

# Use QUARK
OPTION(USE_QUARK "Use Quark Tasking Library" OFF)
IF(USE_QUARK)
    ADD_DEFINITIONS(-DQUARK)
    MESSAGE(STATUS "Enabling Quark")
    ADD_SUBDIRECTORY(external/quark)
    INCLUDE_DIRECTORIES(external/quark)
ENDIF(USE_QUARK)


INCLUDE(CTest)
IF(BUILD_TESTING)
    add_subdirectory(unit_test)
ENDIF(BUILD_TESTING)

# Set documentation target
OPTION(BUILD_DOCS "Build doxygen documentation" OFF)
IF(BUILD_DOCS)
    SET(EXAFMM_VERSION 0.0)
    configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    )
    ADD_CUSTOM_TARGET(docs ALL COMMAND doxygen ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()


# Set/Check compiler flags
INCLUDE(CheckCXXCompilerFlag)
IF(CMAKE_COMPILER_IS_GNUCXX)
    OPTION(ENABLE_WARNINGS "Enable/Disable Warnings in g++" OFF)
    IF(ENABLE_WARNINGS)
        SET(G++_WARNINGS "-Wuninitialized -Wshadow -Weffc++ -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wold-style-cast")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${G++_WARNINGS}")
        check_cxx_compiler_flag("-Wno-variadic-macros" COMPILER_SUPPORT_WNOVARIADICMACRO)
        IF(COMPILER_SUPPORT_WNOVARIADICMACRO)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-variadic-macros")
        ENDIF()
        check_cxx_compiler_flag("-Wextra" COMPILER_SUPPORT_WEXTRA)
        IF(COMPILER_SUPPORT_WEXTRA)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
        ENDIF()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
    ENDIF()
    OPTION(ENABLE_EXTRA_OPT_FLAGS "Enable/Disable Fags in g++" OFF)
    IF(ENABLE_EXTRA_OPT_FLAGS)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -fno-check-new -fno-common -fstrict-aliasing -ffast-math -funroll-loops -fforce-addr")
        check_cxx_compiler_flag("-fPIC" COMPILER_SUPPORT_PIC)
        IF(COMPILER_SUPPORT_PIC)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        ENDIF()
    ENDIF()
    check_cxx_compiler_flag("-fopenmp" CXXCOMPILER_SUPPORT_OPENMP)
    IF(CXXCOMPILER_SUPPORT_OPENMP)
        OPTION(ENABLE_OPENMP_GLOBALY "Enable/Disable GCC OpenMP" OFF)
        IF(ENABLE_OPENMP_GLOBALY)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            MESSAGE(STATUS "Enabling OpenMP")
        ENDIF()
    ENDIF()
ELSEIF(CMAKE_CXX_COMPILER MATCHES "/icpc$" OR CMAKE_CXX_COMPILER MATCHES "/icpc.exe$")
    SET(CMAKE_AR ${XIAR})
    OPTION(ENABLE_OPENMP_GLOBALY "Enable/Disable Intel OpenMP" OFF)
    IF(ENABLE_OPENMP_GLOBALY)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
        IF(USE_CUDA)
            SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -Xcompiler -fopenmp -lgomp)
        ENDIF()
        MESSAGE(STATUS "Enabling Intel OpenMP")
    ENDIF()
    OPTION(ENABLE_EXTRA_OPT_FLAGS "Enable/Disable Fags in icpc" OFF)
    IF(ENABLE_EXTRA_OPT_FLAGS)
        IF(COMPILER_SUPPORT_PIC)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        ENDIF()
        SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -unroll-aggressive -parallel -fomit-frame-pointer -ipo -static -no-prec-div -xHost")
    ENDIF()
ENDIF()

ADD_SUBDIRECTORY(kernel)
ADD_SUBDIRECTORY(fast)
ADD_SUBDIRECTORY(wrapper)

