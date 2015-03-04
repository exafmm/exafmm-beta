# Reject in-place build
IF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  MESSAGE(FATAL_ERROR "Rejecting in-place build. Please use a build directory")
ENDIF()
# Project header
PROJECT(ExaFMM)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# All options are set here
SET(EQUATION Laplace) # {Laplace, Yukawa, Helmholtz, Stokes}
SET(BASIS Cartesian) # {Cartesian, Spherical, Rotation, Planewave}
OPTION(USE_MPI "Use MPI" ON)
OPTION(USE_GPU "Use GPUs" OFF)
OPTION(USE_VTK "Use VTK" OFF)

# MPI
IF(USE_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  MESSAGE(STATUS "Enabling MPI")
  SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
ENDIF()

# CUDA
IF(USE_GPU)
  SET(DEVICE GPU)
  FIND_PACKAGE(CUDA REQUIRED)
  MESSAGE(STATUS "Enabling CUDA ${CUDA_VERSION}")
  IF(CUDA_SDK_ROOT_DIR MATCHES "CUDA_SDK_ROOT_DIR-NOTFOUND")
    MESSAGE(FATAL_ERROR "Could not locate CUDA SDK. Try setting environment variable NVSDKCOMPUTE_ROOT")
  ENDIF()
  SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -Xcompiler -fopenmp --ptxas-options=-v -O3 -use_fast_math -arch=sm_21)
  INCLUDE_DIRECTORIES(${CUDA_SDK_ROOT_DIR}/common/inc)
ELSE()
  SET(DEVICE CPU)
ENDIF()

# VTK
IF(USE_VTK)
  FIND_PACKAGE(VTK REQUIRED)
  MESSAGE(STATUS "Enabling VTK ${VTK_VERSION}")
  ADD_DEFINITIONS(-DVTK)
  INCLUDE(${VTK_USE_FILE})
  LINK_LIBRARIES(vtkCommon vtkFiltering vtkGraphics vtkIO vtkRendering vtkViews vtkWidgets)
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wuninitialized -O3")
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)
ADD_DEFINITIONS(-D${EQUATION})
ADD_DEFINITIONS(-D${BASIS})
ADD_DEFINITIONS(-D${DEVICE})
ADD_SUBDIRECTORY(kernels)

# Unit tests
INCLUDE(CTest)
ADD_SUBDIRECTORY(examples)
