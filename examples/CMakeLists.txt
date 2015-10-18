ADD_EXECUTABLE(serial fmm.cxx)
TARGET_LINK_LIBRARIES(serial Kernels P2PKernels)
ADD_TEST(serial ${CMAKE_CURRENT_BINARY_DIR}/serial)

IF(USE_MPI)
  ADD_EXECUTABLE(parallel fmm_mpi.cxx)
  TARGET_LINK_LIBRARIES(parallel Kernels P2PKernels)
  ADD_TEST(parallel ${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/parallel)
ENDIF()
