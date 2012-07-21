ADD_EXECUTABLE(serial serial.cxx)
TARGET_LINK_LIBRARIES(serial Kernels)
ADD_TEST(serial ${CMAKE_CURRENT_BINARY_DIR}/serial)

IF(USE_MPI)
  ADD_EXECUTABLE(parallel parallel.cxx)
  TARGET_LINK_LIBRARIES(parallel Kernels)
  ADD_TEST(parallel ${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/parallel)
ENDIF()
