# These need mpi
IF(USE_MPI)
    ADD_LIBRARY(coulomb_shr SHARED coulomb.cxx)
    ADD_LIBRARY(coulomb_stc STATIC coulomb.cxx)
    ADD_LIBRARY(coulombVdW_shr SHARED coulombVdW.cxx)
    ADD_LIBRARY(coulombVdW_stc STATIC coulombVdW.cxx)
    # Rename libraries
    SET_TARGET_PROPERTIES(coulomb_shr PROPERTIES OUTPUT_NAME coulomb)
    SET_TARGET_PROPERTIES(coulombVdW_shr PROPERTIES OUTPUT_NAME coulombVdW)
    SET_TARGET_PROPERTIES(coulomb_stc PROPERTIES OUTPUT_NAME coulomb)
    SET_TARGET_PROPERTIES(coulombVdW_stc PROPERTIES OUTPUT_NAME coulombVdW)

    IF(USE_GPU)
        cuda_add_library(mr3 STATIC mr3.cu)
        TARGET_LINK_LIBRARIES(mr3 ${CUDA_LIBRARIES} gpuKernel)
    ELSE()
        ADD_LIBRARY(mr3 STATIC mr3_host.cxx)
        TARGET_LINK_LIBRARIES(mr3 cpuKernel)
    ENDIF()

    ADD_EXECUTABLE(test_coulomb_shr test_coulomb.cxx)
    TARGET_LINK_LIBRARIES(test_coulomb_shr coulomb_shr mr3)

    ADD_EXECUTABLE(test_coulomb_stc test_coulomb.cxx)
    TARGET_LINK_LIBRARIES(test_coulomb_stc coulomb_stc mr3)

    ADD_EXECUTABLE(test_coulombVdW_shr test_coulombVdW.cxx)
    TARGET_LINK_LIBRARIES(test_coulombVdW_shr coulombVdW_shr mr3)

    ADD_EXECUTABLE(test_coulombVdW_stc test_coulombVdW.cxx)
    TARGET_LINK_LIBRARIES(test_coulombVdW_stc coulombVdW_stc mr3)
    ADD_EXECUTABLE(ewald ewald.cxx)
    TARGET_LINK_LIBRARIES(ewald coulombVdW_stc mr3)
    IF(BUILD_TESTING)
        add_test(coulomb_shr ${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_coulomb_shr)
        add_test(coulombVdW_shr ${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_coulombVdW_shr)
        add_test(coulomb_stc ${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_coulomb_stc)
        add_test(coulombVdW_stc ${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_coulombVdW_stc)
        add_test(test_ewald_1 ${CMAKE_CURRENT_BINARY_DIR}/ewald 0 12)
        add_test(test_ewald_2 ${CMAKE_CURRENT_BINARY_DIR}/ewald 1 12)
    ENDIF()
ENDIF()


