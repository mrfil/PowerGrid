 
    find_package(MPI)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    include_directories(../../Support)
    include_directories(../)
    set(mpi_POWERGRID_PcSense
            mpipcSENSE.cpp
            PowerGridPcSenseMPI.cpp
            )

    add_executable(PowerGridPcSenseMPI ${mpi_POWERGRID_PcSense})
    find_library(ISMRMRD_LIBRARIES ismrmrd HINTS /usr/lib/)
    set(LIBS ${LIBS} ${ISMRMRD_LIBRARIES})
    set(mpiLIBS ${LIBS})
    set(mpiLIBS ${mpiLIBS} ${MPI_CXX_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}  ${MPI_CXX_LINK_FLAGS}")
    target_link_libraries(PowerGridPcSenseMPI ${mpiLIBS} PowerGrid)

    install(TARGETS PowerGridPcSenseMPI DESTINATION bin)