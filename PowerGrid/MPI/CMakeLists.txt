    add_definitions(-DPowerGridMPI)

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

    set(mpi_POWERGRID_PcSense_TS
    mpipcSENSETimeSeg.cpp
    PowerGridPcSenseMPI_TS.cpp
    )

    add_executable(PowerGridPcSenseMPI_TS ${mpi_POWERGRID_PcSense_TS})
    find_library(ISMRMRD_LIBRARIES ismrmrd HINTS /usr/lib/)

    set(mpi_POWERGRID_Sense
    PowerGridSenseMPI.cpp
    )

    add_executable(PowerGridSenseMPI ${mpi_POWERGRID_Sense})
    find_library(ISMRMRD_LIBRARIES ismrmrd HINTS /usr/lib/)

    set(LIBS ${LIBS} ${ISMRMRD_LIBRARIES})
    set(mpiLIBS ${LIBS})
    set(mpiLIBS ${mpiLIBS} ${MPI_CXX_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}  ${MPI_CXX_LINK_FLAGS}")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/lib/openmpi/include/openmpi/opal/mca/event/libevent2021/libevent -I/usr/lib/openmpi/include/openmpi/opal/mca/event/libevent2021/libevent/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi")
    #set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS} -Wl,-rpath -Wl,/usr/lib/openmpi/lib -Wl,--enable-new-dtags -L/usr/lib/openmpi/lib -lmpi_cxx -lmpi ")


    target_link_libraries(PowerGridPcSenseMPI ${mpiLIBS} PowerGrid)
    target_link_libraries(PowerGridSenseMPI ${mpiLIBS} PowerGrid)
    target_link_libraries(PowerGridPcSenseMPI_TS ${mpiLIBS} PowerGrid)
    
    install(TARGETS PowerGridPcSenseMPI DESTINATION bin)
    install(TARGETS PowerGridSenseMPI DESTINATION bin)
    install(TARGETS PowerGridPcSenseMPI_TS DESTINATION bin)