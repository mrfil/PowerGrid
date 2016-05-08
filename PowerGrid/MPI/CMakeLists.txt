    add_definitions(-DPowerGridMPI)
    find_package(MPI)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    include_directories(../../Support)
    set(mpi_POWERGRID_DWI3DDFT
            mpiPowerGridDWI3DDFT.cpp
            ../config.cxx
            )
    add_executable(mpiPowerGrid3DDFT ${mpi_POWERGRID_DWI3DDFT})
    set(mpi_POWERGRID_STEAM
            mpiPowerGridSTEAM.cpp
            ../config.cxx
            )
    add_executable(mpiPowerGridSTEAM ${mpi_POWERGRID_STEAM})
    set(mpiLIBS ${LIBS})
    set(mpiLIBS ${mpiLIBS} ${MPI_CXX_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}  ${MPI_CXX_LINK_FLAGS}")
    target_link_libraries(mpiPowerGrid3DDFT ${mpiLIBS})
    target_link_libraries(mpiPowerGridSTEAM ${mpiLIBS})