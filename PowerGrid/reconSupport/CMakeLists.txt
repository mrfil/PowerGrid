    include_directories(../../Support)
    include_directories(../)
    set(POWERGRID_GridCoilImages
            gridCoilImages.cpp
            )

    add_executable(gridCoilImages ${POWERGRID_GridCoilImages})
    find_library(ISMRMRD_LIBRARIES ismrmrd HINTS /usr/lib/)

    set(LIBS ${LIBS} ${ISMRMRD_LIBRARIES})

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS} ")

    target_link_libraries(gridCoilImages ${LIBS} PowerGrid)
    
    install(TARGETS gridCoilImages DESTINATION bin)
