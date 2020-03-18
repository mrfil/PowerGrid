    include_directories(../../Support)
    include_directories(../)
    set(POWERGRID_oscillateRecon
            oscillateRecon.cpp
            R_lowRank.cpp
            LRobj.cpp
            )

    add_executable(oscillateRecon ${POWERGRID_oscillateRecon})
    find_library(ISMRMRD_LIBRARIES ismrmrd HINTS /usr/lib/)

    set(LIBS ${LIBS} ${ISMRMRD_LIBRARIES})

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS} ")

    target_link_libraries(oscillateRecon ${LIBS} PowerGrid)
    
    install(TARGETS oscillateRecon DESTINATION bin)