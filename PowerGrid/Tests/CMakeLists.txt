cmake_minimum_required(VERSION 3.2) #CMake version 3.2 required to properly support CMAKE_CXX_STANDARD

## Flags and settings that change the compiler (i.e. from default to C++11 or C++14)
## need to be set BEFORE calling project command (which does introspection and sets
## variables about the behavior of the compiler).
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_EXTENSIONS ON)

project("PowerGridTests" CXX)


# Prepare "Catch" library for other executables
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../Support/Catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

set(TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pgComplexTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pgColTests.cpp  
    ${CMAKE_CURRENT_SOURCE_DIR}/pgColCplxTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pgMatTests.cpp   
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -ta=tesla -fast")

# Make test executable
add_executable(tests ${TEST_SOURCES})
target_link_libraries(tests Catch)