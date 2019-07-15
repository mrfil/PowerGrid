cmake_minimum_required(VERSION 3.2) #CMake version 3.2 required to properly support CMAKE_CXX_STANDARD

## Flags and settings that change the compiler (i.e. from default to C++11 or C++14)
## need to be set BEFORE calling project command (which does introspection and sets
## variables about the behavior of the compiler).
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_EXTENSIONS ON)

project("PowerGridpgColBench" CXX)

set(TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -Minfo=acc -ta=nvidia:cc60,fastmath,flushz -fast -Minline -fpermissive ")
add_definitions(-DARMA_NO_DEBUG)
# Make test executable
add_executable(pgCol ${TEST_SOURCES})
#target_link_libraries(tests Catchf