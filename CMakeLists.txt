cmake_minimum_required(VERSION 3.16)
project(GEM_QuarkModel VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Collect source files
file(GLOB_RECURSE SOURCES
    ${CMAKE_SOURCE_DIR}/src/*.cpp
)

# Create library
add_library(gem_quark STATIC ${SOURCES})
target_link_libraries(gem_quark PUBLIC Eigen3::Eigen)
target_include_directories(gem_quark PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Main executable
add_executable(gem_calc apps/main.cpp)
target_link_libraries(gem_calc PRIVATE gem_quark)

# Charmonium test
add_executable(test_charmonium apps/test_charmonium.cpp)
target_link_libraries(test_charmonium PRIVATE gem_quark)

# Debug potential test
add_executable(debug_potential apps/debug_potential.cpp)
target_link_libraries(debug_potential PRIVATE gem_quark)

# Debug charmonium test
add_executable(debug_charmonium apps/debug_charmonium.cpp)
target_link_libraries(debug_charmonium PRIVATE gem_quark)

# Debug compare test
add_executable(debug_compare apps/debug_compare.cpp)
target_link_libraries(debug_compare PRIVATE gem_quark)

# Testing
option(BUILD_TESTS "Build unit tests" ON)
if(BUILD_TESTS)
    enable_testing()
    find_package(GTest QUIET)

    if(GTest_FOUND)
        file(GLOB TEST_SOURCES ${CMAKE_SOURCE_DIR}/tests/*.cpp)
        foreach(TEST_FILE ${TEST_SOURCES})
            get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
            add_executable(${TEST_NAME} ${TEST_FILE})
            target_link_libraries(${TEST_NAME} PRIVATE gem_quark GTest::gtest GTest::gtest_main)
            add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
        endforeach()
    else()
        message(STATUS "GoogleTest not found, tests will not be built")
    endif()
endif()

# Installation
install(TARGETS gem_quark gem_calc
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(DIRECTORY include/ DESTINATION include)

add_executable(test_charmed apps/test_charmed.cpp)
target_link_libraries(test_charmed PRIVATE gem_quark)
