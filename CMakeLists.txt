cmake_minimum_required(VERSION 3.2)
project(EEBO)

find_program(CCACHE_FOUND ccache)
if (CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif (CCACHE_FOUND)

find_package(MPI REQUIRED)

# Only debug version for now
set(METHOD dbg)

# libMesh
if (LIBMESH_DIR)
    set(ENV{PKG_CONFIG_PATH} ${LIBMESH_DIR}/lib/pkgconfig)

    find_package(PkgConfig REQUIRED)
    pkg_check_modules(LIBMESH REQUIRED libmesh-${METHOD})

    add_definitions(${LIBMESH_CFLAGS_OTHER})
    link_directories(${LIBMESH_LIBRARY_DIRS})
endif ()

include_directories(
        ${LIBMESH_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/base
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kernels
        ${CMAKE_CURRENT_SOURCE_DIR}/src/model)

set(LIBS ${LIBMESH_LIBRARIES} pthread)

set(SOURCE_FILES
        src/base/Init.cpp
        src/base/App.cpp
        src/base/FEProblem.h
        src/base/SystemBase.cpp
        src/model/HeatTransfer.cpp)

add_library(eeboCommon ${SOURCE_FILES})

add_executable(fbc_driver
        src/driver/fbc/fbc.cpp)
target_link_libraries(fbc_driver eeboCommon ${LIBS})

message(STATUS "")
message(STATUS "CXXFLAGS:  ${CMAKE_CXX_FLAGS}")
message(STATUS "LIBRARIES: ${LIBS}")
message(STATUS "")

if (ENABLE_TESTING)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests/gtest)
    enable_testing()
    include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

    add_executable(test_main
            tests/test_main.cpp
            tests/test_driver.cpp)
    target_link_libraries(test_main gtest eeboCommon ${LIBS})
    add_test(NAME all COMMAND test_main)
endif ()