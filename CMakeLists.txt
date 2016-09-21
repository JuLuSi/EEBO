cmake_minimum_required(VERSION 3.2)
project(EEBO)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unknown-pragmas -Wpedantic -Wdisabled-optimization")

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# Only debug version for now
set(METHOD dbg)

if (METHOD STREQUAL "dbg")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC")
endif ()

# libMesh
if (LIBMESH_DIR)
    set(ENV{PKG_CONFIG_PATH} ${LIBMESH_DIR}/lib/pkgconfig)
endif ()
find_package(PkgConfig REQUIRED)
pkg_check_modules(LIBMESH REQUIRED libmesh-${METHOD})

link_directories(${LIBMESH_LIBRARY_DIRS})

include_directories(
        ${LIBMESH_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/common)

set(LIBS ${LIBMESH_LIBRARIES})

set(SOURCE_FILES
        src/common/EEBOApp.cc)
add_library(eeboCommon ${SOURCE_FILES})

add_executable(fbc_driver
        src/driver/fbc/fbc.cc)
target_link_libraries(fbc_driver eeboCommon ${LIBS})

message(STATUS "")
message(STATUS "CXXFLAGS:  ${CMAKE_CXX_FLAGS}")
message(STATUS "LIBRARIES: ${LIBS}")
message(STATUS "")