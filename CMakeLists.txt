cmake_minimum_required(VERSION 3.2)
project(EEBO)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unknown-pragmas")

find_program(CCACHE_FOUND ccache)
if (CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif (CCACHE_FOUND)

find_package(MPI REQUIRED)

# Only debug version for now
set(METHOD dbg)

if (METHOD STREQUAL "dbg")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -felide-constructors \
    -Wextra -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses -Woverloaded-virtual -Wno-variadic-macros")
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
        ${CMAKE_CURRENT_SOURCE_DIR}/src/base
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kernels
        ${CMAKE_CURRENT_SOURCE_DIR}/src/model)

set(LIBS ${LIBMESH_LIBRARIES})

set(SOURCE_FILES
        src/base/Init.cc
        src/base/App.cc
        src/base/FEProblem.hh
        src/model/HeatTransfer.cc
        )
add_library(eeboCommon ${SOURCE_FILES})

add_executable(fbc_driver
        src/driver/fbc/fbc.cc)
target_link_libraries(fbc_driver eeboCommon ${LIBS})

message(STATUS "")
message(STATUS "CXXFLAGS:  ${CMAKE_CXX_FLAGS}")
message(STATUS "LIBRARIES: ${LIBS}")
message(STATUS "")
