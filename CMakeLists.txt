cmake_minimum_required(VERSION 3.16)
project(cvrpbcp LANGUAGES C CXX)

# cvrpbcp: branch-cut-and-price for CVRP / CVRPTW on bplex.
# bplex is an external dependency, never built from this tree's CMake:
#   cmake -S ../bplex -B build/bplex -DBPLEX_USE_MKL=OFF   (out-of-tree)
# and linked here by path through BPLEX_ROOT (source tree, for the header)
# and BPLEX_BUILD (its build directory, for libbplex.a / libbasiclu.a).
# When bplex ships an install/export target this collapses to find_package.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(BPLEX_ROOT "${CMAKE_SOURCE_DIR}/../bplex" CACHE PATH "bplex source tree (include/bplex.h)")
set(BPLEX_BUILD "${CMAKE_BINARY_DIR}/bplex" CACHE PATH "bplex out-of-tree build directory")
option(CVRPBCP_SANITIZE "Build with ASan/UBSan" OFF)

find_package(Threads REQUIRED)

# ---- bplex (imported, prebuilt) ----
find_library(BPLEX_LIB NAMES bplex PATHS "${BPLEX_BUILD}" NO_DEFAULT_PATH)
find_library(BASICLU_LIB NAMES basiclu PATHS "${BPLEX_BUILD}" NO_DEFAULT_PATH)
if(NOT BPLEX_LIB OR NOT BASICLU_LIB)
  message(FATAL_ERROR "bplex not built: expected libbplex.a and libbasiclu.a in ${BPLEX_BUILD}. "
                      "Run: cmake -S ${BPLEX_ROOT} -B ${BPLEX_BUILD} -DBPLEX_USE_MKL=OFF && cmake --build ${BPLEX_BUILD}")
endif()
add_library(bplex::bplex STATIC IMPORTED)
set_target_properties(bplex::bplex PROPERTIES
  IMPORTED_LOCATION "${BPLEX_LIB}"
  INTERFACE_INCLUDE_DIRECTORIES "${BPLEX_ROOT}/include"
  INTERFACE_LINK_LIBRARIES "${BASICLU_LIB};Threads::Threads;m")

# ---- CVRPSEP (vendored, Lysgaard, CPL 1.0) ----
file(GLOB CVRPSEP_SRC CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/third_party/cvrpsep/*.cpp")
add_library(cvrpsep STATIC ${CVRPSEP_SRC})
target_include_directories(cvrpsep PUBLIC "${CMAKE_SOURCE_DIR}/third_party/cvrpsep")
target_compile_options(cvrpsep PRIVATE -w)   # third-party warnings are not ours

# ---- cvrpbcp library ----
file(GLOB_RECURSE CVRPBCP_SRC CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
list(FILTER CVRPBCP_SRC EXCLUDE REGEX "/src/app/")
add_library(cvrpbcp STATIC ${CVRPBCP_SRC})
target_include_directories(cvrpbcp PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(cvrpbcp PUBLIC bplex::bplex cvrpsep)
target_compile_options(cvrpbcp PRIVATE -Wall -Wextra -Wpedantic)
if(CVRPBCP_SANITIZE)
  target_compile_options(cvrpbcp PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(cvrpbcp PUBLIC -fsanitize=address,undefined)
endif()

# ---- executable ----
add_executable(cvrpbcp_cli src/app/main.cpp)
target_link_libraries(cvrpbcp_cli PRIVATE cvrpbcp)
set_target_properties(cvrpbcp_cli PROPERTIES OUTPUT_NAME cvrpbcp)

# ---- tests ----
enable_testing()
file(GLOB CVRPBCP_TESTS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/tests/test_*.cpp")
foreach(t ${CVRPBCP_TESTS})
  get_filename_component(tn ${t} NAME_WE)
  add_executable(${tn} ${t})
  target_link_libraries(${tn} PRIVATE cvrpbcp)
  add_test(NAME ${tn} COMMAND ${tn} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
