cmake_minimum_required(VERSION 3.16)
project(cvrpbcp LANGUAGES C CXX)

# cvrpbcp: branch-cut-and-price for CVRP / CVRPTW on bplex.
# bplex is an external dependency consumed as an installed CMake package
# (API 0.2+): build it out of tree, install it to a prefix and put that prefix
# on CMAKE_PREFIX_PATH (scripts/build_wsl.sh does all three):
#   cmake -S ../bplex -B build/bplex -DBPLEX_USE_MKL=OFF && cmake --build build/bplex
#   cmake --install build/bplex --prefix build/bplex-install
# Neither project's build reads the other's sources.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(CVRPBCP_SANITIZE "Build with ASan/UBSan" OFF)

find_package(Threads REQUIRED)
find_package(OpenMP)

# ---- bplex (installed CMake package: target bplex::bplex) ----
find_package(bplex 0.2 CONFIG REQUIRED)

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
if(OpenMP_CXX_FOUND)
  target_link_libraries(cvrpbcp PRIVATE OpenMP::OpenMP_CXX)
endif()
option(CVRPBCP_NATIVE "compile cvrpbcp for the build machine's CPU (-march=native)" ON)
if(CVRPBCP_NATIVE)
  target_compile_options(cvrpbcp PRIVATE -march=native)
endif()
if(CVRPBCP_SANITIZE)
  target_compile_options(cvrpbcp PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(cvrpbcp PUBLIC -fsanitize=address,undefined)
endif()

# ---- executable ----
add_executable(cvrpbcp_cli src/app/main.cpp)

# LP trace replay (bplex backend; see tools/lpreplay.c for the CPLEX build)
add_executable(lpreplay tools/lpreplay.c)
target_link_libraries(lpreplay PRIVATE bplex::bplex m)
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
