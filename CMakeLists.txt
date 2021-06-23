# CMakeLists files in this project can
# refer to the root source directory of the project as ${VERIFYCURVES_SOURCE_DIR} and
# to the root binary directory of the project as ${VERIFYCURVES_BINARY_DIR}.
cmake_minimum_required (VERSION 3.1)
project (VERIFYCURVES)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(gflags REQUIRED)
find_package(OpenMP REQUIRED)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
add_compile_options(-std=c++11 -fopenmp -O3 -march=native)

include_directories(include)
add_subdirectory(src)
