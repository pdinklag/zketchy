cmake_minimum_required(VERSION 3.28.0 FATAL_ERROR)

# set build type to Release if nothing was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" FORCE)
endif()

project(zketchy)

# options
option(BENCHMARK "Compile performance benchmarking and statistics output" OFF)
option(BUILD_TOOLS "Build extra tools, mainly used for development and diagnostics" OFF)

# include CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmakemodules")

# OpenMP
find_package(OpenMP REQUIRED)

# options
option(LIBSAIS_USE_OPENMP "Use OpenMP for parallelization" ON)

# set C++ build flags
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

if(${BENCHMARK})
    add_compile_definitions(BENCHMARK)
endif()

# dependencies
include(cmake/CPM.cmake)
CPMAddPackage("gh:martinus/unordered_dense#v4.8.1")

CPMAddPackage("gh:pdinklag/cmdline#v1.0.1")
CPMAddPackage("gh:pdinklag/code#v1.0.1")
CPMAddPackage("gh:pdinklag/fp#v1.0.1")
CPMAddPackage("gh:pdinklag/iopp#v1.0.1")
CPMAddPackage("gh:pdinklag/pm#v1.0.1")

CPMAddPackage("gh:pdinklag/lz77#v1.0.4")
CPMAddPackage("gh:pdinklag/alz#experimental")

# include
include_directories(${CMAKE_SOURCE_DIR}/include)

# subdirectories
add_subdirectory(src)
