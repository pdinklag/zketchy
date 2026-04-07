cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)
project(zketchy)

# options
option(BENCHMARK "Compile performance benchmarking and statistics output" OFF)

# include CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmakemodules")

# OpenMP
find_package(OpenMP REQUIRED)

# TBB
find_package(TBB REQUIRED)

# set C++ build flags
set(CXX_STANDARD c++23)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=gnu++23 ${GCC_WARNINGS} ${OpenMP_CXX_FLAGS} -Wstringop-overflow=0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb -march=native")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

if(${BENCHMARK})
    add_compile_definitions(BENCHMARK)
endif()

# dependencies
include(cmake/CPM.cmake)
CPMAddPackage("gh:martinus/unordered_dense#v4.8.1")
CPMAddPackage(
    URI "gh:TooBiased/growt#master"
    DOWNLOAD_ONLY YES)

CPMAddPackage("gh:pdinklag/cmdline#v1.0.0")
CPMAddPackage("gh:pdinklag/code#v1.0.0")
CPMAddPackage("gh:pdinklag/fp#v1.0.0")
CPMAddPackage("gh:pdinklag/iopp#v1.0.0")
CPMAddPackage("gh:pdinklag/pm#v1.0.0")
CPMAddPackage("gh:pdinklag/word-packing#v1.0.0")

CPMAddPackage("gh:pdinklag/alz#experimental")
CPMAddPackage("gh:pdinklag/lz77#v1.0.0")

add_library(growt INTERFACE)
target_include_directories(growt INTERFACE ${growt_SOURCE_DIR})

# include
include_directories(${CMAKE_SOURCE_DIR}/include)

# subdirectories
add_subdirectory(src)
