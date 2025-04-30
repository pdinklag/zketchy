cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)
project(zketchy)

# include CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmakemodules")

# OpenMP
find_package(OpenMP REQUIRED)

# TBB
find_package(TBB REQUIRED)

# set C++ build flags
set(CXX_STANDARD c++20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=gnu++20 ${GCC_WARNINGS} ${OpenMP_CXX_FLAGS} -Wstringop-overflow=0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb -march=native")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# init git submodules
include(InitGitSubmodules)

# extlib
add_subdirectory(extlib)

# include
include_directories(${CMAKE_SOURCE_DIR}/include)

# subdirectories
add_subdirectory(src)
