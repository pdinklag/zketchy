add_subdirectory(cmdline)
add_subdirectory(code)
add_subdirectory(fp)

add_library(gsaca INTERFACE)
target_include_directories(gsaca INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/gsaca-double-sort)
target_link_libraries(gsaca INTERFACE atomic)

add_subdirectory(iopp)

option(LIBSAIS_USE_OPENMP "Use OpenMP for parallelization" ON)
add_subdirectory(lz77)
target_compile_definitions(lz77 PUBLIC "LIBSAIS_OPENMP")

add_subdirectory(pm)
add_subdirectory(word-packing)

add_library(growt INTERFACE)
target_include_directories(growt INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/growt)

#add_library(libsais ${CMAKE_CURRENT_SOURCE_DIR}/libsais/src/libsais.c ${CMAKE_CURRENT_SOURCE_DIR}/libsais/src/libsais64.c)
#target_include_directories(libsais PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libsais/include)

add_library(unordered_dense INTERFACE)
target_include_directories(unordered_dense INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/unordered_dense/include)
