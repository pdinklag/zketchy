add_subdirectory(fp)
add_subdirectory(iopp)
add_subdirectory(oocmd)
add_subdirectory(pm)
add_subdirectory(word-packing)

add_library(growt INTERFACE)
target_include_directories(growt INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/growt)
