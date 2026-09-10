add_library(zk-topk-lz77 INTERFACE)
target_link_libraries(zk-topk-lz77 INTERFACE alz code fp iopp unordered_dense pm-malloc tbb)

add_library(zk-topk-lz78 INTERFACE)
target_link_libraries(zk-topk-lz78 INTERFACE code iopp pm-malloc tbb)

if(${BUILD_TOOLS})
    add_subdirectory(tools)
endif()

add_executable(zk zk.cpp)
target_compile_definitions(zk PUBLIC ALZ_CARDINALITY_CALLBACK_ENABLED)
target_link_libraries(zk alz zk-topk-lz77 cmdline)
