add_library(zk-cscore INTERFACE)
target_link_libraries(zk-cscore INTERFACE fp iopp unordered_dense pm-malloc)

if(${PSAMPLZ})
    add_library(zk-psamplz INTERFACE)
    target_link_libraries(zk-psamplz INTERFACE fp growt iopp lz77 unordered_dense word-packing pm-malloc tbb)
endif()

add_library(zk-topk-lz77 INTERFACE)
target_link_libraries(zk-topk-lz77 INTERFACE alz code fp iopp unordered_dense pm-malloc tbb)

add_library(zk-topk-lz78 INTERFACE)
target_link_libraries(zk-topk-lz78 INTERFACE code iopp pm-malloc tbb)

add_subdirectory(tools)

add_executable(zk zk.cpp)
target_link_libraries(zk alz zk-topk-lz77 cmdline)
