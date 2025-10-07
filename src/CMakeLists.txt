add_library(zk-cscore INTERFACE)
target_link_libraries(zk-cscore INTERFACE fp iopp unordered_dense pm-malloc)

add_library(zk-psamplz INTERFACE)
target_link_libraries(zk-psamplz INTERFACE ${TBB_LIBRARIES} fp growt iopp lz77 word-packing pm-malloc)

add_library(zk-topk-lz77 INTERFACE)
target_link_libraries(zk-topk-lz77 INTERFACE ${TBB_LIBRARIES} code fp iopp lz77 unordered_dense pm-malloc)

add_library(zk-topk-lz78 INTERFACE)
target_link_libraries(zk-topk-lz78 INTERFACE ${TBB_LIBRARIES} code iopp pm-malloc)

add_subdirectory(tools)

add_executable(zk zk.cpp)
target_link_libraries(zk zk-cscore zk-psamplz zk-topk-lz77 cmdline)
