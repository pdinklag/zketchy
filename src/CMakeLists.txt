add_library(zk-cscore INTERFACE)
target_link_libraries(zk-cscore INTERFACE fp iopp pm-malloc)

add_library(zk-psamplz INTERFACE)
target_link_libraries(zk-psamplz INTERFACE ${TBB_LIBRARIES} fp growt iopp word-packing pm-malloc)

add_library(zk-topk-lz77 INTERFACE)
target_link_libraries(zk-topk-lz77 INTERFACE ${TBB_LIBRARIES} code fp iopp lz77 pm-malloc)

add_subdirectory(tools)
