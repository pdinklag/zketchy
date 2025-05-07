add_library(zk-cscore INTERFACE)
target_link_libraries(zk-cscore INTERFACE fp iopp pm-malloc)

add_library(zk-psamplz INTERFACE)
target_link_libraries(zk-psamplz INTERFACE ${TBB_LIBRARIES} fp growt iopp word-packing pm-malloc)

add_subdirectory(tools)
