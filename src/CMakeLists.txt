add_library(zk-cscore INTERFACE)
target_link_libraries(zk-cscore INTERFACE fp iopp oocmd pm-malloc)

add_subdirectory(tools)
