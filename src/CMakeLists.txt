add_library(zk-cscore STATIC cscore.cpp)
target_link_libraries(zk-cscore fp iopp oocmd pm-malloc)

add_subdirectory(tools)
