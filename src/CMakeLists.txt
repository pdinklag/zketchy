add_library(zk-cscore STATIC cscore.cpp)
target_link_libraries(zk-cscore fp iopp oocmd)

add_subdirectory(tools)
