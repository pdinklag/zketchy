add_executable(cscore-tool cscore_tool.cpp)
target_link_libraries(cscore-tool zk-cscore)
set_target_properties(cscore-tool PROPERTIES OUTPUT_NAME "cscore")
