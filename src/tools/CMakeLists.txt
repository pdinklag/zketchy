add_executable(cscore-tool cscore_tool.cpp)
target_link_libraries(cscore-tool zk-cscore oocmd)
set_target_properties(cscore-tool PROPERTIES OUTPUT_NAME "cscore")

add_executable(psamplz-tool psamplz_tool.cpp)
target_link_libraries(psamplz-tool zk-psamplz oocmd)
set_target_properties(psamplz-tool PROPERTIES OUTPUT_NAME "psamplz")
