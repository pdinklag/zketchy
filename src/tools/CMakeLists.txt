add_executable(cscore-tool cscore_tool.cpp)
target_link_libraries(cscore-tool zk-cscore oocmd)
set_target_properties(cscore-tool PROPERTIES OUTPUT_NAME "cscore")

add_executable(lz77-tool lz77_tool.cpp)
target_link_libraries(lz77-tool code iopp lz77 pm-malloc oocmd)
set_target_properties(lz77-tool PROPERTIES OUTPUT_NAME "lz77")

add_executable(psamplz-tool psamplz_tool.cpp)
target_link_libraries(psamplz-tool zk-psamplz oocmd)
set_target_properties(psamplz-tool PROPERTIES OUTPUT_NAME "psamplz")

add_executable(topk-lz77-tool topk_lz77_tool.cpp)
target_link_libraries(topk-lz77-tool zk-topk-lz77 oocmd)
set_target_properties(topk-lz77-tool PROPERTIES OUTPUT_NAME "topk-lz77")
