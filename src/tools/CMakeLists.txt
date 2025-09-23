add_executable(cscore-tool cscore_tool.cpp)
target_link_libraries(cscore-tool zk-cscore oocmd)
set_target_properties(cscore-tool PROPERTIES OUTPUT_NAME "cscore")

add_executable(lz77-tool lz77_tool.cpp)
target_link_libraries(lz77-tool code iopp lz77 pm-malloc oocmd)
set_target_properties(lz77-tool PROPERTIES OUTPUT_NAME "lz77")

add_executable(psamplz-tool psamplz_tool.cpp)
target_link_libraries(psamplz-tool zk-psamplz oocmd)
set_target_properties(psamplz-tool PROPERTIES OUTPUT_NAME "psamplz")

add_executable(slz77-tool slz77_tool.cpp)
target_link_libraries(slz77-tool ${TBB_LIBRARIES} code fp iopp lz77 pm-malloc oocmd unordered_dense)
set_target_properties(slz77-tool PROPERTIES OUTPUT_NAME "slz77")

add_executable(rearrange-tool rearrange_tool.cpp)
target_link_libraries(rearrange-tool fp iopp pm-malloc oocmd unordered_dense)
set_target_properties(rearrange-tool PROPERTIES OUTPUT_NAME "rearrange")

add_executable(rlz-tool rlz_tool.cpp)
target_link_libraries(rlz-tool cmdline code fp iopp lz77 pm-malloc unordered_dense)
set_target_properties(rlz-tool PROPERTIES OUTPUT_NAME "rlz")

add_executable(topk-lz77-tool topk_lz77_tool.cpp)
target_link_libraries(topk-lz77-tool zk-topk-lz77 oocmd)
set_target_properties(topk-lz77-tool PROPERTIES OUTPUT_NAME "topk-lz77")
