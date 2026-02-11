add_executable(cscore-tool cscore_tool.cpp)
target_link_libraries(cscore-tool zk-cscore cmdline)
set_target_properties(cscore-tool PROPERTIES OUTPUT_NAME "cscore")

add_executable(lz77-tool lz77_tool.cpp)
target_link_libraries(lz77-tool cmdline code iopp lz77 pm-malloc)
set_target_properties(lz77-tool PROPERTIES OUTPUT_NAME "lz77")

add_executable(lz78-tool lz78_tool.cpp)
target_link_libraries(lz78-tool cmdline code iopp pm-malloc)
set_target_properties(lz78-tool PROPERTIES OUTPUT_NAME "lz78")

add_executable(lz78-lru-tool lz78_lru_tool.cpp)
target_link_libraries(lz78-lru-tool cmdline code iopp pm-malloc)
set_target_properties(lz78-lru-tool PROPERTIES OUTPUT_NAME "lz78-lru")

add_executable(psamplz-tool psamplz_tool.cpp)
target_link_libraries(psamplz-tool zk-psamplz cmdline)
set_target_properties(psamplz-tool PROPERTIES OUTPUT_NAME "psamplz")

add_executable(rearrange-tool rearrange_tool.cpp)
target_link_libraries(rearrange-tool cmdline fp iopp pm-malloc unordered_dense)
set_target_properties(rearrange-tool PROPERTIES OUTPUT_NAME "rearrange")

add_executable(rlz-tool rlz_tool.cpp)
target_link_libraries(rlz-tool cmdline code fp iopp lz77 pm-malloc unordered_dense)
set_target_properties(rlz-tool PROPERTIES OUTPUT_NAME "rlz")

add_executable(topk-lz77-tool topk_lz77_tool.cpp)
target_link_libraries(topk-lz77-tool zk-topk-lz77 cmdline)
set_target_properties(topk-lz77-tool PROPERTIES OUTPUT_NAME "topk-lz77")

add_executable(topk-lz78-tool topk_lz78_tool.cpp)
target_link_libraries(topk-lz78-tool zk-topk-lz78 cmdline)
set_target_properties(topk-lz78-tool PROPERTIES OUTPUT_NAME "topk-lz78")

add_executable(topk-substrings-tool topk_substrings_tool.cpp)
target_link_libraries(topk-substrings-tool cmdline fp iopp lz77 pm-malloc)
set_target_properties(topk-substrings-tool PROPERTIES OUTPUT_NAME "topk-substrings")
