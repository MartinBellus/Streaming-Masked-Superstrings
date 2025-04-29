add_executable(hll_test hll_test.cpp)
add_dependencies(hll_test create_data_symlink)

target_link_libraries(hll_test hash io)
