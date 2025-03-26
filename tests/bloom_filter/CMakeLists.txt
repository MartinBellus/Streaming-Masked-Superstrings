add_executable(benchmark benchmark_test.cpp)
target_link_libraries(benchmark PRIVATE bloom_filter)

target_compile_options(benchmark PRIVATE -O3)
