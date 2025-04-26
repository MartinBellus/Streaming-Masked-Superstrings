add_subdirectory(bloom_filter)
add_subdirectory(math)
add_subdirectory(helper)
add_subdirectory(io)

add_executable(${MAIN} main.cpp)
target_link_libraries(${MAIN} io bloom_filter)
