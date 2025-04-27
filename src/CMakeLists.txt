add_subdirectory(bloom_filter)
add_subdirectory(math)
add_subdirectory(helper)
add_subdirectory(io)
add_subdirectory(algorithm)

add_executable(${MAIN} main.cpp)
target_link_libraries(${MAIN} io algorithm bloom_filter)
