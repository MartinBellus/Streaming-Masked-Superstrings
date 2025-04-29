add_library(hash poly_hash.cpp murmur_hash.cpp)
target_link_libraries(hash math helper)
