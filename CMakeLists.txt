cmake_minimum_required(VERSION 3.20)

project(MaskedSuperstring)

include(cmake/StandardProjectSettings.cmake)

add_subdirectory(src)

enable_testing()

add_subdirectory(tests)
