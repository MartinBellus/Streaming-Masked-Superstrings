# Set Release as the default build type
set(default_build_type "Release")

# Set the build type to Release if not already set
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build (Debug, Release, RelWithDebInfo)" FORCE)
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(MAIN streaming-masked-superstrings)
