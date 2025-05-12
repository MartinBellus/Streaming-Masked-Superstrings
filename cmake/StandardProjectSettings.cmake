# Set Release as the default build type
set(default_build_type "Release")

# Set the build type to Release if not already set
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build (Debug, Release, RelWithDebInfo)" FORCE)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.0)
        message(FATAL_ERROR "GCC version 14.0 or higher is required for C++23 support.")
    endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
        message(FATAL_ERROR "Clang version 18.0 or higher is required for C++23 support.")
    endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.34)
        message(FATAL_ERROR "MSVC version 19.34 or higher is required for C++23 support.")
    endif()
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(MAIN streaming-masked-superstrings)
