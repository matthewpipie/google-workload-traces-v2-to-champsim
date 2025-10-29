cmake_minimum_required(VERSION 3.16)
project(converter)

# Set basic build properties
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O2 -g")

# Define DynamoRIO paths
set(champsim_base "${CMAKE_CURRENT_SOURCE_DIR}/ChampSim")
set(dynamorio_base "${CMAKE_CURRENT_SOURCE_DIR}/dynamorio")
set(dynamorio_dir "${dynamorio_base}/clients/drcachesim")
set(CMAKE_PREFIX_PATH "${dynamorio_base}/build/cmake")

# Find required packages
find_package(DynamoRIO REQUIRED)
find_package(Threads REQUIRED)

# Define preprocessor definitions
add_definitions(-DLINUX -DX86_64)

# Create an object library for main.cpp with strict warnings
add_library(main_obj OBJECT src/main.cpp)
target_compile_options(main_obj PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Werror
)

# Then add executable with the object file and other sources
add_executable(converter
    $<TARGET_OBJECTS:main_obj>
    # ${dynamorio_dir}/tools/common/decode_cache.cpp
    # ${dynamorio_dir}/tracer/raw2trace_shared.cpp
)

set_target_properties(converter PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/"
)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(converter PRIVATE -Wno-deprecated-declarations)
endif()

# Configure DynamoRIO for the target
configure_DynamoRIO_standalone(converter)

# Add include directories
target_include_directories(main_obj SYSTEM PRIVATE
    ${PROJECT_SOURCE_DIR}
    "${dynamorio_base}/build/clients/include"
    # "${dynamorio_base}/build/clients/include/drmemtrace"
    ${dynamorio_dir}/simulator
    ${dynamorio_dir}/common
    ${dynamorio_dir}/reader
    ${dynamorio_dir}/tracer
    ${dynamorio_dir}/scheduler
    ${dynamorio_dir}/tools/common
    ${dynamorio_dir}/tools
    "${dynamorio_base}/core/ir/"
    ${dynamorio_dir}
    "${dynamorio_base}/include"
    "${dynamorio_base}/tools/include"
    "${champsim_base}/inc/"
)

# Function to add all subdirectories in a given base directory to include directories
function(add_all_subdirs_include_dirs BASE_DIR)
  file(GLOB SUB_DIRS RELATIVE "${BASE_DIR}" "${BASE_DIR}/*")
  foreach(subdir ${SUB_DIRS})
    set(FULL_PATH "${BASE_DIR}/${subdir}")
    if(IS_DIRECTORY "${FULL_PATH}")
      message(STATUS "Adding include directory: ${FULL_PATH}")
      include_directories("${FULL_PATH}")
    endif()
  endforeach()
endfunction()

# Add all subdirectories under dynamorio_base/ext/ as include directories
add_all_subdirs_include_dirs("${dynamorio_base}/ext/")

# Specify link directories and libraries
target_link_directories(converter PRIVATE "${dynamorio_base}/build/clients/lib64/release")
target_link_libraries(converter
    Threads::Threads
    drmemtrace
    drmemtrace_analyzer
    drcovlib_static
)
