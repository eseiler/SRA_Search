# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (SRA_search CXX)

set (CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake/")
set (SEQAN_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include")
set (SDSL_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include/sdsl-lite/include")
# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_PATH})
include_directories (${SDSL_INCLUDE_PATH})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (build src/build.cpp
                      src/helper.h)
add_executable (count_single src/count.cpp
                      src/helper.h)
add_executable (count src/count_single.cpp
                      src/helper.h)
add_executable (time  src/time.cpp
                      src/helper.h)
target_link_libraries (build ${SEQAN_LIBRARIES})
target_link_libraries (count ${SEQAN_LIBRARIES})
