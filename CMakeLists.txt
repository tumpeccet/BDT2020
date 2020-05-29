cmake_minimum_required(VERSION 3.13)
project(BDT2020)
set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES src/main.cpp)

add_subdirectory(bioparser EXCLUDE_FROM_ALL)
add_subdirectory(spoa EXCLUDE_FROM_ALL)
add_subdirectory(googletest EXCLUDE_FROM_ALL)
add_subdirectory(alignment EXCLUDE_FROM_ALL)
add_subdirectory(KMeans EXCLUDE_FROM_ALL)

add_executable(BDT2020 ${SOURCE_FILES})

target_link_libraries(
        BDT2020 PUBLIC
        bioparser
        spoa
        alignment
        kmeans
        )

target_include_directories(
        BDT2020 PUBLIC
        "${PROJECT_BINARY_DIR}/alignment"
)

target_include_directories(
        BDT2020 PUBLIC
        "${PROJECT_BINARY_DIR}/kmeans"
)

add_executable(alignment_test tests/alignment_test.cpp)

target_link_libraries(alignment_test alignment gtest)
