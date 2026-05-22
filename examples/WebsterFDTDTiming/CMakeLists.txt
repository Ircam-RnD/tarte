cmake_minimum_required(VERSION 3.14)
project(tarte-add LANGUAGES CXX)

include("../../cmake/utils.cmake")
string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" is_top_level)

if(is_top_level)
    find_package(tarte REQUIRED)
endif()

set(sources main.cpp)
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${sources})

add_executable(tarte-websterFDTDTiming)
target_sources(tarte-websterFDTDTiming PRIVATE ${sources})
target_link_libraries(tarte-websterFDTDTiming PRIVATE tarte::tarte)

if(NOT is_top_level)
    win_copy_deps_to_target_dir(tarte-websterFDTDTiming tarte::tarte)
endif()