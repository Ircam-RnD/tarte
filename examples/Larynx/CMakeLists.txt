cmake_minimum_required(VERSION 3.14)
project(tarte-add LANGUAGES CXX)

include("../../cmake/utils.cmake")
string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" is_top_level)

if(is_top_level)
    find_package(tarte REQUIRED)
endif()

set(sources main.cpp)
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${sources})

add_executable(tarte-larynx)
target_sources(tarte-larynx PRIVATE ${sources})
target_link_libraries(tarte-larynx PRIVATE tarte::tarte)

if(NOT is_top_level)
    win_copy_deps_to_target_dir(tarte-larynx tarte::tarte)
endif()