cmake_minimum_required(VERSION 3.30)
project(tarte LANGUAGES CXX)

if(TARTE_USE_HDF5)
    include("../../cmake/utils.cmake")
    string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" is_top_level)

    if(is_top_level)
        find_package(tarte REQUIRED)
    endif()

    set(sources main.cpp)
    source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${sources})

    add_executable(tarte-voiceReadWrite)
    target_sources(tarte-voiceReadWrite PRIVATE ${sources})
    target_link_libraries(tarte-voiceReadWrite PRIVATE tarte::tarte)

    if(NOT is_top_level)
        win_copy_deps_to_target_dir(tarte-voiceReadWrite tarte::tarte)
    endif()
endif()