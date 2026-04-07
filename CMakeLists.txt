cmake_minimum_required(VERSION 3.14)
project(tarte
    VERSION 1.0.0
    DESCRIPTION "Physical audio synthesis engine"
    LANGUAGES CXX)

#----------------------------------------------------------------------------------------------------------------------
# general settings and options
#----------------------------------------------------------------------------------------------------------------------

include(cmake/utils.cmake)
include(GNUInstallDirs)

string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" is_top_level)

# Export compile config
add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON)

# TARTE_SHARED_LIBS option (undefined by default) can be used to force shared/static build
option(TARTE_BUILD_TESTS "Build tarte tests" OFF)
option(TARTE_BUILD_EXAMPLES "Build tarte examples" ON)
option(TARTE_BUILD_DOCS "Build tarte documentation" OFF)
option(TARTE_INSTALL "Generate target for installing tarte" ${is_top_level})
set_if_undefined(TARTE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/tarte" CACHE STRING
    "Install path for tarte package-related CMake files")

if(DEFINED TARTE_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ${TARTE_SHARED_LIBS})
endif()

if(NOT DEFINED CMAKE_BUILD_TYPE AND NOT DEFINED CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set_if_undefined(CMAKE_CXX_VISIBILITY_PRESET hidden)
set_if_undefined(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

add_library(tarte) # initialized below
add_library(tarte::tarte ALIAS tarte)

#----------------------------------------------------------------------------------------------------------------------
# tarte dependencies
#----------------------------------------------------------------------------------------------------------------------
target_compile_features(tarte PUBLIC cxx_std_23)
# Search for your dependencies here
find_package(Eigen3 REQUIRED NO_MODULE)

#----------------------------------------------------------------------------------------------------------------------
# tarte sources
#----------------------------------------------------------------------------------------------------------------------

include(GenerateExportHeader)
set(export_file_name "export_shared.h")

if(NOT BUILD_SHARED_LIBS)
    set(export_file_name "export_static.h")
endif()

generate_export_header(tarte EXPORT_FILE_NAME include/${export_file_name})

file(GLOB_RECURSE sources
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/external/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/utility/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${sources})

#----------------------------------------------------------------------------------------------------------------------
# tarte target
#----------------------------------------------------------------------------------------------------------------------

include(CMakePackageConfigHelpers)

target_sources(tarte PRIVATE ${sources})
target_compile_definitions(tarte PUBLIC "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:TARTE_STATIC_DEFINE>")

target_include_directories(tarte
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>")

set_target_properties(tarte PROPERTIES
    SOVERSION ${PROJECT_VERSION_MAJOR}
    VERSION ${PROJECT_VERSION})

target_link_libraries(tarte PUBLIC Eigen3::Eigen)


if(TARTE_INSTALL AND NOT CMAKE_SKIP_INSTALL_RULES)
    configure_package_config_file(cmake/tarte-config.cmake.in tarte-config.cmake
        INSTALL_DESTINATION "${TARTE_INSTALL_CMAKEDIR}")

    write_basic_package_version_file(tarte-config-version.cmake
        COMPATIBILITY SameMajorVersion)

    install(TARGETS tarte EXPORT tarte_export
        RUNTIME COMPONENT tarte
        LIBRARY COMPONENT tarte NAMELINK_COMPONENT tarte-dev
        ARCHIVE COMPONENT tarte-dev
        INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
    install(DIRECTORY include/
        TYPE INCLUDE
        COMPONENT tarte-dev)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/include/tarte/${export_file_name}"
        COMPONENT tarte-dev
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/tarte")

    set(targets_file "tarte-shared-targets.cmake")

    if(NOT BUILD_SHARED_LIBS)
        set(targets_file "tarte-static-targets.cmake")
    endif()

    install(EXPORT tarte_export
        COMPONENT tarte-dev
        FILE "${targets_file}"
        DESTINATION "${TARTE_INSTALL_CMAKEDIR}"
        NAMESPACE tarte::)

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/tarte-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/tarte-config-version.cmake"
        COMPONENT tarte-dev
        DESTINATION "${TARTE_INSTALL_CMAKEDIR}")

    if(MSVC)
        set(pdb_file "")
        set(pdb_file_destination "")

        if(BUILD_SHARED_LIBS)
            set(pdb_file "$<TARGET_PDB_FILE:tarte>")
            set(pdb_file_destination "${CMAKE_INSTALL_BINDIR}")
        else()
            # TARGET_PDB_FILE does not work for pdb file generated for static library build, determining it manually
            set(pdb_file "$<TARGET_FILE_DIR:tarte>/$<TARGET_FILE_PREFIX:tarte>$<TARGET_FILE_BASE_NAME:tarte>.pdb")
            set(pdb_file_destination "${CMAKE_INSTALL_LIBDIR}")
        endif()

        install(FILES "${pdb_file}"
            COMPONENT tarte-dev
            CONFIGURATIONS Debug RelWithDebInfo
            DESTINATION "${pdb_file_destination}"
            OPTIONAL)
    endif()
endif()

#----------------------------------------------------------------------------------------------------------------------
# other targets
#----------------------------------------------------------------------------------------------------------------------

if(TARTE_BUILD_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif()

if(TARTE_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(TARTE_BUILD_DOCS)
    find_package(Doxygen REQUIRED)
    doxygen_add_docs(docs include)
endif()