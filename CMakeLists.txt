cmake_minimum_required(VERSION 3.21)
project(simforode
        VERSION 0.0.1
        DESCRIPTION "SIMFOR ODE standalone library"
        LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(SIMFOR_BUILD_EXAMPLES "Build SIMFOR examples" OFF)
option(SIMFOR_BUILD_DOCS "Build SIMFOR documentation" OFF)


find_package(Boost 1.83.0 COMPONENTS mpi REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

include_directories(${Boost_INCLUDE_DIR})

add_library(simforode)

add_compile_options(-O3)



set(SIMFOR_INCLUDES
        # Public API includes
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/odu.hpp
        # Private API includes
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/internal/types.hpp
)

set(SIMFOR_SOURCES
        ${SIMFOR_INCLUDES}

        ${CMAKE_CURRENT_SOURCE_DIR}/src/odu.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${SIMFOR_SOURCES})

target_sources(simforode PRIVATE ${SIMFOR_SOURCES})

target_include_directories(simforode PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/simfor>
        $<INSTALL_INTERFACE:include/simfor>
)

target_include_directories(simforode PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries( simforode boost_mpi MPI::MPI_CXX OpenMP::OpenMP_CXX)


if (SIMFOR_BUILD_EXAMPLES)
    include(cmake/utils.cmake)
    set(EXAMPLES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/examples")
    subdirlist(SUBDIRS ${EXAMPLES_DIR})
    foreach (subdir ${SUBDIRS})
        add_subdirectory("${EXAMPLES_DIR}/${subdir}")
    endforeach ()
endif ()

install(TARGETS simforode
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/simfor/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/simforode
)

install(TARGETS simforode
    EXPORT simforode
)

install(EXPORT simforode
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/simforode
    NAMESPACE simforode::
)

configure_package_config_file(
    cmake/simforodeConfig.cmake.in                  # the template file
    ${CMAKE_CURRENT_BINARY_DIR}/simforodeConfig.cmake   # destination
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/simforode
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/simforodeConfigVersion.cmake
  COMPATIBILITY AnyNewerVersion
)


install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/simforodeConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/simforodeConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/simforode)



if (SIMFOR_BUILD_DOCS)
    include(cmake/doxygen.cmake)
endif ()
