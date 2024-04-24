find_package(Doxygen REQUIRED)
if (DOXYGEN_FOUND)
    doxygen_add_docs(txt include)

    set(DOXYGEN_INPUT ${PROJECT_BINARY_DIR}/Doxyfile.txt)
    set(DOXYGEN_OUTPUT /doc/html)

    add_custom_command(
            OUTPUT ${DOXYGEN_OUTPUT}
            COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
            COMMAND ${CMAKE_COMMAND} -E echo "Done."
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            DEPENDS ${DOXYGEN_INPUT}
    )

    add_custom_target(apidoc ALL DEPENDS ${DOXYGEN_OUTPUT})

endif (DOXYGEN_FOUND)
