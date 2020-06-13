# Copyright (c) 2017, Mate Soos
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

include_directories(
    ${PROJECT_SOURCE_DIR}
)

if (NOT WIN32)
    add_cxx_flag_if_supported("-Wno-bitfield-constant-conversion")
    #add_cxx_flag_if_supported("-Wduplicated-cond")
    #add_cxx_flag_if_supported("-Wduplicated-branches")
    add_cxx_flag_if_supported("-Wlogical-op")
    add_cxx_flag_if_supported("-Wrestrict")
    add_cxx_flag_if_supported("-Wnull-dereference")
    add_cxx_flag_if_supported("-Wjump-misses-init")
    add_cxx_flag_if_supported("-Wdouble-promotion")
    add_cxx_flag_if_supported("-Wshadow")
    add_cxx_flag_if_supported("-Wformat=2")
    add_cxx_flag_if_supported("-Wextra-semi")
    add_cxx_flag_if_supported("-pedantic")
    #add_cxx_flag_if_supported("-Wdeprecated")
endif()

macro(unigen_add_public_header LIBTARGET HEADER)
    get_target_property(EXISTING_PUBLIC_HEADERS ${LIBTARGET} PUBLIC_HEADER)
    if(EXISTING_PUBLIC_HEADERS)
        list(APPEND EXISTING_PUBLIC_HEADERS "${HEADER}")
    else()
        # Do not append to empty list
        set(EXISTING_PUBLIC_HEADERS "${HEADER}")
    endif()
    set_target_properties(
        ${LIBTARGET}
        PROPERTIES
        PUBLIC_HEADER "${EXISTING_PUBLIC_HEADERS}"
     )
endmacro()

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS})
include_directories(${APPROXMC_INCLUDE_DIRS})
include_directories(${GMP_INCLUDE_DIRS})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

set(unigen_lib_files
    unigen.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp
)


add_library(unigen
    ${unigen_lib_files}
)

set_target_properties(unigen PROPERTIES
    PUBLIC_HEADER "${unigen_public_headers}"
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

install(TARGETS unigen
    EXPORT ${APPROXMC_EXPORT_NAME}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/unigen"
)

unigen_add_public_header(unigen ${CMAKE_CURRENT_SOURCE_DIR}/unigen.h )
unigen_add_public_header(unigen ${CMAKE_CURRENT_SOURCE_DIR}/unigenconfig.h )

#     set_target_properties(
#       unigen PROPERTIES LINK_FLAGS "-s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"ccall\", \"cwrap\"]'  -s LINKABLE=1 -s EXPORT_ALL=1"
#     )
#
#
#     SET_SOURCE_FILES_PROPERTIES(main_min.cpp
#         PROPERTIES COMPILE_FLAGS "-s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"ccall\", \"cwrap\"]'  -s LINKABLE=1 -s EXPORT_ALL=1"
#     )


# -----------------------------------------------------------------------------
# Copy public headers into build directory include directory.
# The unigenConfig.cmake we generate in the build directory depends on
# this.
# -----------------------------------------------------------------------------
set(HEADER_DEST "${PROJECT_BINARY_DIR}/include/unigen")
add_custom_target(CopyPublicHeaders ALL)
get_target_property(unigen_public_headers unigen PUBLIC_HEADER)
foreach(public_header ${unigen_public_headers})
    get_filename_component(HEADER_NAME ${public_header} NAME)
    add_custom_command(TARGET CopyPublicHeaders PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${HEADER_DEST}"
        COMMAND ${CMAKE_COMMAND} -E echo "Copying ${HEADER_NAME} to ${HEADER_DEST}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${public_header} "${HEADER_DEST}"
    )
endforeach()

add_executable(unigen-bin
    main.cpp
    ${unigen_lib_files}
)

set(unigen_exec_link_libs
    ${Boost_LIBRARIES}
    ${GMP_LIBRARY}
)

IF (ZLIB_FOUND)
    SET(unigen_exec_link_libs ${unigen_exec_link_libs} ${ZLIB_LIBRARY})
ENDIF()

target_link_libraries(unigen-bin
    ${CRYPTOMINISAT5_LIBRARIES}
    ${APPROXMC_LIBRARIES}
    ${unigen_exec_link_libs}
    unigen
)

set_target_properties(unigen-bin PROPERTIES
        OUTPUT_NAME unigen
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS unigen-bin
    EXPORT ${unigen_EXPORT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/unigen"

)
SET(CPACK_PACKAGE_EXECUTABLES "unigen")
