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


if(NOT WIN32)
    add_cxx_flag_if_supported("-Wno-bitfield-constant-conversion")
    #add_cxx_flag_if_supported("-Wduplicated-cond")
    #add_cxx_flag_if_supported("-Wduplicated-branches")
    add_cxx_flag_if_supported("-Wlogical-op")
    add_cxx_flag_if_supported("-Wrestrict")
    add_cxx_flag_if_supported("-Wnull-dereference")
    add_cxx_flag_if_supported("-Wdouble-promotion")
    add_cxx_flag_if_supported("-Wshadow")
    add_cxx_flag_if_supported("-Wformat=2")
    add_cxx_flag_if_supported("-Wextra-semi")
    add_cxx_flag_if_supported("-pedantic")
    #add_cxx_flag_if_supported("-Wdeprecated")
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

set(unigen_lib_files
    unigen.cpp
    sampler.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp
)
add_library(unigen ${unigen_lib_files})

target_include_directories(unigen
  PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(unigen
  PUBLIC
    approxmc
    arjun
)

add_executable(unigen-bin main.cpp ${unigen_lib_files})

target_include_directories(unigen-bin
  PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(unigen-bin
  PRIVATE
    PkgConfig::GMP
    PkgConfig::MPFR
    arjun
    approxmc
    unigen
)

if(ZLIB_FOUND)
    target_include_directories(unigen-bin PRIVATE ${ZLIB_INCLUDE_DIR})
    target_link_libraries(unigen-bin PRIVATE ${ZLIB_LIBRARY})
endif()

set_target_properties(unigen PROPERTIES
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)
# On MinGW both DLL and EXE import libs share the libunigen.dll.a name,
# so the unigen library and unigen-bin executable (OUTPUT_NAME=unigen)
# collide on the import lib. Rename the library's output on Windows to
# avoid the clash.
if(WIN32)
    set_target_properties(unigen PROPERTIES OUTPUT_NAME unigenwin)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
    set_target_properties(unigen-bin PROPERTIES
        OUTPUT_NAME unigen
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )
    target_link_options(unigen-bin PRIVATE
        "SHELL:-s WASM=1"
        "SHELL:-s ALLOW_MEMORY_GROWTH=1"
        "SHELL:-s EXPORTED_RUNTIME_METHODS='[\"callMain\", \"ccall\", \"cwrap\", \"FS\", \"print\"]'"
        "SHELL:-s FORCE_FILESYSTEM=1"
    )
else()
    set_target_properties(unigen-bin PROPERTIES
        OUTPUT_NAME unigen
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )
endif()

unigen_add_public_header(unigen ${CMAKE_CURRENT_SOURCE_DIR}/unigen.h)

# Ensure headers from sub-builds are copied before we compile against them.
foreach(_dep_copy_target
        CopyPublicHeaders_cryptominisat5
        CopyPublicHeaders_arjun
        CopyPublicHeaders_approxmc
        CopyPublicHeaders)
    if(TARGET ${_dep_copy_target})
        add_dependencies(unigen ${_dep_copy_target})
        add_dependencies(unigen-bin ${_dep_copy_target})
    endif()
endforeach()

# Public headers are fully known now; the CopyPublicHeaders target below reads
# this property, so it must be set before that target is created.
get_target_property(unigen_public_headers unigen PUBLIC_HEADER)

# -----------------------------------------------------------------------------
# Copy public headers into build directory include directory.
# The unigenConfig.cmake we generate in the build directory depends on
# this.
# -----------------------------------------------------------------------------
set(HEADER_DEST "${PROJECT_BINARY_DIR}/include/unigen/")
add_custom_target(CopyPublicHeaders_unigen ALL)
foreach(public_header ${unigen_public_headers})
    get_filename_component(HEADER_NAME ${public_header} NAME)
    add_custom_command(TARGET CopyPublicHeaders_unigen PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${HEADER_DEST}"
        COMMAND ${CMAKE_COMMAND} -E echo "Copying ${HEADER_NAME} to ${HEADER_DEST}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${public_header} "${HEADER_DEST}/"
    )
endforeach()

add_dependencies(unigen CopyPublicHeaders_unigen)
add_dependencies(unigen-bin CopyPublicHeaders_unigen)

install(TARGETS unigen
    EXPORT ${UNIGEN_EXPORT_NAME}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/unigen"
)

install(TARGETS unigen-bin
    EXPORT ${UNIGEN_EXPORT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
set(CPACK_PACKAGE_EXECUTABLES "unigen-bin" "unigen")
