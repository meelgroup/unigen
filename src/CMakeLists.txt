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
    add_cxx_flag_if_supported("-Wdouble-promotion")
    add_cxx_flag_if_supported("-Wshadow")
    add_cxx_flag_if_supported("-Wformat=2")
    add_cxx_flag_if_supported("-Wextra-semi")
    add_cxx_flag_if_supported("-pedantic")
    #add_cxx_flag_if_supported("-Wdeprecated")
endif()

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS})
include_directories(${APPROXMC_INCLUDE_DIRS})
include_directories(${GMP_INCLUDE_DIRS})
include_directories(${ARJUN_INCLUDE_DIRS})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

set(unigen_lib_files
    sampler.cpp
    unigen.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp
)

add_library(unigen ${unigen_lib_files})

target_include_directories(unigen PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

set_target_properties(unigen PROPERTIES
    PUBLIC_HEADER "${unigen_public_headers}"
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

target_link_libraries(unigen
    LINK_PUBLIC ${CRYPTOMINISAT5_LIBRARIES}
    LINK_PUBLIC ${ARJUN_LIBRARIES}
    LINK_PUBLIC ${APPROXMC_LIBRARIES}
    LINK_PUBLIC ${MPFR_LIBRARIES}
    LINK_PUBLIC ${GMP_LIBRARY}
)

set_target_properties(unigen PROPERTIES
  PUBLIC_HEADER "${unigen_public_headers}"
  VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
  SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  INSTALL_RPATH_USE_LINK_PATH TRUE
)

set(unigen_exec_link_libs "")
IF (ZLIB_FOUND)
    SET(unigen_exec_link_libs ${unigen_exec_link_libs} ${ZLIB_LIBRARY})
ENDIF()
add_executable(unigen-bin main.cpp)
target_link_libraries(unigen-bin ${unigen_exec_link_libs})

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
    set_target_properties(unigen-bin PROPERTIES
      OUTPUT_NAME unigen
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
      INSTALL_RPATH_USE_LINK_PATH TRUE
      LINK_FLAGS "-s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s EXPORTED_RUNTIME_METHODS='[\"callMain\", \"ccall\", \"cwrap\", \"FS\", \"print\"]' -s FORCE_FILESYSTEM=1"
    )
else()
    set_target_properties(unigen-bin PROPERTIES
      OUTPUT_NAME unigen
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
      INSTALL_RPATH_USE_LINK_PATH TRUE
    )
endif()

unigen_add_public_header(unigen ${CMAKE_CURRENT_SOURCE_DIR}/unigen.h )

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
