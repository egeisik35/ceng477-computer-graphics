cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(OpenGL_GL_PREFERENCE GLVND)

project(CENG477-PlanetRenderer LANGUAGES C CXX)

# ============================= #
#  Some useful Compile Options  #
# ============================= #
add_library(compile_options INTERFACE)
if(UNIX AND NOT APPLE)
    target_compile_options(compile_options INTERFACE
        -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align
        -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion
        -Wnull-dereference -Wdouble-promotion -Wformat=2 -Wimplicit-fallthrough
    )
else()
    target_compile_options(compile_options INTERFACE /permissive-)
endif()

# Your source files
set(SRC_ALL
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utility.h
)

set(CENG_SHADER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/working_dir/shaders)
set(SRC_SHADERS
    ${CENG_SHADER_DIR}/generic.vert
    ${CENG_SHADER_DIR}/debug.frag
)

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm CONFIG REQUIRED)

# Keep only small vendored deps (glad + stb_image)
add_subdirectory(ext/glad)
add_subdirectory(ext/stb_image)

add_executable(PlanetRenderer)
target_sources(PlanetRenderer PRIVATE ${SRC_ALL} ${SRC_SHADERS})

target_link_libraries(PlanetRenderer PRIVATE
    glad
    stb_image
    glfw
    glm::glm
    compile_options
    OpenGL::GL
)

set_target_properties(PlanetRenderer PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/working_dir
)

if(UNIX AND NOT APPLE)
    target_link_libraries(PlanetRenderer PRIVATE X11)
endif()
