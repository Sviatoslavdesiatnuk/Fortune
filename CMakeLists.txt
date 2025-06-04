cmake_minimum_required(VERSION 3.29)
project(Fortune_algo)

set(CMAKE_CXX_STANDARD 20)

# Include the GL directory for header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(Fortune_algo main.cpp
        voronoi.h
        voronoi.cpp)

# Find OpenGL
find_package(OpenGL REQUIRED)

# Link directly with the local freeglut library
target_link_libraries(Fortune_algo ${OPENGL_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/freeglut.lib)

# Copy the DLL to the output directory post-build
add_custom_command(TARGET Fortune_algo POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/freeglut.dll"
        $<TARGET_FILE_DIR:Fortune_algo>)