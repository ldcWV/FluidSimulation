cmake_minimum_required(VERSION 3.17.0)
project(Sim)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(CUDA_FEATS "Set to On to use CUDA features" ON)
if (CUDA_FEATS)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

set(SOURCES
    src/main.cpp
    src/Scene.cpp
    src/SequentialSimulator.cpp
    src/Renderer.cpp
    src/ParallelSimulator.cu
    src/Shader.cpp
    src/Camera.cpp
    src/ReplayManagers.cpp
    src/Timer.cpp
)

find_package( OpenGL REQUIRED )

include_directories( ${OPENGL_INCLUDE_DIRS} )

set( GLFW_BUILD_DOCS OFF CACHE BOOL  "GLFW lib only" )
set( GLFW_INSTALL OFF CACHE BOOL  "GLFW lib only" )

add_subdirectory( glfw )

add_subdirectory( glad )

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries( ${PROJECT_NAME} ${OPENGL_LIBRARIES} glfw glad )

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/lib
)

target_compile_definitions(${PROJECT_NAME} PUBLIC SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}/src/")
target_compile_definitions(${PROJECT_NAME} PUBLIC SCENE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/scenes/")
target_compile_definitions(${PROJECT_NAME} PUBLIC REPLAY_DIR="${CMAKE_CURRENT_SOURCE_DIR}/replays/")
