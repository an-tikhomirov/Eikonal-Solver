cmake_minimum_required(VERSION 3.12)
project(eikonal-solver)

# АБСОЛЮТНЫЕ пути к INMOST
set(INMOST_INCLUDE_DIR "/home/andrey/Projects/INM_RAS/inmost-root/INMOST-0.1/Source/Headers")
set(INMOST_LIBRARY "/home/andrey/Projects/INM_RAS/inmost-root/INMOST-build/libinmost.a")

# Проверка существования
if(NOT EXISTS ${INMOST_INCLUDE_DIR}/inmost.h)
    message(FATAL_ERROR "inmost.h not found at: ${INMOST_INCLUDE_DIR}")
endif()

if(NOT EXISTS ${INMOST_LIBRARY})
    message(FATAL_ERROR "INMOST library not found at: ${INMOST_LIBRARY}")
endif()

message(STATUS "INMOST include: ${INMOST_INCLUDE_DIR}")
message(STATUS "INMOST library: ${INMOST_LIBRARY}")


include_directories(include)

add_executable(eikonal-solver 
    src/main.cpp
    src/Mesh.cpp

    src/local_solver_isotropic.cpp
    src/global_solver_isotropic.cpp

    src/local_solver_anisotropic.cpp
    src/global_solver_anisotropic.cpp
)

target_include_directories(eikonal-solver PRIVATE 
    ${INMOST_INCLUDE_DIR} /home/andrey/Projects/INM_RAS/inmost-root/INMOST-build
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(eikonal-solver ${INMOST_LIBRARY})