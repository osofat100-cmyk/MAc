cmake_minimum_required(VERSION 3.20)
project(tokamak_sim_cpp LANGUAGES CXX)

option(SIM_ENABLE_OPENMP "Enable OpenMP support" ON)
option(SIM_BUILD_TESTS "Build unit tests" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(tokamak_sim_core STATIC)
target_include_directories(
  tokamak_sim_core
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

if(SIM_ENABLE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(tokamak_sim_core PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

find_package(Eigen3 CONFIG QUIET)
if(Eigen3_FOUND)
  target_link_libraries(tokamak_sim_core PUBLIC Eigen3::Eigen)
  target_compile_definitions(tokamak_sim_core PUBLIC SIM_HAVE_EIGEN)
else()
  message(STATUS "Eigen3 not found - expect to add vendor copy or provide include path later.")
endif()

target_sources(
  tokamak_sim_core
  PRIVATE
    src/field_solver.cpp
    src/particle.cpp
    src/simulation.cpp
)

if(SIM_BUILD_TESTS)
  enable_testing()
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
  add_subdirectory(tests)
endif()
