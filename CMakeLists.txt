cmake_minimum_required(VERSION 3.12)

project(fragmentation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Akantu REQUIRED)

# ============================================================================ #
# Cohesive Zone Model

add_akantu_simulation(
  fragmentation_czm
  code/CZM/main.cc
  code/CZM/generate_model.cc
  code/CZM/model_functions.cc
  code/common/generate_mesh.cc
  code/common/diagnostic_writer.cc
)

target_include_directories(
  fragmentation_czm
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/code/CZM
  ${CMAKE_CURRENT_SOURCE_DIR}/code/common
)

# ============================================================================ #
# Phase Field

add_akantu_simulation(
  fragmentation_pf
  code/PF/main.cc
  code/PF/generate_model.cc
  code/PF/model_functions.cc
  code/common/generate_mesh.cc
  code/common/diagnostic_writer.cc
)

target_include_directories(
  fragmentation_pf
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/code/PF
  ${CMAKE_CURRENT_SOURCE_DIR}/code/common
)

# ============================================================================ #
# Peridynamics

add_executable(
  fragmentation_pd
  code/PD/main.cc
  code/PD/parameters.cc
  code/PD/make_bar.cc
  code/PD/random_field.cc
  code/PD/simulate_bar.cc
  code/PD/bond_health.cc
  code/PD/count_fragments.cc
  code/PD/fragment_size_theory.cc
)

target_include_directories(
  fragmentation_pd
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/code/PD
)