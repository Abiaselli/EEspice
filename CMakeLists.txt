# ============================================================================
# EEspice — SPICE-compatible circuit simulator
#
# Modern CMake build
#   * Choose a BLAS/LAPACK backend: MKL | AOCL | OpenBLAS | Accelerate
#   * Three build types: Debug | Benchmark | Release
#
# Quick start with presets:
#     cmake --preset mkl-release
#     cmake --build --preset mkl-release
#
# Manual configure:
#     cmake -S . -B build \
#           -DEESPICE_BLAS_BACKEND=MKL \
#           -DCMAKE_BUILD_TYPE=Release
#     cmake --build build -j
# ============================================================================

cmake_minimum_required(VERSION 3.21)

# ---- Custom "Benchmark" build type ---------------------------------------
# Must be set up before project() so multi-config generators see it.
# We override the standard Debug / Release flag sets too, so the build
# semantics exactly match the legacy tasks.json variants.
set(CMAKE_CONFIGURATION_TYPES "Debug;Benchmark;Release"
    CACHE STRING "Supported configurations" FORCE)

project(EEspice
    VERSION     0.1.0
    DESCRIPTION "SPICE-compatible circuit simulator"
    LANGUAGES   CXX C)   # C enabled for the transitive SuiteSparse/OpenMP probe

# Single-config generators: default to Release, restrict the GUI dropdown.
get_property(_is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT _is_multi_config)
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    endif()
    set_property(CACHE CMAKE_BUILD_TYPE
        PROPERTY STRINGS "Debug" "Benchmark" "Release")
endif()

# ---- Per-configuration flag sets -----------------------------------------
# Release   : fully optimized, no debug info                   (-O3 -DARMA_NO_DEBUG)
# Benchmark : optimized + frame pointers and line info for perf
# Debug     : -O0, frame pointers, Armadillo bounds checks on
set(CMAKE_CXX_FLAGS_RELEASE   "-O3"                             CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_BENCHMARK "-O3 -g -fno-omit-frame-pointer"  CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -fno-omit-frame-pointer"  CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_BENCHMARK "" CACHE STRING "" FORCE)
mark_as_advanced(
    CMAKE_CXX_FLAGS_BENCHMARK
    CMAKE_EXE_LINKER_FLAGS_BENCHMARK)

# ---- Project-wide settings -----------------------------------------------
set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Install prefixes (Armadillo, SuiteSparse, SuperLU, OpenBLAS, ...) are NOT
# hardcoded here — point CMake at them via CMAKE_PREFIX_PATH, either on the
# command line, via the CMAKE_PREFIX_PATH environment variable, or through a
# preset's cacheVariables. See CMakePresets.json for this machine's defaults.

# ---- User-facing options -------------------------------------------------
set(EESPICE_BLAS_BACKEND "MKL" CACHE STRING
    "Which BLAS/LAPACK implementation to link against")
set_property(CACHE EESPICE_BLAS_BACKEND
    PROPERTY STRINGS "MKL" "AOCL" "OpenBLAS" "Accelerate")

# Backend-specific hint variables (user overridable on command line).
# MKL: use MKLROOT only if it points at an install that still contains mkl.h.
# Users often carry a stale MKLROOT from a previously-sourced setvars.sh.
# If MKLROOT is missing or invalid, pass -DMKL_ROOT=... or source setvars.sh.
set(_mkl_root_default "")
if(DEFINED ENV{MKLROOT} AND EXISTS "$ENV{MKLROOT}/include/mkl.h")
    set(_mkl_root_default "$ENV{MKLROOT}")
endif()
set(MKL_ROOT        "${_mkl_root_default}" CACHE PATH "Intel oneAPI MKL root directory")
set(AOCL_BLIS_ROOT  ""                     CACHE PATH "AMD AOCL BLIS root (contains include/LP64 and lib/LP64)")
set(AOCL_FLAME_ROOT ""                     CACHE PATH "AMD AOCL libFLAME root (contains include/LP64 and lib/LP64)")

# ---- Target --------------------------------------------------------------
add_executable(eespice src/eespice.cpp)

target_include_directories(eespice PRIVATE "${CMAKE_SOURCE_DIR}/src")

target_compile_options(eespice PRIVATE
    -fdiagnostics-color=always
    -march=native)

# Armadillo is header-only for this project (ARMA_DONT_USE_WRAPPER);
# only the header search path matters. find_path searches CMAKE_PREFIX_PATH
# and the standard system locations.
find_path(ARMADILLO_INCLUDE_DIR NAMES armadillo REQUIRED)
target_include_directories(eespice PRIVATE "${ARMADILLO_INCLUDE_DIR}")

# ARMA_NO_DEBUG disables bounds checks — keep it on for every non-Debug config.
target_compile_definitions(eespice PRIVATE
    ARMA_DONT_USE_WRAPPER
    ARMA_USE_SUPERLU
    $<$<NOT:$<CONFIG:Debug>>:ARMA_NO_DEBUG>)

# ---- Required dependencies -----------------------------------------------
find_package(OpenMP  REQUIRED)
find_package(Threads REQUIRED)
find_package(PkgConfig REQUIRED)

# SuiteSparse ships modern CMake configs. Finding KLU transitively pulls in
# AMD, BTF, COLAMD and SuiteSparse_config, but we list them explicitly for
# readability and to match the legacy link line exactly.
find_package(KLU                CONFIG REQUIRED)
find_package(AMD                CONFIG REQUIRED)
find_package(BTF                CONFIG REQUIRED)
find_package(COLAMD             CONFIG REQUIRED)
find_package(SuiteSparse_config CONFIG REQUIRED)

# SuperLU has a pkg-config file in ~/.local but no CMake config.
pkg_check_modules(SuperLU REQUIRED IMPORTED_TARGET superlu)

target_link_libraries(eespice PRIVATE
    OpenMP::OpenMP_CXX
    Threads::Threads
    SuiteSparse::KLU
    SuiteSparse::AMD
    SuiteSparse::BTF
    SuiteSparse::COLAMD
    SuiteSparse::SuiteSparseConfig
    PkgConfig::SuperLU
    ${CMAKE_DL_LIBS}   # -ldl on Linux, empty on Apple
    m)

# ---- BLAS/LAPACK backend -------------------------------------------------
if(EESPICE_BLAS_BACKEND STREQUAL "MKL")
    # The MKL CMake package is installed with oneAPI.
    set(MKL_INTERFACE_FULL intel_lp64)
    set(MKL_THREADING      sequential)  # matches the legacy static link line
    set(MKL_LINK           static)
    set(MKL_ARCH           intel64)

    set(_mkl_hint)
    if(MKL_ROOT)
        set(_mkl_hint HINTS "${MKL_ROOT}/lib/cmake/mkl")
    endif()
    find_package(MKL CONFIG REQUIRED ${_mkl_hint})

    target_compile_definitions(eespice PRIVATE ARMA_USE_MKL_ALLOC)
    target_link_libraries(eespice PRIVATE MKL::MKL)

elseif(EESPICE_BLAS_BACKEND STREQUAL "AOCL")
    if(NOT AOCL_BLIS_ROOT OR NOT AOCL_FLAME_ROOT)
        message(FATAL_ERROR
            "EESPICE_BLAS_BACKEND=AOCL requires both AOCL_BLIS_ROOT and "
            "AOCL_FLAME_ROOT to be set (pass -DAOCL_BLIS_ROOT=... and "
            "-DAOCL_FLAME_ROOT=... on the command line, or use a preset).")
    endif()

    find_path(AOCL_BLIS_INCLUDE_DIR  blis.h
        HINTS "${AOCL_BLIS_ROOT}/include/LP64" "${AOCL_BLIS_ROOT}/include"
        REQUIRED NO_DEFAULT_PATH)
    find_library(AOCL_BLIS_LIBRARY   NAMES blis
        HINTS "${AOCL_BLIS_ROOT}/lib/LP64"     "${AOCL_BLIS_ROOT}/lib"
        REQUIRED NO_DEFAULT_PATH)

    find_path(AOCL_FLAME_INCLUDE_DIR FLAME.h
        HINTS "${AOCL_FLAME_ROOT}/include/LP64" "${AOCL_FLAME_ROOT}/include"
        REQUIRED NO_DEFAULT_PATH)
    find_library(AOCL_FLAME_LIBRARY  NAMES flame
        HINTS "${AOCL_FLAME_ROOT}/lib/LP64"     "${AOCL_FLAME_ROOT}/lib"
        REQUIRED NO_DEFAULT_PATH)

    target_compile_definitions(eespice PRIVATE ARMA_USE_BLAS ARMA_USE_LAPACK)
    target_include_directories(eespice PRIVATE
        "${AOCL_BLIS_INCLUDE_DIR}"
        "${AOCL_FLAME_INCLUDE_DIR}")
    # libFLAME calls into BLIS, so order matters: flame before blis.
    target_link_libraries(eespice PRIVATE
        "${AOCL_FLAME_LIBRARY}"
        "${AOCL_BLIS_LIBRARY}")

elseif(EESPICE_BLAS_BACKEND STREQUAL "OpenBLAS")
    pkg_check_modules(OpenBLAS REQUIRED IMPORTED_TARGET openblas)
    target_compile_definitions(eespice PRIVATE ARMA_USE_BLAS ARMA_USE_LAPACK)
    target_link_libraries(eespice PRIVATE PkgConfig::OpenBLAS)

elseif(EESPICE_BLAS_BACKEND STREQUAL "Accelerate")
    if(NOT APPLE)
        message(FATAL_ERROR
            "EESPICE_BLAS_BACKEND=Accelerate is only available on Apple platforms.")
    endif()
    find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
    # Opt into the modern LAPACK API (macOS 13.3+); avoids the deprecated CLAPACK.
    target_compile_definitions(eespice PRIVATE
        ARMA_USE_BLAS
        ARMA_USE_LAPACK
        ACCELERATE_NEW_LAPACK)
    target_link_libraries(eespice PRIVATE "${ACCELERATE_FRAMEWORK}")

else()
    message(FATAL_ERROR
        "Unknown EESPICE_BLAS_BACKEND='${EESPICE_BLAS_BACKEND}'. "
        "Expected one of: MKL, AOCL, OpenBLAS, Accelerate.")
endif()

# ---- Configuration summary -----------------------------------------------
message(STATUS "")
message(STATUS "=== EEspice configuration ===")
message(STATUS "  Build type        : ${CMAKE_BUILD_TYPE}")
message(STATUS "  BLAS backend      : ${EESPICE_BLAS_BACKEND}")
message(STATUS "  C++ compiler      : ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "  Output directory  : ${CMAKE_BINARY_DIR}")
message(STATUS "")
