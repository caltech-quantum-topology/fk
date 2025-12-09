cmake_minimum_required(VERSION 3.20)
project(fkcompute_cpp LANGUAGES CXX)

# ---- scikit-build-core gives this so we know where to install in the wheel ----
message(STATUS "Wheel platlib dir: ${SKBUILD_PLATLIB_DIR}")

# ---- Clean previous builds using Makefile ----
# Run make clean in the cpp directory to ensure clean build
execute_process(
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/cpp make clean
    RESULT_VARIABLE CLEAN_RESULT
    OUTPUT_QUIET
    ERROR_QUIET
)
if(CLEAN_RESULT EQUAL 0)
    message(STATUS "Successfully ran make clean")
else()
    message(WARNING "make clean failed or not available")
endif()

# ---- Collect sources -----------------------------------------------------------
# Define specific source files for fk_main (following Makefile pattern)
set(FK_MAIN_SOURCES
  "cpp/main.cpp"
  "cpp/src/fk_computation.cpp"
  "cpp/src/multivariable_polynomial.cpp"
  "cpp/src/fmpoly.cpp"
  "cpp/src/qpolynomial.cpp"
  "cpp/src/fmpoly_class.cpp"
  "cpp/src/bmpoly.cpp"
  "cpp/src/hmpoly.cpp"
  "cpp/src/zmpoly.cpp"
  "cpp/src/qalg_links.cpp"
  "cpp/src/linalg.cpp"
  "cpp/src/string_to_int.cpp"
  "cpp/src/inequality_solver.cpp"
)

add_executable(fk_main ${FK_MAIN_SOURCES})

# ---- C++ standard & compile options -------------------------------------------
set_target_properties(fk_main PROPERTIES
  CXX_STANDARD 20        # or 17 if that’s what you use
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
)

# ---- Include directories -------------------------------------------------------
target_include_directories(fk_main PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
)

# ---- Linking (Threads, OpenMP, FLINT, etc.) -----------------------------------
find_package(Threads REQUIRED)

# On macOS, help CMake find Homebrew's libomp if installed
if(APPLE)
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE HOMEBREW_LIBOMP
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(HOMEBREW_LIBOMP)
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP}/include")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP}/include")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "${HOMEBREW_LIBOMP}/lib/libomp.dylib")
        message(STATUS "Using Homebrew's libomp from: ${HOMEBREW_LIBOMP}")
    endif()
endif()

find_package(OpenMP)

# Find FLINT library using pkg-config
find_package(PkgConfig)
if(PkgConfig_FOUND)
    pkg_check_modules(FLINT flint)
endif()

# Detect platform and set appropriate library names
set(BLAS_LIBRARIES "")
set(FORTRAN_LIBRARIES "")

if(APPLE)
    # macOS: Use OpenBLAS from Homebrew if available
    find_library(OPENBLAS_LIB openblas)
    if(OPENBLAS_LIB)
        set(BLAS_LIBRARIES ${OPENBLAS_LIB})
        message(STATUS "Found OpenBLAS: ${OPENBLAS_LIB}")
    else()
        # Fall back to Accelerate framework (built-in macOS BLAS)
        find_library(ACCELERATE_FRAMEWORK Accelerate)
        if(ACCELERATE_FRAMEWORK)
            set(BLAS_LIBRARIES ${ACCELERATE_FRAMEWORK})
            message(STATUS "Using Accelerate framework for BLAS")
        else()
            message(WARNING "No BLAS library found on macOS")
        endif()
    endif()

    # macOS may have gfortran from Homebrew or may not need it
    find_library(GFORTRAN_LIB gfortran)
    if(GFORTRAN_LIB)
        set(FORTRAN_LIBRARIES ${GFORTRAN_LIB})
        message(STATUS "Found gfortran: ${GFORTRAN_LIB}")
    else()
        message(STATUS "gfortran not found - will try to link without it")
    endif()
else()
    # Linux: Use standard BLAS/CBLAS/gfortran
    set(BLAS_LIBRARIES cblas blas)
    set(FORTRAN_LIBRARIES gfortran)
endif()

# Link libraries (FLINT requires BLAS)
set(LINK_LIBS Threads::Threads)

if(OpenMP_CXX_FOUND)
    list(APPEND LINK_LIBS OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found - parallel processing enabled")
else()
    message(WARNING "OpenMP not found - parallel processing disabled")
    if(APPLE)
        message(STATUS "  On macOS, install with: brew install libomp")
    endif()
endif()

if(FLINT_FOUND)
    list(APPEND LINK_LIBS ${FLINT_LIBRARIES})
    target_include_directories(fk_main PRIVATE ${FLINT_INCLUDE_DIRS})
    message(STATUS "FLINT found via pkg-config")
else()
    list(APPEND LINK_LIBS flint)
    message(WARNING "FLINT pkg-config not found - trying -lflint directly")
endif()

if(BLAS_LIBRARIES)
    list(APPEND LINK_LIBS ${BLAS_LIBRARIES})
endif()

if(FORTRAN_LIBRARIES)
    list(APPEND LINK_LIBS ${FORTRAN_LIBRARIES})
endif()

target_link_libraries(fk_main PRIVATE ${LINK_LIBS})

# ---- Installation into the Python wheel ---------------------------------------
include(GNUInstallDirs)
install(TARGETS fk_main
  RUNTIME DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"
  BUNDLE  DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"
  LIBRARY DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"
  ARCHIVE DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"
)
