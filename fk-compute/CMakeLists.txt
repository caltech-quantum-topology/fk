cmake_minimum_required(VERSION 3.20)
project(fkcompute_cpp LANGUAGES CXX)

# ---- scikit-build-core gives this so we know where to install in the wheel ----
message(STATUS "Wheel platlib dir: ${SKBUILD_PLATLIB_DIR}")

# ---- Collect sources -----------------------------------------------------------
# Define specific source files for fk_main (following Makefile pattern)
set(FK_MAIN_SOURCES
  "cpp/main.cpp"
  "cpp/src/fk_computation.cpp"
  "cpp/src/multivariable_polynomial.cpp"
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

# ---- Linking (Threads etc.) ----------------------------------------------------
find_package(Threads REQUIRED)
target_link_libraries(fk_main PRIVATE Threads::Threads)

# ---- Installation into the Python wheel ---------------------------------------
include(GNUInstallDirs)
install(TARGETS fk_main
  RUNTIME DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"
  BUNDLE  DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"
  LIBRARY DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"
  ARCHIVE DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"
)
