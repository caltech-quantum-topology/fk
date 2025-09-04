cmake_minimum_required(VERSION 3.20)
project(fkcompute_cpp LANGUAGES CXX)

# Discover Python installation paths from scikit-build-core
# (SKBUILD_* variables are set automatically)
message(STATUS "Wheel platlib dir: ${SKBUILD_PLATLIB_DIR}")

# Adjust if you have multiple executables/libraries.
# Example: one executable called fk_segments_links built from all .cpp files.
file(GLOB FKCOMPUTE_SOURCES CONFIGURE_DEPENDS "cpp/*.cpp")
add_executable(fk_segments_links ${FKCOMPUTE_SOURCES})

# Set your C++ standard & any compile options you need
set_target_properties(fk_segments_links PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED YES
)

# Example: add include dirs if you have headers in cpp/include
# target_include_directories(fk_segments_links PRIVATE cpp/include)

# If you need extra libraries:
# find_package(Threads REQUIRED)
# target_link_libraries(fk_segments_links PRIVATE Threads::Threads)

# Install the executable inside the Python package in the wheel:
# On Windows this will produce fk_segments_links.exe
include(GNUInstallDirs)
install(TARGETS fk_segments_links
  RUNTIME DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"
  BUNDLE  DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"  # macOS app bundles
  LIBRARY DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"  # just in case
  ARCHIVE DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"
)
