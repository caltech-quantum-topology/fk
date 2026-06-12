cmake_minimum_required(VERSION 3.20)
project(fkcompute_cpp LANGUAGES CXX)

message(STATUS "Wheel platlib dir: ${SKBUILD_PLATLIB_DIR}")

# ---- Clean previous build ----
execute_process(
    COMMAND make clean
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cpp
    RESULT_VARIABLE CLEAN_RESULT
    OUTPUT_QUIET
    ERROR_QUIET
)
if(CLEAN_RESULT EQUAL 0)
    message(STATUS "Successfully ran make clean")
else()
    message(WARNING "make clean failed or not available")
endif()

# ---- Build fk_main via Makefile ----
execute_process(
    COMMAND make main
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cpp
    RESULT_VARIABLE BUILD_RESULT
    OUTPUT_VARIABLE BUILD_OUTPUT
    ERROR_VARIABLE BUILD_ERROR
)
if(NOT BUILD_RESULT EQUAL 0)
    message(FATAL_ERROR "make main failed:\n${BUILD_OUTPUT}\n${BUILD_ERROR}")
else()
    message(STATUS "Successfully built fk_main via make main")
endif()

# ---- Copy binary to local _bin for development ----
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/fkcompute/_bin")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/cpp/fk_main"
     DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/src/fkcompute/_bin/")

# ---- Install pre-built binary into the wheel ----
install(PROGRAMS
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/fk_main"
    DESTINATION "${SKBUILD_PLATLIB_DIR}/fkcompute/_bin"
)
