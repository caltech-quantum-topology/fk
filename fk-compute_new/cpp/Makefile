# Makefile for the fkcompute C++ backend.
#
# Primary target: fk_main (used by the Python package via `make main`).
# Examples and tests build by their base name, e.g. `make trefoil_demo`.

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    CXX = clang++
    HOMEBREW_OPENMP := $(shell brew --prefix libomp 2>/dev/null)
    ifneq ($(HOMEBREW_OPENMP),)
        OPENMP_FLAGS = -Xpreprocessor -fopenmp -I$(HOMEBREW_OPENMP)/include
        OPENMP_LIBS = -L$(HOMEBREW_OPENMP)/lib -lomp
    else
        OPENMP_FLAGS = -DNO_OPENMP
        OPENMP_LIBS =
        $(warning OpenMP not found on macOS. Install with: brew install libomp)
    endif
    OS_FLAGS =
else
    CXX = g++
    OPENMP_FLAGS = -fopenmp
    OPENMP_LIBS = -fopenmp
    OS_FLAGS = -march=native
endif

CXXFLAGS = -std=c++17 -Wall -Wextra -O3 $(OS_FLAGS) -DNDEBUG $(OPENMP_FLAGS) -Iinclude \
           $(shell pkg-config --cflags flint 2>/dev/null)
PROFILE_FLAGS = -g -pg -fno-omit-frame-pointer -O2 -DNDEBUG
FLINT_LIBS = $(shell pkg-config --libs flint 2>/dev/null || echo "-lflint") -lgmp

SRC_DIR = src
EXAMPLE_DIR = examples
TEST_DIR = tests
BUILD_DIR = build

# Core sources always linked into fk_main
CORE_SOURCES = $(SRC_DIR)/fk_computation.cpp \
               $(SRC_DIR)/qalg_links.cpp \
               $(SRC_DIR)/linalg.cpp \
               $(SRC_DIR)/string_to_int.cpp

# FLINT-backed polynomial implementation
POLY_SOURCES = $(SRC_DIR)/fmpoly_class.cpp $(SRC_DIR)/qpolynomial.cpp

LIB_SOURCES = $(CORE_SOURCES) $(POLY_SOURCES)
LIB_OBJECTS = $(LIB_SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

EXAMPLE_NAMES = $(notdir $(basename $(wildcard $(EXAMPLE_DIR)/*.cpp)))
TEST_NAMES = $(notdir $(basename $(wildcard $(TEST_DIR)/*.cpp)))

all: fk_main

main: fk_main
	@echo "Built main executable: fk_main"

fk_main: $(BUILD_DIR)/main.o $(LIB_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS) $(OPENMP_LIBS)

# Examples and tests build by base name (e.g. `make trefoil_demo`)
$(EXAMPLE_NAMES): %: $(BUILD_DIR)/$(EXAMPLE_DIR)/%.o $(LIB_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS) $(OPENMP_LIBS)

$(TEST_NAMES): %: $(BUILD_DIR)/$(TEST_DIR)/%.o $(LIB_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS) $(OPENMP_LIBS)

examples: $(EXAMPLE_NAMES)
tests: $(TEST_NAMES)

# Compilation rules with automatic header dependency tracking (-MMD)
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/main.o: main.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/$(EXAMPLE_DIR)/%.o: $(EXAMPLE_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/$(TEST_DIR)/%.o: $(TEST_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(shell find $(BUILD_DIR) -name '*.d' 2>/dev/null)

profile:
	@$(MAKE) clean
	@$(MAKE) fk_main CXXFLAGS="-std=c++17 -Wall -Wextra $(PROFILE_FLAGS) $(OPENMP_FLAGS) -Iinclude $(shell pkg-config --cflags flint 2>/dev/null)"
	@echo "Built profiling executable: fk_main (use perf/gprof/valgrind)"

clean:
	rm -rf $(BUILD_DIR) fk_main $(EXAMPLE_NAMES) $(TEST_NAMES) gmon.out

format:
	@if command -v clang-format >/dev/null 2>&1; then \
		find $(SRC_DIR) include $(TEST_DIR) $(EXAMPLE_DIR) -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i; \
		echo "Code formatted with clang-format"; \
	else \
		echo "clang-format not found, skipping formatting"; \
	fi

help:
	@echo "Targets:"
	@echo "  fk_main / main  - Build the main FK computation executable (default)"
	@echo "  <name>          - Build any example or test by base name, e.g. 'make trefoil_demo'"
	@echo "  examples, tests - Build all examples / all tests"
	@echo "  profile         - Build fk_main with profiling flags"
	@echo "  clean           - Remove build artifacts and executables"
	@echo "  format          - Run clang-format (if available)"

.PHONY: all main examples tests profile clean format help
