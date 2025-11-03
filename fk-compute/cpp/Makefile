# Makefile for FK Computation Project
# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -Iinclude
DEBUG_FLAGS = -g -DDEBUG

# Directories
SRC_DIR = src
INCLUDE_DIR = include
TEST_DIR = tests
EXAMPLE_DIR = examples
BUILD_DIR = build

# Create build directory if it doesn't exist
$(shell mkdir -p $(BUILD_DIR))

# Source files
SOURCES = $(SRC_DIR)/fk_segments_links.cpp \
          $(SRC_DIR)/multivariable_polynomial.cpp \
          $(SRC_DIR)/qalg_links.cpp \
          $(SRC_DIR)/linalg.cpp \
          $(SRC_DIR)/parallel_pool.cpp \
          $(SRC_DIR)/string_to_int.cpp \
          $(SRC_DIR)/inequality_solver.cpp \
          $(SRC_DIR)/fk_computation.cpp

# Object files (placed in build directory)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Header files for dependency tracking
HEADERS = $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp \
          $(INCLUDE_DIR)/fk/qalg_links.hpp \
          $(INCLUDE_DIR)/fk/linalg.hpp \
          $(INCLUDE_DIR)/fk/parallel_pool.hpp \
          $(INCLUDE_DIR)/fk/string_to_int.hpp \
          $(INCLUDE_DIR)/fk/bilvector.hpp \
          $(INCLUDE_DIR)/fk/btree.hpp \
          $(INCLUDE_DIR)/fk/inequality_solver.hpp \
          $(INCLUDE_DIR)/fk/fk_computation.hpp

# Main targets
MAIN_TARGET = fk_segments_links
POLYNOMIAL_EXAMPLE_TARGET = polynomial_example
PARALLEL_DEMO_TARGET = parallel_demo
ARITHMETIC_TEST_TARGET = arithmetic_test
EXAMPLE_TARGET = example
INEQUALITY_EXAMPLE_TARGET = inequality_example
FK_COMPUTATION_TARGET = fk_computation_refactored
FK_TEST_TARGET = fk_computation_test
TREFOIL_EXAMPLE_TARGET = trefoil_fk_example
TREFOIL_SIMPLE_TARGET = trefoil_simple_example
TREFOIL_DEMO_TARGET = trefoil_demo
TREFOIL_COMPLETE_TARGET = trefoil_complete_example
PARSER_TEST_TARGET = parser_comparison_test
TEST_TARGET = testing
FK_MAIN_TARGET = fk_main
CHECKPOINT_DEMO_TARGET = checkpoint_demo
CHECKPOINT_RESUME_DEMO_TARGET = checkpoint_resume_demo

# Default target - now builds main.cpp
all: $(FK_MAIN_TARGET) $(EXAMPLE_TARGET)

# Main target is now fk_main (from main.cpp) - explicit rule to override built-in
main: $(FK_MAIN_TARGET)
	@echo "✓ Built main executable: $(FK_MAIN_TARGET)"

# Prevent Make from using built-in rule for main
.PHONY: main

# Main FK computation executable
$(MAIN_TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Polynomial example executable
$(POLYNOMIAL_EXAMPLE_TARGET): $(BUILD_DIR)/polynomial_example.o $(BUILD_DIR)/multivariable_polynomial.o
	$(CXX) $(CXXFLAGS) -o $@ $^


# Parallel demo executable
$(PARALLEL_DEMO_TARGET): $(BUILD_DIR)/parallel_demo.o $(BUILD_DIR)/parallel_pool.o $(BUILD_DIR)/solution_pool_1a_double_links.o
	$(CXX) $(CXXFLAGS) -pthread -o $@ $^

# Arithmetic test executable
$(ARITHMETIC_TEST_TARGET): $(BUILD_DIR)/arithmetic_test.o $(BUILD_DIR)/multivariable_polynomial.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Example executable
$(EXAMPLE_TARGET): $(BUILD_DIR)/example.o $(BUILD_DIR)/multivariable_polynomial.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Inequality example executable
$(INEQUALITY_EXAMPLE_TARGET): $(BUILD_DIR)/inequality_example.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/inequality_solver.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# FK computation refactored executable
$(FK_COMPUTATION_TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# FK computation test executable
$(FK_TEST_TARGET): $(BUILD_DIR)/fk_computation_test.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Trefoil FK example executable
$(TREFOIL_EXAMPLE_TARGET): $(BUILD_DIR)/trefoil_fk_example.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Trefoil simple example executable
$(TREFOIL_SIMPLE_TARGET): $(BUILD_DIR)/trefoil_simple_example.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Trefoil demo executable
$(TREFOIL_DEMO_TARGET): $(BUILD_DIR)/trefoil_demo.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Trefoil complete example executable
$(TREFOIL_COMPLETE_TARGET): $(BUILD_DIR)/trefoil_complete_example.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# FK Main executable (refactored API)
$(FK_MAIN_TARGET): $(BUILD_DIR)/main.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Parser comparison test executable
$(PARSER_TEST_TARGET): $(BUILD_DIR)/parser_comparison_test.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Testing executable
$(TEST_TARGET): $(BUILD_DIR)/testing.o $(filter-out $(BUILD_DIR)/fk_segments_links.o,$(OBJECTS))
	$(CXX) $(CXXFLAGS) -o $@ $^

# Checkpoint demo executable
$(CHECKPOINT_DEMO_TARGET): $(BUILD_DIR)/checkpoint_demo.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Checkpoint resume demo executable
$(CHECKPOINT_RESUME_DEMO_TARGET): $(BUILD_DIR)/checkpoint_resume_demo.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Debug versions
debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: clean $(MAIN_TARGET)

debug-example: CXXFLAGS += $(DEBUG_FLAGS)
debug-example: clean $(EXAMPLE_TARGET)

# Object file compilation rules
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(EXAMPLE_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule for root directory files (like main.cpp)
$(BUILD_DIR)/%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Specific dependencies for object files
$(BUILD_DIR)/fk_segments_links.o: $(SRC_DIR)/fk_segments_links.cpp $(HEADERS)
$(BUILD_DIR)/multivariable_polynomial.o: $(SRC_DIR)/multivariable_polynomial.cpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp $(INCLUDE_DIR)/fk/bilvector.hpp
$(BUILD_DIR)/parallel_pool.o: $(SRC_DIR)/parallel_pool.cpp $(INCLUDE_DIR)/fk/parallel_pool.hpp $(INCLUDE_DIR)/fk/solution_pool_1a_double_links.hpp
$(BUILD_DIR)/qalg_links.o: $(SRC_DIR)/qalg_links.cpp $(INCLUDE_DIR)/fk/qalg_links.hpp $(INCLUDE_DIR)/fk/linalg.hpp $(INCLUDE_DIR)/fk/bilvector.hpp
$(BUILD_DIR)/linalg.o: $(SRC_DIR)/linalg.cpp $(INCLUDE_DIR)/fk/linalg.hpp $(INCLUDE_DIR)/fk/bilvector.hpp
$(BUILD_DIR)/solution_pool_1a_double_links.o: $(SRC_DIR)/solution_pool_1a_double_links.cpp $(INCLUDE_DIR)/fk/solution_pool_1a_double_links.hpp $(INCLUDE_DIR)/fk/btree.hpp
$(BUILD_DIR)/string_to_int.o: $(SRC_DIR)/string_to_int.cpp $(INCLUDE_DIR)/fk/string_to_int.hpp
$(BUILD_DIR)/polynomial_example.o: $(EXAMPLE_DIR)/polynomial_example.cpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp
$(BUILD_DIR)/parallel_demo.o: $(EXAMPLE_DIR)/parallel_demo.cpp $(INCLUDE_DIR)/fk/parallel_pool.hpp $(INCLUDE_DIR)/fk/solution_pool_1a_double_links.hpp
$(BUILD_DIR)/arithmetic_test.o: $(TEST_DIR)/arithmetic_test.cpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp
$(BUILD_DIR)/example.o: $(EXAMPLE_DIR)/example.cpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp
$(BUILD_DIR)/inequality_example.o: $(EXAMPLE_DIR)/inequality_example.cpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp $(INCLUDE_DIR)/fk/inequality_solver.hpp
$(BUILD_DIR)/inequality_solver.o: $(SRC_DIR)/inequality_solver.cpp $(INCLUDE_DIR)/fk/inequality_solver.hpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp
$(BUILD_DIR)/fk_computation.o: $(SRC_DIR)/fk_computation.cpp $(INCLUDE_DIR)/fk/fk_computation.hpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp
$(BUILD_DIR)/fk_computation_test.o: $(EXAMPLE_DIR)/fk_computation_test.cpp $(INCLUDE_DIR)/fk/fk_computation.hpp
$(BUILD_DIR)/trefoil_fk_example.o: $(EXAMPLE_DIR)/trefoil_fk_example.cpp $(INCLUDE_DIR)/fk/fk_computation.hpp
$(BUILD_DIR)/trefoil_simple_example.o: $(EXAMPLE_DIR)/trefoil_simple_example.cpp $(INCLUDE_DIR)/fk/fk_computation.hpp
$(BUILD_DIR)/trefoil_demo.o: $(EXAMPLE_DIR)/trefoil_demo.cpp $(INCLUDE_DIR)/fk/fk_computation.hpp
$(BUILD_DIR)/trefoil_complete_example.o: $(EXAMPLE_DIR)/trefoil_complete_example.cpp $(INCLUDE_DIR)/fk/fk_computation.hpp
$(BUILD_DIR)/parser_comparison_test.o: $(TEST_DIR)/parser_comparison_test.cpp $(INCLUDE_DIR)/fk/fk_computation.hpp $(INCLUDE_DIR)/fk/string_to_int.hpp
$(BUILD_DIR)/main.o: main.cpp $(INCLUDE_DIR)/fk/fk_computation.hpp
$(BUILD_DIR)/testing.o: $(TEST_DIR)/testing.cpp $(HEADERS)
$(BUILD_DIR)/checkpoint_demo.o: checkpoint_demo.cpp $(INCLUDE_DIR)/fk/fk_computation.hpp
$(BUILD_DIR)/checkpoint_resume_demo.o: checkpoint_resume_demo.cpp $(INCLUDE_DIR)/fk/fk_computation.hpp

# Utility targets
clean:
	rm -rf $(BUILD_DIR) $(MAIN_TARGET) $(EXAMPLE_TARGET) $(PARALLEL_DEMO_TARGET) $(ARITHMETIC_TEST_TARGET) $(INEQUALITY_EXAMPLE_TARGET) $(FK_COMPUTATION_TARGET) $(FK_TEST_TARGET) $(TREFOIL_EXAMPLE_TARGET) $(TREFOIL_SIMPLE_TARGET) $(TREFOIL_DEMO_TARGET) $(TREFOIL_COMPLETE_TARGET) $(PARSER_TEST_TARGET) $(TEST_TARGET) $(FK_MAIN_TARGET) $(CHECKPOINT_DEMO_TARGET) $(CHECKPOINT_RESUME_DEMO_TARGET)

# Clean only object files
clean-obj:
	rm -rf $(BUILD_DIR)

# Run targets
run: $(MAIN_TARGET)
	./$(MAIN_TARGET)

run-example: $(EXAMPLE_TARGET)
	./$(EXAMPLE_TARGET)


run-parallel-demo: $(PARALLEL_DEMO_TARGET)
	./$(PARALLEL_DEMO_TARGET)

run-test: $(TEST_TARGET)
	./$(TEST_TARGET)

run-arithmetic: $(ARITHMETIC_TEST_TARGET)
	./$(ARITHMETIC_TEST_TARGET)

# Installation (optional)
install: $(MAIN_TARGET)
	cp $(MAIN_TARGET) /usr/local/bin/

# Formatting (if clang-format is available)
format:
	@if command -v clang-format >/dev/null 2>&1; then \
		find $(SRC_DIR) $(INCLUDE_DIR) $(TEST_DIR) $(EXAMPLE_DIR) -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i; \
		echo "Code formatted with clang-format"; \
	else \
		echo "clang-format not found, skipping formatting"; \
	fi

# Check compilation without linking
check: $(OBJECTS) $(BUILD_DIR)/polynomial_example.o $(BUILD_DIR)/arithmetic_test.o $(BUILD_DIR)/testing.o
	@echo "All source files compile successfully"

# Show project structure
structure:
	@echo "Project Structure:"
	@echo "├── include/fk/     - Header files"
	@echo "├── src/            - Source files"
	@echo "├── tests/          - Test files"
	@echo "├── examples/       - Example files"
	@echo "├── build/          - Object files (generated)"
	@echo "├── docs/           - Documentation"
	@echo "└── Makefile        - Build configuration"

# Show help
help:
	@echo "Available targets:"
	@echo "  all          - Build main executable and example (default)"
	@echo "  fk_segments_links - Build main FK computation executable"
	@echo "  polynomial_example - Build polynomial example"
	@echo "  parallel_demo   - Build parallel processing demonstration"
	@echo "  arithmetic_test - Build arithmetic operations test"
	@echo "  testing      - Build testing executable"
	@echo "  debug        - Build main executable with debug flags"
	@echo "  debug-example - Build example with debug flags"
	@echo "  clean        - Remove all generated files"
	@echo "  clean-obj    - Remove only object files"
	@echo "  run          - Build and run main executable"
	@echo "  run-example  - Build and run polynomial example"
	@echo "  run-parallel-demo - Build and run parallel demo"
	@echo "  run-arithmetic - Build and run arithmetic test"
	@echo "  run-test     - Build and run testing executable"
	@echo "  check        - Check that all files compile"
	@echo "  format       - Format code with clang-format (if available)"
	@echo "  structure    - Show project directory structure"
	@echo "  install      - Install main executable to /usr/local/bin"
	@echo "  help         - Show this help message"

# Phony targets
.PHONY: all debug debug-example clean clean-obj run run-example run-test run-arithmetic install format check structure help
