# Makefile for FK Computation Project
# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -fopenmp -Iinclude $(shell pkg-config --cflags flint)
DEBUG_FLAGS = -g -DDEBUG
PROFILE_FLAGS = -g -pg -fno-omit-frame-pointer -O2 -DNDEBUG
FLINT_LIBS = $(shell pkg-config --libs flint)
OPENMP_LIBS = -fopenmp

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
          $(SRC_DIR)/fmpoly.cpp \
          $(SRC_DIR)/qpolynomial.cpp \
          $(SRC_DIR)/fmpoly_class.cpp \
          $(SRC_DIR)/bmpoly.cpp \
          $(SRC_DIR)/qalg_links.cpp \
          $(SRC_DIR)/linalg.cpp \
          $(SRC_DIR)/string_to_int.cpp \
          $(SRC_DIR)/inequality_solver.cpp \
          $(SRC_DIR)/fk_computation.cpp

# Object files (placed in build directory)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Header files for dependency tracking
HEADERS = $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp \
          $(INCLUDE_DIR)/fk/fmpoly.hpp \
          $(INCLUDE_DIR)/fk/qpolynomial.hpp \
          $(INCLUDE_DIR)/fk/fmpoly_class.hpp \
          $(INCLUDE_DIR)/fk/bmpoly.hpp \
          $(INCLUDE_DIR)/fk/qalg_links.hpp \
          $(INCLUDE_DIR)/fk/linalg.hpp \
          $(INCLUDE_DIR)/fk/string_to_int.hpp \
          $(INCLUDE_DIR)/fk/bilvector.hpp \
          $(INCLUDE_DIR)/fk/btree.hpp \
          $(INCLUDE_DIR)/fk/inequality_solver.hpp \
          $(INCLUDE_DIR)/fk/fk_computation.hpp

# Main targets
MAIN_TARGET = fk_segments_links
POLYNOMIAL_EXAMPLE_TARGET = polynomial_example
FMPOLY_EXAMPLE_TARGET = fmpoly_example
FMPOLY_COMPATIBILITY_TARGET = fmpoly_compatibility_test
PERFORMANCE_COMPARISON_TARGET = performance_comparison_test
TEST_GETQPOLY_TARGET = test_getQPolynomial
DEBUG_FMPOLY_TARGET = debug_fmpoly
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
BILVECTOR_OPERATIONS = bilvector_operations
NEW_QBINOM = new_qbinomials
XQ_EXMPL = xq_poch_example
RANDOM_FMPOLY_TARGET = random_fmpoly
BMPOLY_TEST_TARGET = bmpoly_test 


# Default target - now builds main.cpp
all: $(FK_MAIN_TARGET) $(EXAMPLE_TARGET)

# Main target is now fk_main (from main.cpp) - explicit rule to override built-in
main: $(FK_MAIN_TARGET)
	@echo "✓ Built main executable: $(FK_MAIN_TARGET)"

# Prevent Make from using built-in rule for main
.PHONY: main

# Main FK computation executable
$(MAIN_TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS) $(OPENMP_LIBS)

# Polynomial example executable
$(POLYNOMIAL_EXAMPLE_TARGET): $(BUILD_DIR)/polynomial_example.o $(BUILD_DIR)/multivariable_polynomial.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS)

# FMPoly example executable
$(FMPOLY_EXAMPLE_TARGET): $(BUILD_DIR)/fmpoly_example.o $(BUILD_DIR)/fmpoly.o $(BUILD_DIR)/qpolynomial.o $(BUILD_DIR)/fmpoly_class.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS)

# FMPoly compatibility test executable
$(FMPOLY_COMPATIBILITY_TARGET): $(BUILD_DIR)/fmpoly_compatibility_test.o $(BUILD_DIR)/fmpoly.o $(BUILD_DIR)/qpolynomial.o $(BUILD_DIR)/fmpoly_class.o $(BUILD_DIR)/multivariable_polynomial.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS)

# Performance comparison test executable
$(PERFORMANCE_COMPARISON_TARGET): $(BUILD_DIR)/performance_comparison_test.o $(BUILD_DIR)/fmpoly.o $(BUILD_DIR)/qpolynomial.o $(BUILD_DIR)/fmpoly_class.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/bmpoly.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS)

# Test getQPolynomial executable
$(TEST_GETQPOLY_TARGET): $(BUILD_DIR)/test_getQPolynomial.o $(BUILD_DIR)/fmpoly.o $(BUILD_DIR)/qpolynomial.o $(BUILD_DIR)/fmpoly_class.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS)

# Debug FMPoly executable
$(DEBUG_FMPOLY_TARGET): $(BUILD_DIR)/debug_fmpoly.o $(BUILD_DIR)/fmpoly.o $(BUILD_DIR)/qpolynomial.o $(BUILD_DIR)/fmpoly_class.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS)

# Arithmetic test executable
$(ARITHMETIC_TEST_TARGET): $(BUILD_DIR)/arithmetic_test.o $(BUILD_DIR)/multivariable_polynomial.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS)

# Example executable
$(EXAMPLE_TARGET): $(BUILD_DIR)/example.o $(BUILD_DIR)/multivariable_polynomial.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS)

# Inequality example executable
$(INEQUALITY_EXAMPLE_TARGET): $(BUILD_DIR)/inequality_example.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/inequality_solver.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# FK computation refactored executable
$(FK_COMPUTATION_TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS) $(OPENMP_LIBS)

# FK computation test executable
$(FK_TEST_TARGET): $(BUILD_DIR)/fk_computation_test.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(OPENMP_LIBS)

# Trefoil FK example executable
$(TREFOIL_EXAMPLE_TARGET): $(BUILD_DIR)/trefoil_fk_example.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(OPENMP_LIBS)

# Trefoil simple example executable
$(TREFOIL_SIMPLE_TARGET): $(BUILD_DIR)/trefoil_simple_example.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(OPENMP_LIBS)

# Trefoil demo executable
$(TREFOIL_DEMO_TARGET): $(BUILD_DIR)/trefoil_demo.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(OPENMP_LIBS)

# Trefoil complete example executable
$(TREFOIL_COMPLETE_TARGET): $(BUILD_DIR)/trefoil_complete_example.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(OPENMP_LIBS)

# FK Main executable (refactored API)
$(FK_MAIN_TARGET): $(BUILD_DIR)/main.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/fmpoly.o $(BUILD_DIR)/qpolynomial.o $(BUILD_DIR)/fmpoly_class.o $(BUILD_DIR)/bmpoly.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS) $(OPENMP_LIBS)

# Parser comparison test executable
$(PARSER_TEST_TARGET): $(BUILD_DIR)/parser_comparison_test.o $(BUILD_DIR)/fk_computation.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/string_to_int.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(OPENMP_LIBS)

# Testing executable
$(TEST_TARGET): $(BUILD_DIR)/testing.o $(filter-out $(BUILD_DIR)/fk_segments_links.o,$(OBJECTS))
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS) $(OPENMP_LIBS)


# Debug versions
debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: clean $(MAIN_TARGET)

debug-example: CXXFLAGS += $(DEBUG_FLAGS)
debug-example: clean $(EXAMPLE_TARGET)

# Profile versions - optimized build with profiling support
profile:
	@$(MAKE) clean
	@mkdir -p $(BUILD_DIR)
	@$(MAKE) $(FK_MAIN_TARGET) CXXFLAGS="-std=c++17 -Wall -Wextra $(PROFILE_FLAGS) -Iinclude"
	@echo "✓ Built profiling executable: $(FK_MAIN_TARGET)"
	@echo "📊 Profiling flags: $(PROFILE_FLAGS)"
	@echo ""
	@echo "Usage instructions:"
	@echo "  1. gprof profiling:   Run './$(FK_MAIN_TARGET) <args>' then 'gprof ./$(FK_MAIN_TARGET) gmon.out'"
	@echo "  2. perf profiling:    Run 'perf record ./$(FK_MAIN_TARGET) <args>' then 'perf report'"
	@echo "  3. Valgrind callgrind: Run 'valgrind --tool=callgrind ./$(FK_MAIN_TARGET) <args>'"
	@echo "  4. Valgrind cachegrind: Run 'valgrind --tool=cachegrind ./$(FK_MAIN_TARGET) <args>'"
	@echo ""

profile-example:
	@$(MAKE) clean
	@mkdir -p $(BUILD_DIR)
	@$(MAKE) $(EXAMPLE_TARGET) CXXFLAGS="-std=c++17 -Wall -Wextra $(PROFILE_FLAGS) -Iinclude"
	@echo "✓ Built profiling example: $(EXAMPLE_TARGET)"

# Object file compilation rules
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(EXAMPLE_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Bilvector / multivariable polynomial operations demo
$(BILVECTOR_OPERATIONS): $(BUILD_DIR)/bilvector_operations.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# New qbinomial tests
$(NEW_QBINOM): $(TEST_DIR)/new_qbinomials.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o
	$(CXX) $(CXXFLAGS) -o $@ $^
	
# XQ poch
$(XQ_EXMPL): $(TEST_DIR)/xq_poch_example.o $(BUILD_DIR)/multivariable_polynomial.o $(BUILD_DIR)/qalg_links.o $(BUILD_DIR)/linalg.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Random FMPoly example
$(RANDOM_FMPOLY_TARGET): $(BUILD_DIR)/random_fmpoly.o $(BUILD_DIR)/fmpoly.o $(BUILD_DIR)/qpolynomial.o $(BUILD_DIR)/fmpoly_class.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FLINT_LIBS)

# BMPoly test
$(BMPOLY_TEST_TARGET): $(BUILD_DIR)/bmpoly_test.o $(BUILD_DIR)/bmpoly.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule for root directory files (like main.cpp)
$(BUILD_DIR)/%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Specific dependencies for object files
$(BUILD_DIR)/fk_segments_links.o: $(SRC_DIR)/fk_segments_links.cpp $(HEADERS)
$(BUILD_DIR)/multivariable_polynomial.o: $(SRC_DIR)/multivariable_polynomial.cpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp $(INCLUDE_DIR)/fk/bilvector.hpp
$(BUILD_DIR)/fmpoly.o: $(SRC_DIR)/fmpoly.cpp $(INCLUDE_DIR)/fk/fmpoly.hpp $(INCLUDE_DIR)/fk/qpolynomial.hpp $(INCLUDE_DIR)/fk/fmpoly_class.hpp
$(BUILD_DIR)/qpolynomial.o: $(SRC_DIR)/qpolynomial.cpp $(INCLUDE_DIR)/fk/qpolynomial.hpp
$(BUILD_DIR)/fmpoly_class.o: $(SRC_DIR)/fmpoly_class.cpp $(INCLUDE_DIR)/fk/fmpoly_class.hpp $(INCLUDE_DIR)/fk/qpolynomial.hpp
$(BUILD_DIR)/bmpoly.o: $(SRC_DIR)/bmpoly.cpp $(INCLUDE_DIR)/fk/bmpoly.hpp $(INCLUDE_DIR)/fk/bilvector.hpp
$(BUILD_DIR)/qalg_links.o: $(SRC_DIR)/qalg_links.cpp $(INCLUDE_DIR)/fk/qalg_links.hpp $(INCLUDE_DIR)/fk/linalg.hpp $(INCLUDE_DIR)/fk/bilvector.hpp
$(BUILD_DIR)/linalg.o: $(SRC_DIR)/linalg.cpp $(INCLUDE_DIR)/fk/linalg.hpp $(INCLUDE_DIR)/fk/bilvector.hpp
$(BUILD_DIR)/solution_pool_1a_double_links.o: $(SRC_DIR)/solution_pool_1a_double_links.cpp $(INCLUDE_DIR)/fk/solution_pool_1a_double_links.hpp $(INCLUDE_DIR)/fk/btree.hpp
$(BUILD_DIR)/string_to_int.o: $(SRC_DIR)/string_to_int.cpp $(INCLUDE_DIR)/fk/string_to_int.hpp
$(BUILD_DIR)/polynomial_example.o: $(EXAMPLE_DIR)/polynomial_example.cpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp
$(BUILD_DIR)/fmpoly_example.o: $(EXAMPLE_DIR)/fmpoly_example.cpp $(INCLUDE_DIR)/fk/fmpoly.hpp
$(BUILD_DIR)/fmpoly_compatibility_test.o: $(TEST_DIR)/fmpoly_compatibility_test.cpp $(INCLUDE_DIR)/fk/fmpoly.hpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp
$(BUILD_DIR)/performance_comparison_test.o: $(TEST_DIR)/performance_comparison_test.cpp $(INCLUDE_DIR)/fk/fmpoly.hpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp $(INCLUDE_DIR)/fk/bmpoly.hpp
$(BUILD_DIR)/test_getQPolynomial.o: $(TEST_DIR)/test_getQPolynomial.cpp $(INCLUDE_DIR)/fk/fmpoly.hpp
$(BUILD_DIR)/debug_fmpoly.o: $(TEST_DIR)/debug_fmpoly.cpp $(INCLUDE_DIR)/fk/fmpoly.hpp
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
$(BUILD_DIR)/random_fmpoly.o: $(EXAMPLE_DIR)/random_fmpoly.cpp $(INCLUDE_DIR)/fk/fmpoly.hpp
$(BUILD_DIR)/bmpoly_test.o: $(TEST_DIR)/bmpoly_test.cpp $(INCLUDE_DIR)/fk/bmpoly.hpp

# Utility targets
clean:
	rm -rf $(BUILD_DIR) $(MAIN_TARGET) $(EXAMPLE_TARGET) $(ARITHMETIC_TEST_TARGET) $(INEQUALITY_EXAMPLE_TARGET) $(FK_COMPUTATION_TARGET) $(FK_TEST_TARGET) $(TREFOIL_EXAMPLE_TARGET) $(TREFOIL_SIMPLE_TARGET) $(TREFOIL_DEMO_TARGET) $(TREFOIL_COMPLETE_TARGET) $(PARSER_TEST_TARGET) $(TEST_TARGET) $(FK_MAIN_TARAGET) $(BILVECTOR_OPERATIONS) $(NEW_QBINOM) $(XQ_EXMPL) $(PERFORMANCE_COMPARISON_TARGET) $(FMPOLY_COMPATIBILITY_TARGET) $(FMPOLY_EXAMPLE_TARGET) $(POLYNOMIAL_EXAMPLE_TARGET) $(RANDOM_FMPOLY_TARGET) $(BMPOLY_TEST_TARGET)

# Clean only object files
clean-obj:
	rm -rf $(BUILD_DIR)

# Run targets
run: $(MAIN_TARGET)
	./$(MAIN_TARGET)

run-example: $(EXAMPLE_TARGET)
	./$(EXAMPLE_TARGET)



run-test: $(TEST_TARGET)
	./$(TEST_TARGET)

run-arithmetic: $(ARITHMETIC_TEST_TARGET)
	./$(ARITHMETIC_TEST_TARGET)

run-random-fmpoly: $(RANDOM_FMPOLY_TARGET)
	./$(RANDOM_FMPOLY_TARGET)

run-bmpoly-test: $(BMPOLY_TEST_TARGET)
	./$(BMPOLY_TEST_TARGET)

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
	@echo "  main         - Build main FK computation executable (fk_main)"
	@echo "  fk_segments_links - Build main FK computation executable"
	@echo "  polynomial_example - Build polynomial example"
	@echo "  arithmetic_test - Build arithmetic operations test"
	@echo "  testing      - Build testing executable"
	@echo "  debug        - Build main executable with debug flags"
	@echo "  debug-example - Build example with debug flags"
	@echo "  profile      - Build main with profiling support (gprof, perf, valgrind)"
	@echo "  profile-example - Build example with profiling support"
	@echo "  clean        - Remove all generated files"
	@echo "  clean-obj    - Remove only object files"
	@echo "  run          - Build and run main executable"
	@echo "  run-example  - Build and run polynomial example"
	@echo "  run-arithmetic - Build and run arithmetic test"
	@echo "  run-test     - Build and run testing executable"
	@echo "  run-bmpoly-test - Build and run BMPoly unit tests"
	@echo "  check        - Check that all files compile"
	@echo "  format       - Format code with clang-format (if available)"
	@echo "  structure    - Show project directory structure"
	@echo "  install      - Install main executable to /usr/local/bin"
	@echo "  help         - Show this help message"

# Phony targets
.PHONY: all main debug debug-example profile profile-example clean clean-obj run run-example run-test run-arithmetic run-bmpoly-test install format check structure help bilvector_operations
