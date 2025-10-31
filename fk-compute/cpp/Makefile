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
          $(SRC_DIR)/solution_pool_1a_double_links.cpp \
          $(SRC_DIR)/string_to_int.cpp

# Object files (placed in build directory)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Header files for dependency tracking
HEADERS = $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp \
          $(INCLUDE_DIR)/fk/qalg_links.hpp \
          $(INCLUDE_DIR)/fk/linalg.hpp \
          $(INCLUDE_DIR)/fk/solution_pool_1a_double_links.hpp \
          $(INCLUDE_DIR)/fk/string_to_int.hpp \
          $(INCLUDE_DIR)/fk/bilvector.hpp \
          $(INCLUDE_DIR)/fk/btree.hpp

# Main targets
MAIN_TARGET = fk_segments_links
EXAMPLE_TARGET = polynomial_example
ARITHMETIC_TEST_TARGET = arithmetic_test
TEST_TARGET = testing

# Default target
all: $(MAIN_TARGET) $(EXAMPLE_TARGET)

# Main FK computation executable
$(MAIN_TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Polynomial example executable
$(EXAMPLE_TARGET): $(BUILD_DIR)/polynomial_example.o $(BUILD_DIR)/multivariable_polynomial.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Arithmetic test executable
$(ARITHMETIC_TEST_TARGET): $(BUILD_DIR)/arithmetic_test.o $(BUILD_DIR)/multivariable_polynomial.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Testing executable
$(TEST_TARGET): $(BUILD_DIR)/testing.o $(filter-out $(BUILD_DIR)/fk_segments_links.o,$(OBJECTS))
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

# Specific dependencies for object files
$(BUILD_DIR)/fk_segments_links.o: $(SRC_DIR)/fk_segments_links.cpp $(HEADERS)
$(BUILD_DIR)/multivariable_polynomial.o: $(SRC_DIR)/multivariable_polynomial.cpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp $(INCLUDE_DIR)/fk/bilvector.hpp
$(BUILD_DIR)/qalg_links.o: $(SRC_DIR)/qalg_links.cpp $(INCLUDE_DIR)/fk/qalg_links.hpp $(INCLUDE_DIR)/fk/linalg.hpp $(INCLUDE_DIR)/fk/bilvector.hpp
$(BUILD_DIR)/linalg.o: $(SRC_DIR)/linalg.cpp $(INCLUDE_DIR)/fk/linalg.hpp $(INCLUDE_DIR)/fk/bilvector.hpp
$(BUILD_DIR)/solution_pool_1a_double_links.o: $(SRC_DIR)/solution_pool_1a_double_links.cpp $(INCLUDE_DIR)/fk/solution_pool_1a_double_links.hpp $(INCLUDE_DIR)/fk/btree.hpp
$(BUILD_DIR)/string_to_int.o: $(SRC_DIR)/string_to_int.cpp $(INCLUDE_DIR)/fk/string_to_int.hpp
$(BUILD_DIR)/polynomial_example.o: $(EXAMPLE_DIR)/polynomial_example.cpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp
$(BUILD_DIR)/arithmetic_test.o: $(TEST_DIR)/arithmetic_test.cpp $(INCLUDE_DIR)/fk/multivariable_polynomial.hpp
$(BUILD_DIR)/testing.o: $(TEST_DIR)/testing.cpp $(HEADERS)

# Utility targets
clean:
	rm -rf $(BUILD_DIR) $(MAIN_TARGET) $(EXAMPLE_TARGET) $(ARITHMETIC_TEST_TARGET) $(TEST_TARGET)

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
	@echo "  arithmetic_test - Build arithmetic operations test"
	@echo "  testing      - Build testing executable"
	@echo "  debug        - Build main executable with debug flags"
	@echo "  debug-example - Build example with debug flags"
	@echo "  clean        - Remove all generated files"
	@echo "  clean-obj    - Remove only object files"
	@echo "  run          - Build and run main executable"
	@echo "  run-example  - Build and run polynomial example"
	@echo "  run-arithmetic - Build and run arithmetic test"
	@echo "  run-test     - Build and run testing executable"
	@echo "  check        - Check that all files compile"
	@echo "  format       - Format code with clang-format (if available)"
	@echo "  structure    - Show project directory structure"
	@echo "  install      - Install main executable to /usr/local/bin"
	@echo "  help         - Show this help message"

# Phony targets
.PHONY: all debug debug-example clean clean-obj run run-example run-test run-arithmetic install format check structure help