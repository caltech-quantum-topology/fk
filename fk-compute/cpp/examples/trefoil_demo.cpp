#include <iostream>
#include <iomanip>
#include "fk/fk_computation.hpp"

void displayConfiguration(const fk::FKConfiguration& config) {
    std::cout << "\n=== Trefoil Configuration ===" << std::endl;
    std::cout << "Degree: " << config.degree << std::endl;
    std::cout << "Components: " << config.components << std::endl;
    std::cout << "Writhe: " << config.writhe << std::endl;
    std::cout << "Crossings: " << config.crossings << std::endl;
    std::cout << "Prefactors: " << config.prefactors << std::endl;

    std::cout << "\nClosed strand components: ";
    for (size_t i = 0; i < config.closed_strand_components.size(); ++i) {
        std::cout << config.closed_strand_components[i];
        if (i < config.closed_strand_components.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    std::cout << "Crossing relation types: ";
    for (size_t i = 0; i < config.crossing_relation_types.size(); ++i) {
        std::cout << config.crossing_relation_types[i];
        if (i < config.crossing_relation_types.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    std::cout << "Top crossing components: ";
    for (size_t i = 0; i < config.top_crossing_components.size(); ++i) {
        std::cout << config.top_crossing_components[i];
        if (i < config.top_crossing_components.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    std::cout << "Bottom crossing components: ";
    for (size_t i = 0; i < config.bottom_crossing_components.size(); ++i) {
        std::cout << config.bottom_crossing_components[i];
        if (i < config.bottom_crossing_components.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    if (!config.criteria.empty()) {
        std::cout << "\nCriteria vectors: " << config.criteria.size() << " total" << std::endl;
        for (size_t i = 0; i < std::min(config.criteria.size(), size_t(3)); ++i) {
            std::cout << "  Criteria[" << i << "]: ";
            for (size_t j = 0; j < config.criteria[i].size(); ++j) {
                std::cout << config.criteria[i][j];
                if (j < config.criteria[i].size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
        if (config.criteria.size() > 3) {
            std::cout << "  ... and " << (config.criteria.size() - 3) << " more" << std::endl;
        }
    }

    if (!config.inequalities.empty()) {
        std::cout << "\nInequality vectors: " << config.inequalities.size() << " total" << std::endl;
        for (size_t i = 0; i < std::min(config.inequalities.size(), size_t(3)); ++i) {
            std::cout << "  Inequality[" << i << "]: ";
            for (size_t j = 0; j < config.inequalities[i].size(); ++j) {
                std::cout << config.inequalities[i][j];
                if (j < config.inequalities[i].size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
        if (config.inequalities.size() > 3) {
            std::cout << "  ... and " << (config.inequalities.size() - 3) << " more" << std::endl;
        }
    }
}

int main() {
    try {
        std::cout << "=== Trefoil FK Computation Demo ===" << std::endl;
        std::cout << "Demonstrating refactored FK architecture with trefoil_ilp.csv" << std::endl;

        // Step 1: Parse input file
        std::cout << "\n[1/5] Parsing input file..." << std::endl;
        fk::FKInputParser parser;
        fk::FKConfiguration config = parser.parseFromFile("examples/trefoil_ilp");
        std::cout << "✓ Successfully parsed trefoil_ilp.csv" << std::endl;

        // Display the parsed configuration
        displayConfiguration(config);

        // Step 2: Validate configuration
        std::cout << "\n[2/5] Validating configuration..." << std::endl;
        if (config.isValid()) {
            std::cout << "✓ Configuration is valid" << std::endl;
        } else {
            std::cout << "✗ Configuration validation failed" << std::endl;
            return 1;
        }

        // Step 3: Initialize computation engine
        std::cout << "\n[3/5] Initializing computation engine..." << std::endl;
        fk::FKComputationEngine engine(config);
        std::cout << "✓ Computation engine initialized successfully" << std::endl;

        // Step 4: Perform actual FK computation
        std::cout << "\n[4/5] Performing FK computation..." << std::endl;
        fk::FKComputation computation;

        // Progress callback to show computation progress
        auto progress_callback = [](double progress) {
            std::cout << "\rProgress: " << std::fixed << std::setprecision(1)
                      << (progress * 100.0) << "%" << std::flush;
        };

        computation.compute(config, "trefoil_demo_output", progress_callback);
        std::cout << "\n✓ FK computation completed successfully" << std::endl;

        // Step 5: Display the computed result
        std::cout << "\n[5/5] Displaying computed FK polynomial..." << std::endl;
        const MultivariablePolynomial& result = computation.getLastResult();

        std::cout << "Computed FK polynomial for trefoil knot:" << std::endl;
        result.print(10);

        std::cout << "\n✓ Result written to trefoil_demo_output.json and trefoil_demo_output.txt" << std::endl;

        std::cout << "\n=== FK Computation Complete! ===" << std::endl;
        std::cout << "Successfully computed the FK polynomial for the trefoil knot using the refactored architecture." << std::endl;
        std::cout << "\nRefactored components used:" << std::endl;
        std::cout << "  ✓ FKConfiguration - Clean data structure" << std::endl;
        std::cout << "  ✓ FKInputParser - CSV parsing with error handling" << std::endl;
        std::cout << "  ✓ FKComputationEngine - Core computation logic" << std::endl;
        std::cout << "  ✓ FKResultWriter - Output formatting" << std::endl;
        std::cout << "  ✓ FKComputation - Main orchestrator with progress tracking" << std::endl;

        std::cout << "\nThe computation demonstrates:" << std::endl;
        std::cout << "  • Complete FK polynomial calculation" << std::endl;
        std::cout << "  • Progress tracking during computation" << std::endl;
        std::cout << "  • Proper error handling" << std::endl;
        std::cout << "  • Clean separation of concerns" << std::endl;
        std::cout << "  • Maintainable architecture" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "\nThis demonstrates proper exception handling in the refactored architecture!" << std::endl;
        return 1;
    }
}
