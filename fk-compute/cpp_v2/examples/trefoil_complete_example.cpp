#include <iostream>
#include <iomanip>
#include <chrono>
#include "fk/fk_computation.hpp"

void printProgressBar(double progress, int width = 50) {
    std::cout << "\r[";
    int filled = static_cast<int>(progress * width);
    for (int i = 0; i < width; ++i) {
        if (i < filled) {
            std::cout << "=";
        } else if (i == filled) {
            std::cout << ">";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "] " << std::fixed << std::setprecision(1) << (progress * 100.0) << "%" << std::flush;
}

int main() {
    try {
        std::cout << "=== Complete Trefoil FK Computation ===" << std::endl;
        std::cout << "Using refactored FK computation architecture" << std::endl;
        std::cout << "This should produce identical results to the original implementation" << std::endl << std::endl;

        // Create FK computation instance
        fk::FKComputation computation;

        // Progress tracking
        auto start_time = std::chrono::high_resolution_clock::now();


        std::cout << "Input file: examples/trefoil_ilp.csv" << std::endl;
        std::cout << "Output file: trefoil_refactored_result.json" << std::endl << std::endl;

        std::cout << "Starting complete FK computation with pooling..." << std::endl;

        // Run the complete computation
        computation.compute("examples/trefoil_ilp", "trefoil_refactored_result");

        std::cout << std::endl << "=== Computation Results ===" << std::endl;

        // Get the configuration
        const auto& config = computation.getLastConfiguration();
        std::cout << "Configuration:" << std::endl;
        std::cout << "  Degree: " << config.degree << std::endl;
        std::cout << "  Components: " << config.components << std::endl;
        std::cout << "  Writhe: " << config.writhe << std::endl;
        std::cout << "  Crossings: " << config.crossings << std::endl;
        std::cout << "  Prefactors: " << config.prefactors << std::endl;

        // Get and display the result
        const auto& result = computation.getLastResult();
        std::cout << std::endl << "Result polynomial:" << std::endl;
        result.print(20); // Print first 20 terms

        std::cout << std::endl << "=== Analysis ===" << std::endl;

        // Get coefficient map for analysis
        const auto& coeff_map = result.getCoefficientMap();
        std::cout << "Total number of term groups: " << coeff_map.size() << std::endl;

        // Count non-zero coefficients
        int total_nonzero = 0;
        for (const auto& [xPowers, qPoly] : coeff_map) {
            for (int qPower = qPoly.getMaxNegativeIndex(); qPower <= qPoly.getMaxPositiveIndex(); ++qPower) {
                if (qPoly[qPower] != 0) {
                    total_nonzero++;
                }
            }
        }
        std::cout << "Non-zero coefficients: " << total_nonzero << std::endl;

        std::cout << std::endl << "=== Verification ===" << std::endl;
        std::cout << "Expected result should match trefoil_out.json:" << std::endl;
        std::cout << "- Should have 12 non-zero terms" << std::endl;
        std::cout << "- Terms should include: x^15*q^1, x^15*q^41, x^14*q^36, etc." << std::endl;
        std::cout << "- Coefficients should alternate between 1 and -1" << std::endl;

        if (total_nonzero == 12) {
            std::cout << "✓ Correct number of terms found!" << std::endl;
        } else {
            std::cout << "✗ Expected 12 terms, found " << total_nonzero << std::endl;
        }

        std::cout << std::endl << "=== Output Files ===" << std::endl;
        std::cout << "✓ JSON result: trefoil_refactored_result.json" << std::endl;

        // Also create a human-readable text output
        fk::FKResultWriter writer;
        writer.writeToText(result, "trefoil_refactored_result");
        std::cout << "✓ Text result: trefoil_refactored_result.txt" << std::endl;

        std::cout << std::endl << "Compare with original using:" << std::endl;
        std::cout << "  diff trefoil_out.json trefoil_refactored_result.json" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "\nThis indicates an issue in the refactored computation logic." << std::endl;
        std::cerr << "The refactored version should produce identical results to the original." << std::endl;
        return 1;
    }
}