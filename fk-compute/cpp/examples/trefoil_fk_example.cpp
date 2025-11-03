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
        std::cout << "=== Trefoil FK Computation Example ===" << std::endl;
        std::cout << "Using refactored FK computation architecture" << std::endl << std::endl;

        // Create FK computation instance
        fk::FKComputation computation;

        // Progress tracking
        auto start_time = std::chrono::high_resolution_clock::now();

        auto progress_callback = [&start_time](double progress) {
            printProgressBar(progress);

            if (progress == 1.0) {
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
                std::cout << "\nComputation completed in " << duration.count() << " ms" << std::endl;
            }
        };

        std::cout << "Input file: examples/trefoil_ilp.csv" << std::endl;
        std::cout << "Output file: trefoil_result.json" << std::endl << std::endl;

        std::cout << "Starting FK computation..." << std::endl;

        // Run the computation with progress tracking
        computation.compute("examples/trefoil_ilp", "trefoil_result", progress_callback);

        std::cout << std::endl << "=== Computation Results ===" << std::endl;

        // Get the configuration that was used
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
        std::cout << "Total number of terms: " << coeff_map.size() << std::endl;

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

        // Analyze degree distribution
        int max_total_degree = 0;
        for (const auto& [xPowers, qPoly] : coeff_map) {
            int total_degree = 0;
            for (int power : xPowers) {
                total_degree += std::abs(power);
            }
            max_total_degree = std::max(max_total_degree, total_degree);
        }
        std::cout << "Maximum total degree in x variables: " << max_total_degree << std::endl;

        std::cout << std::endl << "=== Output Files ===" << std::endl;
        std::cout << "✓ JSON result: trefoil_result.json" << std::endl;

        // Also create a human-readable text output
        fk::FKResultWriter writer;
        writer.writeToText(result, "trefoil_result");
        std::cout << "✓ Text result: trefoil_result.txt" << std::endl;

        std::cout << std::endl << "=== Success! ===" << std::endl;
        std::cout << "Trefoil FK computation completed successfully using refactored architecture." << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "\nPossible issues:" << std::endl;
        std::cerr << "- Check that examples/trefoil_ilp.csv exists and is readable" << std::endl;
        std::cerr << "- Verify the CSV format matches expected structure" << std::endl;
        std::cerr << "- Ensure output directory is writable" << std::endl;
        return 1;
    }
}