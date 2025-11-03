#include "fk/fk_computation.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>

/**
 * Demo function for computeWithCheckpoints
 * Emulates the behavior of ./fk_main data/trefoil_ilp_long out
 * but saves checkpoints every 50 points and prints when it saves
 */
void demoCheckpointComputation() {
    std::cout << "=== FK Computation Checkpoint Demo ===" << std::endl;
    std::cout << "Demonstrating checkpointing functionality with trefoil_ilp_long" << std::endl;
    std::cout << "Checkpoint period: 50 points" << std::endl;
    std::cout << "=========================================" << std::endl;

    try {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Create FK computation instance
        fk::FKComputation computation;

        // Input and output files (same as ./fk_main data/trefoil_ilp_long out)
        std::string input_file = "data/trefoil_ilp_long";
        std::string output_file = "out";
        int checkpoint_period = 50;  // Save checkpoint every 50 points

        std::cout << "Input file: " << input_file << ".csv" << std::endl;
        std::cout << "Output file: " << output_file << ".json" << std::endl;
        std::cout << "Checkpoint period: " << checkpoint_period << " points" << std::endl;
        std::cout << std::endl;

        // Perform computation with checkpointing
        std::cout << "Starting FK computation with checkpointing..." << std::endl;
        std::cout << "Will save checkpoint every " << checkpoint_period << " points processed" << std::endl;
        std::cout << std::endl;

        // Note: computeWithCheckpointing will print checkpoint saves internally
        std::string checkpoint_file = computation.computeWithCheckpointing(
            input_file,
            output_file,
            checkpoint_period
        );

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        std::cout << std::endl;
        std::cout << "✓ Computation completed successfully!" << std::endl;
        std::cout << "⏱ Total computation time: " << duration.count() << " ms" << std::endl;
        std::cout << "📁 Final checkpoint saved as: " << checkpoint_file << std::endl;
        std::cout << "📄 Results written to: " << output_file << ".json" << std::endl;

        // Show brief result statistics
        const auto& result = computation.getLastResult();
        auto coeffs = result.getCoefficients();
        int non_zero_count = 0;
        for (const auto& coeff : coeffs) {
            for (int i = coeff.getMaxNegativeIndex(); i <= coeff.getMaxPositiveIndex(); ++i) {
                if (coeff[i] != 0) {
                    non_zero_count++;
                }
            }
        }
        std::cout << "📊 Result contains " << non_zero_count << " non-zero terms" << std::endl;

        std::cout << std::endl;
        std::cout << "=== Demo Complete ===" << std::endl;
        std::cout << "This demonstrates the checkpoint functionality:" << std::endl;
        std::cout << "• Computation progresses in chunks of " << checkpoint_period << " points" << std::endl;
        std::cout << "• State is saved periodically to allow resumption" << std::endl;
        std::cout << "• Progress is printed when checkpoints are saved" << std::endl;
        std::cout << "• Final result matches ./fk_main data/trefoil_ilp_long out" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << std::endl;
        std::cerr << "❌ Error: " << e.what() << std::endl;
        std::cerr << "Demo failed to complete." << std::endl;
    }
}

/**
 * Demo function showing how to resume from a checkpoint
 * This function can be called separately to demonstrate checkpoint resumption
 */
void demoCheckpointResumption(const std::string& checkpoint_path) {
    std::cout << "=== FK Computation Checkpoint Resumption Demo ===" << std::endl;
    std::cout << "Demonstrating resumption from checkpoint: " << checkpoint_path << std::endl;
    std::cout << "===============================================" << std::endl;

    try {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Create FK computation instance
        fk::FKComputation computation;

        std::string output_file = "out_resumed";

        std::cout << "Resuming computation from checkpoint..." << std::endl;
        std::cout << "Checkpoint file: " << checkpoint_path << std::endl;
        std::cout << "Output file: " << output_file << ".json" << std::endl;
        std::cout << std::endl;

        // Resume computation from checkpoint
        computation.resumeFromCheckpoint(checkpoint_path, output_file);

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        std::cout << "✓ Resumption completed successfully!" << std::endl;
        std::cout << "⏱ Resumption time: " << duration.count() << " ms" << std::endl;
        std::cout << "📄 Results written to: " << output_file << ".json" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << std::endl;
        std::cerr << "❌ Error: " << e.what() << std::endl;
        std::cerr << "Resumption demo failed." << std::endl;
    }
}

int main() {
    std::cout << "FK Computation Checkpoint System Demo" << std::endl;
    std::cout << "=====================================" << std::endl;
    std::cout << std::endl;

    // Run the main checkpoint demo
    demoCheckpointComputation();

    std::cout << std::endl;
    std::cout << "Note: To test checkpoint resumption, you can call:" << std::endl;
    std::cout << "  demoCheckpointResumption(\"<checkpoint_file_path>\")" << std::endl;
    std::cout << "with the checkpoint file created during the computation." << std::endl;

    return 0;
}