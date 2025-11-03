#include "fk/fk_computation.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <filesystem>
#include <fstream>

/**
 * Demo function for resuming computation from checkpoints
 * Shows how to resume from a saved checkpoint file
 */
void demoCheckpointResumption(const std::string& checkpoint_path) {
    std::cout << "=== FK Computation Checkpoint Resumption Demo ===" << std::endl;
    std::cout << "Demonstrating resumption from checkpoint: " << checkpoint_path << std::endl;
    std::cout << "===============================================" << std::endl;

    try {
        // Check if checkpoint file exists
        if (!std::filesystem::exists(checkpoint_path)) {
            std::cerr << "❌ Error: Checkpoint file does not exist: " << checkpoint_path << std::endl;
            return;
        }

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

        std::cout << std::endl;
        std::cout << "✓ Resumption completed successfully!" << std::endl;
        std::cout << "⏱ Resumption time: " << duration.count() << " ms" << std::endl;
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
        std::cout << "=== Resume Demo Complete ===" << std::endl;
        std::cout << "Successfully resumed computation from checkpoint!" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << std::endl;
        std::cerr << "❌ Error: " << e.what() << std::endl;
        std::cerr << "Resumption demo failed." << std::endl;
    }
}

/**
 * Find the most recent checkpoint file in the current directory
 */
std::string findMostRecentCheckpoint() {
    std::string most_recent_checkpoint;
    std::filesystem::file_time_type most_recent_time{};

    for (const auto& entry : std::filesystem::directory_iterator(".")) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();
            // Look for checkpoint files with pattern fk_checkpoint_*.json
            if (filename.find("fk_checkpoint_") == 0 &&
                filename.size() >= 5 && filename.substr(filename.size() - 5) == ".json") {
                auto file_time = entry.last_write_time();
                if (most_recent_checkpoint.empty() || file_time > most_recent_time) {
                    most_recent_checkpoint = filename;
                    most_recent_time = file_time;
                }
            }
        }
    }

    return most_recent_checkpoint;
}

/**
 * List all available checkpoint files
 */
void listAvailableCheckpoints() {
    std::cout << "Available checkpoint files:" << std::endl;
    bool found_any = false;

    for (const auto& entry : std::filesystem::directory_iterator(".")) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();
            if (filename.find("fk_checkpoint_") == 0 &&
                filename.size() >= 5 && filename.substr(filename.size() - 5) == ".json") {
                std::cout << "  • " << filename << std::endl;
                found_any = true;
            }
        }
    }

    if (!found_any) {
        std::cout << "  (No checkpoint files found)" << std::endl;
    }
}

/**
 * Interactive demo that finds and resumes from the most recent checkpoint
 */
void interactiveDemoResumption() {
    std::cout << "=== Interactive Checkpoint Resume Demo ===" << std::endl;
    std::cout << "This demo will find and resume from the most recent checkpoint" << std::endl;
    std::cout << "=========================================" << std::endl;

    // List available checkpoints
    listAvailableCheckpoints();
    std::cout << std::endl;

    // Find the most recent checkpoint
    std::string most_recent = findMostRecentCheckpoint();

    if (most_recent.empty()) {
        std::cout << "❌ No checkpoint files found!" << std::endl;
        std::cout << "Please run the checkpoint demo first to generate checkpoint files." << std::endl;
        std::cout << "Example: ./checkpoint_demo" << std::endl;
        return;
    }

    std::cout << "📁 Most recent checkpoint: " << most_recent << std::endl;
    std::cout << "Resuming computation from this checkpoint..." << std::endl;
    std::cout << std::endl;

    // Resume from the most recent checkpoint
    demoCheckpointResumption(most_recent);
}

/**
 * Demo that uses an existing partial checkpoint and resumes from it
 * This simulates stopping a computation partway through and resuming later
 */
void demoPartialComputationAndResume() {
    std::cout << "=== Partial Computation and Resume Demo ===" << std::endl;
    std::cout << "This demo will:" << std::endl;
    std::cout << "1. Find an existing checkpoint with remaining points" << std::endl;
    std::cout << "2. Resume computation from that checkpoint to complete the work" << std::endl;
    std::cout << "3. Compare with a fresh complete computation" << std::endl;
    std::cout << "================================================================" << std::endl;

    try {
        // Step 1: Find an intermediate checkpoint with remaining points
        std::cout << "[Step 1] Finding intermediate checkpoint with remaining work..." << std::endl;

        std::string intermediate_checkpoint;
        int remaining_count = 0;

        // Look for any checkpoint with remaining points
        for (const auto& entry : std::filesystem::directory_iterator(".")) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                if (filename.find("fk_checkpoint_") == 0 &&
                    filename.size() >= 5 && filename.substr(filename.size() - 5) == ".json") {

                    // Check if this checkpoint has remaining points
                    std::ifstream check_file(filename);
                    std::string line;
                    int points_found = 0;
                    bool in_remaining_points = false;

                    while (std::getline(check_file, line)) {
                        if (line.find("\"remaining_points\"") != std::string::npos) {
                            in_remaining_points = true;
                            continue;
                        }
                        if (in_remaining_points && line.find("[") != std::string::npos) {
                            // Count actual point entries
                            while (std::getline(check_file, line)) {
                                if (line.find("]") != std::string::npos &&
                                    line.find("\"remaining_points\"") == std::string::npos) {
                                    break;
                                }
                                if (line.find("[") != std::string::npos && line.find("]") != std::string::npos) {
                                    points_found++;
                                }
                            }
                            break;
                        }
                    }
                    check_file.close();

                    if (points_found > 0) {
                        intermediate_checkpoint = filename;
                        remaining_count = points_found;
                        break;  // Use the first one we find
                    }
                }
            }
        }

        if (intermediate_checkpoint.empty()) {
            std::cout << "❌ No checkpoint with remaining points found!" << std::endl;
            std::cout << "Please run the regular checkpoint demo first to create checkpoints." << std::endl;
            return;
        }

        std::cout << "✓ Found checkpoint: " << intermediate_checkpoint << std::endl;
        std::cout << "✓ Remaining points to process: " << remaining_count << std::endl;
        std::cout << std::endl;

        // Step 2: Resume from the intermediate checkpoint
        std::cout << "[Step 2] Resuming computation from checkpoint..." << std::endl;
        std::cout << "This will process the remaining " << remaining_count << " points." << std::endl;
        std::cout << std::endl;

        fk::FKComputation resume_computation;
        resume_computation.resumeFromCheckpoint(intermediate_checkpoint, "partial_resumed_out");

        std::cout << "✓ Partial resume computation completed!" << std::endl;
        std::cout << std::endl;

        // Step 3: Run a fresh complete computation for comparison
        std::cout << "[Step 3] Running fresh complete computation for comparison..." << std::endl;

        fk::FKComputation fresh_computation;
        fresh_computation.compute("data/trefoil_ilp_long", "fresh_complete_out");

        std::cout << "✓ Fresh complete computation finished!" << std::endl;
        std::cout << std::endl;

        // Compare results
        const auto& resume_result = resume_computation.getLastResult();
        const auto& fresh_result = fresh_computation.getLastResult();

        std::cout << "=== Result Comparison ===" << std::endl;
        std::cout << "Resumed computation non-zero terms: ";
        auto coeffs1 = resume_result.getCoefficients();
        int count1 = 0;
        for (const auto& coeff : coeffs1) {
            for (int i = coeff.getMaxNegativeIndex(); i <= coeff.getMaxPositiveIndex(); ++i) {
                if (coeff[i] != 0) count1++;
            }
        }
        std::cout << count1 << std::endl;

        std::cout << "Fresh complete computation non-zero terms: ";
        auto coeffs2 = fresh_result.getCoefficients();
        int count2 = 0;
        for (const auto& coeff : coeffs2) {
            for (int i = coeff.getMaxNegativeIndex(); i <= coeff.getMaxPositiveIndex(); ++i) {
                if (coeff[i] != 0) count2++;
            }
        }
        std::cout << count2 << std::endl;

        if (count1 == count2) {
            std::cout << "✓ Results match! Partial checkpoint and resume functionality working correctly." << std::endl;
            std::cout << "✓ This proves that a computation can be stopped partway through and resumed later." << std::endl;
        } else {
            std::cout << "⚠ Results differ. Resumed=" << count1 << ", Fresh=" << count2 << std::endl;
        }

        // Verify the outputs are identical
        std::cout << "\n=== File Comparison ===" << std::endl;
        int diff_result = std::system("diff partial_resumed_out.json fresh_complete_out.json > /dev/null 2>&1");
        if (diff_result == 0) {
            std::cout << "✓ Output files are identical!" << std::endl;
            std::cout << "✓ Partial resume produces exactly the same result as a fresh computation!" << std::endl;
        } else {
            std::cout << "⚠ Output files differ." << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "❌ Error in partial computation demo: " << e.what() << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::cout << "FK Computation Checkpoint Resume Demo" << std::endl;
    std::cout << "====================================" << std::endl;
    std::cout << std::endl;

    if (argc == 1) {
        // No arguments - run interactive demo
        interactiveDemoResumption();
    } else if (argc == 2) {
        // One argument - resume from specific checkpoint file
        std::string checkpoint_path = argv[1];
        demoCheckpointResumption(checkpoint_path);
    } else if (argc == 3 && std::string(argv[1]) == "--partial") {
        // Special mode: demonstrate partial computation and resume
        demoPartialComputationAndResume();
    } else {
        // Show usage
        std::cout << "Usage:" << std::endl;
        std::cout << "  " << argv[0] << "                          # Resume from most recent checkpoint" << std::endl;
        std::cout << "  " << argv[0] << " <checkpoint_file>        # Resume from specific checkpoint" << std::endl;
        std::cout << "  " << argv[0] << " --partial demo           # Demo partial computation and resume" << std::endl;
        std::cout << std::endl;
        std::cout << "Examples:" << std::endl;
        std::cout << "  " << argv[0] << std::endl;
        std::cout << "  " << argv[0] << " fk_checkpoint_1234567890_5678.json" << std::endl;
        std::cout << "  " << argv[0] << " --partial demo" << std::endl;
        return 1;
    }

    return 0;
}
