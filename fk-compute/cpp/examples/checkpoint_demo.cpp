#include "fk/checkpoint.hpp"
#include "fk/solution_pool_1a_double_links.hpp"

#include <iostream>
#include <vector>
#include <csignal>
#include <cstring>

// Global variables for signal handling
volatile sig_atomic_t interrupted = 0;
std::string checkpoint_file_global;

void signalHandler(int signal) {
    std::cout << "\nReceived signal " << signal << " (";
    if (signal == SIGINT) {
        std::cout << "SIGINT";
    } else if (signal == SIGTERM) {
        std::cout << "SIGTERM";
    }
    std::cout << "). Saving checkpoint and exiting..." << std::endl;
    interrupted = 1;
}

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n"
              << "Options:\n"
              << "  --checkpoint FILE   Use FILE for checkpointing\n"
              << "  --resume           Resume from checkpoint if available\n"
              << "  --interval N       Set checkpoint interval (default: 100)\n"
              << "  --help             Show this help message\n"
              << "\nExample:\n"
              << "  " << program_name << " --checkpoint demo.chk --resume --interval 50\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    std::string checkpoint_file;
    bool resume_from_checkpoint = false;
    size_t checkpoint_interval = 100;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--checkpoint") == 0 && i + 1 < argc) {
            checkpoint_file = argv[++i];
            checkpoint_file_global = checkpoint_file;
        } else if (std::strcmp(argv[i], "--resume") == 0) {
            resume_from_checkpoint = true;
        } else if (std::strcmp(argv[i], "--interval") == 0 && i + 1 < argc) {
            checkpoint_interval = std::stoul(argv[++i]);
        } else if (std::strcmp(argv[i], "--help") == 0) {
            printUsage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown option: " << argv[i] << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    // Set up signal handlers for graceful shutdown
    std::signal(SIGINT, signalHandler);
    std::signal(SIGTERM, signalHandler);

    std::cout << "=== Checkpoint Demo ===" << std::endl;

    if (!checkpoint_file.empty()) {
        std::cout << "Checkpoint file: " << checkpoint_file << std::endl;
        std::cout << "Checkpoint interval: " << checkpoint_interval << std::endl;
        std::cout << "Resume from checkpoint: " << (resume_from_checkpoint ? "yes" : "no") << std::endl;
    } else {
        std::cout << "No checkpoint file specified. Running without checkpointing." << std::endl;
    }

    std::cout << std::endl;

    // Example problem setup (small test case)
    std::vector<std::vector<double>> main_inequalities = {
        {10.0, -1.0, -1.0},  // 10 - x1 - x2 >= 0
        {8.0, -1.0, 0.0},    // 8 - x1 >= 0
        {6.0, 0.0, -1.0}     // 6 - x2 >= 0
    };

    std::vector<std::vector<double>> supporting_inequalities = {
        {0.0, 1.0, 0.0},     // x1 >= 0
        {0.0, 0.0, 1.0}      // x2 >= 0
    };

    // Solution counter
    int solution_count = 0;

    // Function to process each found solution
    auto solution_handler = [&](const std::vector<int>& point) {
        solution_count++;
        std::cout << "Solution " << solution_count << ": (";
        for (size_t i = 0; i < point.size(); i++) {
            std::cout << point[i];
            if (i < point.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;

        // Check for interruption during solution processing
        if (interrupted) {
            std::cout << "Interrupted during solution processing." << std::endl;
            return;
        }
    };

    try {
        if (!checkpoint_file.empty()) {
            // Use checkpointing version
            poolingWithCheckpoints(main_inequalities, supporting_inequalities,
                                 solution_handler, checkpoint_file,
                                 resume_from_checkpoint, checkpoint_interval);
        } else {
            // Use regular version
            pooling(main_inequalities, supporting_inequalities, solution_handler);
        }

        std::cout << "\nComputation completed successfully!" << std::endl;
        std::cout << "Total solutions found: " << solution_count << std::endl;

        // Clean up checkpoint file if computation completed successfully
        if (!checkpoint_file.empty()) {
            CheckpointManager manager(checkpoint_file, checkpoint_interval);
            if (manager.checkpointExists()) {
                std::cout << "Cleaning up checkpoint file..." << std::endl;
                manager.deleteCheckpoint();
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error during computation: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}