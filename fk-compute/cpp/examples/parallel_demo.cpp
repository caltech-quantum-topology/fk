#include "fk/parallel_pool.hpp"
#include "fk/solution_pool_1a_double_links.hpp"

#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <vector>

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n"
              << "Options:\n"
              << "  --threads N        Use N threads (default: auto-detect)\n"
              << "  --sequential       Run sequential version for comparison\n"
              << "  --benchmark        Run both sequential and parallel for benchmarking\n"
              << "  --problem SIZE     Set problem size (small/medium/large, default: medium)\n"
              << "  --help             Show this help message\n"
              << "\nExample:\n"
              << "  " << program_name << " --threads 8 --benchmark\n"
              << std::endl;
}

// Create test problems of different sizes
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
createTestProblem(const std::string& size) {
    std::vector<std::vector<double>> main_inequalities;
    std::vector<std::vector<double>> supporting_inequalities;

    if (size == "small") {
        // Small problem: 3 variables
        main_inequalities = {
            {10.0, -1.0, -1.0, -1.0},  // 10 - x1 - x2 - x3 >= 0
            {8.0, -1.0, 0.0, 0.0},     // 8 - x1 >= 0
            {6.0, 0.0, -1.0, 0.0},     // 6 - x2 >= 0
            {4.0, 0.0, 0.0, -1.0}      // 4 - x3 >= 0
        };

        supporting_inequalities = {
            {0.0, 1.0, 0.0, 0.0},      // x1 >= 0
            {0.0, 0.0, 1.0, 0.0},      // x2 >= 0
            {0.0, 0.0, 0.0, 1.0}       // x3 >= 0
        };
    } else if (size == "medium") {
        // Medium problem: 4 variables
        main_inequalities = {
            {20.0, -1.0, -1.0, -1.0, -1.0},  // 20 - x1 - x2 - x3 - x4 >= 0
            {15.0, -2.0, -1.0, 0.0, 0.0},    // 15 - 2*x1 - x2 >= 0
            {12.0, 0.0, -1.0, -2.0, 0.0},    // 12 - x2 - 2*x3 >= 0
            {10.0, 0.0, 0.0, -1.0, -2.0},    // 10 - x3 - 2*x4 >= 0
            {8.0, -1.0, 0.0, 0.0, -1.0}      // 8 - x1 - x4 >= 0
        };

        supporting_inequalities = {
            {0.0, 1.0, 0.0, 0.0, 0.0},       // x1 >= 0
            {0.0, 0.0, 1.0, 0.0, 0.0},       // x2 >= 0
            {0.0, 0.0, 0.0, 1.0, 0.0},       // x3 >= 0
            {0.0, 0.0, 0.0, 0.0, 1.0}        // x4 >= 0
        };
    } else { // large
        // Large problem: 5 variables
        main_inequalities = {
            {30.0, -1.0, -1.0, -1.0, -1.0, -1.0},  // 30 - x1 - x2 - x3 - x4 - x5 >= 0
            {25.0, -2.0, -1.0, -1.0, 0.0, 0.0},    // 25 - 2*x1 - x2 - x3 >= 0
            {20.0, -1.0, -2.0, 0.0, -1.0, 0.0},    // 20 - x1 - 2*x2 - x4 >= 0
            {18.0, 0.0, -1.0, -2.0, 0.0, -1.0},    // 18 - x2 - 2*x3 - x5 >= 0
            {15.0, 0.0, 0.0, -1.0, -2.0, -1.0},    // 15 - x3 - 2*x4 - x5 >= 0
            {12.0, -1.0, 0.0, 0.0, -1.0, -2.0}     // 12 - x1 - x4 - 2*x5 >= 0
        };

        supporting_inequalities = {
            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},        // x1 >= 0
            {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},        // x2 >= 0
            {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},        // x3 >= 0
            {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},        // x4 >= 0
            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}         // x5 >= 0
        };
    }

    return {main_inequalities, supporting_inequalities};
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    int num_threads = -1;  // auto-detect
    bool run_sequential = false;
    bool run_benchmark = false;
    std::string problem_size = "medium";

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--sequential") == 0) {
            run_sequential = true;
        } else if (std::strcmp(argv[i], "--benchmark") == 0) {
            run_benchmark = true;
        } else if (std::strcmp(argv[i], "--problem") == 0 && i + 1 < argc) {
            problem_size = argv[++i];
        } else if (std::strcmp(argv[i], "--help") == 0) {
            printUsage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown option: " << argv[i] << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    std::cout << "=== Parallel Solution Pool Demo ===" << std::endl;

    // Determine optimal thread count
    int optimal_threads = getOptimalThreadCount();
    if (num_threads <= 0) {
        num_threads = optimal_threads;
    }

    std::cout << "Hardware threads detected: " << optimal_threads << std::endl;
    std::cout << "Using threads: " << num_threads << std::endl;
    std::cout << "Problem size: " << problem_size << std::endl;
    std::cout << std::endl;

    // Create test problem
    auto [main_inequalities, supporting_inequalities] = createTestProblem(problem_size);

    std::cout << "Problem configuration:" << std::endl;
    std::cout << "  Variables: " << (main_inequalities[0].size() - 1) << std::endl;
    std::cout << "  Main inequalities: " << main_inequalities.size() << std::endl;
    std::cout << "  Supporting inequalities: " << supporting_inequalities.size() << std::endl;
    std::cout << std::endl;

    // Solution counters and timing
    int sequential_solutions = 0;
    int parallel_solutions = 0;
    std::chrono::milliseconds sequential_time(0);
    std::chrono::milliseconds parallel_time(0);

    // Sequential version
    if (run_sequential || run_benchmark) {
        std::cout << "Running sequential version..." << std::endl;
        sequential_solutions = 0;

        auto sequential_handler = [&](const std::vector<int>& point) {
            sequential_solutions++;
            if (sequential_solutions <= 10 || sequential_solutions % 100 == 0) {
                std::cout << "Sequential solution " << sequential_solutions << ": (";
                for (size_t i = 0; i < point.size(); i++) {
                    std::cout << point[i];
                    if (i < point.size() - 1) std::cout << ", ";
                }
                std::cout << ")" << std::endl;
            }
        };

        auto start_time = std::chrono::steady_clock::now();
        pooling(main_inequalities, supporting_inequalities, sequential_handler);
        auto end_time = std::chrono::steady_clock::now();

        sequential_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        std::cout << "Sequential computation completed:" << std::endl;
        std::cout << "  Solutions found: " << sequential_solutions << std::endl;
        std::cout << "  Time: " << sequential_time.count() << "ms" << std::endl;
        std::cout << std::endl;
    }

    // Parallel version
    if (!run_sequential || run_benchmark) {
        std::cout << "Running parallel version..." << std::endl;
        parallel_solutions = 0;

        auto parallel_handler = [&](const std::vector<int>& point) {
            parallel_solutions++;
            if (parallel_solutions <= 10 || parallel_solutions % 100 == 0) {
                std::cout << "Parallel solution " << parallel_solutions << ": (";
                for (size_t i = 0; i < point.size(); i++) {
                    std::cout << point[i];
                    if (i < point.size() - 1) std::cout << ", ";
                }
                std::cout << ")" << std::endl;
            }
        };

        auto start_time = std::chrono::steady_clock::now();
        parallelPooling(main_inequalities, supporting_inequalities, parallel_handler, num_threads);
        auto end_time = std::chrono::steady_clock::now();

        parallel_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        std::cout << "Parallel computation completed:" << std::endl;
        std::cout << "  Solutions found: " << parallel_solutions << std::endl;
        std::cout << "  Time: " << parallel_time.count() << "ms" << std::endl;
        std::cout << std::endl;
    }

    // Benchmark results
    if (run_benchmark) {
        std::cout << "=== Benchmark Results ===" << std::endl;
        std::cout << "Sequential: " << sequential_solutions << " solutions in "
                  << sequential_time.count() << "ms" << std::endl;
        std::cout << "Parallel (" << num_threads << " threads): " << parallel_solutions
                  << " solutions in " << parallel_time.count() << "ms" << std::endl;

        if (sequential_solutions == parallel_solutions) {
            std::cout << "✓ Solution counts match!" << std::endl;
        } else {
            std::cout << "⚠ Warning: Solution counts differ!" << std::endl;
        }

        if (sequential_time.count() > 0) {
            double speedup = static_cast<double>(sequential_time.count()) / parallel_time.count();
            double efficiency = speedup / num_threads;
            std::cout << "Speedup: " << std::fixed << std::setprecision(2) << speedup << "x" << std::endl;
            std::cout << "Efficiency: " << std::fixed << std::setprecision(1) << (efficiency * 100) << "%" << std::endl;
        }
    }

    return 0;
}