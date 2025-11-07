#include "fk/fk_computation.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " <input_file> <output_file> [options]\n";
    std::cout << "\nArguments:\n";
    std::cout << "  input_file   Input CSV filename (without .csv extension)\n";
    std::cout << "  output_file  Output JSON filename (without .json extension)\n";
    std::cout << "\nOptions:\n";
    std::cout << "  --threads N  Number of computation engines/threads to use (default: 1)\n";
    std::cout << "  --verbose    Show detailed configuration information\n";
    std::cout << "  --help       Show this help message\n";
    std::cout << "\nExamples:\n";
    std::cout << "  " << program_name << " trefoil_ilp trefoil_output\n";
    std::cout << "  " << program_name << " input output --threads 4 --verbose\n";
}


void printConfiguration(const fk::FKConfiguration& config) {
    std::cout << "\n=== FK Configuration ===\n";
    std::cout << "Degree: " << config.degree << "\n";
    std::cout << "Components: " << config.components << "\n";
    std::cout << "Writhe: " << config.writhe << "\n";
    std::cout << "Crossings: " << config.crossings << "\n";
    std::cout << "Prefactors: " << config.prefactors << "\n";

    if (!config.closed_strand_components.empty()) {
        std::cout << "\nClosed strand components: ";
        for (size_t i = 0; i < config.closed_strand_components.size(); ++i) {
            std::cout << config.closed_strand_components[i];
            if (i < config.closed_strand_components.size() - 1) std::cout << ", ";
        }
        std::cout << "\n";
    }

    if (!config.crossing_relation_types.empty()) {
        std::cout << "Crossing relation types: ";
        for (size_t i = 0; i < config.crossing_relation_types.size(); ++i) {
            std::cout << config.crossing_relation_types[i];
            if (i < config.crossing_relation_types.size() - 1) std::cout << ", ";
        }
        std::cout << "\n";
    }

    if (!config.top_crossing_components.empty() && !config.bottom_crossing_components.empty()) {
        std::cout << "Top crossing components: ";
        for (size_t i = 0; i < config.top_crossing_components.size(); ++i) {
            std::cout << config.top_crossing_components[i];
            if (i < config.top_crossing_components.size() - 1) std::cout << ", ";
        }
        std::cout << "\n";
        std::cout << "Bottom crossing components: ";
        for (size_t i = 0; i < config.bottom_crossing_components.size(); ++i) {
            std::cout << config.bottom_crossing_components[i];
            if (i < config.bottom_crossing_components.size() - 1) std::cout << ", ";
        }
        std::cout << "\n";
    }

    std::cout << "\nCriteria vectors: " << config.criteria.size() << " total\n";
    if (!config.criteria.empty()) {
        for (size_t i = 0; i < std::min(config.criteria.size(), size_t(3)); ++i) {
            std::cout << "  Criteria[" << i << "]: ";
            for (size_t j = 0; j < config.criteria[i].size(); ++j) {
                std::cout << config.criteria[i][j];
                if (j < config.criteria[i].size() - 1) std::cout << ", ";
            }
            std::cout << "\n";
        }
        if (config.criteria.size() > 3) {
            std::cout << "  ... and " << (config.criteria.size() - 3) << " more\n";
        }
    }

    std::cout << "\nInequalities vectors: " << config.inequalities.size() << " total\n";
    if (!config.inequalities.empty()) {
        for (size_t i = 0; i < std::min(config.inequalities.size(), size_t(3)); ++i) {
            std::cout << "  Inequality[" << i << "]: ";
            for (size_t j = 0; j < config.inequalities[i].size(); ++j) {
                std::cout << config.inequalities[i][j];
                if (j < config.inequalities[i].size() - 1) std::cout << ", ";
            }
            std::cout << "\n";
        }
        if (config.inequalities.size() > 3) {
            std::cout << "  ... and " << (config.inequalities.size() - 3) << " more\n";
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        printUsage(argv[0]);
        return 1;
    }

    // Parse command line arguments
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    bool verbose = false;
    int num_threads = 1;

    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--verbose") {
            verbose = true;
        } else if (arg == "--threads") {
            if (i + 1 < argc) {
                try {
                    num_threads = std::stoi(argv[++i]);
                    if (num_threads < 1) {
                        std::cerr << "Number of threads must be at least 1\n";
                        return 1;
                    }
                } catch (const std::exception& e) {
                    std::cerr << "Invalid thread count: " << argv[i] << "\n";
                    return 1;
                }
            } else {
                std::cerr << "--threads requires a number\n";
                printUsage(argv[0]);
                return 1;
            }
        } else if (arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    std::cout << "FK Computation Tool\n";
    std::cout << "===================\n";
    std::cout << "Input file: " << input_file << ".csv\n";
    std::cout << "Output file: " << output_file << ".json\n";
    std::cout << "Threads: " << num_threads << "\n\n";

    try {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Create FK computation instance
        fk::FKComputation computation;

        // Perform computation
        std::cout << "Starting FK computation...\n";

        computation.compute(input_file, output_file, num_threads);


        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        std::cout << "✓ Computation completed successfully!\n";
        std::cout << "⏱ Computation time: " << duration.count() << " ms\n";

        // Show configuration details if requested
        if (verbose) {
            const auto& config = computation.getLastConfiguration();
            printConfiguration(config);
        }

        // Show output information
        std::cout << "\n📄 Results written to: " << output_file << ".json\n";

        // Show brief result statistics
        const auto& result = computation.getLastResult();

        std::cout << "📊 Result contains " << result.nTerms()
                  << " non-zero terms\n";

    } catch (const std::exception& e) {
        std::cerr << "\n❌ Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
