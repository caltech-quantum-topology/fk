#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include "fk/fk_computation.hpp"
#include "fk/string_to_int.hpp"

// Original FK parsing logic extracted from fk_segments_links.cpp
struct OriginalFKConfig {
    int degree;
    int components;
    int writhe;
    int prefactors;
    int crossings;

    std::vector<int> closed_strand_components;
    std::vector<std::vector<std::array<int, 2>>> crossing_matrices;
    std::vector<int> crossing_relation_types;
    std::vector<int> top_crossing_components;
    std::vector<int> bottom_crossing_components;

    std::vector<std::vector<double>> criteria;
    std::vector<std::vector<double>> inequalities;
    std::vector<std::vector<std::vector<int>>> variable_assignments;
};

OriginalFKConfig parseWithOriginalLogic(const std::string& filename) {
    OriginalFKConfig config;

    std::ifstream infile;
    infile.open(filename + ".csv");
    if (!infile.is_open()) {
        throw std::runtime_error("Unable to open file '" + filename + ".csv'");
    }

    std::string line;

    // Parse degree (line 1)
    std::getline(infile, line, '\n');
    int index = line.find(",");
    config.degree = parseStringToInteger(line.substr(0, index));

    // Parse components (line 2)
    std::getline(infile, line, '\n');
    index = line.find(",");
    config.components = parseStringToInteger(line.substr(0, index));

    // Parse writhe (line 3)
    std::getline(infile, line, '\n');
    index = line.find(",");
    config.writhe = parseStringToInteger(line.substr(0, index));

    // Parse crossing data (line 4)
    std::getline(infile, line, '\n');
    int height = 0;
    index = line.find(",");
    while (true) {
        if (index == -1) {
            break;
        }

        int c = parseStringToInteger(line.substr(0, index));
        config.crossing_matrices.push_back({{height, c - 1},
                                           {height, c},
                                           {height + 1, c - 1},
                                           {height + 1, c}});
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");

        config.crossing_relation_types.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);

        height++;
        index = line.find(",");
    }

    // Parse closed strand components (line 5)
    std::getline(infile, line, '\n');
    index = line.find(",");
    while (true) {
        if (index == -1) {
            break;
        }
        config.closed_strand_components.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");
    }
    config.prefactors = config.closed_strand_components.size();
    config.crossings = config.crossing_relation_types.size();

    // Parse top/bottom crossing components (line 6)
    std::getline(infile, line, '\n');
    index = line.find(",");
    while (true) {
        if (index == -1) {
            break;
        }
        config.top_crossing_components.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");
        config.bottom_crossing_components.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");
    }

    // Parse remaining sections
    int stage = 0;
    int criteria_index = 0;
    int inequality_index = 0;
    int extension_index = 0;

    config.variable_assignments.resize(config.crossings + 1,
                                     std::vector<std::vector<int>>(config.prefactors + 1));

    while (std::getline(infile, line, '\n')) {
        if (line[0] == '/') {
            stage++;
        } else if (stage == 0) {
            // Criteria section
            config.criteria.push_back({});
            index = line.find(",");
            while (true) {
                if (index == -1) {
                    break;
                }
                config.criteria[criteria_index].push_back(
                    parseStringToDouble(line.substr(0, index)));
                line = line.substr(index + 1, line.size() - index - 1);
                index = line.find(",");
            }
            criteria_index++;
        } else if (stage == 1) {
            // Inequalities section
            config.inequalities.push_back({});
            index = line.find(",");
            while (true) {
                if (index == -1) {
                    break;
                }
                config.inequalities[inequality_index].push_back(
                    parseStringToInteger(line.substr(0, index)));
                line = line.substr(index + 1, line.size() - index - 1);
                index = line.find(",");
            }
            inequality_index++;
        } else if (stage == 2) {
            // Variable assignments section
            index = line.find(",");
            while (true) {
                if (index == -1) {
                    break;
                }
                config.variable_assignments[extension_index % (config.crossings + 1)]
                                          [extension_index / (config.crossings + 1)]
                                              .push_back(parseStringToInteger(
                                                  line.substr(0, index)));
                line = line.substr(index + 1, line.size() - index - 1);
                index = line.find(",");
            }
            extension_index++;
        }
    }

    return config;
}

void printOriginalConfig(const OriginalFKConfig& config) {
    std::cout << "=== Original Parser Configuration ===" << std::endl;
    std::cout << "Basic parameters:" << std::endl;
    std::cout << "  Degree: " << config.degree << std::endl;
    std::cout << "  Components: " << config.components << std::endl;
    std::cout << "  Writhe: " << config.writhe << std::endl;
    std::cout << "  Crossings: " << config.crossings << std::endl;
    std::cout << "  Prefactors: " << config.prefactors << std::endl;

    std::cout << "Closed strand components (" << config.closed_strand_components.size() << "): ";
    for (size_t i = 0; i < config.closed_strand_components.size(); ++i) {
        std::cout << config.closed_strand_components[i];
        if (i < config.closed_strand_components.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    std::cout << "Crossing relation types (" << config.crossing_relation_types.size() << "): ";
    for (size_t i = 0; i < config.crossing_relation_types.size(); ++i) {
        std::cout << config.crossing_relation_types[i];
        if (i < config.crossing_relation_types.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    std::cout << "Top crossing components (" << config.top_crossing_components.size() << "): ";
    for (size_t i = 0; i < config.top_crossing_components.size(); ++i) {
        std::cout << config.top_crossing_components[i];
        if (i < config.top_crossing_components.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    std::cout << "Bottom crossing components (" << config.bottom_crossing_components.size() << "): ";
    for (size_t i = 0; i < config.bottom_crossing_components.size(); ++i) {
        std::cout << config.bottom_crossing_components[i];
        if (i < config.bottom_crossing_components.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    std::cout << "Crossing matrices (" << config.crossing_matrices.size() << "):" << std::endl;
    for (size_t i = 0; i < config.crossing_matrices.size(); ++i) {
        std::cout << "  [" << i << "]: ";
        for (const auto& pair : config.crossing_matrices[i]) {
            std::cout << "{" << pair[0] << "," << pair[1] << "} ";
        }
        std::cout << std::endl;
    }

    std::cout << "Criteria vectors (" << config.criteria.size() << "):" << std::endl;
    for (size_t i = 0; i < std::min(config.criteria.size(), size_t(3)); ++i) {
        std::cout << "  [" << i << "]: ";
        for (size_t j = 0; j < config.criteria[i].size(); ++j) {
            std::cout << config.criteria[i][j];
            if (j < config.criteria[i].size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }
    if (config.criteria.size() > 3) {
        std::cout << "  ... and " << (config.criteria.size() - 3) << " more" << std::endl;
    }

    std::cout << "Inequality vectors (" << config.inequalities.size() << "):" << std::endl;
    for (size_t i = 0; i < std::min(config.inequalities.size(), size_t(3)); ++i) {
        std::cout << "  [" << i << "]: ";
        for (size_t j = 0; j < config.inequalities[i].size(); ++j) {
            std::cout << config.inequalities[i][j];
            if (j < config.inequalities[i].size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }
    if (config.inequalities.size() > 3) {
        std::cout << "  ... and " << (config.inequalities.size() - 3) << " more" << std::endl;
    }

    std::cout << "Variable assignments (" << config.variable_assignments.size() << "x"
              << (config.variable_assignments.empty() ? 0 : config.variable_assignments[0].size()) << "):" << std::endl;
    for (size_t i = 0; i < std::min(config.variable_assignments.size(), size_t(3)); ++i) {
        for (size_t j = 0; j < config.variable_assignments[i].size(); ++j) {
            if (!config.variable_assignments[i][j].empty()) {
                std::cout << "  [" << i << "][" << j << "]: ";
                for (size_t k = 0; k < config.variable_assignments[i][j].size(); ++k) {
                    std::cout << config.variable_assignments[i][j][k];
                    if (k < config.variable_assignments[i][j].size() - 1) std::cout << ", ";
                }
                std::cout << std::endl;
            }
        }
    }
}

void printRefactoredConfig(const fk::FKConfiguration& config) {
    std::cout << "=== Refactored Parser Configuration ===" << std::endl;
    std::cout << "Basic parameters:" << std::endl;
    std::cout << "  Degree: " << config.degree << std::endl;
    std::cout << "  Components: " << config.components << std::endl;
    std::cout << "  Writhe: " << config.writhe << std::endl;
    std::cout << "  Crossings: " << config.crossings << std::endl;
    std::cout << "  Prefactors: " << config.prefactors << std::endl;

    std::cout << "Closed strand components (" << config.closed_strand_components.size() << "): ";
    for (size_t i = 0; i < config.closed_strand_components.size(); ++i) {
        std::cout << config.closed_strand_components[i];
        if (i < config.closed_strand_components.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    std::cout << "Crossing relation types (" << config.crossing_relation_types.size() << "): ";
    for (size_t i = 0; i < config.crossing_relation_types.size(); ++i) {
        std::cout << config.crossing_relation_types[i];
        if (i < config.crossing_relation_types.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    std::cout << "Top crossing components (" << config.top_crossing_components.size() << "): ";
    for (size_t i = 0; i < config.top_crossing_components.size(); ++i) {
        std::cout << config.top_crossing_components[i];
        if (i < config.top_crossing_components.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    std::cout << "Bottom crossing components (" << config.bottom_crossing_components.size() << "): ";
    for (size_t i = 0; i < config.bottom_crossing_components.size(); ++i) {
        std::cout << config.bottom_crossing_components[i];
        if (i < config.bottom_crossing_components.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    std::cout << "Crossing matrices (" << config.crossing_matrices.size() << "):" << std::endl;
    for (size_t i = 0; i < config.crossing_matrices.size(); ++i) {
        std::cout << "  [" << i << "]: ";
        for (const auto& pair : config.crossing_matrices[i]) {
            std::cout << "{" << pair[0] << "," << pair[1] << "} ";
        }
        std::cout << std::endl;
    }

    std::cout << "Criteria vectors (" << config.criteria.size() << "):" << std::endl;
    for (size_t i = 0; i < std::min(config.criteria.size(), size_t(3)); ++i) {
        std::cout << "  [" << i << "]: ";
        for (size_t j = 0; j < config.criteria[i].size(); ++j) {
            std::cout << config.criteria[i][j];
            if (j < config.criteria[i].size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }
    if (config.criteria.size() > 3) {
        std::cout << "  ... and " << (config.criteria.size() - 3) << " more" << std::endl;
    }

    std::cout << "Inequality vectors (" << config.inequalities.size() << "):" << std::endl;
    for (size_t i = 0; i < std::min(config.inequalities.size(), size_t(3)); ++i) {
        std::cout << "  [" << i << "]: ";
        for (size_t j = 0; j < config.inequalities[i].size(); ++j) {
            std::cout << config.inequalities[i][j];
            if (j < config.inequalities[i].size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }
    if (config.inequalities.size() > 3) {
        std::cout << "  ... and " << (config.inequalities.size() - 3) << " more" << std::endl;
    }

    std::cout << "Variable assignments (" << config.variable_assignments.size() << "x"
              << (config.variable_assignments.empty() ? 0 : config.variable_assignments[0].size()) << "):" << std::endl;
    for (size_t i = 0; i < std::min(config.variable_assignments.size(), size_t(3)); ++i) {
        for (size_t j = 0; j < config.variable_assignments[i].size(); ++j) {
            if (!config.variable_assignments[i][j].empty()) {
                std::cout << "  [" << i << "][" << j << "]: ";
                for (size_t k = 0; k < config.variable_assignments[i][j].size(); ++k) {
                    std::cout << config.variable_assignments[i][j][k];
                    if (k < config.variable_assignments[i][j].size() - 1) std::cout << ", ";
                }
                std::cout << std::endl;
            }
        }
    }
}

bool compareConfigurations(const OriginalFKConfig& orig, const fk::FKConfiguration& refact) {
    bool identical = true;

    std::cout << "\n=== Configuration Comparison ===" << std::endl;

    // Basic parameters
    if (orig.degree != refact.degree) {
        std::cout << "✗ Degree mismatch: " << orig.degree << " vs " << refact.degree << std::endl;
        identical = false;
    } else {
        std::cout << "✓ Degree matches: " << orig.degree << std::endl;
    }

    if (orig.components != refact.components) {
        std::cout << "✗ Components mismatch: " << orig.components << " vs " << refact.components << std::endl;
        identical = false;
    } else {
        std::cout << "✓ Components matches: " << orig.components << std::endl;
    }

    if (orig.writhe != refact.writhe) {
        std::cout << "✗ Writhe mismatch: " << orig.writhe << " vs " << refact.writhe << std::endl;
        identical = false;
    } else {
        std::cout << "✓ Writhe matches: " << orig.writhe << std::endl;
    }

    if (orig.crossings != refact.crossings) {
        std::cout << "✗ Crossings mismatch: " << orig.crossings << " vs " << refact.crossings << std::endl;
        identical = false;
    } else {
        std::cout << "✓ Crossings matches: " << orig.crossings << std::endl;
    }

    if (orig.prefactors != refact.prefactors) {
        std::cout << "✗ Prefactors mismatch: " << orig.prefactors << " vs " << refact.prefactors << std::endl;
        identical = false;
    } else {
        std::cout << "✓ Prefactors matches: " << orig.prefactors << std::endl;
    }

    // Vector comparisons
    if (orig.closed_strand_components != refact.closed_strand_components) {
        std::cout << "✗ Closed strand components mismatch" << std::endl;
        identical = false;
    } else {
        std::cout << "✓ Closed strand components match" << std::endl;
    }

    if (orig.crossing_relation_types != refact.crossing_relation_types) {
        std::cout << "✗ Crossing relation types mismatch" << std::endl;
        identical = false;
    } else {
        std::cout << "✓ Crossing relation types match" << std::endl;
    }

    if (orig.top_crossing_components != refact.top_crossing_components) {
        std::cout << "✗ Top crossing components mismatch" << std::endl;
        identical = false;
    } else {
        std::cout << "✓ Top crossing components match" << std::endl;
    }

    if (orig.bottom_crossing_components != refact.bottom_crossing_components) {
        std::cout << "✗ Bottom crossing components mismatch" << std::endl;
        identical = false;
    } else {
        std::cout << "✓ Bottom crossing components match" << std::endl;
    }

    // Crossing matrices comparison
    if (orig.crossing_matrices.size() != refact.crossing_matrices.size()) {
        std::cout << "✗ Crossing matrices size mismatch: " << orig.crossing_matrices.size() << " vs " << refact.crossing_matrices.size() << std::endl;
        identical = false;
    } else {
        bool matrices_match = true;
        for (size_t i = 0; i < orig.crossing_matrices.size(); ++i) {
            for (size_t j = 0; j < 4; ++j) {
                if (orig.crossing_matrices[i][j][0] != refact.crossing_matrices[i][j][0] ||
                    orig.crossing_matrices[i][j][1] != refact.crossing_matrices[i][j][1]) {
                    matrices_match = false;
                    break;
                }
            }
            if (!matrices_match) break;
        }
        if (matrices_match) {
            std::cout << "✓ Crossing matrices match" << std::endl;
        } else {
            std::cout << "✗ Crossing matrices content mismatch" << std::endl;
            identical = false;
        }
    }

    // Criteria comparison (convert to int for comparison)
    if (orig.criteria.size() != refact.criteria.size()) {
        std::cout << "✗ Criteria size mismatch: " << orig.criteria.size() << " vs " << refact.criteria.size() << std::endl;
        identical = false;
    } else {
        bool criteria_match = true;
        for (size_t i = 0; i < orig.criteria.size(); ++i) {
            if (orig.criteria[i].size() != refact.criteria[i].size()) {
                criteria_match = false;
                break;
            }
            for (size_t j = 0; j < orig.criteria[i].size(); ++j) {
                if (std::abs(orig.criteria[i][j] - refact.criteria[i][j]) > 1e-9) {
                    criteria_match = false;
                    break;
                }
            }
            if (!criteria_match) break;
        }
        if (criteria_match) {
            std::cout << "✓ Criteria vectors match" << std::endl;
        } else {
            std::cout << "✗ Criteria vectors mismatch" << std::endl;
            identical = false;
        }
    }

    // Inequalities comparison
    if (orig.inequalities.size() != refact.inequalities.size()) {
        std::cout << "✗ Inequalities size mismatch: " << orig.inequalities.size() << " vs " << refact.inequalities.size() << std::endl;
        identical = false;
    } else {
        bool inequalities_match = true;
        for (size_t i = 0; i < orig.inequalities.size(); ++i) {
            if (orig.inequalities[i].size() != refact.inequalities[i].size()) {
                inequalities_match = false;
                break;
            }
            for (size_t j = 0; j < orig.inequalities[i].size(); ++j) {
                if (std::abs(orig.inequalities[i][j] - refact.inequalities[i][j]) > 1e-9) {
                    inequalities_match = false;
                    break;
                }
            }
            if (!inequalities_match) break;
        }
        if (inequalities_match) {
            std::cout << "✓ Inequality vectors match" << std::endl;
        } else {
            std::cout << "✗ Inequality vectors mismatch" << std::endl;
            identical = false;
        }
    }

    // Variable assignments comparison
    if (orig.variable_assignments.size() != refact.variable_assignments.size()) {
        std::cout << "✗ Variable assignments outer size mismatch" << std::endl;
        identical = false;
    } else {
        bool assignments_match = true;
        for (size_t i = 0; i < orig.variable_assignments.size(); ++i) {
            if (orig.variable_assignments[i].size() != refact.variable_assignments[i].size()) {
                assignments_match = false;
                break;
            }
            for (size_t j = 0; j < orig.variable_assignments[i].size(); ++j) {
                if (orig.variable_assignments[i][j] != refact.variable_assignments[i][j]) {
                    assignments_match = false;
                    break;
                }
            }
            if (!assignments_match) break;
        }
        if (assignments_match) {
            std::cout << "✓ Variable assignments match" << std::endl;
        } else {
            std::cout << "✗ Variable assignments mismatch" << std::endl;
            identical = false;
        }
    }

    return identical;
}

int main() {
    try {
        std::cout << "=== FK Parser Comparison Test ===" << std::endl;
        std::cout << "Testing trefoil_ilp.csv parsing between original and refactored code" << std::endl << std::endl;

        // Parse with original logic
        std::cout << "[1/4] Parsing with original logic..." << std::endl;
        OriginalFKConfig original_config = parseWithOriginalLogic("examples/trefoil_ilp");
        std::cout << "✓ Original parsing completed" << std::endl;

        // Parse with refactored logic
        std::cout << "\n[2/4] Parsing with refactored logic..." << std::endl;
        fk::FKInputParser refactored_parser;
        fk::FKConfiguration refactored_config = refactored_parser.parseFromFile("examples/trefoil_ilp");
        std::cout << "✓ Refactored parsing completed" << std::endl;

        // Display both configurations
        std::cout << "\n[3/4] Displaying configurations..." << std::endl;
        printOriginalConfig(original_config);
        std::cout << std::endl;
        printRefactoredConfig(refactored_config);

        // Compare configurations
        std::cout << "\n[4/4] Comparing configurations..." << std::endl;
        bool identical = compareConfigurations(original_config, refactored_config);

        std::cout << "\n=== Test Results ===" << std::endl;
        if (identical) {
            std::cout << "🎉 SUCCESS: Configurations are identical!" << std::endl;
            std::cout << "✓ The refactored parser produces exactly the same results" << std::endl;
            std::cout << "✓ All data structures match perfectly" << std::endl;
            std::cout << "✓ Parser refactoring is correct" << std::endl;
        } else {
            std::cout << "❌ FAILURE: Configurations differ!" << std::endl;
            std::cout << "✗ The refactored parser has parsing bugs" << std::endl;
            std::cout << "✗ Data structures do not match" << std::endl;
            std::cout << "✗ Parser refactoring needs fixes" << std::endl;
        }

        return identical ? 0 : 1;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}