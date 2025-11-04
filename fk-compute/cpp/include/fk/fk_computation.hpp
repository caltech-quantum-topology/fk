#pragma once

#include "fk/multivariable_polynomial.hpp"
#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <array>
#include <list>
#include <queue>
#include <set>

namespace fk {

/**
 * Configuration data for FK computation
 * Holds all the parameters and structural data needed for computation
 */
struct FKConfiguration {
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

    FKConfiguration() = default;

    bool isValid() const;
    void clear();
};

/**
 * Parser for FK input files
 * Handles reading and parsing CSV input files into FKConfiguration
 */
class FKInputParser {
public:
    FKInputParser() = default;

    /**
     * Parse FK configuration from CSV file
     * @param filename Input CSV filename (without extension)
     * @return Parsed configuration
     * @throws std::runtime_error if file cannot be read or is malformed
     */
    FKConfiguration parseFromFile(const std::string& filename);

private:
    int parseInteger(const std::string& str);
    double parseDouble(const std::string& str);
    std::vector<std::string> splitLine(const std::string& line, char delimiter = ',');
};

/**
 * FK computation engine
 * Handles the core mathematical computation logic
 */
class FKComputationEngine {
public:
    explicit FKComputationEngine(const FKConfiguration& config);

    /**
     * Compute FK polynomial for given angles
     * @param angles Input angle vector
     * @return Computed polynomial result
     */
    MultivariablePolynomial computeForAngles(const std::vector<int>& angles);
    MultivariablePolynomial computeForAngles_new(const std::vector<int>& angles);

    /**
     * Get the current accumulated result
     */
    const MultivariablePolynomial& getResult() const { return result_; }

    /**
     * Reset computation state
     */
    void reset();

private:
    const FKConfiguration& config_;
    MultivariablePolynomial result_;
    std::vector<int> accumulator_block_sizes_;
    std::vector<std::vector<int>> numerical_assignments_;

    void initializeAccumulatorBlockSizes();
    std::vector<std::vector<int>> computeNumericalAssignments(const std::vector<int>& angles);
    void performCrossingComputations(std::vector<bilvector<int>>& polynomial_terms,
                                    const std::vector<int>& max_x_degrees,
                                    const std::vector<int>& block_sizes);
    void applyCrossingRelation(std::vector<bilvector<int>>& polynomial_terms,
                              int crossing_index, int relation_type);
    void accumulateResult(const std::vector<bilvector<int>>& polynomial_terms,
                         const std::vector<int>& x_power_accumulator,
                         int q_power_accumulator,
                         const std::vector<int>& max_x_degrees);
};

/**
 * Result writer for FK computation
 * Handles output formatting and file writing
 */
class FKResultWriter {
public:
    FKResultWriter() = default;

    /**
     * Write polynomial result to JSON file
     * @param result Polynomial to write
     * @param filename Output filename
     */
    void writeToJson(const MultivariablePolynomial& result, const std::string& filename);

    /**
     * Write polynomial result to human-readable format
     * @param result Polynomial to write
     * @param filename Output filename
     */
    void writeToText(const MultivariablePolynomial& result, const std::string& filename);
};



/**
 * Main FK computation orchestrator
 * Coordinates parsing, computation, and output
 */
class FKComputation {
public:
    FKComputation() = default;

    /**
     * Run complete FK computation from input file to output file
     * @param input_filename Input CSV file (without extension)
     * @param output_filename Output file (without extension)
     */
    void compute(const std::string& input_filename,
                const std::string& output_filename);

    /**
     * Run computation with custom configuration
     * @param config Custom configuration
     * @param output_filename Output file
     */
    void compute(const FKConfiguration& config,
                const std::string& output_filename);

    /**
     * Get the last computed result
     */
    const MultivariablePolynomial& getLastResult() const;

    /**
     * Get the configuration used in the last computation
     */
    const FKConfiguration& getLastConfiguration() const { return config_; }


private:
    FKConfiguration config_;
    std::unique_ptr<FKComputationEngine> engine_;
    FKInputParser parser_;
    FKResultWriter writer_;

    void initializeEngine();
    void runPooledComputation();

    // Pooling functionality - moved from solution_pool_1a_double_links.cpp
    struct EnumerationState {
        std::vector<std::vector<double>> criteria;
        std::list<std::array<int, 2>> bounds;
        std::vector<std::vector<double>> supporting_inequalities;
        std::vector<int> point;
        int current_bound_index;
        int current_value;
        int upper_bound;
    };

    struct VariableAssignmentState {
        std::vector<std::vector<double>> new_criteria;
        std::vector<double> degrees;
        std::vector<std::vector<double>> criteria;
        std::list<std::array<int, 2>> first;
        std::list<std::array<int, 2>> bounds;
        std::vector<std::vector<double>> supporting_inequalities;
        std::vector<int> point;
        size_t current_var_index;
        int current_value;
        int max_value;
    };

    struct BoundedVariables {
        std::vector<int> bounded_v;
        int bounded_count;
        std::list<std::array<int, 2>> first;
    };

    struct ValidatedCriteria {
        std::vector<std::vector<double>> criteria;
        std::vector<double> degrees;
        std::list<std::array<int, 2>> first_bounds;
        std::list<std::array<int, 2>> additional_bounds;
        std::vector<int> initial_point;
        bool is_valid;

        ValidatedCriteria() : is_valid(false) {}
    };

    struct AssignmentResult {
        std::vector<std::vector<double>> criteria;
        std::list<std::array<int, 2>> bounds;
        std::vector<std::vector<double>> supporting_inequalities;
        std::vector<int> point;
    };

    // Private pooling methods
    bool satisfiesConstraints(const std::vector<int>& point,
                             const std::vector<std::vector<double>>& constraints);

    std::vector<std::vector<int>> enumeratePoints(const AssignmentResult& assignment);

    std::vector<AssignmentResult> assignVariables(const ValidatedCriteria& valid_criteria);

    BoundedVariables identifyBoundedVariables(const std::vector<std::vector<double>>& inequalities,
                                             int size);

    std::list<std::array<int, 2>> findAdditionalBounds(
        std::vector<int>& bounded_v,
        int& bounded_count,
        int size,
        const std::vector<std::vector<double>>& supporting_inequalities);

    std::vector<double> extractDegrees(const std::vector<std::vector<double>>& inequalities);

    bool validCriteria(const std::vector<std::vector<double>>& criteria,
                       const std::vector<std::vector<double>>& supporting_inequalities,
                       int size);

    ValidatedCriteria findValidCriteria();

    void pooling(std::vector<std::vector<double>> main_inequalities,
                 std::vector<std::vector<double>> supporting_inequalities,
                 const std::function<void(const std::vector<int>&)>& function);
};

} // namespace fk
