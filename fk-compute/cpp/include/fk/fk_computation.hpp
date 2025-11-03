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
 * Checkpoint state for FK computation
 * Holds all necessary data to resume computation from a saved point
 */
struct CheckpointState {
    int resume_count;                              // Number of times computation has been resumed
    FKConfiguration config;                        // Configuration needed for engine recreation
    std::vector<std::vector<int>> remaining_points; // Points not yet processed
    MultivariablePolynomial accumulated_result;    // Current accumulated result
    int total_points_count;                        // Total points for progress tracking

    CheckpointState() : resume_count(0), config(), remaining_points(),
                       accumulated_result(1, 10), total_points_count(0) {}
    CheckpointState(const FKConfiguration& cfg, int total_count)
        : resume_count(0), config(cfg), remaining_points(),
          accumulated_result(cfg.components, cfg.degree), total_points_count(total_count) {}
};

/**
 * Checkpoint manager for saving and loading computation state
 * Handles file I/O and state serialization
 */
class CheckpointManager {
public:
    /**
     * Constructor for creating a new checkpoint
     * @param state Initial checkpoint state to save
     */
    explicit CheckpointManager(const CheckpointState& state);

    /**
     * Constructor for loading existing checkpoint
     * @param checkpoint_path Path to existing checkpoint file
     */
    explicit CheckpointManager(const std::string& checkpoint_path);

    /**
     * Save current checkpoint state
     * @param state State to save
     * @return true if successful, false otherwise
     */
    bool save(const CheckpointState& state);

    /**
     * Load checkpoint state
     * @return Loaded checkpoint state
     */
    CheckpointState load();

    /**
     * Get the checkpoint filename
     * @return Checkpoint filename
     */
    const std::string& getFilename() const { return checkpoint_filename_; }

private:
    std::string checkpoint_filename_;
    bool is_existing_checkpoint_;

    std::string generateUniqueFilename();
    void saveToJson(const CheckpointState& state, const std::string& filename);
    CheckpointState loadFromJson(const std::string& filename);
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

    /**
     * Run FK computation with periodic checkpointing
     * @param input_filename Input CSV file (without extension)
     * @param output_filename Output file (without extension)
     * @param checkpoint_period How often to save checkpoints (number of points processed)
     * @return Path to checkpoint file for potential resumption
     */
    std::string computeWithCheckpointing(const std::string& input_filename,
                                       const std::string& output_filename,
                                       int checkpoint_period = 1000);

    /**
     * Resume FK computation from checkpoint
     * @param checkpoint_path Path to checkpoint file
     * @param output_filename Output file (without extension)
     */
    void resumeFromCheckpoint(const std::string& checkpoint_path,
                            const std::string& output_filename);

private:
    FKConfiguration config_;
    std::unique_ptr<FKComputationEngine> engine_;
    FKInputParser parser_;
    FKResultWriter writer_;

    void initializeEngine();
    void runPooledComputation();

    // Checkpoint helper methods
    CheckpointState saveCheckpoint(const std::vector<std::vector<int>>& remaining_points,
                                 int total_points) const;
    void loadCheckpoint(const CheckpointState& state);

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
