#pragma once

#include "fk/polynomial_config.hpp"
#include <array>
#include <functional>
#include <list>
#include <memory>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <vector>
#include <unordered_map>
#include <shared_mutex>

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
  FKConfiguration parseFromFile(const std::string &filename);

private:
  int parseInteger(const std::string &str);
  double parseDouble(const std::string &str);
  std::vector<std::string> splitLine(const std::string &line,
                                     char delimiter = ',');
};

/**
 * FK computation engine
 * Handles the core mathematical computation logic
 */
class FKComputationEngine {
public:
  explicit FKComputationEngine(const FKConfiguration &config);

  /**
   * Compute FK polynomial for given angles
   * @param angles Input angle vector
   * @return Computed polynomial result
   */
  PolynomialType computeForAngles(const std::vector<int> &angles);

  /**
   * Get the current accumulated result
   */
  const PolynomialType &getResult() const { return result_; }
  PolynomialType &getResult() { return result_; }

  /**
   * Reset computation state
   */
  void reset();

private:
  const FKConfiguration &config_;
  PolynomialType result_;
  std::vector<int> accumulator_block_sizes_;
  std::vector<std::vector<int>> numerical_assignments_;

  // Cache for crossingFactor results
  struct CrossingFactorKey {
    std::vector<std::vector<int>> numerical_assignments;
    std::vector<int> max_x_degrees;

    bool operator==(const CrossingFactorKey& other) const {
      return numerical_assignments == other.numerical_assignments &&
             max_x_degrees == other.max_x_degrees;
    }
  };

  struct CrossingFactorKeyHash {
    std::size_t operator()(const CrossingFactorKey& k) const {
      std::size_t h1 = hash_vector(k.max_x_degrees);
      std::size_t h2 = 0;
      for (const auto& vec : k.numerical_assignments) {
        h2 ^= hash_vector(vec) + 0x9e3779b9 + (h2 << 6) + (h2 >> 2);
      }
      return h1 ^ (h2 << 1);
    }

  private:
    std::size_t hash_vector(const std::vector<int>& vec) const {
      std::size_t seed = vec.size();
      for(auto& i : vec) {
        seed ^= std::hash<int>()(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
  };

  mutable std::unordered_map<CrossingFactorKey, PolynomialType, CrossingFactorKeyHash> crossing_factor_cache_;
  mutable std::shared_mutex crossing_factor_mutex_;

  void initializeAccumulatorBlockSizes();
  void
  computeNumericalAssignments(const std::vector<int> &angles);

  PolynomialType
  crossingFactor(const std::vector<int> &max_x_degrees); 

  void accumulateResultPoly(const PolynomialType &poly,
                            const std::vector<int> &x_power_accumulator,
                            int q_power_accumulator);
  void performOffsetAdditionPoly(const PolynomialType &source_poly,
                                 const std::vector<int> &x_offset, int q_offset,
                                 int sign_multiplier);
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
  void writeToJson(const PolynomialType &result,
                   const std::string &filename);

  /**
   * Write polynomial result to human-readable format
   * @param result Polynomial to write
   * @param filename Output filename
   */
  void writeToText(const PolynomialType &result,
                   const std::string &filename);
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
   * @param num_threads Number of threads/engines to use (default: 1)
   */
  void compute(const std::string &input_filename,
               const std::string &output_filename,
               int num_threads = 1);

  /**
   * Run computation with custom configuration
   * @param config Custom configuration
   * @param output_filename Output file
   * @param num_threads Number of threads/engines to use (default: 1)
   */
  void compute(const FKConfiguration &config,
               const std::string &output_filename,
               int num_threads = 1);

  /**
   * Get the last computed result
   */
  const PolynomialType &getLastResult() const;

  /**
   * Get the configuration used in the last computation
   */
  const FKConfiguration &getLastConfiguration() const { return config_; }

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

private:
  FKConfiguration config_;
  std::vector<std::unique_ptr<FKComputationEngine>> engines_;
  FKInputParser parser_;
  FKResultWriter writer_;

  void initializeEngine(int num_threads);
  void setupWorkStealingComputation(const std::vector<std::vector<int>> &all_points);
  void combineEngineResults();
  void performFinalOffsetComputation();

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

    bool operator<(const VariableAssignmentState& other) const {
      if (point != other.point) return point < other.point;
      if (current_var_index != other.current_var_index) return current_var_index < other.current_var_index;
      return current_value < other.current_value;
    }

    bool operator==(const VariableAssignmentState& other) const {
      return point == other.point &&
             current_var_index == other.current_var_index &&
             current_value == other.current_value;
    }
  };

  // Private pooling methods
  bool
  satisfiesConstraints(const std::vector<int> &point,
                       const std::vector<std::vector<double>> &constraints);

  std::vector<std::vector<int>>
  enumeratePoints(const AssignmentResult &assignment);

  std::vector<std::vector<int>>
  enumeratePointsFromValue(const AssignmentResult &assignment,
                          const std::vector<std::array<int, 2>> &bounds_vec,
                          int first_value,
                          int first_index);

  std::vector<AssignmentResult>
  assignVariables(const ValidatedCriteria &valid_criteria);

  BoundedVariables
  identifyBoundedVariables(const std::vector<std::vector<double>> &inequalities,
                           int size);

  std::list<std::array<int, 2>> findAdditionalBounds(
      std::vector<int> &bounded_v, int &bounded_count, int size,
      const std::vector<std::vector<double>> &supporting_inequalities);

  std::vector<double>
  extractDegrees(const std::vector<std::vector<double>> &inequalities);

  bool
  validCriteria(const std::vector<std::vector<double>> &criteria,
                const std::vector<std::vector<double>> &supporting_inequalities,
                int size);

  ValidatedCriteria findValidCriteria();

  void pooling(std::vector<std::vector<double>> main_inequalities,
               std::vector<std::vector<double>> supporting_inequalities,
               const std::function<void(const std::vector<int> &)> &function);

  // Helper functions for assignVariables refactoring
  AssignmentResult createSingleAssignment(const ValidatedCriteria &valid_criteria);

  std::vector<std::array<int, 2>> convertBoundsToVector(const std::list<std::array<int, 2>> &bounds);

  VariableAssignmentState createInitialAssignmentState(const ValidatedCriteria &valid_criteria,
                                                       const std::vector<std::array<int, 2>> &bounds_vector);

  std::stack<VariableAssignmentState> initializeAssignmentStack(const ValidatedCriteria &valid_criteria,
                                                               const std::vector<std::array<int, 2>> &bounds_vector);

  int calculateMaxValue(const std::vector<std::vector<double>> &criteria,
                       const std::vector<double> &degrees,
                       const std::array<int, 2> &bound);

  bool isStateExhausted(const VariableAssignmentState &state);

  void processCurrentVariable(VariableAssignmentState &state,
                             const std::vector<std::array<int, 2>> &bounds_vector);

  std::vector<double> calculateUpdatedDegrees(const VariableAssignmentState &state,
                                            const std::vector<std::array<int, 2>> &bounds_vector);

  bool isLastVariable(const VariableAssignmentState &state,
                     const std::vector<std::array<int, 2>> &bounds_vector);

  AssignmentResult createAssignmentResult(const VariableAssignmentState &state);

  VariableAssignmentState createNextState(const VariableAssignmentState &current_state,
                                         std::vector<double> &&updated_degrees,
                                         const std::vector<std::array<int, 2>> &bounds_vector);

  // Constants for criteria validation
  static constexpr double COMBINATION_FACTOR = 0.5;


  // Core helper functions for findValidCriteria
  ValidatedCriteria tryInitialCriteria(int variable_count);
  ValidatedCriteria searchForValidCriteria(int variable_count);
  ValidatedCriteria buildValidatedCriteriaFromValid(const std::vector<std::vector<double>>& criteria, int variable_count);

  // Criteria exploration helpers
  std::vector<std::vector<std::vector<double>>> generateCriteriaNeighbors(
    const std::vector<std::vector<double>>& current_criteria, int variable_count);

  // Criteria combination and validation helpers
  std::vector<std::vector<double>> combineInequalitiesAndCriteria() const;
  std::vector<std::vector<double>> createCombinedCriteria(const std::vector<std::vector<double>>& base_criteria,
                                                         const std::vector<double>& inequality,
                                                         int criterion_index) const;
  bool isPotentiallyBeneficial(const std::vector<double>& criterion,
                              const std::vector<double>& inequality) const;
  bool shouldExploreCombination(const std::vector<std::vector<double>>& current_criteria,
                               const std::vector<double>& inequality,
                               int criterion_index) const;
};

} // namespace fk
