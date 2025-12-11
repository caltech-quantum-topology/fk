#include "fk/fk_computation.hpp"
#include "fk/linalg.hpp"
#include "fk/qalg_links.hpp"
#include "fk/string_to_int.hpp"
#include "graph_search.h"

#include <atomic>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <mutex>
#include <set>
#include <shared_mutex>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fk {
// FKConfiguration implementation
bool FKConfiguration::isValid() const {
  if (degree <= 0 || components <= 0 || crossings < 0 || prefactors < 0) {
    return false;
  }

  if (crossing_matrices.size() != static_cast<size_t>(crossings) ||
      crossing_relation_types.size() != static_cast<size_t>(crossings) ||
      top_crossing_components.size() != static_cast<size_t>(crossings) ||
      bottom_crossing_components.size() != static_cast<size_t>(crossings)) {
    return false;
  }

  if (closed_strand_components.size() != static_cast<size_t>(prefactors)) {
    return false;
  }

  return true;
}

void FKConfiguration::clear() {
  degree = components = writhe = prefactors = crossings = 0;
  closed_strand_components.clear();
  crossing_matrices.clear();
  crossing_relation_types.clear();
  top_crossing_components.clear();
  bottom_crossing_components.clear();
  criteria.clear();
  inequalities.clear();
  variable_assignments.clear();
}

// FKInputParser implementation
FKConfiguration FKInputParser::parseFromFile(const std::string &filename) {
  FKConfiguration config;
  std::ifstream infile(filename + ".csv");

  if (!infile.is_open()) {
    throw std::runtime_error("Unable to open file '" + filename + ".csv'");
  }

  std::string line;

  // Parse degree
  if (!std::getline(infile, line)) {
    throw std::runtime_error("Cannot read degree from file");
  }
  auto parts = splitLine(line);
  if (parts.empty()) {
    throw std::runtime_error("Invalid degree line");
  }
  config.degree = parseInteger(parts[0]);

  // Parse components
  if (!std::getline(infile, line)) {
    throw std::runtime_error("Cannot read components from file");
  }
  parts = splitLine(line);
  if (parts.empty()) {
    throw std::runtime_error("Invalid components line");
  }
  config.components = parseInteger(parts[0]);

  // Parse writhe
  if (!std::getline(infile, line)) {
    throw std::runtime_error("Cannot read writhe from file");
  }
  parts = splitLine(line);
  if (parts.empty()) {
    throw std::runtime_error("Invalid writhe line");
  }
  config.writhe = parseInteger(parts[0]);

  // Parse crossing data
  if (!std::getline(infile, line)) {
    throw std::runtime_error("Cannot read crossing data from file");
  }
  parts = splitLine(line);
  int height = 0;
  for (size_t i = 0; i + 1 < parts.size(); i += 2) {
    int c = parseInteger(parts[i]);
    int relation_type = parseInteger(parts[i + 1]);

    config.crossing_matrices.push_back(
        {{height, c - 1}, {height, c}, {height + 1, c - 1}, {height + 1, c}});
    config.crossing_relation_types.push_back(relation_type);
    height++;
  }
  config.crossings = config.crossing_relation_types.size();

  // Parse closed strand components
  if (!std::getline(infile, line)) {
    throw std::runtime_error("Cannot read closed strand components from file");
  }
  parts = splitLine(line);
  for (const auto &part : parts) {
    if (!part.empty()) {
      config.closed_strand_components.push_back(parseInteger(part));
    }
  }

  config.prefactors = config.closed_strand_components.size();

  // Parse crossing components
  if (!std::getline(infile, line)) {
    throw std::runtime_error("Cannot read crossing components from file");
  }
  parts = splitLine(line);
  for (size_t i = 0; i + 1 < parts.size(); i += 2) {
    config.top_crossing_components.push_back(parseInteger(parts[i]));
    config.bottom_crossing_components.push_back(parseInteger(parts[i + 1]));
  }

  // Parse remaining sections
  int stage = 0;
  while (std::getline(infile, line)) {
    if (line.empty())
      continue;

    if (line[0] == '/') {
      stage++;
      continue;
    }

    parts = splitLine(line);
    if (parts.empty())
      continue;

    if (stage == 0) {
      // Criteria section
      std::vector<double> criteria_row;
      for (const auto &part : parts) {
        if (!part.empty()) {
          criteria_row.push_back(parseDouble(part));
        }
      }
      config.criteria.push_back(criteria_row);
    } else if (stage == 1) {
      // Inequalities section
      std::vector<double> inequality_row;
      for (const auto &part : parts) {
        if (!part.empty()) {
          inequality_row.push_back(parseDouble(part));
        }
      }
      config.inequalities.push_back(inequality_row);
    } else if (stage == 2) {
      // Variable assignments section
      static int extension_index = 0;
      if (config.variable_assignments.empty()) {
        config.variable_assignments.resize(
            config.crossings + 1,
            std::vector<std::vector<int>>(config.prefactors + 1));
      }

      for (const auto &part : parts) {
        if (!part.empty()) {
          int row = extension_index % (config.crossings + 1);
          int col = extension_index / (config.crossings + 1);
          config.variable_assignments[row][col].push_back(parseInteger(part));
        }
      }
      extension_index++;
    }
  }

  if (!config.isValid()) {
    throw std::runtime_error("Parsed configuration is invalid");
  }

  return config;
}

int FKInputParser::parseInteger(const std::string &str) {
  return parseStringToInteger(str);
}

double FKInputParser::parseDouble(const std::string &str) {
  return parseStringToDouble(str);
}

std::vector<std::string> FKInputParser::splitLine(const std::string &line,
                                                  char delimiter) {
  std::vector<std::string> parts;
  std::stringstream ss(line);
  std::string part;

  while (std::getline(ss, part, delimiter)) {
    parts.push_back(part);
  }

  return parts;
}

// FKComputationEngine implementation
FKComputationEngine::FKComputationEngine(const FKConfiguration &config)
    : config_(config), result_(config.components, config.degree) {
  initializeAccumulatorBlockSizes();
  numerical_assignments_.resize(config_.crossings + 1,
                                std::vector<int>(config_.prefactors + 1));
}

void FKComputationEngine::initializeAccumulatorBlockSizes() {
  accumulator_block_sizes_.clear();
  accumulator_block_sizes_.push_back(1);

  for (int i = 1; i < config_.components; i++) {
    accumulator_block_sizes_.push_back(accumulator_block_sizes_[i - 1] *
                                       (config_.degree + 1));
  }
}

PolynomialType
FKComputationEngine::computeForAngles(const std::vector<int> &angles) {
  computeNumericalAssignments(angles);
  //
  // Calculate power accumulators
  // Writhe is the writhe of the link.
  // Prefactors is sum_{closed_strands} sign(strand_closure)*1/2
  double q_power_accumulator_double =
      (config_.writhe / 2. - config_.prefactors) / 2.0;

  std::vector<double> x_power_accumulator_double(config_.components, 0);
  int initial_coefficient = 1;

  // Apply prefactor adjustments
  for (int i = 0; i < config_.prefactors; i++) {
    q_power_accumulator_double -= numerical_assignments_[0][i + 1];
    x_power_accumulator_double[config_.closed_strand_components[i]] -= 0.5;
  }
  // q_power_accumulator_double -= config_.writhe/4.0;
  x_power_accumulator_double[0] -= 0.5;

  // Apply crossing adjustments
  for (int crossing_index = 0; crossing_index < config_.crossings;
       crossing_index++) {
    const auto &crossing_matrix = config_.crossing_matrices[crossing_index];

    int param_i =
        numerical_assignments_[crossing_matrix[0][0]][crossing_matrix[0][1]];
    int param_j =
        numerical_assignments_[crossing_matrix[1][0]][crossing_matrix[1][1]];
    int param_ip =
        numerical_assignments_[crossing_matrix[2][0]][crossing_matrix[2][1]];
    int param_jp =
        numerical_assignments_[crossing_matrix[3][0]][crossing_matrix[3][1]];

    int top_component = config_.top_crossing_components[crossing_index];
    int bottom_component = config_.bottom_crossing_components[crossing_index];
    int relation_type = config_.crossing_relation_types[crossing_index];

    if (relation_type == 1 || relation_type == 2) {
      q_power_accumulator_double +=
          (param_j + param_jp + 0.5) / 2.0 + param_j * param_jp;
      x_power_accumulator_double[top_component] +=
          ((param_j + param_ip + 1) / 4.0);
      x_power_accumulator_double[bottom_component] +=
          ((3 * param_jp - param_i + 1) / 4.0);
    } else if (relation_type == 3 || relation_type == 4) {
      // Canonical contribution from R-matrix
      q_power_accumulator_double +=
          -(param_i + param_ip + 0.5) / 2.0 - param_i * param_ip;
      x_power_accumulator_double[top_component] +=
          -((3 * param_ip - param_j + 1) / 4.0);
      x_power_accumulator_double[bottom_component] +=
          -((param_jp + param_i + 1) / 4.0);

      // Contribution from q binomial reversal
      q_power_accumulator_double += -param_ip * (param_j - param_ip);

      // Contribution from q-pochhammer x-variable reversal
      if (((param_j - param_ip) % 2 == 0 && relation_type == 3) ||
          ((param_j - param_ip) % 2 == 1 && relation_type == 4)) {
        initial_coefficient *= -1;
      }
      x_power_accumulator_double[top_component] += -(param_j - param_ip);
      q_power_accumulator_double +=
          -param_i * (param_j - param_ip) -
          (param_j - param_ip) * (param_j - param_ip + 1) / 2.0;
    }
  }

  // Convert to integer power accumulators
  int q_power_accumulator =
      static_cast<int>(std::floor(q_power_accumulator_double));
  std::vector<int> x_power_accumulator(config_.components);
  std::vector<int> max_x_degrees(config_.components);

  for (int n = 0; n < config_.components; n++) {
    x_power_accumulator[n] =
        static_cast<int>(std::floor(x_power_accumulator_double[n]));
    max_x_degrees[n] = config_.degree - x_power_accumulator[n];
  }

  // Accumulate result
  PolynomialType poly(config_.components, 0);
  poly.setCoefficient(q_power_accumulator, x_power_accumulator,
                      initial_coefficient);
  poly *= crossingFactor(max_x_degrees);
  result_ += poly;

  // Periodic truncation to prevent unbounded growth (every 1000 points)
  static thread_local int points_processed = 0;
  points_processed++;
  if (points_processed % 1000 == 0) {
    result_ = result_.truncate(config_.degree * 4);  // Very conservative bound
  }

  return result_;
}

void FKComputationEngine::computeNumericalAssignments(
    const std::vector<int> &angles) {
  // assume numerical_assignments_ is already sized correctly in ctor/reset
  for (int i = 0; i < config_.crossings + 1; i++) {
    for (int j = 0; j < config_.prefactors + 1; j++) {
      numerical_assignments_[i][j] =
          computeDotProduct(config_.variable_assignments[i][j], angles);
    }
  }
}

PolynomialType
FKComputationEngine::crossingFactor(const std::vector<int> &max_x_degrees) {

  // Create cache key
  CrossingFactorKey cache_key{numerical_assignments_, max_x_degrees};

  // Check cache (read lock)
  {
    std::shared_lock<std::shared_mutex> lock(crossing_factor_mutex_);
    auto it = crossing_factor_cache_.find(cache_key);
    if (it != crossing_factor_cache_.end()) {
      return it->second;  // Cache hit
    }
  }

  PolynomialType result(config_.components, 0);
  result.setCoefficient(0, std::vector<int>(config_.components, 0), 1);

  // Single pass: process each crossing once with both binomial and Pochhammer
  for (int crossing_index = 0; crossing_index < config_.crossings;
       ++crossing_index) {
    const auto &crossing_matrix = config_.crossing_matrices[crossing_index];
    const int relation_type = config_.crossing_relation_types[crossing_index];

    // Extract parameters once
    const int param_i =
        numerical_assignments_[crossing_matrix[0][0]][crossing_matrix[0][1]];
    const int param_j =
        numerical_assignments_[crossing_matrix[1][0]][crossing_matrix[1][1]];
    const int param_ip =
        numerical_assignments_[crossing_matrix[2][0]][crossing_matrix[2][1]];
    const int param_jp =
        numerical_assignments_[crossing_matrix[3][0]][crossing_matrix[3][1]];

    const int top_comp = config_.top_crossing_components[crossing_index];
    const int bottom_comp = config_.bottom_crossing_components[crossing_index];
    // Apply binomial and Pochhammer operations based on relation type
    PolynomialType factor(config_.components, 0);
    switch (relation_type) {
    case 1: {
      // Binomial part
      const QPolynomialType binomial = QBinomial(param_i, param_i - param_jp);

      // Pochhammer part
      const PolynomialType poch(
          qpochhammer_xq_q(param_i - param_jp, param_j + 1), config_.components,
          bottom_comp);
      factor = poch * binomial;
      break;
    }
    case 2: {
      // Binomial part
      auto binomial = QBinomial(param_i, param_jp);

      // Pochhammer part - limit q-power to reduce polynomial size
      const int max_q_power = config_.degree * 2;  // More aggressive bound
      const PolynomialType poch(
          inverse_qpochhammer_xq_q(param_jp - param_i,
                                   param_j - param_jp + param_i + 1,
                                   max_x_degrees[bottom_comp],
                                   max_q_power),
          config_.components, bottom_comp);
      factor = poch * binomial;
      break;
    }
    case 3: {
      // Binomial part
      auto binomial = QBinomial(param_j, param_ip);

      // Pochhammer part - limit q-power to reduce polynomial size
      const int max_q_power = config_.degree * 2;  // More aggressive bound
      const PolynomialType poch(
          inverse_qpochhammer_xq_q(param_ip - param_j,
                                   param_i - param_ip + param_j + 1,
                                   max_x_degrees[top_comp],
                                   max_q_power),
          config_.components, top_comp);
      factor = poch * binomial;
      break;
    }
    case 4: {
      // Binomial part
      const QPolynomialType binomial = QBinomial(param_j, param_j - param_ip);

      // Pochhammer part
      const PolynomialType poch(
          qpochhammer_xq_q(param_j - param_ip, param_i + 1), config_.components,
          top_comp);
      factor = poch * binomial;
      break;
    }
    }
    // Multiply and truncate in one pass to avoid intermediate polynomial growth
    result = result.multiplyAndTruncate(factor, max_x_degrees);
  }

  // Store in cache (write lock) with size limit
  {
    std::unique_lock<std::shared_mutex> lock(crossing_factor_mutex_);

    // Limit cache size to prevent unbounded memory growth for large degrees
    // If cache is full, clear it (simple eviction strategy)
    const size_t MAX_CACHE_SIZE = 10000;
    if (crossing_factor_cache_.size() >= MAX_CACHE_SIZE) {
      crossing_factor_cache_.clear();
    }

    // Use insert instead of [] to avoid needing a default constructor
    crossing_factor_cache_.insert({cache_key, result});
  }

  return result;
}

void FKComputationEngine::reset() {
  result_ = PolynomialType(config_.components, config_.degree);
  for (auto &row : numerical_assignments_) {
    std::fill(row.begin(), row.end(), 0);
  }

  // Clear the crossingFactor cache
  {
    std::unique_lock<std::shared_mutex> lock(crossing_factor_mutex_);
    crossing_factor_cache_.clear();
  }
}

// FKResultWriter implementation
void FKResultWriter::writeToJson(const PolynomialType &result,
                                 const std::string &filename) {
  result.exportToJson(filename);
}

void FKResultWriter::writeToText(const PolynomialType &result,
                                 const std::string &filename) {
  std::ofstream outfile(filename + ".txt");
  if (!outfile.is_open()) {
    throw std::runtime_error("Cannot open output file: " + filename + ".txt");
  }

  outfile << "FK Computation Result\n";
  outfile << "====================\n\n";

  // Redirect cout to capture print output
  std::streambuf *orig = std::cout.rdbuf();
  std::cout.rdbuf(outfile.rdbuf());

  result.print(50); // Print up to 50 terms

  std::cout.rdbuf(orig);
  outfile.close();
}

// FKComputation implementation
void FKComputation::compute(const std::string &input_filename,
                            const std::string &output_filename,
                            int num_threads) {

  FKConfiguration config = parser_.parseFromFile(input_filename);

  // Call the config-based compute method
  compute(config, output_filename, num_threads);
}

void FKComputation::compute(const FKConfiguration &config,
                            const std::string &output_filename,
                            int num_threads) {

  if (!config.isValid()) {
    throw std::runtime_error("Invalid configuration provided");
  }

  if (num_threads < 1) {
    throw std::runtime_error("Number of threads must be at least 1");
  }

  config_ = config;
  std::cout << "Computing to degree: " << config_.degree << std::endl;
  initializeEngine(num_threads);

  // Find valid criteria
  ValidatedCriteria valid_criteria = findValidCriteria();
  if (!valid_criteria.is_valid) {
    throw std::runtime_error("No valid criteria found");
  }

  // Assign variables to get list of variable assignments
  std::vector<AssignmentResult> assignments = assignVariables(valid_criteria);
  std::cout << assignments.size() << " assignments found" << std::endl;

  // Process each assignment's points with callback-based streaming
  // This eliminates point vector storage entirely
  for (size_t assign_idx = 0; assign_idx < assignments.size(); ++assign_idx) {
    const AssignmentResult &assignment = assignments[assign_idx];
    std::cout << "Processing assignment " << (assign_idx + 1) << "/" << assignments.size() << std::endl;

    std::atomic<int> points_processed{0};

    // Callback processes each point immediately as it's enumerated
    auto process_point = [this, &points_processed](const std::vector<int>& point) {
#ifdef _OPENMP
      int thread_id = omp_get_thread_num();
#else
      int thread_id = 0;
#endif
      engines_[thread_id]->computeForAngles(point);

      // Progress tracking
      int count = ++points_processed;
      if (count % 1000 == 0 && thread_id == 0) {
        std::cout << "  Processed " << count << " points" << std::endl;
      }
    };

    enumeratePointsWithCallback(assignment, process_point);
    std::cout << "  Completed " << points_processed.load() << " points" << std::endl;
  }

  // Combine results and perform final computations
  PolynomialType result(config_.components, 0);
  int num_engines = engines_.size();
  for (int engine_idx = 0; engine_idx < num_engines; ++engine_idx) {
    result += engines_[engine_idx]->getResult();
  }

  PolynomialType offset(config_.components, 0);
  std::vector<int> xPowers(config_.components, 0);
  offset.setCoefficient(0, xPowers, -1);
  xPowers[0] += 1;
  offset.setCoefficient(0, xPowers, 1);
  result *= offset;
  result = result.truncate(config_.degree - 1);
  writer_.writeToJson(result, output_filename);
}

const PolynomialType &FKComputation::getLastResult() const {
  if (engines_.empty()) {
    throw std::runtime_error("No computation has been performed yet");
  }
  return engines_[0]->getResult();
}

void FKComputation::initializeEngine(int num_threads) {
  engines_.clear();
  engines_.reserve(num_threads);

  for (int i = 0; i < num_threads; ++i) {
    engines_.push_back(std::make_unique<FKComputationEngine>(config_));
  }
}

// Pooling functionality implementations (moved from
// solution_pool_1a_double_links.cpp)

bool FKComputation::satisfiesConstraints(
    const std::vector<int> &point,
    const std::vector<std::vector<double>> &constraints) {
  for (const auto &constraint : constraints) {
    int acc = static_cast<int>(constraint[0]);
    for (size_t j = 0; j < point.size(); j++) {
      acc += point[j] * static_cast<int>(constraint[1 + j]);
    }
    if (acc < 0) {
      return false;
    }
  }
  return true;
}

std::vector<std::vector<int>>
FKComputation::enumeratePoints(const AssignmentResult &assignment) {

  if (assignment.bounds.empty()) {
    // Base case: check constraints and add point if valid
    if (satisfiesConstraints(assignment.point,
                             assignment.supporting_inequalities) &&
        satisfiesConstraints(assignment.point, assignment.criteria)) {
      return {assignment.point};
    }
    return {};
  }

  // Convert bounds list to vector for easier iteration
  std::vector<std::array<int, 2>> bounds_vec(assignment.bounds.begin(),
                                             assignment.bounds.end());

  // Calculate upper bound for first variable to determine parallelization range
  int first_index = bounds_vec[0][0];
  int first_inequality = bounds_vec[0][1];
  int first_upper =
      static_cast<int>(assignment.supporting_inequalities[first_inequality][0]);
  for (size_t i = 0; i < assignment.point.size(); i++) {
    if (static_cast<int>(i) != first_index) {
      first_upper +=
          static_cast<int>(
              assignment.supporting_inequalities[first_inequality][1 + i]) *
          assignment.point[i];
    }
  }
  first_upper /= -static_cast<int>(
      assignment.supporting_inequalities[first_inequality][1 + first_index]);

  // Parallel region: process each value of the first variable in parallel
  std::vector<std::vector<int>> valid_points;

#ifdef _OPENMP
  std::mutex points_mutex;

#pragma omp parallel for schedule(dynamic)
  for (int first_value = 0; first_value <= first_upper; ++first_value) {
    std::vector<std::vector<int>> local_points = enumeratePointsFromValue(
        assignment, bounds_vec, first_value, first_index);

    // Combine local results into global result
    std::lock_guard<std::mutex> lock(points_mutex);
    valid_points.insert(valid_points.end(), local_points.begin(),
                        local_points.end());
  }
#else
  // Sequential fallback
  for (int first_value = 0; first_value <= first_upper; ++first_value) {
    std::vector<std::vector<int>> local_points = enumeratePointsFromValue(
        assignment, bounds_vec, first_value, first_index);
    valid_points.insert(valid_points.end(), local_points.begin(),
                        local_points.end());
  }
#endif

  return valid_points;
}

std::vector<std::vector<int>> FKComputation::enumeratePointsFromValue(
    const AssignmentResult &assignment,
    const std::vector<std::array<int, 2>> &bounds_vec, int first_value,
    int first_index) {

  std::vector<std::vector<int>> local_points;

  // Stack to manage iteration state
  struct IterationFrame {
    std::vector<int> point;
    size_t bound_index;
    int current_value;
    int upper_bound;
  };

  std::stack<IterationFrame> stack;

  // Initialize first frame with the given first_value
  IterationFrame initial_frame;
  initial_frame.point = assignment.point;
  initial_frame.point[first_index] = first_value;
  initial_frame.bound_index =
      1; // Start from second variable since first is set
  initial_frame.current_value = 0;

  // If there's only one variable, check and return
  if (bounds_vec.size() == 1) {
    if (satisfiesConstraints(initial_frame.point,
                             assignment.supporting_inequalities) &&
        satisfiesConstraints(initial_frame.point, assignment.criteria)) {
      local_points.push_back(initial_frame.point);
    }
    return local_points;
  }

  // Calculate upper bound for second variable
  int second_index = bounds_vec[1][0];
  int second_inequality = bounds_vec[1][1];
  int second_upper = static_cast<int>(
      assignment.supporting_inequalities[second_inequality][0]);
  for (size_t i = 0; i < initial_frame.point.size(); i++) {
    if (static_cast<int>(i) != second_index) {
      second_upper +=
          static_cast<int>(
              assignment.supporting_inequalities[second_inequality][1 + i]) *
          initial_frame.point[i];
    }
  }
  second_upper /= -static_cast<int>(
      assignment.supporting_inequalities[second_inequality][1 + second_index]);
  initial_frame.upper_bound = second_upper;

  stack.push(initial_frame);

  while (!stack.empty()) {
    IterationFrame &current = stack.top();

    if (current.current_value > current.upper_bound) {
      stack.pop();
      continue;
    }

    // Set current variable value
    current.point[bounds_vec[current.bound_index][0]] = current.current_value;

    if (current.bound_index == bounds_vec.size() - 1) {
      // Last variable - check constraints and add point if valid
      if (satisfiesConstraints(current.point,
                               assignment.supporting_inequalities) &&
          satisfiesConstraints(current.point, assignment.criteria)) {
        local_points.push_back(current.point);
      }
      current.current_value++;
    } else {
      // More variables to process
      IterationFrame next_frame;
      next_frame.point = current.point;
      next_frame.bound_index = current.bound_index + 1;
      next_frame.current_value = 0;

      // Calculate upper bound for next variable
      int next_index = bounds_vec[next_frame.bound_index][0];
      int next_inequality = bounds_vec[next_frame.bound_index][1];
      int next_upper = static_cast<int>(
          assignment.supporting_inequalities[next_inequality][0]);
      for (size_t i = 0; i < next_frame.point.size(); i++) {
        if (static_cast<int>(i) != next_index) {
          next_upper +=
              static_cast<int>(
                  assignment.supporting_inequalities[next_inequality][1 + i]) *
              next_frame.point[i];
        }
      }
      next_upper /= -static_cast<int>(
          assignment.supporting_inequalities[next_inequality][1 + next_index]);
      next_frame.upper_bound = next_upper;

      current.current_value++;
      stack.push(next_frame);
    }
  }

  return local_points;
}

void FKComputation::enumeratePointsFromValueWithCallback(
    const AssignmentResult &assignment,
    const std::vector<std::array<int, 2>> &bounds_vec, int first_value,
    int first_index,
    const std::function<void(const std::vector<int>&)>& callback) {

  // Stack to manage iteration state
  struct IterationFrame {
    std::vector<int> point;
    size_t bound_index;
    int current_value;
    int upper_bound;
  };

  std::stack<IterationFrame> stack;

  // Initialize first frame with the given first_value
  IterationFrame initial_frame;
  initial_frame.point = assignment.point;
  initial_frame.point[first_index] = first_value;
  initial_frame.bound_index =
      1; // Start from second variable since first is set
  initial_frame.current_value = 0;

  // If there's only one variable, check and invoke callback
  if (bounds_vec.size() == 1) {
    if (satisfiesConstraints(initial_frame.point,
                             assignment.supporting_inequalities) &&
        satisfiesConstraints(initial_frame.point, assignment.criteria)) {
      callback(initial_frame.point);
    }
    return;
  }

  // Calculate upper bound for second variable
  int second_index = bounds_vec[1][0];
  int second_inequality = bounds_vec[1][1];
  int second_upper = static_cast<int>(
      assignment.supporting_inequalities[second_inequality][0]);
  for (size_t i = 0; i < initial_frame.point.size(); i++) {
    if (static_cast<int>(i) != second_index) {
      second_upper +=
          static_cast<int>(
              assignment.supporting_inequalities[second_inequality][1 + i]) *
          initial_frame.point[i];
    }
  }
  second_upper /= -static_cast<int>(
      assignment.supporting_inequalities[second_inequality][1 + second_index]);
  initial_frame.upper_bound = second_upper;

  stack.push(initial_frame);

  while (!stack.empty()) {
    IterationFrame &current = stack.top();

    if (current.current_value > current.upper_bound) {
      stack.pop();
      continue;
    }

    // Set current variable value
    current.point[bounds_vec[current.bound_index][0]] = current.current_value;

    if (current.bound_index == bounds_vec.size() - 1) {
      // Last variable - check constraints and invoke callback if valid
      if (satisfiesConstraints(current.point,
                               assignment.supporting_inequalities) &&
          satisfiesConstraints(current.point, assignment.criteria)) {
        callback(current.point);
      }
      current.current_value++;
    } else {
      // More variables to process
      IterationFrame next_frame;
      next_frame.point = current.point;
      next_frame.bound_index = current.bound_index + 1;
      next_frame.current_value = 0;

      // Calculate upper bound for next variable
      int next_index = bounds_vec[next_frame.bound_index][0];
      int next_inequality = bounds_vec[next_frame.bound_index][1];
      int next_upper = static_cast<int>(
          assignment.supporting_inequalities[next_inequality][0]);
      for (size_t i = 0; i < next_frame.point.size(); i++) {
        if (static_cast<int>(i) != next_index) {
          next_upper +=
              static_cast<int>(
                  assignment.supporting_inequalities[next_inequality][1 + i]) *
              next_frame.point[i];
        }
      }
      next_upper /= -static_cast<int>(
          assignment.supporting_inequalities[next_inequality][1 + next_index]);
      next_frame.upper_bound = next_upper;

      current.current_value++;
      stack.push(next_frame);
    }
  }
}

void FKComputation::enumeratePointsWithCallback(
    const AssignmentResult &assignment,
    const std::function<void(const std::vector<int>&)>& callback) {

  if (assignment.bounds.empty()) {
    // Base case: check constraints and invoke callback if valid
    if (satisfiesConstraints(assignment.point,
                             assignment.supporting_inequalities) &&
        satisfiesConstraints(assignment.point, assignment.criteria)) {
      callback(assignment.point);
    }
    return;
  }

  // Convert bounds list to vector for easier iteration
  std::vector<std::array<int, 2>> bounds_vec(assignment.bounds.begin(),
                                             assignment.bounds.end());

  // Calculate upper bound for first variable to determine parallelization range
  int first_index = bounds_vec[0][0];
  int first_inequality = bounds_vec[0][1];
  int first_upper =
      static_cast<int>(assignment.supporting_inequalities[first_inequality][0]);
  for (size_t i = 0; i < assignment.point.size(); i++) {
    if (static_cast<int>(i) != first_index) {
      first_upper +=
          static_cast<int>(
              assignment.supporting_inequalities[first_inequality][1 + i]) *
          assignment.point[i];
    }
  }
  first_upper /= -static_cast<int>(
      assignment.supporting_inequalities[first_inequality][1 + first_index]);

  // Parallel region: process each value of the first variable in parallel
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
  for (int first_value = 0; first_value <= first_upper; ++first_value) {
    enumeratePointsFromValueWithCallback(assignment, bounds_vec, first_value,
                                        first_index, callback);
  }
#else
  // Sequential fallback
  for (int first_value = 0; first_value <= first_upper; ++first_value) {
    enumeratePointsFromValueWithCallback(assignment, bounds_vec, first_value,
                                        first_index, callback);
  }
#endif
}

std::vector<FKComputation::AssignmentResult>
FKComputation::assignVariables(const ValidatedCriteria &valid_criteria) {

  std::vector<AssignmentResult> assignments;

  if (valid_criteria.first_bounds.empty()) {
    return {createSingleAssignment(valid_criteria)};
  }

  const auto bounds_vector = convertBoundsToVector(valid_criteria.first_bounds);

  using AssignmentNode = VariableAssignmentState;

  graph_search::DepthFirstSearchWithVisited<AssignmentNode>::Config dfs_config;

  // We no longer use "goal nodes" to emit assignments. All work happens
  // in get_neighbors, so we disable goal handling here.
  dfs_config.is_valid_goal = [](const AssignmentNode &) {
    return false;
  };

  dfs_config.get_neighbors =
      [this, &bounds_vector, &assignments](const AssignmentNode &current_state) {
        std::vector<AssignmentNode> neighbors;

        // If this state is exhausted, it has no neighbors,
        // mirroring the old "if (isStateExhausted) { continue; }" behavior.
        if (isStateExhausted(current_state)) {
          return neighbors;
        }

        // Enumerate all still-available values for the current variable.
        for (int value = current_state.current_value;
             value <= current_state.max_value; ++value) {
          auto state_copy = current_state;
          state_copy.current_value = value;

          // This matches the old loop body:
          //   processCurrentVariable(...)
          //   updated_degrees = calculateUpdatedDegrees(...)
          processCurrentVariable(state_copy, bounds_vector);
          const auto updated_degrees =
              calculateUpdatedDegrees(state_copy, bounds_vector);

          if (!isLastVariable(state_copy, bounds_vector)) {
            // Non-last variable: create the next state (next variable),
            // as in the old createNextState(...) + stack.push(next_state).
            auto next_state =
                createNextState(state_copy, updated_degrees, bounds_vector);
            neighbors.push_back(next_state);
          } else {
            // Last variable: instead of pushing a neighbor, we emit the
            // assignment directly, just like the old
            //   if (isLastVariable) { assignments.push_back(...); }
            assignments.push_back(createAssignmentResult(state_copy));
          }
        }

        return neighbors;
      };

  // No-op: all work is done in get_neighbors.
  dfs_config.process_node =
      [](const AssignmentNode &) {};

  graph_search::DepthFirstSearchWithVisited<AssignmentNode> dfs_search(
      dfs_config);

  auto initial_state =
      createInitialAssignmentState(valid_criteria, bounds_vector);
  dfs_search.search_all(initial_state);

  return assignments;
}


/*
std::vector<FKComputation::AssignmentResult>
FKComputation::assignVariables(const ValidatedCriteria &valid_criteria) {

  std::vector<AssignmentResult> assignments;

  if (valid_criteria.first_bounds.empty()) {
    return {createSingleAssignment(valid_criteria)};
  }

  const auto bounds_vector = convertBoundsToVector(valid_criteria.first_bounds);

  using AssignmentNode = VariableAssignmentState;

  graph_search::DepthFirstSearchWithVisited<AssignmentNode>::Config dfs_config;

  dfs_config.is_valid_goal = [this,
                              &bounds_vector](const AssignmentNode &state) {
    return isLastVariable(state, bounds_vector) && !isStateExhausted(state);
  };

  dfs_config.get_neighbors =
      [this, &bounds_vector](const AssignmentNode &current_state) {
        std::vector<AssignmentNode> neighbors;

        for (int value = current_state.current_value;
             value <= current_state.max_value; ++value) {
          auto state_copy = current_state;
          state_copy.current_value = value;
          processCurrentVariable(state_copy, bounds_vector);
          const auto updated_degrees =
              calculateUpdatedDegrees(state_copy, bounds_vector);

          if (!isLastVariable(state_copy, bounds_vector)) {
            auto next_state =
                createNextState(state_copy, updated_degrees, bounds_vector);
            neighbors.push_back(next_state);
          } else {
            neighbors.push_back(state_copy);
          }
        }

        return neighbors;
      };

  dfs_config.process_node = [&assignments,
                             this](const AssignmentNode &final_state) {
    assignments.push_back(createAssignmentResult(final_state));
  };

  graph_search::DepthFirstSearchWithVisited<AssignmentNode> dfs_search(
      dfs_config);

  auto initial_state =
      createInitialAssignmentState(valid_criteria, bounds_vector);
  dfs_search.search_all(initial_state);

  return assignments;
}
*/

FKComputation::AssignmentResult
FKComputation::createSingleAssignment(const ValidatedCriteria &valid_criteria) {
  AssignmentResult result;
  result.criteria = valid_criteria.criteria;
  result.bounds = valid_criteria.additional_bounds;
  result.supporting_inequalities = config_.inequalities;
  result.point = valid_criteria.initial_point;
  return result;
}

std::vector<std::array<int, 2>> FKComputation::convertBoundsToVector(
    const std::list<std::array<int, 2>> &bounds) {
  return std::vector<std::array<int, 2>>(bounds.begin(), bounds.end());
}

FKComputation::VariableAssignmentState
FKComputation::createInitialAssignmentState(
    const ValidatedCriteria &valid_criteria,
    const std::vector<std::array<int, 2>> &bounds_vector) {

  VariableAssignmentState initial_state;
  initial_state.new_criteria = valid_criteria.criteria;
  initial_state.degrees = valid_criteria.degrees;
  initial_state.criteria = valid_criteria.criteria;
  initial_state.bounds = valid_criteria.additional_bounds;
  auto all_original = config_.inequalities;
  all_original.insert(all_original.end(), config_.criteria.begin(),
                      config_.criteria.end());
  initial_state.supporting_inequalities = all_original;
  initial_state.point = valid_criteria.initial_point;
  initial_state.current_var_index = 0;
  initial_state.current_value = 0;
  initial_state.max_value = calculateMaxValue(
      valid_criteria.criteria, valid_criteria.degrees, bounds_vector[0]);

  return initial_state;
}

std::stack<FKComputation::VariableAssignmentState>
FKComputation::initializeAssignmentStack(
    const ValidatedCriteria &valid_criteria,
    const std::vector<std::array<int, 2>> &bounds_vector) {
  std::stack<VariableAssignmentState> stack;

  VariableAssignmentState initial_state;
  initial_state.new_criteria = valid_criteria.criteria;
  initial_state.degrees = valid_criteria.degrees;
  initial_state.criteria = valid_criteria.criteria;
  initial_state.bounds = valid_criteria.additional_bounds;
  auto all_original = config_.inequalities;
  all_original.insert(all_original.end(), config_.criteria.begin(),
                      config_.criteria.end());
  initial_state.supporting_inequalities = all_original;
  initial_state.point = valid_criteria.initial_point;
  initial_state.current_var_index = 0;
  initial_state.current_value = 0;
  initial_state.max_value = calculateMaxValue(
      valid_criteria.criteria, valid_criteria.degrees, bounds_vector[0]);

  stack.push(initial_state);
  return stack;
}

int FKComputation::calculateMaxValue(
    const std::vector<std::vector<double>> &criteria,
    const std::vector<double> &degrees, const std::array<int, 2> &bound) {
  const int var_index = bound[0];
  const int constraint_index = bound[1];
  const double slope = -criteria[constraint_index][var_index];
  return static_cast<int>(degrees[constraint_index] / slope);
}

bool FKComputation::isStateExhausted(const VariableAssignmentState &state) {
  return state.current_value > state.max_value;
}

void FKComputation::processCurrentVariable(
    VariableAssignmentState &state,
    const std::vector<std::array<int, 2>> &bounds_vector) {
  const int var_index = bounds_vector[state.current_var_index][0];
  state.point[var_index - 1] = state.current_value;
}

std::vector<double> FKComputation::calculateUpdatedDegrees(
    const VariableAssignmentState &state,
    const std::vector<std::array<int, 2>> &bounds_vector) {
  const int var_index = bounds_vector[state.current_var_index][0];
  const int constraint_index = bounds_vector[state.current_var_index][1];
  const double slope = -state.new_criteria[constraint_index][var_index];

  auto updated_degrees = state.degrees;
  updated_degrees[constraint_index] =
      state.degrees[constraint_index] - state.current_value * slope;

  return updated_degrees;
}

bool FKComputation::isLastVariable(
    const VariableAssignmentState &state,
    const std::vector<std::array<int, 2>> &bounds_vector) {
  return state.current_var_index == bounds_vector.size() - 1;
}

FKComputation::AssignmentResult
FKComputation::createAssignmentResult(const VariableAssignmentState &state) {
  AssignmentResult result;
  result.criteria = state.criteria;
  result.bounds = state.bounds;
  result.supporting_inequalities = state.supporting_inequalities;
  result.point = state.point;
  return result;
}

FKComputation::VariableAssignmentState FKComputation::createNextState(
    const VariableAssignmentState &current_state,
    const std::vector<double> &updated_degrees,
    const std::vector<std::array<int, 2>> &bounds_vector) {
  VariableAssignmentState next_state;
  next_state.new_criteria = current_state.new_criteria;
  next_state.degrees = updated_degrees;
  next_state.criteria = current_state.criteria;
  next_state.bounds = current_state.bounds;
  next_state.supporting_inequalities = current_state.supporting_inequalities;
  next_state.point = current_state.point;
  next_state.current_var_index = current_state.current_var_index + 1;
  next_state.current_value = 0;
  next_state.max_value =
      calculateMaxValue(next_state.new_criteria, updated_degrees,
                        bounds_vector[next_state.current_var_index]);

  return next_state;
}

FKComputation::BoundedVariables FKComputation::identifyBoundedVariables(
    const std::vector<std::vector<double>> &inequalities, int size) {
  BoundedVariables result;
  result.bounded_v.resize(size - 1, 0);
  result.bounded_count = 0;

  int mains = inequalities.size();
  for (int i = 0; i < mains; i++) {
    bool condition = true;
    std::vector<bool> locally_bounded(size - 1, false);

    for (int k = 1; k < size; k++) {
      if (inequalities[i][k] > 0) {
        condition = false;
        break;
      } else if (inequalities[i][k] < 0) {
        locally_bounded[k - 1] = true;
      }
    }

    if (condition) {
      for (int v = 0; v < size - 1; v++) {
        if (locally_bounded[v] && !result.bounded_v[v]) {
          result.bounded_v[v] = true;
          result.first.push_back({v + 1, i});
          result.bounded_count++;
        }
      }
    }
  }

  return result;
}

std::list<std::array<int, 2>> FKComputation::findAdditionalBounds(
    std::vector<int> &bounded_v, int &bounded_count, int size,
    const std::vector<std::vector<double>> &supporting_inequalities) {

  std::list<std::array<int, 2>> bounds;
  int support = supporting_inequalities.size();

  int index = 0;
  while (index < size - 1) {
    if (!bounded_v[index]) {
      for (int l = 0; l < support; l++) {
        if (supporting_inequalities[l][1 + index] < 0) {
          bool useful = true;
          for (int n = 0; n < size - 1; n++) {
            if (n != index && supporting_inequalities[l][1 + n] > 0 &&
                !bounded_v[n]) {
              useful = false;
              break;
            }
          }
          if (useful) {
            bounds.push_back({index, l});
            bounded_v[index] = true;
            bounded_count++;
            if (bounded_count == size - 1) {
              return bounds;
            }
            index = -1;
            break;
          }
        }
      }
    }
    index++;
  }
  return bounds;
}

std::vector<double> FKComputation::extractDegrees(
    const std::vector<std::vector<double>> &inequalities) {
  std::vector<double> degrees;
  for (const auto &x : inequalities) {
    degrees.push_back(x[0]);
  }
  return degrees;
}

bool FKComputation::validCriteria(
    const std::vector<std::vector<double>> &criteria,
    const std::vector<std::vector<double>> &supporting_inequalities, int size) {
  auto bounded_info = identifyBoundedVariables(criteria, size);
  if (bounded_info.bounded_count == 0) {
    return false; // Not valid - no bounded variables
  }

  if (bounded_info.bounded_count == size - 1) {
    return true; // We have enough bounded variables
  }

  // Try to find additional bounds
  auto additional_bounds =
      findAdditionalBounds(bounded_info.bounded_v, bounded_info.bounded_count,
                           size, supporting_inequalities);
  return (bounded_info.bounded_count ==
          size - 1); // Valid if we now have enough
}

FKComputation::ValidatedCriteria FKComputation::findValidCriteria() {
  if (config_.criteria.empty()) {
    return ValidatedCriteria();
  }

  const int variable_count = config_.criteria[0].size();

  // Try initial criteria first (fast path)
  auto initial_result = tryInitialCriteria(variable_count);
  if (initial_result.is_valid) {
    return initial_result;
  }

  // Search for valid criteria using breadth-first exploration
  return searchForValidCriteria(variable_count);
}

FKComputation::ValidatedCriteria
FKComputation::tryInitialCriteria(int variable_count) {
  if (validCriteria(config_.criteria, config_.inequalities, variable_count)) {
    return buildValidatedCriteriaFromValid(config_.criteria, variable_count);
  }
  return ValidatedCriteria();
}

FKComputation::ValidatedCriteria FKComputation::buildValidatedCriteriaFromValid(
    const std::vector<std::vector<double>> &criteria, int variable_count) {
  ValidatedCriteria result;
  const auto bounded_info = identifyBoundedVariables(criteria, variable_count);

  result.criteria = criteria;
  result.degrees = extractDegrees(config_.criteria);
  result.first_bounds = bounded_info.first;
  result.initial_point = std::vector<int>(variable_count - 1, 0);
  result.is_valid = true;

  // Find additional bounds if needed
  if (bounded_info.bounded_count < variable_count - 1) {
    auto all_original = config_.inequalities;
    all_original.insert(all_original.end(), config_.criteria.begin(),
                        config_.criteria.end());

    auto bounded_v_copy = bounded_info.bounded_v;
    auto bounded_count_copy = bounded_info.bounded_count;
    result.additional_bounds = findAdditionalBounds(
        bounded_v_copy, bounded_count_copy, variable_count, all_original);
  }

  return result;
}

FKComputation::ValidatedCriteria
FKComputation::searchForValidCriteria(int variable_count) {
  using CriteriaNode = std::vector<std::vector<double>>;

  // Set to track visited criteria combinations
  std::set<CriteriaNode> visited;

  // Configure BFS search
  graph_search::BreadthFirstSearch<CriteriaNode>::Config bfs_config;

  bfs_config.is_valid_goal = [this,
                              variable_count](const CriteriaNode &criteria) {
    return validCriteria(criteria, combineInequalitiesAndCriteria(),
                         variable_count);
  };

  bfs_config.get_neighbors =
      [this, variable_count](const CriteriaNode &current_criteria) {
        return generateCriteriaNeighbors(current_criteria, variable_count);
      };

  bfs_config.has_been_visited = [&visited](const CriteriaNode &criteria) {
    return visited.find(criteria) != visited.end();
  };

  bfs_config.mark_visited = [&visited](const CriteriaNode &criteria) {
    visited.insert(criteria);
  };

  graph_search::BreadthFirstSearch<CriteriaNode> bfs_search(bfs_config);

  auto result = bfs_search.search(config_.criteria);
  if (result.has_value()) {
    return buildValidatedCriteriaFromValid(result.value(), variable_count);
  }

  return ValidatedCriteria();
}

std::vector<std::vector<std::vector<double>>>
FKComputation::generateCriteriaNeighbors(
    const std::vector<std::vector<double>> &current_criteria,
    int /* variable_count */) {
  std::vector<std::vector<std::vector<double>>> neighbors;

  const int criterion_count = config_.criteria.size();
  const auto combined_inequalities = combineInequalitiesAndCriteria();

  for (int criterion_index = 0; criterion_index < criterion_count;
       ++criterion_index) {
    for (const auto &inequality : combined_inequalities) {
      if (shouldExploreCombination(current_criteria, inequality,
                                   criterion_index)) {
        const auto new_criteria = createCombinedCriteria(
            current_criteria, inequality, criterion_index);
        neighbors.push_back(new_criteria);
      }
    }
  }

  return neighbors;
}

std::vector<std::vector<double>>
FKComputation::combineInequalitiesAndCriteria() const {
  std::vector<std::vector<double>> combined = config_.inequalities;
  combined.insert(combined.end(), config_.criteria.begin(),
                  config_.criteria.end());
  return combined;
}

bool FKComputation::isPotentiallyBeneficial(
    const std::vector<double> &criterion,
    const std::vector<double> &inequality) const {
  // Check if this combination could be beneficial (has opposing signs)
  for (size_t var_index = 1; var_index < criterion.size(); ++var_index) {
    if (criterion[var_index] > 0 && inequality[var_index] < 0) {
      return true;
    }
  }
  return false;
}

std::vector<std::vector<double>> FKComputation::createCombinedCriteria(
    const std::vector<std::vector<double>> &base_criteria,
    const std::vector<double> &inequality, int criterion_index) const {

  auto combined_criteria = base_criteria;
  for (size_t var_index = 0; var_index < inequality.size(); ++var_index) {
    combined_criteria[criterion_index][var_index] +=
        inequality[var_index] * COMBINATION_FACTOR;
  }

  return combined_criteria;
}

bool FKComputation::shouldExploreCombination(
    const std::vector<std::vector<double>> &current_criteria,
    const std::vector<double> &inequality, int criterion_index) const {

  // Check if this combination could be beneficial
  return isPotentiallyBeneficial(current_criteria[criterion_index], inequality);
}

void FKComputation::setupWorkStealingComputation(
    const std::vector<std::vector<int>> &all_points) {
  int total_points = all_points.size();

#ifdef _OPENMP
  int num_engines = engines_.size();
  omp_set_num_threads(num_engines);
  std::cout << "Processing " << total_points << " points with " << num_engines
            << " threads using OpenMP" << std::endl;
#else
  std::cout << "Processing " << total_points
            << " points sequentially (OpenMP not available)" << std::endl;
#endif

  // Use OpenMP parallel for to distribute work across threads
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < total_points; ++i) {
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
    engines_[thread_id]->computeForAngles(all_points[i]);

    // Print progress occasionally from thread 0 only to avoid race conditions
    if (thread_id == 0 && i % 100 == 0) {
      std::cout << "Processing point " << i << "/" << total_points << std::endl;
    }
#else
    engines_[0]->computeForAngles(all_points[i]);

    // Print progress for sequential execution
    if (i % 100 == 0) {
      std::cout << "Processing point " << i << "/" << total_points << std::endl;
    }
#endif
  }

  std::cout << "All points processed!" << std::endl;
}

void FKComputation::combineEngineResults() {
  int num_engines = engines_.size();
  std::cout << "Combining results from " << num_engines << " engines..."
            << std::endl;

  for (int engine_idx = 1; engine_idx < num_engines; ++engine_idx) {
    engines_[0]->getResult() += engines_[engine_idx]->getResult();
  }
}

void FKComputation::performFinalOffsetComputation() {
  // Final offset addition: multiply by (-1 + x₀)
  // This computes the discrete derivative in the x₀ direction

  PolynomialType offset(config_.components, 0);
  std::vector<int> xPowers(config_.components, 0);

  // Coefficient for x₀⁰: -1
  offset.setCoefficient(0, xPowers, -1);

  // Coefficient for x₀¹: +1
  xPowers[0] = 1;
  offset.setCoefficient(0, xPowers, 1);

  // Multiply result by offset polynomial
  const_cast<PolynomialType &>(engines_[0]->getResult()) *= offset;

  // Truncate to final degree
  std::vector<int> maxima(config_.components, config_.degree - 1);
  const_cast<PolynomialType &>(engines_[0]->getResult()) =
      engines_[0]->getResult().truncate(maxima);
}

} // namespace fk
