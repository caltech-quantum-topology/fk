#include "fk/fk_computation.hpp"
#include "fk/linalg.hpp"
#include "fk/polynomial_config.hpp"
#include "fk/qalg_links.hpp"
#include "fk/string_to_int.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

void printCrossingTerm(QPolynomialType binom, PolynomialType poch) {
  std::cout << "Binomial: ";
  binom.print();
  std::cout << "XQPochhammer: ";
  poch.print();
  auto ct = binom * poch;
  std::cout << "Contribution: ";
  (binom * poch).print(100000);
  if (ct.isZero()) {
    std::cout << "\t\t\t\t\t\t" << "Contribution is zero!!!!" << std::endl;
  }
}

void printTerms(const std::vector<Term> &terms) {
  std::cout << "=== Terms (" << terms.size() << ") ===\n";

  for (const auto &[xPowers, bilvec] : terms) {
    // Print the x-powers
    std::cout << "x-powers: [";
    for (size_t i = 0; i < xPowers.size(); ++i) {
      std::cout << xPowers[i];
      if (i + 1 < xPowers.size())
        std::cout << ", ";
    }
    std::cout << "]\n";

    // Delegate to bilvector's own print function
    std::cout << "  q-bilvector: ";
    bilvec.print();
    std::cout << "\n";
  }
}

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
  numerical_assignments_ = computeNumericalAssignments(angles);

  // Calculate power accumulators
  double q_power_accumulator_double =
      (config_.writhe - config_.prefactors) / 2.0;
  std::vector<double> x_power_accumulator_double(config_.components, 0);
  int initial_coefficient = 1;

  // Apply prefactor adjustments
  for (int i = 0; i < config_.prefactors; i++) {
    q_power_accumulator_double -= numerical_assignments_[0][i + 1];
    x_power_accumulator_double[config_.closed_strand_components[i]] -= 0.5;
  }
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
          (param_j + param_jp) / 2.0 + param_j * param_jp;
      x_power_accumulator_double[top_component] +=
          ((param_j + param_ip + 1) / 4.0);
      x_power_accumulator_double[bottom_component] +=
          ((3 * param_jp - param_i + 1) / 4.0);
    } else if (relation_type == 3 || relation_type == 4) {
      q_power_accumulator_double +=
          -(param_i + param_ip) / 2.0 - (param_j - param_ip) * param_i -
          (param_j - param_ip) * (param_j - param_ip + 1) - param_i * param_i;
      /*
        (param_i + param_ip + param_jp * (param_jp + 1) -
         param_i * (param_i + 1)) /
            2.0 +
        param_i * param_ip;
      */
      x_power_accumulator_double[top_component] -=
          ((3 * param_j - param_ip + 1) / 4.0);
      x_power_accumulator_double[bottom_component] -=
          ((param_i + param_jp + 1) / 4.0);

      if ((relation_type == 4) || (relation_type == 3)) {
        initial_coefficient *= -1;
      }
      /*
      if ((relation_type == 4 && (param_j - param_ip) % 2 == 0) ||
          (relation_type == 3 && (param_ip - param_j) % 2 == 1)) {
        initial_coefficient *= -1;
      }
      */
    }
  }

  // Convert to integer power accumulators
  int q_power_accumulator =
      static_cast<int>(std::floor(q_power_accumulator_double));
  std::vector<int> x_power_accumulator(config_.components);
  std::vector<int> max_x_degrees(config_.components);
  std::vector<int> block_sizes(config_.components);

  for (int n = 0; n < config_.components; n++) {
    x_power_accumulator[n] =
        static_cast<int>(std::floor(x_power_accumulator_double[n]));
    max_x_degrees[n] = config_.degree - x_power_accumulator[n];
  }

  PolynomialType poly(config_.components, 0, max_x_degrees);
  poly.setCoefficient(0, std::vector<int>(config_.components, 0),
                      initial_coefficient);

  // Perform crossing computations

  std::cout << "========================" << std::endl;
  std::cout << "========================" << std::endl;
  std::cout << "Angle: ";
  for (auto it = angles.begin(); it != angles.end(); ++it) {
    std::cout << *it;
    if (std::next(it) != angles.end()) {
      std::cout << ", ";
    } else {
      std::cout << std::endl;
    }
  }

  std::cout << "Max x degrees: ";
  for (auto it = max_x_degrees.begin(); it != max_x_degrees.end(); ++it) {
    std::cout << *it;
    if (std::next(it) != max_x_degrees.end()) {
      std::cout << ", ";
    } else {
      std::cout << std::endl;
    }
  }

  std::cout << "========================" << std::endl;

  poly *= crossingFactor(max_x_degrees);

  // Accumulate result
  // performOffsetAdditionPoly(poly, x_power_accumulator, q_power_accumulator,
  // 1);

  PolynomialType offset(config_.components);
  offset.setCoefficient(q_power_accumulator, x_power_accumulator, 1);
  std::cout << "Offset: ";
  offset.print();
  poly *= offset;

  std::cout << "Result: ";
  poly.print(50);
  std::cout << "========================" << std::endl;
  std::cout << "========================" << std::endl;
  std::cout << "\n\n\n";

  result_ += poly;
  return result_;
}

std::vector<std::vector<int>> FKComputationEngine::computeNumericalAssignments(
    const std::vector<int> &angles) {

  std::vector<std::vector<int>> assignments(
      config_.crossings + 1, std::vector<int>(config_.prefactors + 1));
  for (int i = 0; i < config_.crossings + 1; i++) {
    for (int j = 0; j < config_.prefactors + 1; j++) {
      assignments[i][j] =
          computeDotProduct(config_.variable_assignments[i][j], angles);
    }
  }
  return assignments;
}

PolynomialType
FKComputationEngine::crossingFactor(const std::vector<int> &max_x_degrees) {
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
    std::cout << "Relation Type: " << relation_type << std::endl;
    std::cout << "i: " << param_i << "\tj: " << param_j << std::endl;
    std::cout << "ip: " << param_ip << "\tjp: " << param_jp << std::endl;
    switch (relation_type) {
    case 1: {
      // Binomial part
      const QPolynomialType binomial =
          (param_i > 0) ? QBinomialPositive(param_i, param_i - param_jp)
                        : QBinomialNegative(param_i, param_i - param_jp);
      result *= binomial;

      // Pochhammer part
      const PolynomialType poch(
          qpochhammer_xq_q(param_i - param_jp, param_j+1, 1), config_.components,
          bottom_comp);
      /*
      // Binomial part
      const QPolynomialType binomial =
          (param_i > 0) ? QBinomialPositive(param_i, param_i - param_jp)
                        : QBinomialNegative(param_i, param_i - param_jp);
      result *= binomial;

      // Pochhammer part
      const PolynomialType poch(
          qpochhammer_xq_q(param_ip - param_j, param_j + 1, 1),
          config_.components, bottom_comp);
      */
      printCrossingTerm(binomial, poch);
      result *= poch;
      break;
    }
    case 2: {
      /*
      // Binomial part
      const QPolynomialType binomial = QBinomialNegative(param_i, param_jp);
      result *= binomial;

      // Pochhammer part
      const PolynomialType poch(
          inverse_qpochhammer_xq_q(param_jp - param_i, param_j,
                                   max_x_degrees[bottom_comp] + 1, -1, true),
          config_.components, bottom_comp);
      */
      // Binomial part
      const QPolynomialType binomial = QBinomialNegative(param_i, param_jp);
      result *= binomial;

      // Pochhammer part
      const PolynomialType poch(
          inverse_qpochhammer_xq_q(param_jp-param_i,param_i+param_j-param_jp+1,
                                   max_x_degrees[bottom_comp] + 1, 1, false),
          config_.components, bottom_comp);
      printCrossingTerm(binomial, poch);
      result *= poch;
      break;
    }
    case 3: {
      // Binomial part
      const QPolynomialType binomial =
          QBinomialNegative(param_j, param_ip).invertExponents();
      result *= binomial;

      // Pochhammer part
      const PolynomialType poch(
          inverse_qpochhammer_xq_q(param_ip - param_j,
                                   param_i - param_ip + param_j + 1,
                                   max_x_degrees[top_comp] + 1, 1, false),
          config_.components, top_comp);
      /*
      // Binomial part
      const QPolynomialType binomial =
          QBinomialNegative(param_j, param_ip).invertExponents();
      result *= binomial;

      // Pochhammer part
      const PolynomialType poch(
          inverse_qpochhammer_xq_q(param_i - param_jp, param_jp + 1,
                                   max_x_degrees[top_comp] + 1, 1, true),
          config_.components, top_comp);

      */
      printCrossingTerm(binomial, poch);
      result *= poch;
      break;
    }
    case 4: {
      /*
      // Binomial part
      const QPolynomialType binomial =
          (param_j > 0)
              ? QBinomialPositive(param_j, param_j - param_ip)
              : QBinomialNegative(param_j, param_j - param_ip);

      result *= binomial;

      // Pochhammer part
      const PolynomialType poch(
          qpochhammer_xq_q(param_j - param_ip, param_j, -1), config_.components,
          bottom_comp);
      */
      // Binomial part
      const QPolynomialType binomial =
          (param_j > 0)
              ? QBinomialPositive(param_j, param_j - param_ip).invertExponents()
              : QBinomialNegative(param_j, param_j - param_ip)
                    .invertExponents();

      result *= binomial;

      // Pochhammer part
      const PolynomialType poch(
          qpochhammer_xq_q(param_j - param_ip, param_i + 1, 1),
          config_.components, top_comp);
      printCrossingTerm(binomial, poch);
      result *= poch;
      break;
    }
    }
    std::cout << "========================" << std::endl;
  }
  result = result.truncate(max_x_degrees);
  std::cout << "QAlg res: ";
  result.print(100000);
  return result;
}

void FKComputationEngine::performOffsetAdditionPoly(
    const PolynomialType &source_poly, const std::vector<int> &x_offset,
    int q_offset, int sign_multiplier) {

  // Iterate through all coefficients in the source polynomial
  const auto coeffs = source_poly.getCoefficients();

  for (const auto &[x_powers, q_poly] : coeffs) {
    // Calculate the offset x-powers
    std::vector<int> offset_x_powers = x_powers;
    for (size_t i = 0; i < x_powers.size(); ++i) {
      offset_x_powers[i] += x_offset[i];
    }

    // Iterate through all q-powers in the bilvector
    for (int q_power = q_poly.getMaxNegativeIndex();
         q_power <= q_poly.getMaxPositiveIndex(); ++q_power) {

      int coeff = q_poly[q_power];
      if (coeff != 0) {
        // Apply q-offset and add to result polynomial
        result_.addToCoefficient(q_power + q_offset, offset_x_powers,
                                 sign_multiplier * coeff);
      }
    }
  }
}

void FKComputationEngine::accumulateResultPoly(
    const PolynomialType &poly, const std::vector<int> &x_power_accumulator,
    int q_power_accumulator) {

  performOffsetAdditionPoly(poly, x_power_accumulator, q_power_accumulator, 1);
}

void FKComputationEngine::reset() {
  result_ = PolynomialType(config_.components, config_.degree);
  for (auto &row : numerical_assignments_) {
    std::fill(row.begin(), row.end(), 0);
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
  initializeEngine(num_threads);

  // Find valid criteria
  ValidatedCriteria valid_criteria = findValidCriteria();
  if (!valid_criteria.is_valid) {
    throw std::runtime_error("No valid criteria found");
  }

  // Assign variables to get list of variable assignments
  std::vector<AssignmentResult> assignments = assignVariables(valid_criteria);

  // Collect all points from each variable assignment
  std::vector<std::vector<int>> all_points;
  for (const AssignmentResult &assignment : assignments) {
    auto points = enumeratePoints(assignment);
    all_points.insert(all_points.end(), points.begin(), points.end());
  }

  // Execute work stealing computation
  setupWorkStealingComputation(all_points);

  // Combine results and perform final computations
  combineEngineResults();
  performFinalOffsetComputation();

  writer_.writeToJson(engines_[0]->getResult(), output_filename);
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

  std::vector<std::vector<int>> valid_points;

  if (assignment.bounds.empty()) {
    // Base case: check constraints and add point if valid
    if (satisfiesConstraints(assignment.point,
                             assignment.supporting_inequalities) &&
        satisfiesConstraints(assignment.point, assignment.criteria)) {
      valid_points.push_back(assignment.point);
    }
    return valid_points;
  }

  // Convert bounds list to vector for easier iteration
  std::vector<std::array<int, 2>> bounds_vec(assignment.bounds.begin(),
                                             assignment.bounds.end());

  // Stack to manage iteration state
  struct IterationFrame {
    std::vector<int> point;
    size_t bound_index;
    int current_value;
    int upper_bound;
  };

  std::stack<IterationFrame> stack;

  // Initialize first frame
  IterationFrame initial_frame;
  initial_frame.point = assignment.point;
  initial_frame.bound_index = 0;
  initial_frame.current_value = 0;

  // Calculate upper bound for first variable
  int index = bounds_vec[0][0];
  int inequality = bounds_vec[0][1];
  int upper =
      static_cast<int>(assignment.supporting_inequalities[inequality][0]);
  for (size_t i = 0; i < assignment.point.size(); i++) {
    if (static_cast<int>(i) != index) {
      upper += static_cast<int>(
                   assignment.supporting_inequalities[inequality][1 + i]) *
               assignment.point[i];
    }
  }
  upper /= -static_cast<int>(
      assignment.supporting_inequalities[inequality][1 + index]);
  initial_frame.upper_bound = upper;

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
        valid_points.push_back(current.point);
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
          next_upper =
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

  return valid_points;
}

std::vector<FKComputation::AssignmentResult>
FKComputation::assignVariables(const ValidatedCriteria &valid_criteria) {

  std::vector<AssignmentResult> assignments;

  if (valid_criteria.first_bounds.empty()) {
    AssignmentResult result;
    result.criteria = valid_criteria.criteria;
    result.bounds = valid_criteria.additional_bounds;
    result.supporting_inequalities = config_.inequalities;
    result.point = valid_criteria.initial_point;
    assignments.push_back(result);
    return assignments;
  }

  // Convert first list to vector for easier iteration
  std::vector<std::array<int, 2>> first_vec(valid_criteria.first_bounds.begin(),
                                            valid_criteria.first_bounds.end());

  std::stack<VariableAssignmentState> stack;

  // Initialize first state
  VariableAssignmentState initial_state;
  initial_state.new_criteria = valid_criteria.criteria;
  initial_state.degrees = valid_criteria.degrees;
  initial_state.criteria = valid_criteria.criteria;
  initial_state.bounds = valid_criteria.additional_bounds;
  initial_state.supporting_inequalities = config_.inequalities;
  initial_state.point = valid_criteria.initial_point;
  initial_state.current_var_index = 0;
  initial_state.current_value = 0;

  // Calculate max value for first variable
  int var_index = first_vec[0][0];
  int main_index = first_vec[0][1];
  double slope = -valid_criteria.criteria[main_index][var_index];
  initial_state.max_value =
      static_cast<int>(valid_criteria.degrees[main_index] / slope);

  stack.push(initial_state);

  while (!stack.empty()) {
    VariableAssignmentState &current = stack.top();

    if (current.current_value > current.max_value) {
      stack.pop();
      continue;
    }

    // Set current variable value
    int var_idx = first_vec[current.current_var_index][0];
    int main_idx = first_vec[current.current_var_index][1];
    current.point[var_idx - 1] = current.current_value;

    // Update degrees
    double slope = -current.new_criteria[main_idx][var_idx];
    std::vector<double> new_degrees = current.degrees;
    new_degrees[main_idx] =
        current.degrees[main_idx] - current.current_value * slope;

    if (current.current_var_index == first_vec.size() - 1) {
      // Last variable - create assignment result
      AssignmentResult result;
      result.criteria = current.criteria;
      result.bounds = current.bounds;
      result.supporting_inequalities = current.supporting_inequalities;
      result.point = current.point;
      assignments.push_back(result);
      current.current_value++;
    } else {
      // More variables to assign
      VariableAssignmentState next_state;
      next_state.new_criteria = current.new_criteria;
      next_state.degrees = new_degrees;
      next_state.criteria = current.criteria;
      next_state.bounds = current.bounds;
      next_state.supporting_inequalities = current.supporting_inequalities;
      next_state.point = current.point;
      next_state.current_var_index = current.current_var_index + 1;
      next_state.current_value = 0;

      // Calculate max value for next variable
      int next_var_idx = first_vec[next_state.current_var_index][0];
      int next_main_idx = first_vec[next_state.current_var_index][1];
      double next_slope = -next_state.new_criteria[next_main_idx][next_var_idx];
      next_state.max_value =
          static_cast<int>(new_degrees[next_main_idx] / next_slope);

      current.current_value++;
      stack.push(next_state);
    }
  }

  return assignments;
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
    return ValidatedCriteria(); // Return invalid criteria
  }

  const int mains = config_.criteria.size();
  const int size = config_.criteria[0].size();
  const double COMBINATION_FACTOR = 0.5; // Was hardcoded /2.0

  // Helper function to build ValidatedCriteria from valid criteria
  auto buildValidatedCriteria =
      [&](const std::vector<std::vector<double>> &criteria)
      -> ValidatedCriteria {
    ValidatedCriteria result;
    auto bounded_info = identifyBoundedVariables(criteria, size);

    result.criteria = criteria;
    result.degrees = extractDegrees(config_.criteria);
    result.first_bounds = bounded_info.first;
    result.initial_point = std::vector<int>(size - 1, 0);
    result.is_valid = true;

    // Add additional bounds if needed
    if (bounded_info.bounded_count < size - 1) {
      result.additional_bounds = findAdditionalBounds(
          bounded_info.bounded_v, bounded_info.bounded_count, size,
          config_.inequalities);
    }
    return result;
  };

  // Combine all inequalities for processing
  std::vector<std::vector<double>> all_inequalities = config_.inequalities;
  all_inequalities.insert(all_inequalities.end(), config_.criteria.begin(),
                          config_.criteria.end());

  // Use a single set to track visited criterion configurations more efficiently
  std::set<std::vector<std::vector<double>>> visited_criteria;

  // Queue for breadth-first exploration of criterion space
  std::queue<std::vector<std::vector<double>>> criteria_queue;

  // Start with the main inequalities
  visited_criteria.insert(config_.criteria);
  criteria_queue.push(config_.criteria);

  // Check initial criteria first
  if (validCriteria(config_.criteria, config_.inequalities, size)) {
    return buildValidatedCriteria(config_.criteria);
  }

  // Search for satisfactory criteria
  while (!criteria_queue.empty()) {
    auto current_criteria = std::move(criteria_queue.front());
    criteria_queue.pop();

    // Try combining each criterion with each inequality
    for (int criterion_idx = 0; criterion_idx < mains; ++criterion_idx) {
      for (const auto &inequality : all_inequalities) {

        // Check if this combination could be beneficial (has opposing signs)
        bool potentially_beneficial = false;
        for (int var_idx = 1; var_idx < size; ++var_idx) {
          if (current_criteria[criterion_idx][var_idx] > 0 &&
              inequality[var_idx] < 0) {
            potentially_beneficial = true;
            break;
          }
        }

        if (!potentially_beneficial)
          continue;

        // Create new criterion by linear combination
        auto new_criteria = current_criteria;
        for (int var_idx = 0; var_idx < size; ++var_idx) {
          new_criteria[criterion_idx][var_idx] +=
              inequality[var_idx] * COMBINATION_FACTOR;
        }

        // Skip if we've seen this configuration before
        if (visited_criteria.find(new_criteria) != visited_criteria.end()) {
          continue;
        }

        // Mark as visited
        visited_criteria.insert(new_criteria);

        // Check if these criteria are valid
        if (validCriteria(new_criteria, config_.inequalities, size)) {
          return buildValidatedCriteria(new_criteria); // Found valid criteria
        }

        // Add to queue for further exploration
        criteria_queue.push(std::move(new_criteria));
      }
    }
  }

  // No valid criteria found
  return ValidatedCriteria();
}

void FKComputation::setupWorkStealingComputation(
    const std::vector<std::vector<int>> &all_points) {
  int total_points = all_points.size();
  int num_engines = engines_.size();

#ifdef _OPENMP
  omp_set_num_threads(num_engines);
  std::cout << "Processing " << total_points << " points with " << num_engines
            << " threads using OpenMP" << std::endl;
#else
  std::cout << "Processing " << total_points
            << " points sequentially (OpenMP not available)" << std::endl;
#endif

  // Use OpenMP parallel for to distribute work across threads
#pragma omp parallel for schedule(dynamic)
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
    auto engine_coeffs = engines_[engine_idx]->getResult().getCoefficients();

    // Add this engine's results to the first engine's results
    for (const auto &[x_powers, q_poly] : engine_coeffs) {
      for (int q_power = q_poly.getMaxNegativeIndex();
           q_power <= q_poly.getMaxPositiveIndex(); ++q_power) {
        int coeff = q_poly[q_power];
        if (coeff != 0) {
          const_cast<PolynomialType &>(engines_[0]->getResult())
              .addToCoefficient(q_power, x_powers, coeff);
        }
      }
    }
  }
}

void FKComputation::performFinalOffsetComputation() {
  // Final offset addition
  std::vector<int> increment_offset(config_.components, 0);
  increment_offset[0] = 1;
  std::vector<int> maxima(config_.components, config_.degree - 1);

  auto combined_coeffs1 = engines_[0]->getResult().getCoefficients();
  auto combined_coeffs2 = engines_[0]->getResult().getCoefficients();

  std::vector<int> accumulator_block_sizes;
  accumulator_block_sizes.push_back(1);
  for (int i = 1; i < config_.components; i++) {
    accumulator_block_sizes.push_back(accumulator_block_sizes[i - 1] *
                                      (config_.degree + 1));
  }

  int components_copy = config_.components; // Make a mutable copy

  performOffsetAddition(combined_coeffs1, combined_coeffs2, increment_offset, 0,
                        -1,
                        components_copy, // dimensions
                        maxima           // arrayLengths
  );

  const_cast<PolynomialType &>(engines_[0]->getResult())
      .syncFromSparseVector(combined_coeffs1);
}

} // namespace fk
