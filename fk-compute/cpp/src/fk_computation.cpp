#include "fk/fk_computation.hpp"
#include "fk/linalg.hpp"
#include "fk/qalg_links.hpp"
#include "fk/solution_pool_1a_double_links.hpp"
#include "fk/string_to_int.hpp"

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

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

MultivariablePolynomial
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
    int param_k =
        numerical_assignments_[crossing_matrix[2][0]][crossing_matrix[2][1]];
    int param_m =
        numerical_assignments_[crossing_matrix[3][0]][crossing_matrix[3][1]];

    int top_component = config_.top_crossing_components[crossing_index];
    int bottom_component = config_.bottom_crossing_components[crossing_index];
    int relation_type = config_.crossing_relation_types[crossing_index];

    if (relation_type == 1 || relation_type == 2) {
      q_power_accumulator_double +=
          (param_j + param_m) / 2.0 + param_j * param_m;
      x_power_accumulator_double[top_component] +=
          ((param_j + param_k + 1) / 4.0);
      x_power_accumulator_double[bottom_component] +=
          ((3 * param_m - param_i + 1) / 4.0);
    } else if (relation_type == 3 || relation_type == 4) {
      q_power_accumulator_double -=
          (param_i + param_k + param_m * (param_m + 1) -
           param_i * (param_i + 1)) /
              2.0 +
          param_i * param_k;
      x_power_accumulator_double[top_component] -=
          ((3 * param_j - param_k + 1) / 4.0);
      x_power_accumulator_double[bottom_component] -=
          ((param_i + param_m + 1) / 4.0);

      if ((relation_type == 4 && (param_j - param_k) % 2 == 0) ||
          (relation_type == 3 && (param_k - param_j) % 2 == 1)) {
        initial_coefficient *= -1;
      }
    }
  }

  // Convert to integer power accumulators
  int q_power_accumulator =
      static_cast<int>(std::floor(q_power_accumulator_double));
  std::vector<int> x_power_accumulator(config_.components);
  std::vector<int> max_x_degrees(config_.components);
  std::vector<int> block_sizes(config_.components);

  block_sizes[0] = 1;
  for (int n = 0; n < config_.components; n++) {
    x_power_accumulator[n] = static_cast<int>(x_power_accumulator_double[n]);
    max_x_degrees[n] = config_.degree - x_power_accumulator[n];
    if (n != 0) {
      block_sizes[n] = (max_x_degrees[n - 1] + 1) * block_sizes[n - 1];
    }
  }

  int total_product_size = block_sizes[config_.components - 1] *
                           (max_x_degrees[config_.components - 1] + 1);

  // Initialize polynomial terms
  std::vector<bilvector<int>> polynomial_terms(total_product_size,
                                               bilvector<int>(0, 1, 20, 0));
  polynomial_terms[0][0] = initial_coefficient;

  // Perform crossing computations
  performCrossingComputations(polynomial_terms, max_x_degrees, block_sizes);

  // Accumulate result
  accumulateResult(polynomial_terms, x_power_accumulator, q_power_accumulator, max_x_degrees);

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

void FKComputationEngine::performCrossingComputations(
    std::vector<bilvector<int>> &polynomial_terms,
    const std::vector<int>& max_x_degrees,
    const std::vector<int>& block_sizes) {
  // First pass: binomial computations
  for (int crossing_index = 0; crossing_index < config_.crossings;
       crossing_index++) {
    const auto &crossing_matrix = config_.crossing_matrices[crossing_index];
    int relation_type = config_.crossing_relation_types[crossing_index];

    int param_i =
        numerical_assignments_[crossing_matrix[0][0]][crossing_matrix[0][1]];
    int param_j =
        numerical_assignments_[crossing_matrix[1][0]][crossing_matrix[1][1]];
    int param_k =
        numerical_assignments_[crossing_matrix[2][0]][crossing_matrix[2][1]];
    int param_m =
        numerical_assignments_[crossing_matrix[3][0]][crossing_matrix[3][1]];

    if (relation_type == 1) {
      if (param_i > 0) {
        computePositiveQBinomial(polynomial_terms, param_i, param_i - param_m,
                                 false);
      } else {
        computeNegativeQBinomial(polynomial_terms, param_i, param_i - param_m,
                                 false);
      }
    } else if (relation_type == 2) {
      computeNegativeQBinomial(polynomial_terms, param_i, param_m, false);
    } else if (relation_type == 3) {
      computeNegativeQBinomial(polynomial_terms, param_j, param_k, true);
    } else if (relation_type == 4) {
      if (param_j > 0) {
        computePositiveQBinomial(polynomial_terms, param_j, param_j - param_k,
                                 true);
      } else {
        computeNegativeQBinomial(polynomial_terms, param_j, param_j - param_k,
                                 true);
      }
    }
  }

  // Second pass: Pochhammer computations

  for (int crossing_index = 0; crossing_index < config_.crossings;
       crossing_index++) {
    const auto &crossing_matrix = config_.crossing_matrices[crossing_index];
    int relation_type = config_.crossing_relation_types[crossing_index];

    int param_i =
        numerical_assignments_[crossing_matrix[0][0]][crossing_matrix[0][1]];
    int param_j =
        numerical_assignments_[crossing_matrix[1][0]][crossing_matrix[1][1]];
    int param_k =
        numerical_assignments_[crossing_matrix[2][0]][crossing_matrix[2][1]];
    int param_m =
        numerical_assignments_[crossing_matrix[3][0]][crossing_matrix[3][1]];

    int top_comp = config_.top_crossing_components[crossing_index];
    int bottom_comp = config_.bottom_crossing_components[crossing_index];

    if (relation_type == 1) {
      computeXQPochhammer(polynomial_terms, param_k, param_j + 1, bottom_comp,
                          config_.components, max_x_degrees, block_sizes);
    } else if (relation_type == 2) {
      computeXQInversePochhammer(polynomial_terms, param_j, param_k + 1,
                                 bottom_comp, config_.components, max_x_degrees,
                                 block_sizes);
    } else if (relation_type == 3) {
      computeXQInversePochhammer(polynomial_terms, param_i, param_m + 1,
                                 top_comp, config_.components, max_x_degrees,
                                 block_sizes);
    } else if (relation_type == 4) {
      computeXQPochhammer(polynomial_terms, param_m, param_i + 1, top_comp,
                          config_.components, max_x_degrees, block_sizes);
    }
  }
}

void FKComputationEngine::accumulateResult(
    const std::vector<bilvector<int>> &polynomial_terms,
    const std::vector<int> &x_power_accumulator, int q_power_accumulator,
    const std::vector<int> &max_x_degrees) {
  auto result_coeffs = result_.getCoefficients();
  std::vector<int> x_power_accumulator_copy =
      x_power_accumulator;                  // Make a mutable copy
  int components_copy = config_.components; // Make a mutable copy
  std::vector<bilvector<int>> polynomial_terms_copy =
      polynomial_terms; // Make a mutable copy

  performOffsetAddition(result_coeffs, polynomial_terms_copy,
                        x_power_accumulator_copy, q_power_accumulator,
                        components_copy, max_x_degrees, 1,
                        accumulator_block_sizes_, accumulator_block_sizes_);

  result_.syncFromDenseVector(result_coeffs);
}

void FKComputationEngine::reset() {
  result_ = MultivariablePolynomial(config_.components, config_.degree);
  for (auto &row : numerical_assignments_) {
    std::fill(row.begin(), row.end(), 0);
  }
}

// FKResultWriter implementation
void FKResultWriter::writeToJson(const MultivariablePolynomial &result,
                                 const std::string &filename) {
  result.exportToJson(filename);
}

void FKResultWriter::writeToText(const MultivariablePolynomial &result,
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
                            const std::string &output_filename) {
  FKConfiguration config = parser_.parseFromFile(input_filename);

  // Call the config-based compute method
  compute(config, output_filename);
}

void FKComputation::compute(const FKConfiguration &config,
                            const std::string &output_filename) {
  if (!config.isValid()) {
    throw std::runtime_error("Invalid configuration provided");
  }

  config_ = config;

  initializeEngine();

  std::function<void(const std::vector<int> &)> computation_function =
      [this](const std::vector<int> &angles) {
        engine_->computeForAngles(angles);
      };

  pooling(config_.criteria, config_.inequalities, computation_function);

  // Final offset addition (from original implementation)
  std::vector<int> increment_offset(config_.components, 0);
  increment_offset[0] = 1;
  std::vector<int> maxima(config_.components, config_.degree - 1);

  auto result_coeffs1 = engine_->getResult().getCoefficients();
  auto result_coeffs2 = engine_->getResult().getCoefficients();

  std::vector<int> accumulator_block_sizes;
  accumulator_block_sizes.push_back(1);
  for (int i = 1; i < config_.components; i++) {
    accumulator_block_sizes.push_back(accumulator_block_sizes[i - 1] *
                                      (config_.degree + 1));
  }

  int components_copy = config_.components; // Make a mutable copy

  performOffsetAddition(result_coeffs1, result_coeffs2, increment_offset, 0,
                        components_copy, maxima, -1, accumulator_block_sizes,
                        accumulator_block_sizes);

  const_cast<MultivariablePolynomial &>(engine_->getResult())
      .syncFromDenseVector(result_coeffs1);

  writer_.writeToJson(engine_->getResult(), output_filename);
}

const MultivariablePolynomial &FKComputation::getLastResult() const {
  if (!engine_) {
    throw std::runtime_error("No computation has been performed yet");
  }
  return engine_->getResult();
}

void FKComputation::initializeEngine() {
  engine_ = std::make_unique<FKComputationEngine>(config_);
}
} // namespace fk
