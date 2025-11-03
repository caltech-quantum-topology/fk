#include "fk/fk_computation.hpp"
#include "fk/linalg.hpp"
#include "fk/qalg_links.hpp"
#include "fk/string_to_int.hpp"

#include <cmath>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <list>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stack>
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
  std::cout<<"Computing numerical assignments"<<std::endl;
  numerical_assignments_ = computeNumericalAssignments(angles);
  std::cout<<"Numerical assignments computed"<<std::endl;

  // Calculate power accumulators
  double q_power_accumulator_double =
      (config_.writhe - config_.prefactors) / 2.0;
  std::vector<double> x_power_accumulator_double(config_.components, 0);
  int initial_coefficient = 1;
  std::cout<<"Power accumulators computed"<<std::endl;

  // Apply prefactor adjustments
  for (int i = 0; i < config_.prefactors; i++) {
    q_power_accumulator_double -= numerical_assignments_[0][i + 1];
    x_power_accumulator_double[config_.closed_strand_components[i]] -= 0.5;
  }
  x_power_accumulator_double[0] -= 0.5;

  std::cout<<"Prefactor adjustments made"<<std::endl;
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

  std::cout<<"Crossing adjustments made"<<std::endl;

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

MultivariablePolynomial
FKComputationEngine::computeForAngles_new(const std::vector<int> &angles) {
  std::cout<<"Computing numerical assignments"<<std::endl;
  numerical_assignments_ = computeNumericalAssignments(angles);
  std::cout<<"Numerical assignments computed"<<std::endl;

  // Calculate power accumulators
  double q_power_accumulator_double =
      (config_.writhe - config_.prefactors) / 2.0;
  std::vector<double> x_power_accumulator_double(config_.components, 0);
  int initial_coefficient = 1;
  std::cout<<"Power accumulators computed"<<std::endl;

  // Apply prefactor adjustments
  for (int i = 0; i < config_.prefactors; i++) {
    q_power_accumulator_double -= numerical_assignments_[0][i + 1];
    x_power_accumulator_double[config_.closed_strand_components[i]] -= 0.5;
  }
  x_power_accumulator_double[0] -= 0.5;

  std::cout<<"Prefactor adjustments made"<<std::endl;
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

  std::cout<<"Crossing adjustments made"<<std::endl;

  // Convert to integer power accumulators
  int q_power_accumulator =
      static_cast<int>(std::floor(q_power_accumulator_double));
  std::vector<int> x_power_accumulator(config_.components);
  std::vector<int> max_x_degrees(config_.components);

  for (int n = 0; n < config_.components; n++) {
    x_power_accumulator[n] = static_cast<int>(x_power_accumulator_double[n]);
    max_x_degrees[n] = config_.degree - x_power_accumulator[n];
  }

  MultivariablePolynomial polynomial(config_.components, config_.degree, max_x_degrees);
  std::vector<int> zeroDeg(config_.components, 0);
  polynomial.setCoefficient(0,zeroDeg,1);
  polynomial.print();
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

// CheckpointManager implementation
CheckpointManager::CheckpointManager(const CheckpointState& state)
    : is_existing_checkpoint_(false) {
  checkpoint_filename_ = generateUniqueFilename();
  if (!save(state)) {
    throw std::runtime_error("Failed to create initial checkpoint: " + checkpoint_filename_);
  }
}

CheckpointManager::CheckpointManager(const std::string& checkpoint_path)
    : checkpoint_filename_(checkpoint_path), is_existing_checkpoint_(true) {
  // Verify file exists
  std::ifstream file(checkpoint_path);
  if (!file.is_open()) {
    throw std::runtime_error("Checkpoint file not found: " + checkpoint_path);
  }
  file.close();
}

bool CheckpointManager::save(const CheckpointState& state) {
  try {
    saveToJson(state, checkpoint_filename_);
    return true;
  } catch (const std::exception& e) {
    std::cerr << "Error saving checkpoint: " << e.what() << std::endl;
    return false;
  }
}

CheckpointState CheckpointManager::load() {
  return loadFromJson(checkpoint_filename_);
}

std::string CheckpointManager::generateUniqueFilename() {
  // Generate unique filename using timestamp and random component
  auto now = std::chrono::system_clock::now();
  auto time_t = std::chrono::system_clock::to_time_t(now);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(1000, 9999);

  std::stringstream ss;
  ss << "fk_checkpoint_" << time_t << "_" << dis(gen) << ".json";
  return ss.str();
}

void CheckpointManager::saveToJson(const CheckpointState& state, const std::string& filename) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot create checkpoint file: " + filename);
  }

  file << "{\n";
  file << "  \"resume_count\": " << state.resume_count << ",\n";
  file << "  \"total_points_count\": " << state.total_points_count << ",\n";

  // Save configuration
  file << "  \"config\": {\n";
  file << "    \"degree\": " << state.config.degree << ",\n";
  file << "    \"components\": " << state.config.components << ",\n";
  file << "    \"writhe\": " << state.config.writhe << ",\n";
  file << "    \"prefactors\": " << state.config.prefactors << ",\n";
  file << "    \"crossings\": " << state.config.crossings << "\n";
  file << "  },\n";

  // Save remaining points
  file << "  \"remaining_points\": [\n";
  for (size_t i = 0; i < state.remaining_points.size(); ++i) {
    file << "    [";
    for (size_t j = 0; j < state.remaining_points[i].size(); ++j) {
      file << state.remaining_points[i][j];
      if (j < state.remaining_points[i].size() - 1) file << ", ";
    }
    file << "]";
    if (i < state.remaining_points.size() - 1) file << ",";
    file << "\n";
  }
  file << "  ],\n";

  // Save polynomial result to temporary file and reference it
  std::string poly_filename = filename + "_polynomial.json";
  state.accumulated_result.exportToJson(poly_filename.substr(0, poly_filename.find_last_of('.')));
  file << "  \"polynomial_file\": \"" << poly_filename << "\"\n";

  file << "}\n";
  file.close();
}

CheckpointState CheckpointManager::loadFromJson(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open checkpoint file: " + filename);
  }

  CheckpointState state;
  std::string line;
  bool in_config = false;
  bool in_remaining_points = false;
  std::string polynomial_file;

  // Simple JSON parsing - this is a basic implementation
  // In production, would use a proper JSON library
  while (std::getline(file, line)) {
    // Trim whitespace
    line.erase(0, line.find_first_not_of(" \t"));
    line.erase(line.find_last_not_of(" \t") + 1);

    if (line.find("\"resume_count\"") != std::string::npos) {
      size_t colon_pos = line.find(':');
      if (colon_pos != std::string::npos) {
        std::string value = line.substr(colon_pos + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t,") + 1);
        state.resume_count = std::stoi(value) + 1; // Increment on load
      }
    } else if (line.find("\"total_points_count\"") != std::string::npos) {
      size_t colon_pos = line.find(':');
      if (colon_pos != std::string::npos) {
        std::string value = line.substr(colon_pos + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t,") + 1);
        state.total_points_count = std::stoi(value);
      }
    } else if (line.find("\"config\"") != std::string::npos) {
      in_config = true;
    } else if (in_config && line.find("\"degree\"") != std::string::npos) {
      size_t colon_pos = line.find(':');
      if (colon_pos != std::string::npos) {
        std::string value = line.substr(colon_pos + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t,") + 1);
        state.config.degree = std::stoi(value);
      }
    } else if (in_config && line.find("\"components\"") != std::string::npos) {
      size_t colon_pos = line.find(':');
      if (colon_pos != std::string::npos) {
        std::string value = line.substr(colon_pos + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t,") + 1);
        state.config.components = std::stoi(value);
      }
    } else if (in_config && line.find("\"writhe\"") != std::string::npos) {
      size_t colon_pos = line.find(':');
      if (colon_pos != std::string::npos) {
        std::string value = line.substr(colon_pos + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t,") + 1);
        state.config.writhe = std::stoi(value);
      }
    } else if (in_config && line.find("\"prefactors\"") != std::string::npos) {
      size_t colon_pos = line.find(':');
      if (colon_pos != std::string::npos) {
        std::string value = line.substr(colon_pos + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t,") + 1);
        state.config.prefactors = std::stoi(value);
      }
    } else if (in_config && line.find("\"crossings\"") != std::string::npos) {
      size_t colon_pos = line.find(':');
      if (colon_pos != std::string::npos) {
        std::string value = line.substr(colon_pos + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t,") + 1);
        state.config.crossings = std::stoi(value);
      }
    } else if (line.find("}") != std::string::npos && in_config) {
      in_config = false;
    } else if (line.find("\"remaining_points\"") != std::string::npos) {
      in_remaining_points = true;
    } else if (in_remaining_points && line.find("[") != std::string::npos && line.find("]") != std::string::npos) {
      // Parse point array like [1, 2, 3]
      std::vector<int> point;
      size_t start = line.find('[');
      size_t end = line.find(']');
      if (start != std::string::npos && end != std::string::npos && end > start) {
        std::string point_str = line.substr(start + 1, end - start - 1);
        std::stringstream ss(point_str);
        std::string num;
        while (std::getline(ss, num, ',')) {
          num.erase(0, num.find_first_not_of(" \t"));
          num.erase(num.find_last_not_of(" \t") + 1);
          if (!num.empty()) {
            point.push_back(std::stoi(num));
          }
        }
        if (!point.empty()) {
          state.remaining_points.push_back(point);
        }
      }
    } else if (line.find("]") != std::string::npos && in_remaining_points &&
               line.find("\"remaining_points\"") == std::string::npos) {
      in_remaining_points = false;
    } else if (line.find("\"polynomial_file\"") != std::string::npos) {
      size_t colon_pos = line.find(':');
      if (colon_pos != std::string::npos) {
        std::string value = line.substr(colon_pos + 1);
        value.erase(0, value.find_first_not_of(" \t\""));
        value.erase(value.find_last_not_of(" \t\",") + 1);
        polynomial_file = value;
      }
    }
  }

  file.close();

  // Initialize accumulated_result with proper dimensions
  state.accumulated_result = MultivariablePolynomial(state.config.components, state.config.degree);

  // Load polynomial from polynomial_file
  if (!polynomial_file.empty()) {
    try {
      // Parse the polynomial JSON file manually
      std::ifstream poly_file(polynomial_file);
      if (poly_file.is_open()) {
        std::string line;
        bool in_terms = false;
        int terms_loaded = 0;

        while (std::getline(poly_file, line)) {
          // Trim whitespace
          line.erase(0, line.find_first_not_of(" \t\n\r\f\v"));
          line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);

          if (line.find("\"terms\"") != std::string::npos) {
            in_terms = true;
            continue;
          }

          if (in_terms && line.find("\"x\":") != std::string::npos &&
              line.find("\"q\":") != std::string::npos &&
              line.find("\"c\":") != std::string::npos) {
            // Parse term like: {"x": [150], "q": 2, "c": -1}
            std::vector<int> x_vals;
            int q_val = 0;
            int c_val = 0;

            // Extract x values
            size_t x_start = line.find("[");
            size_t x_end = line.find("]");
            if (x_start != std::string::npos && x_end != std::string::npos) {
              std::string x_part = line.substr(x_start + 1, x_end - x_start - 1);
              std::istringstream x_stream(x_part);
              std::string x_num;
              while (std::getline(x_stream, x_num, ',')) {
                x_num.erase(0, x_num.find_first_not_of(" \t"));
                x_num.erase(x_num.find_last_not_of(" \t") + 1);
                if (!x_num.empty()) {
                  x_vals.push_back(std::stoi(x_num));
                }
              }
            }

            // Extract q value
            size_t q_pos = line.find("\"q\":");
            if (q_pos != std::string::npos) {
              size_t q_start = q_pos + 4;
              size_t q_end = line.find(",", q_start);
              if (q_end == std::string::npos) q_end = line.find("}", q_start);
              if (q_end != std::string::npos) {
                std::string q_str = line.substr(q_start, q_end - q_start);
                q_str.erase(0, q_str.find_first_not_of(" \t"));
                q_str.erase(q_str.find_last_not_of(" \t") + 1);
                q_val = std::stoi(q_str);
              }
            }

            // Extract c value
            size_t c_pos = line.find("\"c\":");
            if (c_pos != std::string::npos) {
              size_t c_start = c_pos + 4;
              size_t c_end = line.find("}", c_start);
              if (c_end != std::string::npos) {
                std::string c_str = line.substr(c_start, c_end - c_start);
                c_str.erase(0, c_str.find_first_not_of(" \t"));
                c_str.erase(c_str.find_last_not_of(" \t") + 1);
                c_val = std::stoi(c_str);
              }
            }

            // Add the term to the polynomial
            if (!x_vals.empty() && c_val != 0) {
              state.accumulated_result.setCoefficient(q_val, x_vals, c_val);
              terms_loaded++;
            }
          }

          // Check for end of terms array (line with just "]" or "\t]")
          if (in_terms && (line == "]" || line == "\t]")) {
            in_terms = false;
          }
        }
        poly_file.close();
        // Uncomment for debugging: std::cout << "Loaded " << terms_loaded << " terms from checkpoint polynomial" << std::endl;
      }
    } catch (const std::exception& e) {
      std::cerr << "Warning: Failed to load polynomial from " << polynomial_file
                << ": " << e.what() << std::endl;
      // Continue with empty polynomial if loading fails
    }
  }

  return state;
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

  // Find valid criteria
  auto valid_criteria = findValidCriteria();
  if (!valid_criteria.is_valid) {
    throw std::runtime_error("No valid criteria found");
  }

  // Assign variables to get list of variable assignments
  auto assignments = assignVariables(valid_criteria);

  // Collect all points from each variable assignment
  std::vector<std::vector<int>> all_points;
  for (const auto& assignment : assignments) {
    auto points = enumeratePoints(assignment);
    all_points.insert(all_points.end(), points.begin(), points.end());
  }

  // Run the function on all collected points
  for (const auto& point : all_points) {
    engine_->computeForAngles_new(point);
  }

  // Final offset addition
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

// Pooling functionality implementations (moved from solution_pool_1a_double_links.cpp)

bool FKComputation::satisfiesConstraints(const std::vector<int>& point,
                         const std::vector<std::vector<double>>& constraints) {
  for (const auto& constraint : constraints) {
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

std::vector<std::vector<int>> FKComputation::enumeratePoints(const AssignmentResult& assignment) {

  std::vector<std::vector<int>> valid_points;

  if (assignment.bounds.empty()) {
    // Base case: check constraints and add point if valid
    if (satisfiesConstraints(assignment.point, assignment.supporting_inequalities) &&
        satisfiesConstraints(assignment.point, assignment.criteria)) {
      valid_points.push_back(assignment.point);
    }
    return valid_points;
  }

  // Convert bounds list to vector for easier iteration
  std::vector<std::array<int, 2>> bounds_vec(assignment.bounds.begin(), assignment.bounds.end());

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
  int upper = static_cast<int>(assignment.supporting_inequalities[inequality][0]);
  for (size_t i = 0; i < assignment.point.size(); i++) {
    if (static_cast<int>(i) != index) {
      upper += static_cast<int>(assignment.supporting_inequalities[inequality][1 + i]) * assignment.point[i];
    }
  }
  upper /= -static_cast<int>(assignment.supporting_inequalities[inequality][1 + index]);
  initial_frame.upper_bound = upper;

  stack.push(initial_frame);

  while (!stack.empty()) {
    IterationFrame& current = stack.top();

    if (current.current_value > current.upper_bound) {
      stack.pop();
      continue;
    }

    // Set current variable value
    current.point[bounds_vec[current.bound_index][0]] = current.current_value;

    if (current.bound_index == bounds_vec.size() - 1) {
      // Last variable - check constraints and add point if valid
      if (satisfiesConstraints(current.point, assignment.supporting_inequalities) &&
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
      int next_upper = static_cast<int>(assignment.supporting_inequalities[next_inequality][0]);
      for (size_t i = 0; i < next_frame.point.size(); i++) {
        if (static_cast<int>(i) != next_index) {
          next_upper += static_cast<int>(assignment.supporting_inequalities[next_inequality][1 + i]) * next_frame.point[i];
        }
      }
      next_upper /= -static_cast<int>(assignment.supporting_inequalities[next_inequality][1 + next_index]);
      next_frame.upper_bound = next_upper;

      current.current_value++;
      stack.push(next_frame);
    }
  }

  return valid_points;
}

std::vector<FKComputation::AssignmentResult> FKComputation::assignVariables(const ValidatedCriteria& valid_criteria) {

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
  std::vector<std::array<int, 2>> first_vec(valid_criteria.first_bounds.begin(), valid_criteria.first_bounds.end());

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
  initial_state.max_value = static_cast<int>(valid_criteria.degrees[main_index] / slope);

  stack.push(initial_state);

  while (!stack.empty()) {
    VariableAssignmentState& current = stack.top();

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
    new_degrees[main_idx] = current.degrees[main_idx] - current.current_value * slope;

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
      next_state.max_value = static_cast<int>(new_degrees[next_main_idx] / next_slope);

      current.current_value++;
      stack.push(next_state);
    }
  }

  return assignments;
}

FKComputation::BoundedVariables FKComputation::identifyBoundedVariables(const std::vector<std::vector<double>>& inequalities,
                                         int size) {
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
    std::vector<int>& bounded_v,
    int& bounded_count,
    int size,
    const std::vector<std::vector<double>>& supporting_inequalities) {

  std::list<std::array<int, 2>> bounds;
  int support = supporting_inequalities.size();

  int index = 0;
  while (index < size - 1) {
    if (!bounded_v[index]) {
      for (int l = 0; l < support; l++) {
        if (supporting_inequalities[l][1 + index] < 0) {
          bool useful = true;
          for (int n = 0; n < size - 1; n++) {
            if (n != index && supporting_inequalities[l][1 + n] > 0 && !bounded_v[n]) {
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

std::vector<double> FKComputation::extractDegrees(const std::vector<std::vector<double>>& inequalities) {
  std::vector<double> degrees;
  for (const auto& x : inequalities) {
    degrees.push_back(x[0]);
  }
  return degrees;
}

bool FKComputation::validCriteria(const std::vector<std::vector<double>>& criteria,
                   const std::vector<std::vector<double>>& supporting_inequalities,
                   int size) {

  auto bounded_info = identifyBoundedVariables(criteria, size);

  if (bounded_info.bounded_count == 0) {
    return false; // Not valid - no bounded variables
  }

  if (bounded_info.bounded_count == size - 1) {
    return true; // We have enough bounded variables
  }

  // Try to find additional bounds
  auto additional_bounds = findAdditionalBounds(bounded_info.bounded_v, bounded_info.bounded_count,
                                               size, supporting_inequalities);

  return (bounded_info.bounded_count == size - 1); // Valid if we now have enough
}

FKComputation::ValidatedCriteria FKComputation::findValidCriteria() {

  if (config_.criteria.empty()) {
    return ValidatedCriteria(); // Return invalid criteria
  }

  const int mains = config_.criteria.size();
  const int size = config_.criteria[0].size();
  const double COMBINATION_FACTOR = 0.5;  // Was hardcoded /2.0

  // Helper function to build ValidatedCriteria from valid criteria
  auto buildValidatedCriteria = [&](const std::vector<std::vector<double>>& criteria) -> ValidatedCriteria {
    ValidatedCriteria result;
    auto bounded_info = identifyBoundedVariables(criteria, size);

    result.criteria = criteria;
    result.degrees = extractDegrees(config_.criteria);
    result.first_bounds = bounded_info.first;
    result.initial_point = std::vector<int>(size - 1, 0);
    result.is_valid = true;

    // Add additional bounds if needed
    if (bounded_info.bounded_count < size - 1) {
      result.additional_bounds = findAdditionalBounds(bounded_info.bounded_v, bounded_info.bounded_count,
                                                     size, config_.inequalities);
    }
    return result;
  };

  // Combine all inequalities for processing
  std::vector<std::vector<double>> all_inequalities = config_.inequalities;
  all_inequalities.insert(all_inequalities.end(), config_.criteria.begin(), config_.criteria.end());

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
      for (const auto& inequality : all_inequalities) {

        // Check if this combination could be beneficial (has opposing signs)
        bool potentially_beneficial = false;
        for (int var_idx = 1; var_idx < size; ++var_idx) {
          if (current_criteria[criterion_idx][var_idx] > 0 && inequality[var_idx] < 0) {
            potentially_beneficial = true;
            break;
          }
        }

        if (!potentially_beneficial) continue;

        // Create new criterion by linear combination
        auto new_criteria = current_criteria;
        for (int var_idx = 0; var_idx < size; ++var_idx) {
          new_criteria[criterion_idx][var_idx] += inequality[var_idx] * COMBINATION_FACTOR;
        }

        // Skip if we've seen this configuration before
        if (visited_criteria.find(new_criteria) != visited_criteria.end()) {
          continue;
        }

        // Mark as visited
        visited_criteria.insert(new_criteria);

        // Check if these criteria are valid
        if (validCriteria(new_criteria, config_.inequalities, size)) {
          return buildValidatedCriteria(new_criteria);  // Found valid criteria
        }

        // Add to queue for further exploration
        criteria_queue.push(std::move(new_criteria));
      }
    }
  }

  // No valid criteria found
  return ValidatedCriteria();
}

// Checkpoint methods implementation
std::string FKComputation::computeWithCheckpointing(const std::string& input_filename,
                                                   const std::string& output_filename,
                                                   int checkpoint_period) {
  FKConfiguration config = parser_.parseFromFile(input_filename);

  if (!config.isValid()) {
    throw std::runtime_error("Invalid configuration provided");
  }

  config_ = config;
  initializeEngine();

  // Find valid criteria
  auto valid_criteria = findValidCriteria();
  if (!valid_criteria.is_valid) {
    throw std::runtime_error("No valid criteria found");
  }

  // Assign variables to get list of variable assignments
  auto assignments = assignVariables(valid_criteria);

  // Collect all points from each variable assignment
  std::vector<std::vector<int>> all_points;
  for (const auto& assignment : assignments) {
    auto points = enumeratePoints(assignment);
    all_points.insert(all_points.end(), points.begin(), points.end());
  }

  // Create initial checkpoint
  CheckpointState initial_state(config_, all_points.size());
  initial_state.remaining_points = all_points;
  initial_state.accumulated_result = MultivariablePolynomial(config_.components, config_.degree);

  CheckpointManager checkpoint_manager(initial_state);
  std::cout << "Created checkpoint: " << checkpoint_manager.getFilename() << std::endl;

  // Process points with periodic checkpointing
  int points_processed = 0;
  std::vector<std::vector<int>> remaining_points = all_points;

  for (size_t i = 0; i < all_points.size(); ++i) {
    engine_->computeForAngles(all_points[i]);
    points_processed++;

    // Save checkpoint periodically
    if (checkpoint_period > 0 && points_processed % checkpoint_period == 0) {
      // Create remaining points list
      std::vector<std::vector<int>> current_remaining(
        all_points.begin() + i + 1, all_points.end());

      CheckpointState current_state = saveCheckpoint(current_remaining, all_points.size());
      checkpoint_manager.save(current_state);

      std::cout << "Checkpoint saved after " << points_processed
                << " points. Remaining: " << current_remaining.size() << std::endl;
    }
  }

  // Final offset addition (same as original compute method)
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

  int components_copy = config_.components;

  performOffsetAddition(result_coeffs1, result_coeffs2, increment_offset, 0,
                       components_copy, maxima, -1, accumulator_block_sizes,
                       accumulator_block_sizes);

  const_cast<MultivariablePolynomial &>(engine_->getResult())
     .syncFromDenseVector(result_coeffs1);

  writer_.writeToJson(engine_->getResult(), output_filename);

  // Clean up checkpoint after successful completion
  std::cout << "Computation completed successfully." << std::endl;

  return checkpoint_manager.getFilename();
}

void FKComputation::resumeFromCheckpoint(const std::string& checkpoint_path,
                                        const std::string& output_filename) {
  std::cout<<"Started resume from checkpoint"<<std::endl;
  CheckpointManager checkpoint_manager(checkpoint_path);
  std::cout<<"Created checkpoint manager"<<std::endl;
  CheckpointState state = checkpoint_manager.load();
  std::cout<<"Loaded state"<<std::endl;

  std::cout << "Resuming computation from checkpoint (resume count: "
            << state.resume_count << ")" << std::endl;
  std::cout << "Remaining points: " << state.remaining_points.size()
            << " of " << state.total_points_count << std::endl;

  // Restore state
  config_ = state.config;
  initializeEngine();
  std::cout<<"Engine initialized"<<std::endl;

  // Restore accumulated result
  loadCheckpoint(state);
  std::cout<<"State loaded"<<std::endl;

  // Continue processing remaining points
  for (const auto& point : state.remaining_points) {
    std::cout<<"Point"<<std::endl;
    for (const auto& c : point){
      std::cout<<c<<std::endl;
    }
    engine_->computeForAngles(point);
  }

  // Always perform final offset addition
  // This is needed because checkpoints save the intermediate engine state before final processing
  {
    // Final offset addition (same as original compute method)
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

    int components_copy = config_.components;

    performOffsetAddition(result_coeffs1, result_coeffs2, increment_offset, 0,
                         components_copy, maxima, -1, accumulator_block_sizes,
                         accumulator_block_sizes);

    const_cast<MultivariablePolynomial &>(engine_->getResult())
       .syncFromDenseVector(result_coeffs1);
  }

  writer_.writeToJson(engine_->getResult(), output_filename);

  std::cout << "Resumed computation completed successfully." << std::endl;
}

CheckpointState FKComputation::saveCheckpoint(const std::vector<std::vector<int>>& remaining_points,
                                             int total_points) const {
  CheckpointState state(config_, total_points);
  state.remaining_points = remaining_points;
  state.accumulated_result = engine_->getResult();
  return state;
}

void FKComputation::loadCheckpoint(const CheckpointState& state) {
  // Restore the accumulated result in the engine
  // Note: This is a simplified approach - in practice, we'd need to
  // ensure the engine's internal state is properly restored
  auto result_coeffs = state.accumulated_result.getCoefficients();
  const_cast<MultivariablePolynomial &>(engine_->getResult())
     .syncFromDenseVector(result_coeffs);
}

} // namespace fk
