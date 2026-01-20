// implement multithreading at top level of recursion

#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <vector>

#include "fk/bilvector.hpp"
#include "fk/linalg.hpp"
#include "fk/polynomial_config.hpp"
#include "fk/qalg_links.hpp"
#include "fk/string_to_int.hpp"

// State structure for iterative enumeration
struct EnumerationState {
  std::vector<std::vector<double>> criteria;
  std::list<std::array<int, 2>> bounds;
  std::vector<std::vector<double>> supporting_inequalities;
  std::vector<int> point;
  int current_bound_index;
  int current_value;
  int upper_bound;
};

// State structure for iterative variable assignment
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

// Helper function to identify bounded variables
struct BoundedVariables {
  std::vector<int> bounded_v;
  int bounded_count;
  std::list<std::array<int, 2>> first;
};

// Structure to hold validated criteria ready for variable assignment
struct ValidatedCriteria {
  std::vector<std::vector<double>> criteria;
  std::vector<double> degrees;
  std::list<std::array<int, 2>> first_bounds;
  std::list<std::array<int, 2>> additional_bounds;
  std::vector<int> initial_point;
  bool is_valid;

  ValidatedCriteria() : is_valid(false) {}
};

// Helper function to check if a point satisfies all constraints
bool satisfiesConstraints(const std::vector<int>& point,
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

BoundedVariables identifyBoundedVariables(const std::vector<std::vector<double>>& inequalities,
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

// Helper function to find additional bounds
std::list<std::array<int, 2>> findAdditionalBounds(
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

// Helper function to extract degrees from inequalities
std::vector<double> extractDegrees(const std::vector<std::vector<double>>& inequalities) {
  std::vector<double> degrees;
  for (const auto& x : inequalities) {
    degrees.push_back(x[0]);
  }
  return degrees;
}

// Simple function to check if criteria are valid (just returns bool)
bool validCriteria(const std::vector<std::vector<double>>& criteria,
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

// Find valid criteria through breadth-first search
ValidatedCriteria findValidCriteria(const std::vector<std::vector<double>>& main_inequalities,
                                   const std::vector<std::vector<double>>& supporting_inequalities) {

  if (main_inequalities.empty()) {
    return ValidatedCriteria(); // Return invalid criteria
  }

  const int mains = main_inequalities.size();
  const int size = main_inequalities[0].size();
  const double COMBINATION_FACTOR = 0.5;  // Was hardcoded /2.0

  // Helper function to build ValidatedCriteria from valid criteria
  auto buildValidatedCriteria = [&](const std::vector<std::vector<double>>& criteria) -> ValidatedCriteria {
    ValidatedCriteria result;
    auto bounded_info = identifyBoundedVariables(criteria, size);

    result.criteria = criteria;
    result.degrees = extractDegrees(main_inequalities);
    result.first_bounds = bounded_info.first;
    result.initial_point = std::vector<int>(size - 1, 0);
    result.is_valid = true;

    // Add additional bounds if needed
    if (bounded_info.bounded_count < size - 1) {
      result.additional_bounds = findAdditionalBounds(bounded_info.bounded_v, bounded_info.bounded_count,
                                                     size, supporting_inequalities);
    }
    return result;
  };

  // Combine all inequalities for processing
  std::vector<std::vector<double>> all_inequalities = supporting_inequalities;
  all_inequalities.insert(all_inequalities.end(), main_inequalities.begin(), main_inequalities.end());

  // Use a single set to track visited criterion configurations more efficiently
  std::set<std::vector<std::vector<double>>> visited_criteria;

  // Queue for breadth-first exploration of criterion space
  std::queue<std::vector<std::vector<double>>> criteria_queue;

  // Start with the main inequalities
  visited_criteria.insert(main_inequalities);
  criteria_queue.push(main_inequalities);

  // Check initial criteria first
  if (validCriteria(main_inequalities, supporting_inequalities, size)) {
    return buildValidatedCriteria(main_inequalities);
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
        if (validCriteria(new_criteria, supporting_inequalities, size)) {
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

// Iterative version of recurse_2
void enumeratePoints(std::vector<std::vector<double>>& criteria,
                    std::list<std::array<int, 2>> bounds,
                    std::vector<std::vector<double>> supporting_inequalities,
                    std::vector<int> point,
                    const std::function<void(const std::vector<int>&)>& function) {

  if (bounds.empty()) {
    // Base case: check constraints and call function if valid
    if (satisfiesConstraints(point, supporting_inequalities) &&
        satisfiesConstraints(point, criteria)) {
      function(point);
    }
    return;
  }

  // Convert bounds list to vector for easier iteration
  std::vector<std::array<int, 2>> bounds_vec(bounds.begin(), bounds.end());

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
  initial_frame.point = point;
  initial_frame.bound_index = 0;
  initial_frame.current_value = 0;

  // Calculate upper bound for first variable
  int index = bounds_vec[0][0];
  int inequality = bounds_vec[0][1];
  int upper = static_cast<int>(supporting_inequalities[inequality][0]);
  for (size_t i = 0; i < point.size(); i++) {
    if (static_cast<int>(i) != index) {
      upper += static_cast<int>(supporting_inequalities[inequality][1 + i]) * point[i];
    }
  }
  upper /= -static_cast<int>(supporting_inequalities[inequality][1 + index]);
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
      // Last variable - check constraints and call function
      if (satisfiesConstraints(current.point, supporting_inequalities) &&
          satisfiesConstraints(current.point, criteria)) {
        function(current.point);
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
      int next_upper = static_cast<int>(supporting_inequalities[next_inequality][0]);
      for (size_t i = 0; i < next_frame.point.size(); i++) {
        if (static_cast<int>(i) != next_index) {
          next_upper += static_cast<int>(supporting_inequalities[next_inequality][1 + i]) * next_frame.point[i];
        }
      }
      next_upper /= -static_cast<int>(supporting_inequalities[next_inequality][1 + next_index]);
      next_frame.upper_bound = next_upper;

      current.current_value++;
      stack.push(next_frame);
    }
  }
}

// Iterative version of recurse_1
void assignVariables(std::vector<std::vector<double>>& new_criteria,
                    std::vector<double> degrees,
                    std::vector<std::vector<double>>& criteria,
                    std::list<std::array<int, 2>> first,
                    std::list<std::array<int, 2>> bounds,
                    std::vector<std::vector<double>> supporting_inequalities,
                    std::vector<int> point,
                    const std::function<void(const std::vector<int>&)>& function) {

  if (first.empty()) {
    enumeratePoints(criteria, bounds, supporting_inequalities, point, function);
    return;
  }

  // Convert first list to vector for easier iteration
  std::vector<std::array<int, 2>> first_vec(first.begin(), first.end());

  std::stack<VariableAssignmentState> stack;

  // Initialize first state
  VariableAssignmentState initial_state;
  initial_state.new_criteria = new_criteria;
  initial_state.degrees = degrees;
  initial_state.criteria = criteria;
  initial_state.bounds = bounds;
  initial_state.supporting_inequalities = supporting_inequalities;
  initial_state.point = point;
  initial_state.current_var_index = 0;
  initial_state.current_value = 0;

  // Calculate max value for first variable
  int var_index = first_vec[0][0];
  int main_index = first_vec[0][1];
  double slope = -new_criteria[main_index][var_index];
  initial_state.max_value = static_cast<int>(degrees[main_index] / slope);

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
      // Last variable - proceed to enumeration
      enumeratePoints(current.criteria, current.bounds,
                     current.supporting_inequalities, current.point, function);
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
}

void pooling(std::vector<std::vector<double>> main_inequalities,
             std::vector<std::vector<double>> supporting_inequalities,
             const std::function<void(const std::vector<int>&)>& function) {

  // Find valid criteria
  auto valid_criteria = findValidCriteria(main_inequalities, supporting_inequalities);
  if (!valid_criteria.is_valid) {
    return;
  }

  // Assign variables and enumerate points
  auto criteria_copy = valid_criteria.criteria;
  assignVariables(criteria_copy, valid_criteria.degrees, criteria_copy,
                 valid_criteria.first_bounds, valid_criteria.additional_bounds,
                 supporting_inequalities, valid_criteria.initial_point, function);
}

class FK {
private:
  std::vector<int> angles;
  std::vector<int> angle_signs;
  std::vector<int> braid;
  std::vector<int> signs;
  std::vector<int> segments;
  std::vector<int> accumulatorBlockSizes = {1};
  std::vector<int> closed_strand_components = {};
  std::vector<std::vector<std::array<int, 2>>> crossingMatrices = {};
  std::vector<int> crossingRelationTypes = {};
  std::vector<std::vector<int>> transitions;
  std::vector<bool> trivial_angles_;
  std::vector<int> nontrivial_map;
  std::vector<int> inversion_data;
  PolynomialType result{
      1, 1}; // Initialize with dummy values, will be reassigned
  std::vector<std::vector<std::vector<int>>> variableAssignments;
  std::vector<std::vector<int>> numericalAssignments;
  std::vector<int> top_crossing_components = {};
  std::vector<int> bottom_crossing_components = {};
  int components;
  int writhe = 0;
  int prefactors;
  int crossings;
  int degree;
  int points_found = 0;
  void computeNumericalAssignment(const std::vector<int> &angles) {
    points_found++;
    std::cout<<"Computing for angles: "<<std::endl;
    for (auto angle : angles){   std::cout<<angle<<std::endl;  }

    for (int i = 0; i < crossings + 1; i++) {
      for (int j = 0; j < prefactors + 1; j++) {
        numericalAssignments[i][j] =
            computeDotProduct(variableAssignments[i][j], angles);
      }
    }
    std::cout << "Numerical assignments:" << std::endl;
    for (size_t i = 0; i < numericalAssignments.size(); i++) {
        for (size_t j = 0; j < numericalAssignments[i].size(); j++) {
            std::cout << "assignments[" << i << "][" << j << "] = " << numericalAssignments[i][j] << std::endl;
        }
    }

    double qPowerAccumulatorDouble = (writhe - prefactors) / 2.0;
    std::vector<double> xPowerAccumulatorDouble(components, 0);
    int initialCoefficient = 1;
    for (int i = 0; i < prefactors; i++) {
      qPowerAccumulatorDouble -= numericalAssignments[0][i + 1];
      xPowerAccumulatorDouble[closed_strand_components[i]] -= 0.5;
    }
    xPowerAccumulatorDouble[0] -= 0.5;

    for (int crossingIndex = 0; crossingIndex < crossings; crossingIndex++) {
      int matrixParam_i =
          numericalAssignments[crossingMatrices[crossingIndex][0][0]]
                              [crossingMatrices[crossingIndex][0][1]];
      int matrixParam_k =
          numericalAssignments[crossingMatrices[crossingIndex][2][0]]
                              [crossingMatrices[crossingIndex][2][1]];
      int matrixParam_j =
          numericalAssignments[crossingMatrices[crossingIndex][1][0]]
                              [crossingMatrices[crossingIndex][1][1]];
      int matrixParam_m =
          numericalAssignments[crossingMatrices[crossingIndex][3][0]]
                              [crossingMatrices[crossingIndex][3][1]];
      int topComponent = top_crossing_components[crossingIndex];
      int bottomComponent = bottom_crossing_components[crossingIndex];
      if (crossingRelationTypes[crossingIndex] == 1 ||
          crossingRelationTypes[crossingIndex] == 2) {
        qPowerAccumulatorDouble += (matrixParam_j + matrixParam_m) / 2.0 +
                                   matrixParam_j * matrixParam_m;
        xPowerAccumulatorDouble[topComponent] +=
            ((matrixParam_j + matrixParam_k + 1) / 4.0);
        xPowerAccumulatorDouble[bottomComponent] +=
            ((3 * matrixParam_m - matrixParam_i + 1) / 4.0);
      } else if (crossingRelationTypes[crossingIndex] == 4) {
        qPowerAccumulatorDouble -= (matrixParam_i + matrixParam_k +
                                    matrixParam_m * (matrixParam_m + 1) -
                                    matrixParam_i * (matrixParam_i + 1)) /
                                       2.0 +
                                   matrixParam_i * matrixParam_k;
        xPowerAccumulatorDouble[topComponent] -=
            ((3 * matrixParam_j - matrixParam_k + 1) / 4.0);
        xPowerAccumulatorDouble[bottomComponent] -=
            ((matrixParam_i + matrixParam_m + 1) / 4.0);
        if ((matrixParam_j - matrixParam_k) % 2 == 0) {
          initialCoefficient *= -1;
        }
      } else if (crossingRelationTypes[crossingIndex] == 3) {
        qPowerAccumulatorDouble -= (matrixParam_i + matrixParam_k +
                                    matrixParam_m * (matrixParam_m + 1) -
                                    matrixParam_i * (matrixParam_i + 1)) /
                                       2.0 +
                                   matrixParam_i * matrixParam_k;
        xPowerAccumulatorDouble[topComponent] -=
            ((3 * matrixParam_j - matrixParam_k + 1) / 4.0);
        xPowerAccumulatorDouble[bottomComponent] -=
            ((matrixParam_i + matrixParam_m + 1) / 4.0);
        if ((matrixParam_k - matrixParam_j) % 2 == 1) {
          initialCoefficient *= -1;
        }
      }
    }

    std::cout << "qPowerAccumulatorDouble : " << qPowerAccumulatorDouble
              << "\n";
    int qPowerAccumulator = static_cast<int>(std::floor(
        qPowerAccumulatorDouble)); // currently, we are losing the actual q
                                   // powers because rounding down; later,
                                   // implement output with the exact q powers
    std::vector<int> xPowerAccumulator(components);
    std::vector<int> maxXDegrees(components);
    std::vector<int> blockSizes(components);
    blockSizes[0] = 1;
    for (int n = 0; n < components; n++) {
      xPowerAccumulator[n] = xPowerAccumulatorDouble[n];
      maxXDegrees[n] = degree - xPowerAccumulator[n];
      if (n != 0) {
        blockSizes[n] = (maxXDegrees[n - 1] + 1) * blockSizes[n - 1];
      }
    }
    int totalProductSize =
        blockSizes[components - 1] * (maxXDegrees[components - 1] + 1);
    // std::cout << "x_acc: " << xPowerAccumulatorDouble[0] << "\n\n";
    std::cout << "x_acc: " << xPowerAccumulatorDouble[0] << " "
              << xPowerAccumulatorDouble[1]
              << "\n\n"; // modifying: x_acc's are coming out to half-integers
                         // in general
    // exit(0);
    // if (xPowerAccumulatorDouble[0] > 15) {
    //     std::cout << "exiting due to x_acc large enough\n";
    //     exit(0);
    // }
    std::vector<bilvector<int>> polynomialTerms(
        totalProductSize,
        bilvector<int>(0, 1, 20, 0)); // error is the degree issue; modifying
    //

    polynomialTerms[0][0] = initialCoefficient;
    for (int crossingIndex = 0; crossingIndex < crossings; crossingIndex++) {
        int param_i =
            numericalAssignments[crossingMatrices[crossingIndex][0][0]]
                                [crossingMatrices[crossingIndex][0][1]];
        int param_m =
            numericalAssignments[crossingMatrices[crossingIndex][3][0]]
                                [crossingMatrices[crossingIndex][3][1]];
        int param_j =
            numericalAssignments[crossingMatrices[crossingIndex][1][0]]
                                [crossingMatrices[crossingIndex][1][1]];
        int param_k =
            numericalAssignments[crossingMatrices[crossingIndex][2][0]]
                                [crossingMatrices[crossingIndex][2][1]];
        int relation_type = crossingRelationTypes[crossingIndex];
       std::cout<< "relation_type "<<relation_type<<std::endl;
       std::cout<< "param_i "<<param_i<<std::endl;
       std::cout<< "param_j "<<param_j<<std::endl;
       std::cout<< "param_k "<<param_k<<std::endl;
       std::cout<< "param_m "<<param_m<<std::endl;

      if (relation_type == 1) {
        if (param_i > 0) {
          computePositiveQBinomial(polynomialTerms, param_i, param_i - param_m,
                                   false);
        } else {
          computeNegativeQBinomial(polynomialTerms, param_i, param_i - param_m,
                                   false);
        }
      } else if (relation_type == 2) {
      computeNegativeQBinomial(polynomialTerms, param_i, param_m, false);
      } else if (relation_type == 3) {
        computeNegativeQBinomial(polynomialTerms, param_j, param_k, true);
      } else {
        if (param_j > 0) {
          computePositiveQBinomial(polynomialTerms, param_j, param_j - param_k,
                                   true);
        } else {
          computeNegativeQBinomial(polynomialTerms, param_j, param_j - param_k,
                                   true);
        }
      }
    }
    for (size_t i = 0; i < polynomialTerms.size(); i++) {
      int max_pos = polynomialTerms[i].getMaxPositiveIndex();
      int max_neg = polynomialTerms[i].getMaxNegativeIndex();
      std::cout << "polynomialTerms[" << i << "]: ";
      for (int j = max_neg; j <= max_pos; j++) {
        std::cout << polynomialTerms[i][j] << " ";
      }
      std::cout << std::endl;
    }

    for (int crossingIndex = 0; crossingIndex < crossings; crossingIndex++) {
        int relation_type = crossingRelationTypes[crossingIndex];
        int param_j =
            numericalAssignments[crossingMatrices[crossingIndex][1][0]]
                                [crossingMatrices[crossingIndex][1][1]];
        int param_k =
            numericalAssignments[crossingMatrices[crossingIndex][2][0]]
                                [crossingMatrices[crossingIndex][2][1]];
        int param_i =
            numericalAssignments[crossingMatrices[crossingIndex][0][0]]
                                [crossingMatrices[crossingIndex][0][1]];
        int param_m =
            numericalAssignments[crossingMatrices[crossingIndex][3][0]]
                                [crossingMatrices[crossingIndex][3][1]];
        int bottomComp = bottom_crossing_components[crossingIndex];
        int topComp = top_crossing_components[crossingIndex];

       std::cout<<"================================================="<<std::endl;
       std::cout<< "relation_type "<<relation_type<<std::endl;
       std::cout<< "param_i "<<param_i<<std::endl;
       std::cout<< "param_j "<<param_j<<std::endl;
       std::cout<< "param_k "<<param_k<<std::endl;
       std::cout<< "param_m "<<param_m<<std::endl;
       std::cout<< "topComp "<<topComp<<std::endl;
       std::cout<< "bottomComp "<<bottomComp<<std::endl;
       std::cout<< "maxXDegrees: "<<std::endl;
       for (auto d: maxXDegrees){
            std::cout<<"\t"<<d<<std::endl;
       }
       std::cout<< "blockSizes: "<<std::endl;
       for (auto d: blockSizes){
            std::cout<<"\t"<<d<<std::endl;
       }
      if (crossingRelationTypes[crossingIndex] == 1) {
        computeXQPochhammer(polynomialTerms, param_k, param_j + 1, bottomComp,
                            components, maxXDegrees, blockSizes);
      } else if (crossingRelationTypes[crossingIndex] == 2) {
        computeXQInversePochhammer(polynomialTerms, param_j, param_k + 1,
                                   bottomComp, components, maxXDegrees,
                                   blockSizes);
      } else if (crossingRelationTypes[crossingIndex] == 3) {
        computeXQInversePochhammer(polynomialTerms, param_i, param_m + 1,
                                   topComp, components, maxXDegrees,
                                   blockSizes);
      } else {
        computeXQPochhammer(polynomialTerms, param_m, param_i + 1, topComp,
                            components, maxXDegrees, blockSizes);
      }
    }

    // Add polynomialTerms to result with offset
    // Refactored to use direct addToCoefficient instead of get/sync pattern
    for (size_t i = 0; i < polynomialTerms.size(); ++i) {
        const auto &qPoly = polynomialTerms[i];

        // Build x-powers vector (simplified - assumes single variable for now)
        std::vector<int> xPowers(components, 0);
        xPowers[0] = static_cast<int>(i);

        // Apply offset to x-powers
        for (int j = 0; j < components; ++j) {
            xPowers[j] += xPowerAccumulator[j];
        }

        // Add each q-coefficient with offset
        for (int q = qPoly.getMaxNegativeIndex(); q <= qPoly.getMaxPositiveIndex(); ++q) {
            int coeff = qPoly[q];
            if (coeff != 0) {
                result.addToCoefficient(q + qPowerAccumulator, xPowers, coeff);
            }
        }
    }
  }
  void writeResultsToJson(std::string fileName) {
    result.exportToJson(fileName);
  }

public:
  std::vector<std::vector<double>> inequalities;
  std::vector<std::vector<double>> criteria;
  std::vector<std::vector<int>> extensions;
  std::string metadata;
  FK(std::string infile_, std::string outfile_) {

    std::ifstream infile;
    infile.open(infile_ + ".csv");
    if (infile.is_open()) {
      std::string line;

      std::getline(infile, line, '\n');
      int index = line.find(",");
      degree = parseStringToInteger(line.substr(0, index));

      std::getline(infile, line, '\n');
      index = line.find(",");
      components = parseStringToInteger(line.substr(0, index));

      std::getline(infile, line, '\n');
      index = line.find(",");
      writhe = parseStringToInteger(line.substr(0, index));

      // Initialize the polynomial after we know the degree and components
      result = PolynomialType(components, degree);

      std::getline(infile, line, '\n');
      int height = 0;
      index = line.find(",");
      bool done = false;
      while (true) {
        if (index == -1) {
          break;
        }

        int c = parseStringToInteger(line.substr(0, index));
        crossingMatrices.push_back({{height, c - 1},
                                    {height, c},
                                    {height + 1, c - 1},
                                    {height + 1, c}});
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");

        crossingRelationTypes.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);

        height++;
        index = line.find(",");
      }
      std::getline(infile, line, '\n');
      index = line.find(",");
      done = false;
      while (true) {
        if (index == -1) {
          break;
        }
        closed_strand_components.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");
      }
      prefactors = closed_strand_components.size();
      crossings = crossingRelationTypes.size();
      variableAssignments.resize(crossings + 1,
                                 std::vector<std::vector<int>>(prefactors + 1));
      numericalAssignments.resize(crossings + 1,
                                  std::vector<int>(prefactors + 1));
      std::getline(infile, line, '\n');
      index = line.find(",");
      done = false;
      while (true) {
        if (index == -1) {
          break;
        }
        top_crossing_components.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");
        bottom_crossing_components.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");
      }
      int stage = 0;
      int criteria_index = 0;
      int inequality_index = 0;
      int extension_index = 0;
      while (std::getline(infile, line, '\n')) {
        if (line[0] == '/') {
          stage++;
        } else if (stage == 0) {
          criteria.push_back({});
          done = false;
          index = line.find(",");
          while (true) {
            if (index == -1) {
              break;
            }
            criteria[criteria_index].push_back(
                parseStringToDouble(line.substr(0, index)));
            line = line.substr(index + 1, line.size() - index - 1);
            index = line.find(",");
          }
          criteria_index++;
        } else if (stage == 1) {
          inequalities.push_back({});
          done = false;
          index = line.find(",");
          while (true) {
            if (index == -1) {
              break;
            }
            inequalities[inequality_index].push_back(
                parseStringToInteger(line.substr(0, index)));
            line = line.substr(index + 1, line.size() - index - 1);
            index = line.find(",");
          }
          inequality_index++;
        } else if (stage == 2) {
          done = false;
          index = line.find(",");
          while (true) {
            if (index == -1) {
              break;
            }
            variableAssignments[extension_index % (crossings + 1)]
                               [extension_index / (crossings + 1)]
                                   .push_back(parseStringToInteger(
                                       line.substr(0, index)));
            line = line.substr(index + 1, line.size() - index - 1);
            index = line.find(",");
          }
          extension_index++;
        }
      }
    } else {
      std::cout << "ERROR: Unable to open file '" + infile_ + ".csv'!";
      exit(0);
    }
    for (int i = 1; i < components; i++) {
      accumulatorBlockSizes.push_back(accumulatorBlockSizes[i - 1] *
                                      (degree + 1));
    }
    std::function<void(const std::vector<int> &)> f_wrapper =
        [this](const std::vector<int> &v) { 
        std::cout<<"Launching wrapper"<<std::endl;
        computeNumericalAssignment(v); 
      };

    pooling(criteria, inequalities, f_wrapper);
    std::vector<int> increment_offset(components);
    increment_offset[0] = 1;
    std::vector<int> maxima(components, degree - 1);

    // for (int w = 0; w < degree + 1; w++) {
    //     for (int a = 0; a < degree + 1; a++) {
    //         for (int j = acc[w * (degree + 1) + a].get_max_nindex(); j <=
    //         acc[w * (degree + 1) + a].get_max_pindex(); j++) {
    //             std::cout << acc[w * (degree + 1) + a][j] << " ";
    //         }
    //         std::cout << "\t";
    //     }
    //     std::cout << "\n";
    // }
    // std::cout << std::endl;
    // std::cout << accumulatorBlockSizes.size() << " " <<
    // increment_offset.size() << " " << components << " " << maxima.size() <<
    // "\n";
    // Apply final offset: multiply by (-1 + x₀) and truncate
    // Refactored to use direct polynomial operations instead of get/sync pattern
    PolynomialType offset(components, 0);
    std::vector<int> xPowers(components, 0);
    offset.setCoefficient(0, xPowers, -1);  // Coefficient for x₀⁰: -1
    xPowers[0] = 1;
    offset.setCoefficient(0, xPowers, 1);   // Coefficient for x₀¹: +1

    result *= offset;
    result = result.truncate(maxima);
    std::cout<<"Points found: "<<points_found<<std::endl;
    writeResultsToJson(outfile_);

    // for (int w = 0; w < degree + 1; w++) {
    //     for (int a = 0; a < degree + 1; a++) {
    //         for (int j = acc[w * (degree + 1) + a].get_max_nindex(); j <=
    //         acc[w * (degree + 1) + a].get_max_pindex(); j++) {
    //             std::cout << acc[w * (degree + 1) + a][j] << " ";
    //         }
    //         std::cout << "\t";
    //     }
    //     std::cout << "\n";
    // }
    std::cout << std::endl;
  }
};

int main(int argc, char *argv[]) {
  std::cout << argc << std::endl;
  /*
  if (argc < 4){
    std::cerr << "Usage error: input and output file expected in argv" <<
  std::endl; return 1;
  }
  */

  const std::string in = argv[1];
  const std::string out = argv[2];
  std::cout << in << std::endl;
  std::cout << out << std::endl;
  try {
    FK(in, out);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

// implement multithreading at top level of recursion
