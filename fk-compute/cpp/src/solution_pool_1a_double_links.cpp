#include "fk/solution_pool_1a_double_links.hpp"

#include <queue>
#include <stack>

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
        std::cout<<"Satisfies Constraints"<<std::endl;
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

// Helper function to identify bounded variables
struct BoundedVariables {
  std::vector<int> bounded_v;
  int bounded_count;
  std::list<std::array<int, 2>> first;
};

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

// Helper function to process a set of criteria
bool processCriteria(const std::vector<std::vector<double>>& criteria,
                    const std::vector<std::vector<double>>& main_inequalities,
                    const std::vector<std::vector<double>>& supporting_inequalities,
                    const std::function<void(const std::vector<int>&)>& function,
                    int size) {

  auto bounded_info = identifyBoundedVariables(criteria, size);

  if (bounded_info.bounded_count > 0) {
    std::list<std::array<int, 2>> bounds;

    if (bounded_info.bounded_count == size - 1) {
      std::vector<double> degrees = extractDegrees(main_inequalities);
      std::vector<int> point(size - 1, 0);

      // Make a copy of criteria that can be modified
      std::vector<std::vector<double>> criteria_copy = criteria;
      assignVariables(criteria_copy, degrees, criteria_copy, bounded_info.first, bounds,
                     supporting_inequalities, point, function);
      return true;
    }

    bounds = findAdditionalBounds(bounded_info.bounded_v, bounded_info.bounded_count,
                                 size, supporting_inequalities);

    if (bounded_info.bounded_count == size - 1) {
      std::vector<double> degrees = extractDegrees(main_inequalities);
      std::vector<int> point(size - 1, 0);

      // Make a copy of criteria that can be modified
      std::vector<std::vector<double>> criteria_copy = criteria;
      assignVariables(criteria_copy, degrees, criteria_copy, bounded_info.first, bounds,
                     supporting_inequalities, point, function);
      return true;
    }
  }

  return false;
}

void pooling(std::vector<std::vector<double>> main_inequalities,
             std::vector<std::vector<double>> supporting_inequalities,
             const std::function<void(const std::vector<int>&)>& function) {

  int mains = main_inequalities.size();
  int size = main_inequalities[0].size();

  // Add main inequalities to supporting inequalities
  for (const auto& x : main_inequalities) {
    supporting_inequalities.push_back(x);
  }
  int support = supporting_inequalities.size();

  // Initialize visited tracking
  std::vector<btree<double>> visited(mains);
  for (int i = 0; i < mains; i++) {
    visited[i].insertVector(main_inequalities[i]);
  }

  // Try to process the initial criteria
  if (processCriteria(main_inequalities, main_inequalities, supporting_inequalities,
                     function, size)) {
    return;
  }

  // Initialize criteria and queue for iterative processing
  std::vector<std::vector<double>> criteria(mains, std::vector<double>(size));
  std::vector<std::vector<double>> new_criteria(mains, std::vector<double>(size));
  std::queue<std::vector<std::vector<double>>> processing_queue;
  processing_queue.push(main_inequalities);

  // Main processing loop
  while (!processing_queue.empty()) {
    criteria = processing_queue.front();
    processing_queue.pop();

    // Try different supporting inequalities
    for (int i = 0; i < support; i++) {
      bool found_improvement = false;

      // Check each main criterion
      for (int q = 0; q < mains && !found_improvement; q++) {
        // Check each variable coefficient
        for (int j = 1; j < size; j++) {
          if (criteria[q][j] > 0 && supporting_inequalities[i][j] < 0) {
            // Create new criteria by combining with supporting inequality
            for (int k = 0; k < size; k++) {
              new_criteria[q][k] = criteria[q][k] + supporting_inequalities[i][k] / 2.0;
            }

            // Check if this new criterion set has been visited
            bool is_new = false;
            for (int visdex = 0; visdex < mains; visdex++) {
              if (!visited[visdex].containsVector(new_criteria[visdex])) {
                is_new = true;
                break;
              }
            }

            if (is_new) {
              // Mark as visited
              for (int s = 0; s < mains; s++) {
                visited[s].insertVector(new_criteria[s]);
              }

              // Try to process new criteria
              if (processCriteria(new_criteria, main_inequalities, supporting_inequalities,
                                function, size)) {
                return;
              }

              // Add to queue for further processing
              processing_queue.push(new_criteria);
              found_improvement = true;
            }
            break;
          }
        }
      }
    }
  }
}
