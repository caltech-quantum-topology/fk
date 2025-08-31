/**
 * @file solution_pool_1a_double_links.hpp
 * @brief Specialized Integer Linear Programming solver for Gukov-Manolescu invariant computation
 * 
 * This file implements a custom ILP solver that finds ALL feasible integer solutions
 * to systems of linear inequalities arising from braid segment sign constraints.
 * 
 * The solver uses a breadth-first search strategy with degree-bounding criteria:
 * 1. Uses degree-bounding inequalities (criteria) as the basis for BFS
 * 2. Adds supporting inequalities to create relaxations where all free variable 
 *    coefficients in at least one criterion are negative
 * 3. Bounds variables using negative coefficient criteria and inequality constraints
 * 4. Recursively enumerates all integer points in the bounded polytope
 * 5. Verifies solutions against the original unrelaxed ILP
 * 6. Calls provided function for each valid solution
 * 
 * Key features:
 * - Avoids revisiting criteria combinations using btree hashing
 * - Handles variable bounding through negative coefficient detection
 * - Uses direct recursive enumeration for bounded subproblems
 */

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <list>
#include <functional>

#include "btree.hpp"

// TODO: Handle criterion-bounded variable overlap more intelligently
// Currently degrees are decremented for single criteria only, but when multiple
// criteria bound the same variable, this leads to unnecessary computations.
// Can be optimized by decrementing degrees of other criteria before calling recurse_1.

/**
 * @brief Detect variables that can be bounded by criteria with all-negative coefficients
 * @param criteria_set The criteria to analyze
 * @param bounded_v Output array of bounded variable flags
 * @param first Output list of [variable_index, criterion_index] pairs (1-indexed variables)
 * @return Number of variables that could be bounded
 */
int detect_bounded_variables(
    const std::vector<std::vector<double>>& criteria_set,
    std::vector<int>& bounded_v,
    std::list<std::array<int, 2>>& first
) {
    int size = criteria_set[0].size();
    int mains = criteria_set.size();
    int bounded = 0;
    
    bounded_v.assign(size - 1, false);
    first.clear();
    
    // Check each criterion for variables with all-negative coefficients
    for (int i = 0; i < mains; i++) {
        bool condition = true;  // True if all variable coefficients are negative
        std::vector<bool> locally_bounded(size - 1, false);
        
        // Check each variable coefficient in this criterion
        for (int k = 1; k < size; k++) {
            if (criteria_set[i][k] > 0) {
                condition = false;  // Positive coefficient disqualifies criterion
                break; // Can break early since criterion won't be useful
            }
            else if (criteria_set[i][k] < 0) {
                locally_bounded[k - 1] = true;  // Negative coeff can bound this variable
            }
        }
        
        // If criterion has all negative coefficients, mark bounded variables
        if (condition) {
            for (int v = 0; v < size - 1; v++) {
                if (locally_bounded[v] && !bounded_v[v]) {
                    bounded_v[v] = true;
                    first.push_back({v + 1, i});  // Note: 1-indexed for recurse_1
                    bounded++;
                }
            }
        }
    }
    
    return bounded;
}

/**
 * @brief Try to bound remaining unbounded variables using supporting inequalities
 * @param bounded_v Input/output array of bounded variable flags
 * @param bounds Output list of [variable_index, inequality_index] pairs (0-indexed variables)
 * @param supporting_inequalities All inequality constraints
 * @return Updated count of bounded variables
 */
int try_bound_with_inequalities(
    std::vector<int>& bounded_v,
    std::list<std::array<int, 2>>& bounds,
    const std::vector<std::vector<double>>& supporting_inequalities
) {
    int size = bounded_v.size() + 1;  // +1 for constant term
    int support = supporting_inequalities.size();
    int bounded = std::count(bounded_v.begin(), bounded_v.end(), true);
    
    bounds.clear();
    
    // Try to bound remaining unbounded variables using supporting inequalities
    int index = 0;
    while (index < size - 1) {
        if (!bounded_v[index]) { // Found an unbounded variable
            // Search through all supporting inequalities for potential bounds
            for (int l = 0; l < support; l++) {
                // Check if this inequality has negative coefficient for current variable (can bound it)
                if (supporting_inequalities[l][1 + index] < 0) {
                    bool useful = true;
                    // Check if inequality is "useful" - no positive coefficients for other unbounded variables
                    for (int n = 0; n < size - 1; n++) {
                        if (n != index && supporting_inequalities[l][1 + n] > 0 && !bounded_v[n]) {
                            useful = false;  // Positive coeff on unbounded var makes inequality not useful for bounding
                        }
                    }
                    if (useful) {
                        bounds.push_back({index, l});  // Add this variable-inequality pair to bounds list (0-indexed)
                        bounded_v[index] = true;
                        bounded++;
                        if (bounded == size - 1) {  // All variables now bounded!
                            return bounded;
                        }
                        index = -1;  // Restart search from beginning
                        break;
                    }
                }
            }
        }
        index++;
    }
    
    return bounded;
}

/**
 * @brief Solve the ILP when all variables are bounded
 * @param criteria Current criteria set
 * @param main_inequalities Original degree-bounding criteria
 * @param first Variables bounded by criteria
 * @param bounds Variables bounded by supporting inequalities  
 * @param supporting_inequalities All inequality constraints
 * @param function Callback function for feasible solutions
 */
void solve_bounded_problem(
    const std::vector<std::vector<double>>& criteria,
    const std::vector<std::vector<double>>& main_inequalities,
    const std::list<std::array<int, 2>>& first,
    const std::list<std::array<int, 2>>& bounds,
    const std::vector<std::vector<double>>& supporting_inequalities,
    const std::function<void(const std::vector<int>&)>& function
) {
    int size = criteria[0].size();
    
    // Extract degree bounds from criteria
    std::vector<double> degrees = {};
    for (auto x : criteria) {
        degrees.push_back(x[0]);
    }
    
    // Initialize solution point
    std::vector<int> point(size - 1, 0);
    
    // Call recursive solver
    recurse_1(
        const_cast<std::vector<std::vector<double>>&>(criteria), 
        degrees, 
        const_cast<std::vector<std::vector<double>>&>(main_inequalities),
        first, 
        bounds, 
        const_cast<std::vector<std::vector<double>>&>(supporting_inequalities), 
        point, 
        function
    );
}

/**
 * @brief Recursive function to enumerate integer solutions for remaining unbounded variables
 * @param criteria Original degree-bounding criteria constraints  
 * @param bounds List of [variable_index, inequality_index] pairs for variables bounded by supporting inequalities
 * @param supporting_inequalities All inequality constraints (criteria + additional inequalities)
 * @param point Current partial integer solution being constructed
 * @param function Callback function to execute for each feasible solution
 * 
 * This function handles variables that couldn't be bounded by criteria (negative coefficients)
 * and instead uses supporting inequalities to determine bounds. For each bounded variable,
 * it computes the upper bound and recursively tries all integer values from 0 to upper bound.
 * 
 * IMPORTANT: Variable indices in this function are 0-indexed, but inequality coefficients
 * start at index 1, so we use (1 + index) to access variable coefficients in inequalities.
 */
void recurse_2(
    std::vector<std::vector<double>>& criteria,
    std::list<std::array<int, 2>> bounds, 
    std::vector<std::vector<double>> supporting_inequalities,
    std::vector<int> point,
    const std::function<void(const std::vector<int>&)>& function
) 
{
    if (bounds.size()) {
        // Process next variable bounded by supporting inequality
        int index = bounds.front()[0];        // Variable index (0-indexed)
        int inequality = bounds.front()[1];   // Supporting inequality index
        bounds.pop_front();
        
        // Compute upper bound by mandating: inequality's constant + sum(coeff*var) >= 0
        // Rearranging: var[index] <= (constant + sum(other_terms)) / (-coeff[index])
        int upper = supporting_inequalities[inequality][0];  // Start with constant term
        for (int i = 0; i < point.size(); i++) {
            if (i != index) {
                // Note: variables are 0-indexed but coefficients start at index 1
                upper += supporting_inequalities[inequality][1 + i] * point[i];
            }
        }
        // Divide by negative coefficient (variables 0-indexed, coefficients at 1+index)
        upper /= -supporting_inequalities[inequality][1 + index];
        
        // Recursively try all integer values from 0 to upper bound
        for (int i = 0; i <= upper; i++) {
            point[index] = i;  
            recurse_2(
                criteria,
                bounds, 
                supporting_inequalities, 
                point, 
                function
            );
        }
    }
    else {
        // All variables assigned - verify feasibility against original unrelaxed ILP
        bool cond = true;
        
        // Check feasibility against all supporting inequalities
        for (int i = 0; i < supporting_inequalities.size(); i++) {
            int acc = supporting_inequalities[i][0];  // Start with constant term
            for (int j = 0; j < point.size(); j++){
                // Variables 0-indexed, coefficients at 1+j  
                acc += point[j] * supporting_inequalities[i][1 + j];
            }
            if (acc < 0) {  // Inequality violated
                cond = false;
                break;
            }
        }
        
        // If supporting inequalities satisfied, check original criteria
        if (cond) {
            for (auto x : criteria) {
                int acc = x[0];  // Start with constant term
                for (int j = 0; j < point.size(); j++){
                    // Variables 0-indexed, coefficients at 1+j
                    acc += point[j] * x[1 + j];
                }
                if (acc < 0) {  // Criterion violated
                    cond = false;
                    break;
                }
            }
        }
        
        // If all constraints satisfied, execute callback with feasible solution
        if (cond) {
            function(point);
        }
    }
}

/**
 * @brief Recursive function to handle variables bounded by criteria (negative coefficients)
 * @param new_criteria Modified criteria after adding supporting inequalities
 * @param degrees Current degree bounds for each criterion  
 * @param criteria Original degree-bounding criteria
 * @param first List of [variable_index, criterion_index] pairs for variables bounded by criteria
 * @param bounds List of [variable_index, inequality_index] pairs for variables bounded by supporting inequalities
 * @param supporting_inequalities All inequality constraints  
 * @param point Current partial integer solution being constructed
 * @param function Callback function to execute for each feasible solution
 * 
 * This function processes variables that can be bounded by criteria (those with negative
 * coefficients). It uses the slope of each criterion to determine bounds and recursively
 * enumerates integer values, updating degree bounds as it proceeds.
 * 
 * IMPORTANT: Variable indices in this function are 1-indexed (used directly as inequality
 * coefficient indices), unlike recurse_2 which uses 0-indexed variables.
 */
void recurse_1(
    std::vector<std::vector<double>>& new_criteria, 
    std::vector<double> degrees,
    std::vector<std::vector<double>>& criteria,
    std::list<std::array<int, 2>> first, 
    std::list<std::array<int, 2>> bounds, 
    std::vector<std::vector<double>> supporting_inequalities,
    std::vector<int> point,
    const std::function<void(const std::vector<int>&)>& function
) 
{
    if (first.size()) {
        // Process next variable bounded by a criterion with negative coefficient
        int var_index = first.front()[0];   // Variable index (1-indexed)
        int main_index = first.front()[1];  // Criterion index
        double slope = -new_criteria[main_index][var_index];  // Negative of coefficient (positive slope)
        first.pop_front();
        
        // Copy degree bounds to update as we recurse
        std::vector<double> new_degrees = degrees;
        
        // Try all integer values from 0 up to degree bound divided by slope
        for (int i = 0;  i <= degrees[main_index] / slope; i++) {
            point[var_index - 1] = i;  // Convert 1-indexed var_index to 0-indexed for point array
            recurse_1(
                new_criteria, 
                new_degrees, 
                criteria,
                first, 
                bounds, 
                supporting_inequalities, 
                point, 
                function
            );
            new_degrees[main_index] -= slope;  // Decrement degree bound by slope for next iteration
        }
    }
    else {
        // All criterion-bounded variables processed, handle remaining variables via recurse_2
        recurse_2(
            criteria, 
            bounds, 
            supporting_inequalities, 
            point, 
            function
        );
    }
}

/**
 * @brief Main ILP solver entry point using breadth-first search with degree-bounding criteria
 * @param main_inequalities Degree-bounding criteria constraints (coefficients: [constant, var1, var2, ...])
 * @param supporting_inequalities Additional inequality constraints in same format
 * @param function Callback function to execute for each feasible integer solution
 * 
 * Algorithm:
 * 1. Performs BFS by adding supporting inequalities to criteria to create relaxations
 * 2. Seeks criteria combinations where at least one criterion has all negative free variable coefficients
 * 3. Uses variables bounded by one or more criteria to extend boundedness to all variables
 * 4. If all variables can be bounded, uses recursive enumeration to find integer solutions
 * 5. Verifies each solution against original unrelaxed ILP before calling function
 * 
 * Uses btree memoization to avoid revisiting the same criteria combinations.
 */
void pooling(
    std::vector<std::vector<double>> main_inequalities, 
    std::vector<std::vector<double>> supporting_inequalities,
    const std::function<void(const std::vector<int>&)>& function
) 
{

    // Initialize solver data structures
    int mains = main_inequalities.size();      // Number of degree-bounding criteria
    int size = main_inequalities[0].size();    // Number of variables + 1 (for constant term)
    
    // Combine criteria and supporting inequalities into single constraint set
    for (auto x : main_inequalities) {
        supporting_inequalities.push_back(x);
    }
    int support = supporting_inequalities.size();  // Total number of constraints
    
    // Initialize visited criteria tracking using btree hashing
    std::vector<btree<double>> visited(mains);
    for (int i = 0; i < mains; i++) {
        visited[i].update(main_inequalities[i]);  // Mark original criteria as visited
    }
    // Track which variables are bounded and how
    std::vector<int> bounded_v(size - 1);       // Flags for bounded variables
    int bounded = 0;                            // Count of bounded variables
    std::list<std::array<int, 2>> first = {};   // Variables bounded by criteria [var_index, criterion_index]
    
    // Detect variables that can be bounded by criteria with all-negative coefficients
    bounded = detect_bounded_variables(main_inequalities, bounded_v, first);
    if (bounded > 0) {
        std::list<std::array<int, 2>> bounds = {};
        if (bounded == size - 1) {
            // All variables bounded by criteria - solve directly
            solve_bounded_problem(main_inequalities, main_inequalities, first, bounds, supporting_inequalities, function);
            return;
        }
        // Try to bound remaining variables with supporting inequalities
        bounded = try_bound_with_inequalities(bounded_v, bounds, supporting_inequalities);
        if (bounded == size - 1) {
            // All variables now bounded - solve directly
            solve_bounded_problem(main_inequalities, main_inequalities, first, bounds, supporting_inequalities, function);
            return;
        }
    }
    // Begin BFS: Create relaxations by adding supporting inequalities to criteria
    std::vector<std::vector<double>> criteria(mains, std::vector<double>(size));
    std::vector<std::vector<double>> new_criteria(mains, std::vector<double>(size)); 
    std::list<std::vector<std::vector<double>>> queue = {main_inequalities};  // BFS queue initialized with original criteria
    
    while (true) {
        criteria = queue.front();  // Get next criteria combination to explore
        queue.pop_front();
        // Try adding each supporting inequality to create new relaxations
        for (int i = 0; i < support; i++) {
            bool cond = false;
            int q = 0;
            while (q < mains) {
                // Look for positive coefficients in criteria that can be "fixed" by negative coefficients in supporting inequality
                for (int j = 1 ; j < size; j++) {
                    if (criteria[q][j] > 0 && supporting_inequalities[i][j] < 0) {
                        // Create new criterion by adding half of supporting inequality to current criterion
                        for (int k = 0; k < size; k++) {
                            new_criteria[q][k] = criteria[q][k] + supporting_inequalities[i][k] / 2.0;  
                        }
                        // Check if this new criteria combination has been visited before
                        for (int visdex = 0; visdex < mains; visdex++) {
                            if (!visited[visdex].contains(new_criteria[visdex])) {
                                // Mark new criteria as visited to avoid revisiting
                                for (int s = 0; s < mains; s++) {
                                    visited[s].update(new_criteria[s]);
                                }
                                queue.push_back(new_criteria);  // Add new relaxation to BFS queue
                                // Test new relaxation: check if it provides better variable bounding
                                std::vector<int> bounded_v(size - 1);
                                std::list<std::array<int, 2>> first = {};
                                int bounded = detect_bounded_variables(new_criteria, bounded_v, first);
                                if (bounded > 0) {
                                    std::list<std::array<int, 2>> bounds = {};
                                    if (bounded == size - 1) {
                                        // All variables bounded by new criteria - solve directly
                                        solve_bounded_problem(new_criteria, main_inequalities, first, bounds, supporting_inequalities, function);
                                        return;
                                    }
                                    // Try to bound remaining variables with supporting inequalities
                                    bounded = try_bound_with_inequalities(bounded_v, bounds, supporting_inequalities);
                                    if (bounded == size - 1) {
                                        // All variables now bounded - solve with new criteria
                                        solve_bounded_problem(new_criteria, main_inequalities, first, bounds, supporting_inequalities, function);
                                        return;
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
                q++;
            }
        }       
    }
}
