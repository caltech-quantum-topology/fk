#pragma once

#include <array>
#include <functional>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "fk/btree.hpp"

// Structure to hold assignment results
struct AssignmentResult {
  std::vector<std::vector<double>> criteria;
  std::list<std::array<int, 2>> bounds;
  std::vector<std::vector<double>> supporting_inequalities;
  std::vector<int> point;
  bool is_ready;

  AssignmentResult() : is_ready(false) {}
};

// Modern iterative versions of the algorithms
void enumeratePoints(std::vector<std::vector<double>>& criteria,
                    std::list<std::array<int, 2>> bounds,
                    std::vector<std::vector<double>> supporting_inequalities,
                    std::vector<int> point,
                    const std::function<void(const std::vector<int>&)>& function);

void assignVariables(std::vector<std::vector<double>>& new_criteria,
                    std::vector<double> degrees,
                    std::vector<std::vector<double>>& criteria,
                    std::list<std::array<int, 2>> first,
                    std::list<std::array<int, 2>> bounds,
                    std::vector<std::vector<double>> supporting_inequalities,
                    std::vector<int> point,
                    const std::function<void(const std::vector<int>&)>& function);

// Helper function to process a set of criteria
bool processCriteria(const std::vector<std::vector<double>>& criteria,
                    const std::vector<std::vector<double>>& main_inequalities,
                    const std::vector<std::vector<double>>& supporting_inequalities,
                    const std::function<void(const std::vector<int>&)>& function,
                    int size);

// NEED TO HANDLE CASE WHEN CRITERION-BOUNDED VARIABLES OVERLAP, LEADING TO
// INCONSISTENCIES BETWEEN CRITERIA
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void pooling(std::vector<std::vector<double>> main_inequalities,
             std::vector<std::vector<double>> supporting_inequalities,
             const std::function<void(const std::vector<int> &)> &function);

// when leisurely, add minimum (of individual variables) bounding

// does the "inequalities cannot be cyclic if they're not contradictory" rule
// apply to general ILP problems?

// check stl for multi-threaded manipulations on iterables (comparisons in
// particular; think python's "or" on sets (is this even multi-threaded?))

// !!! make sure to add scaling and translations

// !!! add inequality redundancy reducers

// ! extract python version of ILP feasible region solver (and any other
// potentially generally useful algorithms) out of relations.py=