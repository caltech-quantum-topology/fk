#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <list>
#include <functional>

#include "fk/btree.hpp"

void recurse_2(
    std::vector<std::vector<double>>& criteria,
    std::list<std::array<int, 2>> bounds,
    std::vector<std::vector<double>> supporting_inequalities,
    std::vector<int> point,
    const std::function<void(const std::vector<int>&)>& function
);

void recurse_1(
    std::vector<std::vector<double>>& new_criteria,
    std::vector<double> degrees,
    std::vector<std::vector<double>>& criteria,
    std::list<std::array<int, 2>> first,
    std::list<std::array<int, 2>> bounds,
    std::vector<std::vector<double>> supporting_inequalities,
    std::vector<int> point,
    const std::function<void(const std::vector<int>&)>& function
);

// NEED TO HANDLE CASE WHEN CRITERION-BOUNDED VARIABLES OVERLAP, LEADING TO INCONSISTENCIES BETWEEN CRITERIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void pooling(
    std::vector<std::vector<double>> main_inequalities,
    std::vector<std::vector<double>> supporting_inequalities,
    const std::function<void(const std::vector<int>&)>& function
);

// when leisurely, add minimum (of individual variables) bounding

// does the "inequalities cannot be cyclic if they're not contradictory" rule apply to general ILP problems?

// check stl for multi-threaded manipulations on iterables (comparisons in particular; think python's "or" on sets (is this even multi-threaded?))

// !!! make sure to add scaling and translations

// !!! add inequality redundancy reducers

// ! extract python version of ILP feasible region solver (and any other potentially generally useful algorithms) out of relations.py=