#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <list>
#include <functional>

#include "btree.hpp"

void recurse_2(
    std::vector<std::vector<double>>& criteria,
    std::list<std::array<int, 2>> bounds, 
    std::vector<std::vector<double>> supporting_inequalities,
    std::vector<int> point,
    const std::function<void(const std::vector<int>&)>& function
) 
{ // make "bounds" a list! // code might be displaying weird behavior because you forgot to take into account the unbounded variables that can still be used if they don't enter pathologically into the inequalities. Make sure, separately, that you cannot have loops! It should be okay if two variables can be bounded by a single inequality sequentially, but maybe there are some sus cases. 
    if (bounds.size()) {
        int index = bounds.front()[0];
        int inequality = bounds.front()[1];
        bounds.pop_front();
        int upper = supporting_inequalities[inequality][0];
        for (int i = 0; i < point.size(); i++) {
            if (i != index) {
                upper += supporting_inequalities[inequality][1 + i] * point[i];
            }
        }
        upper /= -supporting_inequalities[inequality][1 + index];
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
        bool cond = true;
        for (int i = 0; i < supporting_inequalities.size(); i++) {
            int acc = supporting_inequalities[i][0];
            for (int j = 0; j < point.size(); j++){
                acc += point[j] * supporting_inequalities[i][1 + j];
            }
            if (acc < 0) {
                cond = false;
                break;
            }
        }
        if (cond) {
            for (auto x : criteria) {
                int acc = x[0];
                for (int j = 0; j < point.size(); j++){
                    acc += point[j] * x[1 + j];
                }
                if (acc < 0) {
                    cond = false;
                    break;
                }
            }
        }
        if (cond) {
            function(point);
        }
    }
}

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
        int var_index = first.front()[0];
        int main_index = first.front()[1];
        double slope = -new_criteria[main_index][var_index];
        first.pop_front();
        std::vector<double> new_degrees = degrees;
        for (int i = 0;  i <= degrees[main_index] / slope; i++) {
            point[var_index - 1] = i;
            new_degrees[main_index] = degrees[main_index] - i * slope; // eventually change this to constant decrements by "slope", to avoid re-accessing the "degrees" vector
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
        }
    }
    else {
        recurse_2(
            criteria, 
            bounds, 
            supporting_inequalities, 
            point, 
            function
        );
    }
}

// NEED TO HANDLE CASE WHEN CRITERION-BOUNDED VARIABLES OVERLAP, LEADING TO INCONSISTENCIES BETWEEN CRITERIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void pooling(
    std::vector<std::vector<double>> main_inequalities, 
    std::vector<std::vector<double>> supporting_inequalities,
    const std::function<void(const std::vector<int>&)>& function
) 
{
    int mains = main_inequalities.size();
    int size = main_inequalities[0].size();
    for (auto x : main_inequalities) {
        supporting_inequalities.push_back(x);
    }
    int support = supporting_inequalities.size();
    std::vector<btree<double>> visited(mains);
    for (int i = 0; i < mains; i++) {
        visited[i].update(main_inequalities[i]);
    }
    std::vector<int> bounded_v(size - 1);
    int bounded = 0;
    std::list<std::array<int, 2>> first = {};
    for (int i = 0; i < mains; i++) {
        bool condition = true;
        std::vector<bool> locally_bounded(size - 1, false);
        for (int k = 1; k < size; k++) {
            if (main_inequalities[i][k] > 0) {
                condition = false;
                break; // can break here, but cannot below, as we need to save the created criteria
            }
            else if (main_inequalities[i][k] < 0) {
                locally_bounded[k - 1] = true;
            }
        }
        if (condition) {
            for (int v = 0; v < size - 1; v++) {
                if (locally_bounded[v] && !bounded_v[v]) {
                    bounded_v[v] = true;
                    first.push_back({v + 1, i});
                    bounded++;
                }
            }
        }
    }
    if (bounded > 0) {
        std::list<std::array<int, 2>> bounds = {};
        if (bounded == size - 1) {
            std::vector<double> degrees = {};
            for (auto x : main_inequalities) {
                degrees.push_back(x[0]);
            }
            std::vector<int> point(size - 1, 0);
            recurse_1(
                main_inequalities, 
                degrees, 
                main_inequalities,
                first, 
                bounds, 
                supporting_inequalities, 
                point, 
                function
            );
            return;
        }
        int index = 0; // issue: you're only searching through supporting inequalities here, but you also want to search through criteria to bound the variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        while (index < size - 1) {
            if (!bounded_v[index]) { // we find a vanishing variable
                for (int l = 0; l < support; l++) {  // we search for inequalities to bound vanishing variable
                    if (supporting_inequalities[l][1 + index] < 0) {  // we can possibly bound vanishing variable with inequalities
                        bool useful = true;
                        for (int n = 0; n < size - 1; n++) {
                            if (n != index && supporting_inequalities[l][1 + n] > 0 && !bounded_v[n]) {
                                useful = false;
                            }
                        }
                        if (useful) {
                            bounds.push_back({index, l});
                            bounded_v[index] = true;
                            bounded++;
                            if (bounded == size - 1) {
                                std::vector<double> degrees = {};
                                for (auto x : main_inequalities) {
                                    degrees.push_back(x[0]);
                                }
                                std::vector<int> point(size - 1, 0);
                                recurse_1(
                                    main_inequalities, 
                                    degrees, 
                                    main_inequalities,
                                    first,
                                    bounds, 
                                    supporting_inequalities, 
                                    point, 
                                    function
                                );
                                return;
                            }
                            index = -1;
                            break;
                        }
                    }
                }
            }
            index++;
        }
    }
    std::vector<std::vector<double>> criteria(mains, std::vector<double>(size));
    std::vector<std::vector<double>> new_criteria(mains, std::vector<double>(size));
    std::list<std::vector<std::vector<double>>> queue = {main_inequalities};
    while (true) {
        criteria = queue.front();
        queue.pop_front();
        for (int i = 0; i < support; i++) {
            bool cond = false;
            int q = 0;
            while (q < mains) {
                for (int j = 1 ; j < size; j++) {
                    if (criteria[q][j] > 0 && supporting_inequalities[i][j] < 0) {
                        for (int k = 0; k < size; k++) {
                            new_criteria[q][k] = criteria[q][k] + supporting_inequalities[i][k] / 2.0;  
                        }
                        for (int visdex = 0; visdex < mains; visdex++) {
                            if (!visited[visdex].contains(new_criteria[visdex])) {
                                for (int s = 0; s < mains; s++) {
                                    visited[s].update(new_criteria[s]);
                                }
                                queue.push_back(new_criteria);
                                std::vector<int> bounded_v(size - 1);
                                int bounded = 0;
                                std::list<std::array<int, 2>> first = {};
                                for (int l = 0; l < mains; l++) {
                                    bool condition = true;
                                    std::vector<bool> locally_bounded(size - 1, false);
                                    for (int k = 1; k < size; k++) {
                                        if (new_criteria[l][k] > 0) {
                                            condition = false;
                                            break; // can break here, but cannot below, as we need to save the created criteria
                                        }
                                        else if (new_criteria[l][k] < 0) {
                                            locally_bounded[k - 1] = true;
                                        }
                                    }
                                    if (condition) {
                                        for (int v = 0; v < size - 1; v++) {
                                            if (locally_bounded[v] && !bounded_v[v]) {
                                                bounded_v[v] = true;
                                                first.push_back({v + 1, l});
                                                bounded++;
                                            }
                                        }
                                    }
                                }
                                if (bounded > 0) {
                                    std::list<std::array<int, 2>> bounds = {};
                                    if (bounded == size - 1) {
                                        std::vector<double> degrees = {};
                                        for (auto x : new_criteria) {
                                            degrees.push_back(x[0]);
                                        }
                                        std::vector<int> point(size - 1, 0);
                                        recurse_1(
                                            new_criteria, 
                                            degrees, 
                                            criteria,
                                            first, 
                                            bounds, 
                                            supporting_inequalities, 
                                            point, 
                                            function
                                        );
                                        return;
                                    }
                                    int index = 0;
                                    while (index < size - 1) {
                                        if (!bounded_v[index]) {
                                            for (int l = 0; l < support; l++) {  // we search for inequalities to bound vanishing variable
                                                if (supporting_inequalities[l][1 + index] < 0) {  // we can possibly bound vanishing variable with inequalities
                                                    bool useful = true;
                                                    for (int n = 0; n < size - 1; n++) {
                                                        if (n != index && supporting_inequalities[l][1 + n] > 0 && !bounded_v[n]) {
                                                            useful = false;
                                                        }
                                                    }
                                                    if (useful) {
                                                        bounds.push_back({index, l});
                                                        bounded_v[index] = true;
                                                        bounded++;
                                                        if (bounded == size - 1) {
                                                            std::vector<double> degrees = {};
                                                            for (auto x : main_inequalities) {
                                                                degrees.push_back(x[0]);
                                                            }
                                                            std::vector<int> point(size - 1, 0);
                                                            recurse_1(
                                                                new_criteria, 
                                                                degrees, 
                                                                criteria,
                                                                first,
                                                                bounds, 
                                                                supporting_inequalities, 
                                                                point, 
                                                                function
                                                            );
                                                            return;
                                                        }
                                                        index = -1;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                        index++;
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

// when leisurely, add minimum (of individual variables) bounding

// does the "inequalities cannot be cyclic if they're not contradictory" rule apply to general ILP problems?

// check stl for multi-threaded manipulations on iterables (comparisons in particular; think python's "or" on sets (is this even multi-threaded?))

// !!! make sure to add scaling and translations

// !!! add inequality redundancy reducers

// ! extract python version of ILP feasible region solver (and any other potentially generally useful algorithms) out of relations.py=
