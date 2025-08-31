/**
 * @file qalg_links.hpp
 * @brief Quantum algebraic operations for Gukov-Manolescu invariant computation
 * 
 * This file implements quantum binomial coefficients and q-Pochhammer symbols
 * used in the Gukov-Manolescu invariant calculation. These functions apply
 * quantum algebraic operations to polynomial terms according to Kim's R-matrix
 * prescription for crossing contributions.
 * 
 * Key quantum algebraic objects:
 * - q-binomial coefficients: [n choose k]_q for positive-positive and negative-positive cases
 * - q-Pochhammer symbols: (x*q)_n for regular and inverse cases
 * - Matrix-as-vector operations on polynomial terms with multi-component indexing
 */

#pragma once

#include <vector>
#include <functional>

#include "bilvector.hpp"
#include "linalg.hpp"

// TODO: Consider saving multiplicands as separate variables before multiplying by term

/**
 * @brief Helper function to compute positive-positive q-binomial coefficient [i choose j]_q
 * @param binom Output array to accumulate binomial coefficient terms  
 * @param i Upper parameter of binomial coefficient
 * @param j Lower parameter of binomial coefficient
 * @param shift q-degree offset for coefficient placement
 * 
 * Recursively computes q-binomial coefficient using the recurrence:
 * [i choose j]_q = q^j * [i-1 choose j]_q + [i-1 choose j-1]_q
 * with base cases [i choose i]_q = 1 and [i choose 0]_q = 1
 */
void pp_q_binom_ (std::vector<int>& binom, int i, int j, int shift) {
    if (i == j) {
        binom[shift] += 1;  // Base case: [i choose i]_q = 1
    }
    else if (j == 0) {
        binom[shift] += 1;  // Base case: [i choose 0]_q = 1
    }
    else {
        // Recurrence: [i choose j]_q = q^j * [i-1 choose j]_q + [i-1 choose j-1]_q
        pp_q_binom_(binom, i - 1, j, shift + j);      // q^j term
        pp_q_binom_(binom, i - 1, j - 1, shift);      // Constant term
    }
}

/**
 * @brief Apply positive-positive q-binomial coefficient to polynomial term
 * @param term Input/output polynomial term (modified in place)
 * @param i Upper parameter of q-binomial coefficient  
 * @param j Lower parameter of q-binomial coefficient
 * @param neg If true, multiply by q^(-k), if false multiply by q^k
 * 
 * Computes [i choose j]_q and multiplies it with the polynomial term.
 * The neg parameter controls whether q-powers are added (false) or subtracted (true).
 */
void pp_q_binom (std::vector<bilvector<int>>& term, int i, int j, bool neg) {
    int Q_MAX_DEGREE = j * (i - j);  // Maximum q-degree in the binomial expansion
    std::vector<int> binom(Q_MAX_DEGREE + 1, 0);
    
    // Compute q-binomial coefficient
    if (i == j) {
        binom[0] = 1;  // [i choose i]_q = 1
    }
    else if (j == 0) {
        binom[0] = 1;  // [i choose 0]_q = 1
    }
    else {
        // Use recursive helper to compute coefficient
        pp_q_binom_(binom, i - 1, j, j);      // q^j * [i-1 choose j]_q  
        pp_q_binom_(binom, i - 1, j - 1, 0);  // [i-1 choose j-1]_q
    }
    // Apply q-binomial coefficient to polynomial term
    if (neg) {
        // Negative case: multiply term by q^(-k) * binom[k] (subtract q-powers)
        bilvector<int> term_(term[0].nsize(), term[0].psize(), term[0].get_component_size(), 0);
        
        // Copy original term and clear
        for (int j = term[0].get_max_nindex(); j <= term[0].get_max_pindex(); j++) {
            term_[j] = term[0][j];
            term[0][j] = 0;
        }
        
        // Apply binomial: sum over k of binom[k] * coeff[j] * q^(j-k)
        for (int j = term_.get_max_nindex(); j <= term_.get_max_pindex(); j++) {
            if (term_[j] != 0) {
                for (int k = 0; k < Q_MAX_DEGREE + 1; k++) {
                    term[0][j - k] += binom[k] * term_[j];  // Subtract k from q-power
                }
            }
        }
    }
    else {
        // Positive case: multiply term by q^k * binom[k] (add q-powers)
        bilvector<int> term_(term[0].nsize(), term[0].psize(), term[0].get_component_size(), 0);
        
        // Copy original term and clear
        for (int j = term[0].get_max_nindex(); j <= term[0].get_max_pindex(); j++) {
            term_[j] = term[0][j];
            term[0][j] = 0;
        }
        
        // Apply binomial: sum over k of binom[k] * coeff[j] * q^(j+k)
        for (int j = term_.get_max_nindex(); j <= term_.get_max_pindex(); j++) {
            if (term_[j] != 0) {
                for (int k = 0; k < Q_MAX_DEGREE + 1; k++) {
                    term[0][j + k] += binom[k] * term_[j];  // Add k to q-power
                }
            }
        }
    }
}

/**
 * @brief Helper function to compute negative-positive q-binomial coefficient [i choose j]_q  
 * @param binom Output array to accumulate binomial coefficient terms
 * @param i Upper parameter (can be negative)
 * @param j Lower parameter (positive)
 * @param shift q-degree offset for coefficient placement
 * @param neg Sign flag for current term
 * 
 * Recursively computes q-binomial coefficient for cases where i can be negative,
 * using specialized recurrence relations and sign tracking.
 */
void np_q_binom_ (std::vector<int>& binom, int i, int j, int shift, bool neg) {
    if (j == 0) {
        // Base case: [i choose 0]_q = 1 (with sign)
        if (neg) {
            binom[shift] += -1;
        }
        else {
            binom[shift] += 1;
        }
    }
    else if (i == -1) {
        // Special case for i = -1: use recurrence with sign flip
        np_q_binom_(binom, -1, j - 1, shift - j, !neg);
    }
    else {
        // General recurrence for negative-positive case with sign tracking
        np_q_binom_(binom, i, j - 1, shift + 1 + i - j, !neg);
        np_q_binom_(binom, i + 1, j, shift, neg);
    }
}

/**
 * @brief Apply negative-positive q-binomial coefficient to polynomial term
 * @param term Input/output polynomial term (modified in place)
 * @param i Upper parameter (can be negative)
 * @param j Lower parameter (positive)
 * @param neg If true, multiply by q^(-k), if false multiply by q^k
 * 
 * Computes [i choose j]_q for negative-positive case and multiplies it with polynomial term.
 * Handles negative upper parameter using specialized degree calculations and offsets.
 */
void np_q_binom (std::vector<bilvector<int>>& term, int i, int j, bool neg) {
    int Q_DEGREE_DELTA = -(1 + i) * j;      // Degree range for negative case
    int Q_MAX_DEGREE = -j * (j + 1) / 2;    // Maximum negative degree
    std::vector<int> binom(Q_DEGREE_DELTA + 1, 0);
    
    // Compute negative-positive q-binomial coefficient
    if (j == 0) {
        binom[0] = 1;  // [i choose 0]_q = 1
    }
    else if (i == -1) {
        // Special case: i = -1
        np_q_binom_(binom, -1, j - 1, Q_DEGREE_DELTA - Q_MAX_DEGREE - j, true);
    }
    else {
        // General case: use recurrence with proper offset calculations
        np_q_binom_(binom, i, j - 1, Q_DEGREE_DELTA - Q_MAX_DEGREE + 1 + i - j, true);
        np_q_binom_(binom, i + 1, j, Q_DEGREE_DELTA - Q_MAX_DEGREE, false);
    }
    // Apply negative-positive q-binomial coefficient to polynomial term
    if (neg) {
        // Negative case: subtract q-powers with offset correction
        bilvector<int> term_(term[0].nsize(), term[0].psize(), term[0].get_component_size(), 0);
        
        // Copy original term and clear
        for (int j = term[0].get_max_nindex(); j <= term[0].get_max_pindex(); j++) {
            term_[j] = term[0][j];
            term[0][j] = 0;
        }
        
        // Apply binomial with negative case offset
        for (int j = term_.get_max_nindex(); j <= term_.get_max_pindex(); j++) {
            if (term_[j] != 0) {
                for (int k = 0; k < Q_DEGREE_DELTA + 1; k++) {
                    term[0][j - k + Q_DEGREE_DELTA - Q_MAX_DEGREE] += binom[k] * term_[j];
                }
            }
        }
    }
    else {
        // Positive case: add q-powers with offset correction
        bilvector<int> term_(term[0].nsize(), term[0].psize(), term[0].get_component_size(), 0);
        
        // Copy original term and clear
        for (int j = term[0].get_max_nindex(); j <= term[0].get_max_pindex(); j++) {
            term_[j] = term[0][j];
            term[0][j] = 0;
        }
        
        // Apply binomial with positive case offset
        for (int j = term_.get_max_nindex(); j <= term_.get_max_pindex(); j++) {
            if (term_[j] != 0) {
                for (int k = 0; k < Q_DEGREE_DELTA + 1; k++) {
                    term[0][j + k - Q_DEGREE_DELTA + Q_MAX_DEGREE] += binom[k] * term_[j];
                }
            }  
        }
    }
}

/**
 * @brief Apply q-Pochhammer symbol (x*q)_n to polynomial term
 * @param term Input/output polynomial term (modified in place)
 * @param up Upper limit of q-power range
 * @param low Lower limit of q-power range  
 * @param component Matrix dimension index to operate on -- controls the variable in which the Pochhammer is expanded
 * @param components Total number of matrix dimensions -- is the number of variables in the polynomial
 * @param lengths Array of component lengths for matrix indexing
 * @param blocks Block structure for matrix-as-vector operations
 * 
 * Applies the q-Pochhammer symbol (x*q)_n = (1-x*q)(1-x*q^2)...(1-x*q^n)
 * to the polynomial term. Each factor contributes (1 - x_ij * q^z) terms
 * where x_ij are matrix elements in the specified component.
 */
void x_q_pochhammer (std::vector<bilvector<int>>& term, int up, int low, int component, int components, std::vector<int> lengths, std::vector<int> blocks) {
    for (int z = low; z <= up; z++) {
        for (int w = lengths[component]; w > 0; w--) {
            matrix_index_column(components, lengths, component, w - 1, term, 1, z, -1, blocks);
        }
    }
}

/**
 * @brief Apply inverse q-Pochhammer symbol 1/(x*q)_n to polynomial term
 * @param term Input/output polynomial term (modified in place)
 * @param up Upper limit of q-power range
 * @param low Lower limit of q-power range
 * @param component Matrix component index to operate on -- controls the variable in which the Pochhammer is expanded
 * @param components Total number of matrix components -- is the number of variables in the polynomial
 * @param lengths Array of component lengths for matrix indexing
 * @param blocks Block structure for matrix-as-vector operations
 * 
 * Applies the inverse q-Pochhammer symbol 1/(x*q)_n to the polynomial term.
 * This involves expanding 1/((1-x*q)(1-x*q^2)...(1-x*q^n)) using geometric
 * series, contributing multiple terms with different powers of matrix elements
 * and q-powers. The triple loop structure generates all necessary coefficient
 * contributions for the inverse Pochhammer expansion.
 */
void x_q_inv_pochhammer (std::vector<bilvector<int>>& term, int up, int low, int component, int components, std::vector<int> lengths, std::vector<int> blocks) {
    for (int z = low; z <= up; z++) {
        for (int w = lengths[component]; w > 0; w--) {
            for (int r = 1; r <= w; r++) {
                matrix_index_column(components, lengths, component, w - r, term, r, z, 1, blocks);
            }
        }
    }
}
