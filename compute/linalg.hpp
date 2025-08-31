/**
 * @file linalg.hpp
 * @brief Linear algebra operations for matrix-polynomial computations in Gukov-Manolescu invariants
 * 
 * This file implements linear algebra utilities for the Gukov-Manolescu invariant computation,
 * including matrix-vector operations, dot products, and specialized multi-dimensional indexing
 * operations for polynomial terms stored as bilvector arrays.
 * 
 * Key mathematical utilities:
 * - Change-of-basis transformations between angle and segment variables
 * - Matrix-vector products for integer and double vectors  
 * - Multi-dimensional matrix indexing with block structure
 * - Offset addition operations for polynomial term manipulation
 * 
 * Coordinate system transformations:
 * - Original inequality: x @ p <= xdeg (in angle variables)
 * - Transformed inequality: (M^-1)^T p @ M(angle x-criterion) <= xdeg (in segment variables)
 * - Relationship: p = M^T (M^-1)^T p (between integer hull points)
 * Note: M^-1 is never explicitly computed, only M^T and (M^-1)^T are used
 * 
 * When both ILPs are in standard form, the nontrivial inequalities of one side
 * become the trivial identity matrix expressing variable signs in the other side.
 *
 * TODO: matrix_index_column is equivalent to offset_addition with bilvector_offset = r * z,
 * but r and z are specified separately in matrix_index_column and can be combined into 
 * a single bilvector offset for efficiency.
 */

#pragma once

#include <vector>
#include "bilvector.hpp"

/**
 * @brief Matrix-vector multiplication for integer matrices and vectors
 * @param matrix Square integer matrix (reference for efficiency)
 * @param vector Integer vector to multiply (reference for efficiency)
 * @return Result of matrix * vector as integer vector
 * 
 * Performs standard matrix-vector multiplication: out[i] = sum_j(matrix[i][j] * vector[j])
 * Used for coordinate transformations in the ILP solver.
 */
std::vector<int> mult(std::vector<std::vector<int>>& matrix, std::vector<int>& vector) {
    std::vector<int> out(vector.size());
    for (int i = 0; i < vector.size(); i++) {
        int acc = 0; 
        for (int j = 0; j < vector.size(); j++) {
            acc += matrix[i][j] * vector[j];
        }
        out[i] = acc;
    }
    return out;
}

/**
 * @brief Matrix-vector multiplication with matrix transpose and 1-based indexing offset
 * @param matrix Integer matrix (reference for efficiency)  
 * @param vector Integer vector to multiply (const reference)
 * @return Result of matrix^T * vector with 1-based indexing offset
 * 
 * Performs transposed matrix-vector multiplication with indexing offset:
 * out[i] = sum_j(matrix[1 + j][1 + i] * vector[j])
 * 
 * The 1-based indexing offset (matrix[1 + j][1 + i]) suggests this function
 * works with matrices that have a header row/column or use 1-based indexing
 * in the underlying data structure.
 */
std::vector<int> multT(std::vector<std::vector<int>>& matrix, const std::vector<int>& vector) {
    std::vector<int> out(vector.size());
    for (int i = 0; i < vector.size(); i++) {
        int acc = 0; 
        for (int j = 0; j < vector.size(); j++) {
            acc += matrix[1 + j][1 + i] * vector[j];
        }
        out[i] = acc;
    }
    return out;
}

/**
 * @brief Matrix-vector multiplication for integer matrix and double vector
 * @param matrix Integer matrix (reference for efficiency)
 * @param vector Double vector to multiply (reference for efficiency)  
 * @return Result of matrix * vector as double vector
 * 
 * Performs matrix-vector multiplication with mixed integer/double precision:
 * out[i] = sum_j(matrix[i][j] * vector[j])
 */
std::vector<double> mult(std::vector<std::vector<int>>& matrix, std::vector<double>& vector) {
    std::vector<double> out(vector.size());
    for (int i = 0; i < vector.size(); i++) {
        double acc = 0.0; 
        for (int j = 0; j < vector.size(); j++) {
            acc += matrix[i][j] * vector[j];
        }
        out[i] = acc;
    }
    return out;
}

/**
 * @brief Dot product with offset indexing for constraint evaluation
 * @param a Integer vector with offset element (includes constant term at index 0)
 * @param b Integer vector for dot product computation
 * @return Dot product result: a[0] + sum_z(a[z+1] * b[z])
 * 
 * Computes dot product with special indexing where a[0] acts as a constant term
 * and the actual dot product starts from a[1]. This pattern is common in
 * constraint evaluation where a[0] represents the constant bound and a[1:]
 * represents the constraint coefficients.
 */
int dot(const std::vector<int>& a, const std::vector<int>& b) {
    int acc = a[0];  // Start with constant term
    for (int z = 0; z < b.size(); z++) {
        acc += a[z + 1] * b[z];  // Dot product with offset indexing
    }
    return acc;
}

/**
 * @brief Recursive helper for multi-dimensional matrix indexing and polynomial coefficient manipulation
 * @param dim Number of matrix dimensions
 * @param lengths Array of dimension lengths for each matrix component
 * @param slice Target dimension to apply fixed value (matrix column selection)
 * @param value Fixed value to use in the slice dimension  
 * @param acc Accumulated linear index for current position
 * @param index Current dimension being processed
 * @param term Array of polynomial terms (bilvectors) to modify
 * @param r Coefficient multiplier for the operation
 * @param z Q-power offset for polynomial degrees
 * @param sign Sign multiplier (+1 or -1) for the coefficient contribution
 * @param blocks Block sizes for converting multi-dimensional indices to linear indices
 * 
 * This recursive function traverses a multi-dimensional matrix structure stored as
 * a linear array, applying polynomial operations at specific matrix positions.
 * When slice == index, it fixes that dimension to 'value'; otherwise it iterates
 * through all possible values in that dimension.
 */
void matrix_index_column_recurse(int& dim, std::vector<int> lengths, int& slice, int value, int acc, int index, std::vector<bilvector<int>>& term, int r, int& z, int sign, std::vector<int> blocks) {
    index++;
    int old_acc = acc;
    
    if (slice == index) {
        // Fixed value case: use specified value in this dimension
        acc = old_acc + value * blocks[index];
        if (dim > index + 1) {
            matrix_index_column_recurse(dim, lengths, slice, value, acc, index, term, r, z, sign, blocks);
        }
        else {
            // Terminal case: apply polynomial coefficient transformation
            int d = acc + r * blocks[slice];
            for (int k = term[acc].get_max_nindex(); k <= term[acc].get_max_pindex(); k++) {
                term[d][k + r * z] += sign * term[acc][k];
            }
        }
    }
    else {
        // Iterate through all values in this dimension
        for (int i = 0; i < lengths[index] + 1; i++) {
            acc = old_acc + i * blocks[index];
            if (dim > index + 1) {
                matrix_index_column_recurse(dim, lengths, slice, value, acc, index, term, r, z, sign, blocks);
            }
            else {
                // Terminal case: apply polynomial coefficient transformation
                int d = acc + r * blocks[slice];
                for (int k = term[acc].get_max_nindex(); k <= term[acc].get_max_pindex(); k++) {
                    term[d][k + r * z] += sign * term[acc][k];
                }
            }
        }
    }
}

/**
 * @brief Apply matrix column operations to polynomial terms with multi-dimensional indexing
 * @param dim Number of matrix dimensions
 * @param lengths Array of dimension lengths for each matrix component
 * @param slice Target dimension for column operation (which matrix variable to modify)
 * @param value Fixed value to apply in the slice dimension
 * @param term Array of polynomial terms (bilvectors) to modify in place
 * @param r Coefficient multiplier for the polynomial transformation
 * @param z Q-power offset for polynomial degree shifts
 * @param sign Sign multiplier (+1 or -1) for coefficient contributions
 * @param blocks Block sizes for multi-dimensional to linear index conversion
 * 
 * This function implements matrix column operations on polynomial terms stored as
 * multi-dimensional arrays. It applies transformations of the form:
 * term[d][k + r*z] += sign * term[source][k]
 * where d and source indices are computed using multi-dimensional indexing.
 * 
 * Used in quantum algebraic operations for Pochhammer symbol expansions.
 */
void matrix_index_column(int& dim, std::vector<int> lengths, int& slice, int value, std::vector<bilvector<int>>& term, int r, int& z, int sign, std::vector<int> blocks) {
    if (slice == 0) {
        // Special case: operating on first dimension
        if (dim > 1) {
            matrix_index_column_recurse(dim, lengths, slice, value, value, 0, term, r, z, sign, blocks);
        }
        else {
            // 1D case: direct polynomial coefficient transformation
            int d = value + r;
            for (int k = term[value].get_max_nindex(); k <= term[value].get_max_pindex(); k++) {
                term[d][k + r * z] += sign * term[value][k];
            }
        }
    }
    else {
        // General case: iterate through first dimension, recurse for others
        for (int i = 0; i < lengths[0] + 1; i++) {
            if (dim > 1) {
                matrix_index_column_recurse(dim, lengths, slice, value, i, 0, term, r, z, sign, blocks);
            }
            else {
                // 1D case: direct polynomial coefficient transformation
                int d = i + r;
                for (int k = term[i].get_max_nindex(); k <= term[i].get_max_pindex(); k++) {
                    term[d][k + r * z] += sign * term[i][k];
                }
            }
        }
    }
}

/**
 * @brief Recursive helper for offset addition between multi-dimensional polynomial arrays
 * @param a Destination polynomial array (modified in place)
 * @param b Source polynomial array to add from
 * @param offsets Array of index offsets for each dimension
 * @param bilvector_offset Q-power offset for polynomial degree adjustment
 * @param dim Number of dimensions in the arrays
 * @param lengths Array of dimension sizes
 * @param index Current dimension being processed
 * @param acc Linear index accumulator for array b
 * @param acc2 Linear index accumulator for array a (with offset corrections)
 * @param sign Sign multiplier (+1 or -1) for the addition
 * @param a_blocks Block sizes for array a indexing
 * @param b_blocks Block sizes for array b indexing
 * 
 * Recursively traverses multi-dimensional polynomial arrays, adding coefficients
 * from array b to array a with specified offsets in each dimension. The offset
 * mechanism allows shifting polynomial terms between different matrix positions.
 */
void offset_addition_recurse (std::vector<bilvector<int>>& a, std::vector<bilvector<int>>& b, std::vector<int>& offsets, int bilvector_offset, int& dim, std::vector<int> lengths, int index, int acc, int acc2, int sign, std::vector<int> a_blocks, std::vector<int> b_blocks) {
    index++;
    int old_acc = acc;
    int old_acc2 = acc2 + offsets[index] * a_blocks[index];
    
    // Iterate through valid range considering offset constraints
    for (int i = std::max(0, -offsets[index]); i < lengths[index] + 1; i++) {
        acc = old_acc + i * b_blocks[index];
        acc2 = old_acc2 + i * a_blocks[index];
        if (dim > index + 1) {
            offset_addition_recurse (a, b, offsets, bilvector_offset, dim, lengths, index, acc, acc2, sign, a_blocks, b_blocks);
        }
        else {
            // Terminal case: add polynomial coefficients with offset
            for (int q = b[acc].get_max_nindex(); q <= b[acc].get_max_pindex(); q++) {
                a[acc2][q + bilvector_offset] += sign * b[acc][q];
            }
        } 
    }
}

/**
 * @brief Add multi-dimensional polynomial arrays with index and degree offsets
 * @param a Destination polynomial array (modified in place)
 * @param b Source polynomial array to add from (copied by value)
 * @param offsets Array of index offsets for each dimension
 * @param bilvector_offset Q-power offset for polynomial degree adjustment
 * @param dim Number of dimensions in the arrays
 * @param lengths Array of dimension sizes  
 * @param sign Sign multiplier (+1 or -1) for the addition
 * @param a_blocks Block sizes for array a multi-dimensional indexing
 * @param b_blocks Block sizes for array b multi-dimensional indexing
 * 
 * Performs offset addition: a[i + offsets] += sign * b[i] with q-power shift.
 * This operation is fundamental for polynomial manipulations in quantum algebraic
 * computations, allowing terms to be moved between different matrix positions
 * with appropriate degree adjustments.
 * 
 * The operation respects boundary constraints: only processes indices where
 * both source and destination positions are valid.
 */
void offset_addition(std::vector<bilvector<int>>& a, std::vector<bilvector<int>> b, std::vector<int>& offsets, int bilvector_offset, int& dim, std::vector<int> lengths, int sign, std::vector<int> a_blocks, std::vector<int> b_blocks) {
    // Iterate through valid range in first dimension considering offset constraints
    for (int i = std::max(0, -offsets[0]); i < lengths[0] + 1; i++) {
        if (dim > 1) {
            offset_addition_recurse (a, b, offsets, bilvector_offset, dim, lengths, 0, i, i + offsets[0], sign, a_blocks, b_blocks);
        }
        else {
            // 1D case: direct coefficient addition with offset
            for (int q = b[i].get_max_nindex(); q <= b[i].get_max_pindex(); q++) {
                a[i + offsets[0]][q + bilvector_offset] += sign * b[i][q];
            }
        }
    } 
}
