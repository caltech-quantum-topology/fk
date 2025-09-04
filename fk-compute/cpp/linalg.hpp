#pragma once

#include <vector>
#include "bilvector.hpp"

// x @ p <= xdeg <------- original inequality in angle variables; use a change-of-basis matrix M to examine the segment side
// (M^-1)^T p @ M(angle x-criterion)  <= xdeg <------- same inequality in segment variables; solve for (M^-1)^T using the transformed criterion
// p = M^T (M^-1)^T p <------- relation between respective points of the integer hulls; interestingly, M^-1 isn't ever invoked explicitly

// when both ILP's are put in standard form, the nontrivial inequalities of a given side just becomes the trivial identity matrix expressing the sign of the variables in the other side

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

// std::vector<int> mult(std::vector<std::vector<int>>& matrix, const std::vector<int>& vector) {
//     std::vector<int> out(vector.size());
//     for (int i = 0; i < vector.size(); i++) {
//         int acc = 0; 
//         for (int j = 0; j < vector.size(); j++) {
//             acc += matrix[1 + i][1 + j] * vector[j];
//         }
//         out[i] = acc;
//     }
//     return out;
// }

// with transpose
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

std::vector<double> mult(std::vector<std::vector<int>>& matrix, std::vector<double>& vector) {
    std::vector<double> out(vector.size());
    for (int i = 0; i < vector.size(); i++) {
        int acc = 0; 
        for (int j = 0; j < vector.size(); j++) {
            acc += matrix[i][j] * vector[j];
        }
        out[i] = acc;
    }
    return out;
}

int dot(const std::vector<int>& a, const std::vector<int>& b) {
    int acc = a[0];
    for (int z = 0; z < b.size(); z++) {
        acc += a[z + 1] * b[z];
    }
    return acc;
}

void matrix_index_column_recurse(int& dim, std::vector<int> lengths, int& slice, int value, int acc, int index, std::vector<bilvector<int>>& term, int r, int& z, int sign, std::vector<int> blocks) {
    index++;
    int old_acc = acc;
    if (slice == index) {
        acc = old_acc + value * blocks[index];
        if (dim > index + 1) {
            matrix_index_column_recurse(dim, lengths, slice, value, acc, index, term, r, z, sign, blocks);
        }
        else {
            int d = acc + r * blocks[slice];
            for (int k = term[acc].get_max_nindex(); k <= term[acc].get_max_pindex(); k++) {
                term[d][k + r * z] += sign * term[acc][k];
            }
        }
    }
    else {
        for (int i = 0; i < lengths[index] + 1; i++) {
            acc = old_acc + i * blocks[index];
            if (dim > index + 1) {
                matrix_index_column_recurse(dim, lengths, slice, value, acc, index, term, r, z, sign, blocks);
            }
            else {
                int d = acc + r * blocks[slice];
                for (int k = term[acc].get_max_nindex(); k <= term[acc].get_max_pindex(); k++) {
                    term[d][k + r * z] += sign * term[acc][k];
                }
            }
        }
    }
}

void matrix_index_column(int& dim, std::vector<int> lengths, int& slice, int value, std::vector<bilvector<int>>& term, int r, int& z, int sign, std::vector<int> blocks) {
    if (slice == 0) {
        if (dim > 1) {
            matrix_index_column_recurse(dim, lengths, slice, value, value, 0, term, r, z, sign, blocks);
        }
        else {
            int d = value + r;
            for (int k = term[value].get_max_nindex(); k <= term[value].get_max_pindex(); k++) {
                term[d][k + r * z] += sign * term[value][k];
            }
        }
    }
    else {
        for (int i = 0; i < lengths[0] + 1; i++) {
            if (dim > 1) {
                matrix_index_column_recurse(dim, lengths, slice, value, i, 0, term, r, z, sign, blocks);
            }
            else {
                int d = i + r;
                for (int k = term[i].get_max_nindex(); k <= term[i].get_max_pindex(); k++) {
                    term[d][k + r * z] += sign * term[i][k];
                }
            }
        }
    }
}

void offset_addition_recurse (std::vector<bilvector<int>>& a, std::vector<bilvector<int>>& b, std::vector<int>& offsets, int bilvector_offset, int& dim, std::vector<int> lengths, int index, int acc, int acc2, int sign, std::vector<int> a_blocks, std::vector<int> b_blocks) {
    index++;
    int old_acc = acc;
    int old_acc2 = acc2 + offsets[index] * a_blocks[index];
    for (int i = std::max(0, -offsets[index]); i < lengths[index] + 1; i++) {
        acc = old_acc + i * b_blocks[index];
        acc2 = old_acc2 + i * a_blocks[index];
        if (dim > index + 1) {
            offset_addition_recurse (a, b, offsets, bilvector_offset, dim, lengths, index, acc, acc2, sign, a_blocks, b_blocks);
        }
        else {
            for (int q = b[acc].get_max_nindex(); q <= b[acc].get_max_pindex(); q++) {
                a[acc2][q + bilvector_offset] += sign * b[acc][q];
            }
        } 
    }
}

void offset_addition(std::vector<bilvector<int>>& a, std::vector<bilvector<int>> b, std::vector<int>& offsets, int bilvector_offset, int& dim, std::vector<int> lengths, int sign, std::vector<int> a_blocks, std::vector<int> b_blocks) {
    for (int i = std::max(0, -offsets[0]); i < lengths[0] + 1; i++) {
        if (dim > 1) {
            offset_addition_recurse (a, b, offsets, bilvector_offset, dim, lengths, 0, i, i + offsets[0], sign, a_blocks, b_blocks);
        }
        else {
            for (int q = b[i].get_max_nindex(); q <= b[i].get_max_pindex(); q++) {
                a[i + offsets[0]][q + bilvector_offset] += sign * b[i][q];
            }
        }
    } 
}