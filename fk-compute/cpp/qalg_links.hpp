#pragma once

#include <vector>
#include <functional>

#include "bilvector.hpp"
#include "linalg.hpp"

// consider saving the multiplicands as separate variables before multiplying by term

void pp_q_binom_ (std::vector<int>& binom, int i, int j, int shift) {
    if (i == j) {
        binom[shift] += 1;
    }
    else if (j == 0) {
        binom[shift] += 1;
    }
    else {
        pp_q_binom_(binom, i - 1, j, shift + j);
        pp_q_binom_(binom, i - 1, j - 1, shift);
    }
}

void pp_q_binom (std::vector<bilvector<int>>& term, int i, int j, bool neg) {
    int Q_MAX_DEGREE = j * (i - j);
    std::vector<int> binom(Q_MAX_DEGREE + 1, 0);
    if (i == j) {
        binom[0] = 1;
    }
    else if (j == 0) {
        binom[0] = 1;
    }
    else {
        pp_q_binom_(binom, i - 1, j, j);
        pp_q_binom_(binom, i - 1, j - 1, 0);
    }
    if (neg) {
        bilvector<int> term_(term[0].nsize(), term[0].psize(), term[0].get_component_size(), 0);
        for (int j = term[0].get_max_nindex(); j <= term[0].get_max_pindex(); j++) {
            term_[j] = term[0][j];
            term[0][j] = 0;
        }
        for (int j = term_.get_max_nindex(); j <= term_.get_max_pindex(); j++) {
            if (term_[j] != 0) {
                for (int k = 0; k < Q_MAX_DEGREE + 1; k++) {
                    term[0][j - k] += binom[k] * term_[j];
                }
            }
        }
    }
    else {
        bilvector<int> term_(term[0].nsize(), term[0].psize(), term[0].get_component_size(), 0);
        for (int j = term[0].get_max_nindex(); j <= term[0].get_max_pindex(); j++) {
            term_[j] = term[0][j];
            term[0][j] = 0;
        }
        for (int j = term_.get_max_nindex(); j <= term_.get_max_pindex(); j++) {
            if (term_[j] != 0) {
                for (int k = 0; k < Q_MAX_DEGREE + 1; k++) {
                    term[0][j + k] += binom[k] * term_[j];
                }
            }
        }
    }
}

void np_q_binom_ (std::vector<int>& binom, int i, int j, int shift, bool neg) {
    if (j == 0) {
        if (neg) {
            binom[shift] += -1;
        }
        else {
            binom[shift] += 1;
        }
    }
    else if (i == -1) {
        np_q_binom_(binom, -1, j - 1, shift - j, !neg);
    }
    else {
        np_q_binom_(binom, i, j - 1, shift + 1 + i - j, !neg);
        np_q_binom_(binom, i + 1, j, shift, neg);
    }
}

void np_q_binom (std::vector<bilvector<int>>& term, int i, int j, bool neg) {
    int Q_DEGREE_DELTA = -(1 + i) * j;
    int Q_MAX_DEGREE = -j * (j + 1) / 2;
    std::vector<int> binom(Q_DEGREE_DELTA + 1, 0);
    if (j == 0) {
        binom[0] = 1;
    }
    else if (i == -1) {
        np_q_binom_(binom, -1, j - 1, Q_DEGREE_DELTA - Q_MAX_DEGREE - j, true);
    }
    else {
        np_q_binom_(binom, i, j - 1, Q_DEGREE_DELTA - Q_MAX_DEGREE + 1 + i - j, true);
        np_q_binom_(binom, i + 1, j, Q_DEGREE_DELTA - Q_MAX_DEGREE, false);
    }
    if (neg) {
        bilvector<int> term_(term[0].nsize(), term[0].psize(), term[0].get_component_size(), 0);
        for (int j = term[0].get_max_nindex(); j <= term[0].get_max_pindex(); j++) {
            term_[j] = term[0][j];
            term[0][j] = 0;
        }
        for (int j = term_.get_max_nindex(); j <= term_.get_max_pindex(); j++) {
            if (term_[j] != 0) {
                for (int k = 0; k < Q_DEGREE_DELTA + 1; k++) {
                    term[0][j - k + Q_DEGREE_DELTA - Q_MAX_DEGREE] += binom[k] * term_[j];
                }
            }
        }
    }
    else {
        bilvector<int> term_(term[0].nsize(), term[0].psize(), term[0].get_component_size(), 0);
        for (int j = term[0].get_max_nindex(); j <= term[0].get_max_pindex(); j++) {
            term_[j] = term[0][j];
            term[0][j] = 0;
        }
        for (int j = term_.get_max_nindex(); j <= term_.get_max_pindex(); j++) {
            if (term_[j] != 0) {
                for (int k = 0; k < Q_DEGREE_DELTA + 1; k++) {
                    term[0][j + k - Q_DEGREE_DELTA + Q_MAX_DEGREE] += binom[k] * term_[j];
                }
            }  
        }
    }
}

void x_q_pochhammer (std::vector<bilvector<int>>& term, int up, int low, int component, int components, std::vector<int> lengths, std::vector<int> blocks) {
    for (int z = low; z <= up; z++) {
        for (int w = lengths[component]; w > 0; w--) {
            matrix_index_column(components, lengths, component, w - 1, term, 1, z, -1, blocks);
        }
    }
}

void x_q_inv_pochhammer (std::vector<bilvector<int>>& term, int up, int low, int component, int components, std::vector<int> lengths, std::vector<int> blocks) {
    for (int z = low; z <= up; z++) {
        for (int w = lengths[component]; w > 0; w--) {
            for (int r = 1; r <= w; r++) {
                matrix_index_column(components, lengths, component, w - r, term, r, z, 1, blocks);
            }
        }
    }
}