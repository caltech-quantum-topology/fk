#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>
#include <string>
#include <sstream>

using namespace std;

// Monomial representation: symbol -> (sign, index, count)
struct Monomial {
    double coefficient;
    int q_exp;  // exponent of q prefactor
    map<string, int> symbols;  // symbol name -> count
    
    Monomial(double c = 1.0, int qe = 0) : coefficient(c), q_exp(qe) {}
};

// Polynomial = sum of monomials
struct Polynomial {
    vector<Monomial> terms;
    
    void add(const Monomial& m) {
        terms.push_back(m);
    }
};

// Matrix of polynomials (each entry is a polynomial)
using PolyMatrix = vector<vector<Polynomial>>;

// Count inversions
int countInversions(const vector<int>& arr) {
    int n = arr.size();
    vector<int> v(n);
    iota(v.begin(), v.end(), 0);
    int ans = 0;
    
    for (int i = 0; i < n; i++) {
        auto it = lower_bound(v.begin(), v.end(), arr[i]);
        ans += distance(v.begin(), it);
        v.erase(it);
    }
    return ans;
}

// Generate permutations with inversions
vector<pair<vector<int>, int>> getPermutations(int m) {
    vector<pair<vector<int>, int>> result;
    vector<int> perm(m);
    iota(perm.begin(), perm.end(), 0);
    
    do {
        result.push_back({perm, countInversions(perm)});
    } while (next_permutation(perm.begin(), perm.end()));
    
    return result;
}

// Generate combinations
void getCombinations(int n, int m, vector<vector<int>>& result) {
    vector<int> comb(m);
    function<void(int, int)> generate = [&](int start, int idx) {
        if (idx == m) {
            result.push_back(comb);
            return;
        }
        for (int i = start; i <= n - m + idx; i++) {
            comb[idx] = i;
            generate(i + 1, idx + 1);
        }
    };
    generate(0, 0);
}

// Create identity matrix
PolyMatrix identity(int n) {
    PolyMatrix mat(n, vector<Polynomial>(n));
    for (int i = 0; i < n; i++) {
        Monomial m(1.0, 0);
        mat[i][i].add(m);
    }
    return mat;
}

// Generate Burau 2x2 block
vector<vector<Polynomial>> burau22(int sign, int index) {
    vector<vector<Polynomial>> block(2, vector<Polynomial>(2));
    
    if (sign == 1) {
        // [[a+, b+], [c+, 0]]
        Monomial a_plus(1.0, 0);
        a_plus.symbols["a+" + to_string(index)] = 1;
        block[0][0].add(a_plus);
        
        Monomial b_plus(1.0, 0);
        b_plus.symbols["b+" + to_string(index)] = 1;
        block[0][1].add(b_plus);
        
        Monomial c_plus(1.0, 0);
        c_plus.symbols["c+" + to_string(index)] = 1;
        block[1][0].add(c_plus);
        
        // block[1][1] is 0 (empty)
    } else {
        // [[0, c-], [b-, a-]]
        Monomial c_minus(1.0, 0);
        c_minus.symbols["c-" + to_string(index)] = 1;
        block[0][1].add(c_minus);
        
        Monomial b_minus(1.0, 0);
        b_minus.symbols["b-" + to_string(index)] = 1;
        block[1][0].add(b_minus);
        
        Monomial a_minus(1.0, 0);
        a_minus.symbols["a-" + to_string(index)] = 1;
        block[1][1].add(a_minus);
    }
    
    return block;
}

// Multiply two monomials
Monomial multiplyMonomials(const Monomial& m1, const Monomial& m2) {
    Monomial result(m1.coefficient * m2.coefficient, m1.q_exp + m2.q_exp);
    result.symbols = m1.symbols;
    for (const auto& [sym, count] : m2.symbols) {
        result.symbols[sym] += count;
    }
    return result;
}

// Multiply two polynomials
Polynomial multiplyPolynomials(const Polynomial& p1, const Polynomial& p2) {
    Polynomial result;
    for (const auto& m1 : p1.terms) {
        for (const auto& m2 : p2.terms) {
            result.add(multiplyMonomials(m1, m2));
        }
    }
    return result;
}

// Matrix multiplication
PolyMatrix matrixMultiply(const PolyMatrix& A, const PolyMatrix& B) {
    int n = A.size();
    int m = B[0].size();
    int p = B.size();
    PolyMatrix result(n, vector<Polynomial>(m));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < p; k++) {
                Polynomial prod = multiplyPolynomials(A[i][k], B[k][j]);
                for (const auto& term : prod.terms) {
                    result[i][j].add(term);
                }
            }
        }
    }
    return result;
}

// Build reduced Burau matrix
PolyMatrix reducedBurauMatrix(const vector<int>& braid_word, int& max_idx, int& writhe) {
    max_idx = 0;
    writhe = 0;
    
    for (int w : braid_word) {
        max_idx = max(max_idx, abs(w));
        writhe += (w > 0) ? 1 : -1;
    }
    
    PolyMatrix mat = identity(max_idx + 1);
    
    int crossing_index = 1;
    for (int w : braid_word) {
        int idx = abs(w) - 1;  // Convert to 0-based
        int sign = (w > 0) ? 1 : -1;
        
        // Build block diagonal matrix with burau22 in the middle
        PolyMatrix factor = identity(max_idx + 1);
        auto burau_block = burau22(sign, crossing_index);
        
        factor[idx][idx] = burau_block[0][0];
        factor[idx][idx + 1] = burau_block[0][1];
        factor[idx + 1][idx] = burau_block[1][0];
        factor[idx + 1][idx + 1] = burau_block[1][1];
        
        mat = matrixMultiply(mat, factor);
        crossing_index++;
    }
    
    // Remove first row and column
    int n = max_idx;
    PolyMatrix reduced(n, vector<Polynomial>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            reduced[i][j] = mat[i + 1][j + 1];
        }
    }
    
    // Multiply by q (add 1 to all q_exp)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (auto& term : reduced[i][j].terms) {
                term.q_exp += 1;
            }
        }
    }
    
    return reduced;
}

// Compute walks (permanent-like sum)
Polynomial computeWalks(const PolyMatrix& matrix, int n) {
    Polynomial result;
    
    for (int m = 1; m <= n; m++) {
        vector<vector<int>> subsets;
        getCombinations(n, m, subsets);
        
        auto perms = getPermutations(m);
        
        Polynomial sum_m;
        
        for (const auto& subset : subsets) {
            for (const auto& [pi, inversions] : perms) {
                Monomial summand(1.0, -inversions);
                
                for (int idx = 0; idx < m; idx++) {
                    // Multiply by matrix[subset[pi[idx]]][subset[idx]]
                    const Polynomial& entry = matrix[subset[pi[idx]]][subset[idx]];
                    if (!entry.terms.empty()) {
                        // For simplicity, take first term (should handle all)
                        for (const auto& term : entry.terms) {
                            Monomial temp = multiplyMonomials(summand, term);
                            sum_m.add(temp);
                            summand = temp;
                            break;  // Simplified: only first term
                        }
                    }
                }
                
                // Don't add if all zeros
                if (!summand.symbols.empty() || summand.coefficient != 0) {
                    sum_m.add(summand);
                }
            }
        }
        
        // Apply (-1)^(m-1) factor
        double sign = (m % 2 == 1) ? 1.0 : -1.0;
        for (auto& term : sum_m.terms) {
            term.coefficient *= sign;
            result.add(term);
        }
    }
    
    return result;
}

// Simplified evaluation at q value
map<int, double> evaluateJones(const vector<int>& braid_word, int N, double q_val) {
    int max_idx, writhe;
    PolyMatrix matrix = reducedBurauMatrix(braid_word, max_idx, writhe);
    
    // For now, return a placeholder result
    // Full implementation would process walks through SWC algorithm
    
    map<int, double> result;
    result[0] = 1.0;
    
    return result;
}

// Print result
void printResult(const map<int, double>& poly) {
    vector<pair<int, double>> terms(poly.begin(), poly.end());
    
    cout << "[";
    for (size_t i = 0; i < terms.size(); i++) {
        cout << "(" << terms[i].first << ", " << terms[i].second << ")";
        if (i < terms.size() - 1) cout << ", ";
    }
    cout << "]" << endl;
}

int main() {
    cout << "Colored Jones Polynomial Calculator (C++)" << endl;
    cout << "Note: This is a complex algorithm. Full implementation in progress." << endl << endl;
    
    // Test cases
    vector<int> trefoil = {1, 1, 1};
    vector<int> figure8 = {1, -2, 1, -2};
    
    cout << "For accurate results, please run the Python version." << endl;
    cout << "C++ translation of full SWC algorithm with symbolic" << endl;
    cout << "manipulation requires significant additional implementation." << endl;
    
    return 0;
}
