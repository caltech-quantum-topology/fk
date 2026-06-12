#include "fk/qalg_links.hpp"
#include "fk/fmpz_wrapper.hpp"
#include <map>
#include <unordered_map>
#include <shared_mutex>
#include <mutex>
#include <functional>

namespace {
  // Hash functions for cache keys
  struct QBinomialKey {
    int upper, lower;
    bool operator==(const QBinomialKey& other) const {
      return upper == other.upper && lower == other.lower;
    }
  };

  struct QBinomialKeyHash {
    std::size_t operator()(const QBinomialKey& k) const {
      return std::hash<int>()(k.upper) ^ (std::hash<int>()(k.lower) << 1);
    }
  };

  struct PochhammerKey {
    int n, qpow;
    bool operator==(const PochhammerKey& other) const {
      return n == other.n && qpow == other.qpow;
    }
  };

  struct PochhammerKeyHash {
    std::size_t operator()(const PochhammerKey& k) const {
      return std::hash<int>()(k.n) ^ (std::hash<int>()(k.qpow) << 1);
    }
  };

  struct InversePochhammerKey {
    int n, qpow, xMax;
    bool operator==(const InversePochhammerKey& other) const {
      return n == other.n && qpow == other.qpow && xMax == other.xMax;
    }
  };

  struct InversePochhammerKeyHash {
    std::size_t operator()(const InversePochhammerKey& k) const {
      return std::hash<int>()(k.n) ^ (std::hash<int>()(k.qpow) << 1) ^ (std::hash<int>()(k.xMax) << 2);
    }
  };

  // Thread-safe caches with shared_mutex for read-heavy access
  std::unordered_map<QBinomialKey, QPolynomialType, QBinomialKeyHash> qbinomial_positive_cache;
  std::shared_mutex qbinomial_positive_mutex;

  std::unordered_map<QBinomialKey, QPolynomialType, QBinomialKeyHash> qbinomial_negative_cache;
  std::shared_mutex qbinomial_negative_mutex;

  std::unordered_map<PochhammerKey, PolynomialType, PochhammerKeyHash> pochhammer_cache;
  std::shared_mutex pochhammer_mutex;

  std::unordered_map<InversePochhammerKey, PolynomialType, InversePochhammerKeyHash> inverse_pochhammer_cache;
  std::shared_mutex inverse_pochhammer_mutex;
}

QPolynomialType QBinomialPositive(int upperLimit, int lowerLimit) {
  QBinomialKey key{upperLimit, lowerLimit};

  // Try to read from cache first
  {
    std::shared_lock<std::shared_mutex> lock(qbinomial_positive_mutex);
    auto it = qbinomial_positive_cache.find(key);
    if (it != qbinomial_positive_cache.end()) {
      return it->second;
    }
  }

  int n = upperLimit;
  int k = lowerLimit;

  // Zero polynomial if out of range
  if (k < 0 || n < 0 || k > n) {
    QPolynomialType result; // zero polynomial

    // Cache the result
    {
      std::unique_lock<std::shared_mutex> lock(qbinomial_positive_mutex);
      qbinomial_positive_cache.emplace(key, result);
    }
    return result;
  }

  auto makeConstOne = []() {
    QPolynomialType p;
    p.setCoefficient(0, 1); // p = 1
    return p;
  };

  // C[i][j] = [i choose j]_q
  std::vector<std::vector<QPolynomialType>> C(
      n + 1, std::vector<QPolynomialType>(k + 1, QPolynomialType()));

  C[0][0] = makeConstOne();

  for (int i = 1; i <= n; ++i) {
    C[i][0] = makeConstOne(); // [i choose 0]_q = 1

    int jMax = std::min(i, k);
    for (int j = 1; j <= jMax; ++j) {
      if (j == i) {
        C[i][j] = makeConstOne(); // [i choose i]_q = 1
      } else {
        // [i choose j]_q = [i-1 choose j]_q + q^(i-j) [i-1 choose j-1]_q
        QPolynomialType shifted = multiplyByQPower(C[i - 1][j - 1], i - j);
        C[i][j] = C[i - 1][j] + shifted;
      }
    }
  }

  QPolynomialType result = C[n][k];

  // Cache the computed result
  {
    std::unique_lock<std::shared_mutex> lock(qbinomial_positive_mutex);
    qbinomial_positive_cache.emplace(key, result);
  }

  return result;
}

QPolynomialType QBinomialNegative(int upperLimit, int lowerLimit) {
  QBinomialKey key{upperLimit, lowerLimit};

  // Try to read from cache first
  {
    std::shared_lock<std::shared_mutex> lock(qbinomial_negative_mutex);
    auto it = qbinomial_negative_cache.find(key);
    if (it != qbinomial_negative_cache.end()) {
      return it->second;
    }
  }

  int k = lowerLimit;
  int u = upperLimit;

  // Handle nonsense k
  if (k < 0) {
    QPolynomialType result; // zero polynomial

    // Cache the result
    {
      std::unique_lock<std::shared_mutex> lock(qbinomial_negative_mutex);
      qbinomial_negative_cache.emplace(key, result);
    }
    return result;
  }

  // If upperLimit >= 0, just reuse the positive version
  if (u >= 0) {
    return QBinomialPositive(u, k);
  }

  // u < 0: write u = -n with n > 0
  int n = -u;

  // Base positive q-binomial: [n + k - 1 choose k]_q
  QPolynomialType base = QBinomialPositive(n + k - 1, k);

  // Shift exponent so that we match computeNegativeQBinomial.
  // Desired shift: u * k - k*(k-1)/2
  // (this gives exponents [- (1+u)k - k(k+1)/2, -k(k+1)/2],
  //  which is exactly what your old code produces)
  int shift = u * k - (k * (k - 1)) / 2;

  QPolynomialType result = multiplyByQPower(base, shift);

  // Apply (-1)^k factor
  if (k % 2 != 0) { // k odd → multiply by -1
    int minExp = result.getMaxNegativeIndex();
    int maxExp = result.getMaxPositiveIndex();
    for (int e = minExp; e <= maxExp; ++e) {
      int coeff = result.getCoefficient(e);
      result.setCoefficient(e, -coeff);
    }
  }

  // Cache the computed result
  {
    std::unique_lock<std::shared_mutex> lock(qbinomial_negative_mutex);
    qbinomial_negative_cache.emplace(key, result);
  }

  return result;
}

QPolynomialType QBinomial(int upperLimit, int lowerLimit) {
  return (upperLimit > 0) ? QBinomialPositive(upperLimit, lowerLimit)
                          : QBinomialNegative(upperLimit, lowerLimit);
}

// Compute (x q; q)_n as a MultivariablePolynomial in one x-variable
// P(q, x) = ∏_{k=1}^n (1 - x q^{qpow + k})
// Optimized using direct coefficient computation
PolynomialType qpochhammer_xq_q(int n, int qpow) {
  PochhammerKey key{n, qpow};

  // Try to read from cache first
  {
    std::shared_lock<std::shared_mutex> lock(pochhammer_mutex);
    auto it = pochhammer_cache.find(key);
    if (it != pochhammer_cache.end()) {
      return it->second;
    }
  }

  const int numXVars = 1;
  const int maxXDegree = n;

  // Coefficients map: coeffs[x_degree][q_power] = coefficient
  std::map<int, std::map<int, fmpz_wrapper>> coeffs;
  coeffs[0][0] = 1; // Initialize with 1

  // For each factor (1 - x q^{qpow + k})
  for (int k = 0; k < n; ++k) {
    const int q_factor = qpow + k;
    std::map<int, std::map<int, fmpz_wrapper>> new_coeffs;

    // Multiply current polynomial by (1 - x q^q_factor)
    for (const auto &[x_deg, q_map] : coeffs) {
      for (const auto &[q_pow, coeff] : q_map) {
        if (!coeff.is_zero()) {
          // "1" term
          new_coeffs[x_deg][q_pow].add(coeff);

          // "-x q^q_factor" term
          if (x_deg + 1 <= maxXDegree) {
            new_coeffs[x_deg + 1][q_pow + q_factor].sub(coeff);
          }
        }
      }
    }
    coeffs = std::move(new_coeffs);
  }

  // Build result polynomial
  PolynomialType result(numXVars, maxXDegree);
  for (const auto &[x_deg, q_map] : coeffs) {
    if (x_deg <= maxXDegree) {
      for (const auto &[q_pow, coeff] : q_map) {
        if (!coeff.is_zero()) {
          result.addToCoefficientFmpz(q_pow, {x_deg}, coeff.val);
        }
      }
    }
  }

  // Cache the computed result
  {
    std::unique_lock<std::shared_mutex> lock(pochhammer_mutex);
    pochhammer_cache.emplace(key, result);
  }

  return result;
}

// Compute 1/(x q^qpow; q)_n as a MultivariablePolynomial in one x-variable
// P(q, x) = ∏_{l=0}^{n-1} ∑_{m=0}^{xMax+1} x^m q^{(l+qpow)*m}
// Optimized using direct coefficient computation


PolynomialType inverse_qpochhammer_xq_q(int n, int qpow, int xMax) {
  InversePochhammerKey key{n, qpow, xMax};

  // Try to read from cache first
  {
    std::shared_lock<std::shared_mutex> lock(inverse_pochhammer_mutex);
    auto it = inverse_pochhammer_cache.find(key);
    if (it != inverse_pochhammer_cache.end()) {
      return it->second;
    }
  }

  const int numXVars = 1;

  // Coefficients map: coeffs[x_degree][q_power] = coefficient
  std::map<int, std::map<int, fmpz_wrapper>> coeffs;
  coeffs[0][0] = 1; // Initialize with 1

  // For each factor (geometric series)
  for (int l = 0; l < n; ++l) {
    const int q_base = l + qpow;
    std::map<int, std::map<int, fmpz_wrapper>> new_coeffs;

    // Multiply current polynomial by the l-th geometric series
    for (const auto &[x_deg, q_map] : coeffs) {
      for (const auto &[q_pow, coeff] : q_map) {
        if (!coeff.is_zero()) {
          // Add terms from geometric series: 1 + x*q^q_base + x^2*q^(2*q_base) + ...
          for (int m = 0; x_deg + m <= xMax; ++m) {
            const int new_x_deg = x_deg + m;
            const int new_q_pow = q_pow + m * q_base;
            new_coeffs[new_x_deg][new_q_pow].add(coeff);
          }
        }
      }
    }
    coeffs = std::move(new_coeffs);
  }

  // Build result polynomial
  PolynomialType result(numXVars, 0);
  for (const auto &[x_deg, q_map] : coeffs) {
    if (x_deg <= xMax) {
      for (const auto &[q_pow, coeff] : q_map) {
        if (!coeff.is_zero()) {
          result.addToCoefficientFmpz(q_pow, {x_deg}, coeff.val);
        }
      }
    }
  }

  // Cache the computed result
  {
    std::unique_lock<std::shared_mutex> lock(inverse_pochhammer_mutex);
    inverse_pochhammer_cache.emplace(key, result);
  }

  return result;
}
