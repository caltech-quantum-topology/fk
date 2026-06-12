// "bi list vectors"
// description:

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>

template <typename T> struct bilvector {
private:
  int componentSize;
  int negativeVectorCount = 0;
  int positiveVectorCount = 0;
  int maxNegativeIndex = 0;
  int maxPositiveIndex = 0;
  T defaultValue;
  std::vector<std::vector<T>> negativeVectors = {};
  std::vector<std::vector<T>> positiveVectors = {};

public:
  /**
   * Default constructor - creates a minimal bilvector
   * Compatible with QPolynomial's default constructor
   */
  bilvector() : bilvector(0, 1, 1, T{}) {}

  bilvector(int initialNegativeVectorCount, int initialPositiveVectorCount,
            int componentSizeParam, T defaultValueParam) {
    componentSize = componentSizeParam;
    defaultValue = defaultValueParam;
    negativeVectorCount = initialNegativeVectorCount;
    negativeVectors.resize(initialNegativeVectorCount, std::vector<T>(componentSize, defaultValue));
    positiveVectorCount = initialPositiveVectorCount;
    positiveVectors.resize(initialPositiveVectorCount, std::vector<T>(componentSize, defaultValue));
  }
  bilvector<T> invertExponents() const;
  void print(const std::string &varName = "q") const;
  int getNegativeSize() { return negativeVectorCount * componentSize; }
  int getNegativeSize() const { return negativeVectorCount * componentSize; }
  int getPositiveSize() { return positiveVectorCount * componentSize; }
  int getPositiveSize() const { return positiveVectorCount * componentSize; }
  T &operator[](int accessIndex) {
    if (componentSize <= 0) {
      throw std::runtime_error("bilvector componentSize is zero or negative");
    }
    if (accessIndex > maxPositiveIndex) {
      maxPositiveIndex = accessIndex;
    } else if (accessIndex < maxNegativeIndex) {
      maxNegativeIndex = accessIndex;
    }
    if (accessIndex >= 0) {
      if (accessIndex >= (*this).getPositiveSize()) {
        int x = (accessIndex - (*this).getPositiveSize()) / componentSize;
        positiveVectorCount += x + 1;
        positiveVectors.resize(positiveVectorCount, std::vector<T>(componentSize, defaultValue));
      }
      int vectorIndex = accessIndex / componentSize;
      return positiveVectors[vectorIndex][accessIndex - vectorIndex * componentSize];
    } else {
      accessIndex = -1 - accessIndex;
      if (accessIndex >= (*this).getNegativeSize()) {
        int x = (accessIndex - (*this).getNegativeSize()) / componentSize;
        negativeVectorCount += x + 1;
        negativeVectors.resize(negativeVectorCount, std::vector<T>(componentSize, defaultValue));
      }
      int vectorIndex = accessIndex / componentSize;
      return negativeVectors[vectorIndex][accessIndex - vectorIndex * componentSize];
    }
  }
  int getNegativeVectorCount() { return negativeVectorCount; }
  int getPositiveVectorCount() { return positiveVectorCount; }
  int getMaxNegativeIndex() { return maxNegativeIndex; }
  int getMaxNegativeIndex() const { return maxNegativeIndex; }
  int getMaxPositiveIndex() { return maxPositiveIndex; }
  int getMaxPositiveIndex() const { return maxPositiveIndex; }
  int getComponentSize() { return componentSize; }
  int getComponentSize() const { return componentSize; }
  int nTerms() const {
      int count = 0;

      // Negative side
      {
          int idx = maxNegativeIndex;
          for (const auto &block : negativeVectors) {
              for (int k = 0; k < componentSize && idx < 0; ++k, ++idx) {
                  if (block[k] != defaultValue) {
                      ++count;
                  }
              }
          }
      }

      // Positive side
      {
          int idx = 0;
          for (const auto &block : positiveVectors) {
              for (int k = 0; k < componentSize && idx <= maxPositiveIndex; ++k, ++idx) {
                  if (block[k] != defaultValue) {
                      ++count;
                  }
              }
          }
      }

      return count;
  }

  // Check if bilvector is zero (all coefficients are default value)
  bool isZero() const {
    // Check all accessible negative indices
    for (int i = maxNegativeIndex; i < 0; ++i) {
      if ((*this)[i] != defaultValue) {
        return false;
      }
    }
    // Check all accessible positive indices
    for (int i = 0; i <= maxPositiveIndex; ++i) {
      if ((*this)[i] != defaultValue) {
        return false;
      }
    }
    return true;
  }

  // Const version of operator[]
  const T &operator[](int accessIndex) const {
    if (componentSize <= 0) {
      throw std::runtime_error("bilvector componentSize is zero or negative");
    }
    // For const access, we can't modify the structure, so we need to handle
    // this carefully
    if (accessIndex >= 0) {
      if (accessIndex >= this->getPositiveSize()) {
        // Return default value for out-of-bounds const access
        static T defaultVal = defaultValue;
        return defaultVal;
      }
      int vectorIndex = accessIndex / componentSize;
      return positiveVectors[vectorIndex][accessIndex - vectorIndex * componentSize];
    } else {
      accessIndex = -1 - accessIndex;
      if (accessIndex >= this->getNegativeSize()) {
        // Return default value for out-of-bounds const access
        static T defaultVal = defaultValue;
        return defaultVal;
      }
      int vectorIndex = accessIndex / componentSize;
      return negativeVectors[vectorIndex][accessIndex - vectorIndex * componentSize];
    }
  }

  /**
   * QPolynomial-compatible interface methods
   * These methods provide the same interface as QPolynomial for generic programming
   */

  /**
   * Get coefficient at a specific power (const)
   * Compatible with QPolynomial::getCoefficient()
   */
  T getCoefficient(int power) const {
    return (*this)[power];
  }

  /**
   * Set coefficient at a specific power
   * Compatible with QPolynomial::setCoefficient()
   */
  void setCoefficient(int power, T value) {
    (*this)[power] = value;
  }

  /**
   * Add to coefficient at a specific power
   * Compatible with QPolynomial::addToCoefficient()
   */
  void addToCoefficient(int power, T value) {
    (*this)[power] += value;
  }

  /**
   * Clear all coefficients (set to default value)
   * Compatible with QPolynomial::clear()
   */
  void clear() {
    for (auto& vec : negativeVectors) {
      std::fill(vec.begin(), vec.end(), defaultValue);
    }
    for (auto& vec : positiveVectors) {
      std::fill(vec.begin(), vec.end(), defaultValue);
    }
  }
};

template <typename T>
bilvector<T> makeLaurentPolynomial(int minExponent, int maxExponent,
                                   T defaultValue = T{}) {
  int componentSize = 1;

  int negativeCount = 0;
  if (minExponent < 0) {
    // exponents covered: -negativeCount, ..., -1
    negativeCount = -minExponent;
  }

  int positiveCount = 0;
  if (maxExponent >= 0) {
    // exponents covered: 0, 1, ..., positiveCount - 1
    positiveCount = maxExponent + 1;
  }

  return bilvector<T>(negativeCount, positiveCount, componentSize,
                      defaultValue);
}

template <typename T>
bilvector<T> makeZeroLike(const bilvector<T> &proto, int minExponent,
                          int maxExponent) {
  int componentSize = proto.getComponentSize();

  int negativeCount = 0;
  if (minExponent < 0) {
    // exponents covered on negative side: -negativeCount, ..., -1
    negativeCount = -minExponent;
  }

  int positiveCount = 0;
  if (maxExponent >= 0) {
    // exponents covered on positive side: 0, 1, ..., positiveCount - 1
    positiveCount = maxExponent + 1;
  }

  // We use T{} as the default value (typically 0 for arithmetic types)
  return bilvector<T>(negativeCount, positiveCount, componentSize, T{});
}

// ========================= Addition of bilvectors ==========================

template <typename T>
bilvector<T> operator+(const bilvector<T> &lhs, const bilvector<T> &rhs) {
  // We assume same component size; if you want, you can add a runtime check.
  // int csL = lhs.getComponentSize();
  // int csR = rhs.getComponentSize();
  // if (csL != csR) { throw std::runtime_error("componentSize mismatch"); }

  int minExp = std::min(lhs.getMaxNegativeIndex(), rhs.getMaxNegativeIndex());
  int maxExp = std::max(lhs.getMaxPositiveIndex(), rhs.getMaxPositiveIndex());

  bilvector<T> result = makeZeroLike(lhs, minExp, maxExp);

  for (int e = minExp; e <= maxExp; ++e) {
    result[e] = lhs[e] + rhs[e];
  }

  return result;
}

template <typename T>
bilvector<T> &operator+=(bilvector<T> &lhs, const bilvector<T> &rhs) {
  lhs = lhs + rhs;
  return lhs;
}

// ====================== Multiplication of bilvectors =======================

template <typename T>
bilvector<T> operator*(const bilvector<T> &lhs, const bilvector<T> &rhs) {
  // Again, we assume same component size; optional check as above.

  // Exponent ranges
  int lhsMin = lhs.getMaxNegativeIndex();
  int lhsMax = lhs.getMaxPositiveIndex();
  int rhsMin = rhs.getMaxNegativeIndex();
  int rhsMax = rhs.getMaxPositiveIndex();

  // Result exponent range: all sums i + j
  int resMin = lhsMin + rhsMin;
  int resMax = lhsMax + rhsMax;

  bilvector<T> result = makeZeroLike(lhs, resMin, resMax);

  for (int i = lhsMin; i <= lhsMax; ++i) {
    T a = lhs[i];
    // If you only use numeric T, you can optionally skip zeros:
    if (a == T{})
      continue;

    for (int j = rhsMin; j <= rhsMax; ++j) {
      T b = rhs[j];
      if (b == T{})
        continue;
      result[i + j] += a * b;
    }
  }

  return result;
}

template <typename T>
bilvector<T> &operator*=(bilvector<T> &lhs, const bilvector<T> &rhs) {
  lhs = lhs * rhs;
  return lhs;
}

// Multiply by q^power: shift exponents
template <typename T>
bilvector<T> multiplyByQPower(const bilvector<T> &poly, int power) {
  if (power == 0) {
    return poly;
  }

  int inMin = poly.getMaxNegativeIndex();
  int inMax = poly.getMaxPositiveIndex();
  int outMin = inMin + power;
  int outMax = inMax + power;

  int componentSize = poly.getComponentSize();

  int negativeCount = 0;
  if (outMin < 0) {
    negativeCount = -outMin;
  }

  int positiveCount = 0;
  if (outMax >= 0) {
    positiveCount = outMax + 1;
  }

  bilvector<T> result(negativeCount, positiveCount, componentSize, T{});

  for (int e = inMin; e <= inMax; ++e) {
    T c = poly[e];
    if (c == T{})
      continue;
    result[e + power] += c;
  }

  return result;
}

// Invert exponents of q in the Laurent polynomial
template <typename T> bilvector<T> bilvector<T>::invertExponents() const {
  bilvector<T> result(this->negativeVectorCount, this->positiveVectorCount,
                      this->componentSize, this->defaultValue);

  // Loop over all valid indices and assign to the mirrored exponent
  for (int j = this->getMaxNegativeIndex(); j <= this->getMaxPositiveIndex();
       ++j) {
    T coeff = (*this)[j];
    if (coeff != this->defaultValue) {
      result[-j] = coeff;
    }
  }

  return result;
}

// Print the Laurent polynomial in human-readable form
template <typename T>
void bilvector<T>::print(const std::string &varName) const {
  bool first = true;
  for (int e = this->getMaxNegativeIndex(); e <= this->getMaxPositiveIndex();
       ++e) {
    T coeff = (*this)[e];
    if (coeff == this->defaultValue || coeff == 0)
      continue;

    if (!first) {
      std::cout << " + ";
    } else {
      first = false;
    }

    std::cout << coeff;

    if (e != 0) {
      std::cout << "*" << varName << "^" << e;
    }
  }

  if (first) {
    std::cout << "0"; // all coefficients were zero
  }
  std::cout << std::endl;
}

void print_pterms(std::vector<bilvector<int>> polynomial_terms);
