// "bi list vectors"
// description:

#pragma once

#include <list>
#include <vector>

template <typename T> struct bilvector {
private:
  int componentSize;
  int negativeVectorCount = 0;
  int positiveVectorCount = 0;
  int maxNegativeIndex = 0;
  int maxPositiveIndex = 0;
  T defaultValue;
  std::list<std::vector<T>> negativeVectors = {};
  std::list<std::vector<T>> positiveVectors = {};

public:
  bilvector(int initialNegativeVectorCount, int initialPositiveVectorCount,
            int componentSizeParam, T defaultValueParam) {
    componentSize = componentSizeParam;
    defaultValue = defaultValueParam;
    negativeVectorCount = initialNegativeVectorCount;
    for (int i = 0; i < initialNegativeVectorCount; i++) {
      negativeVectors.push_back(std::vector<T>(componentSize, defaultValue));
    }
    positiveVectorCount = initialPositiveVectorCount;
    for (int i = 0; i < initialPositiveVectorCount; i++) {
      positiveVectors.push_back(std::vector<T>(componentSize, defaultValue));
    }
  }
  int getNegativeSize() { return negativeVectorCount * componentSize; }
  int getNegativeSize() const { return negativeVectorCount * componentSize; }
  int getPositiveSize() { return positiveVectorCount * componentSize; }
  int getPositiveSize() const { return positiveVectorCount * componentSize; }
  T &operator[](int accessIndex) {
    if (accessIndex > maxPositiveIndex) {
      maxPositiveIndex = accessIndex;
    } else if (accessIndex < maxNegativeIndex) {
      maxNegativeIndex = accessIndex;
    }
    if (accessIndex >= 0) {
      if (accessIndex >= (*this).getPositiveSize()) {
        int x = (accessIndex - (*this).getPositiveSize()) / componentSize;
        positiveVectorCount += x + 1;
        for (int i = 0; i <= x; i++) {
          positiveVectors.push_back(
              std::vector<T>(componentSize, defaultValue));
        }
      }
      auto it = positiveVectors.begin();
      int j;
      for (j = 0; j < accessIndex / componentSize; j++) {
        ++it;
      }
      return (*it)[accessIndex - j * componentSize];
    } else {
      accessIndex = -1 - accessIndex;
      if (accessIndex >= (*this).getNegativeSize()) {
        int x = (accessIndex - (*this).getNegativeSize()) / componentSize;
        negativeVectorCount += x + 1;
        for (int i = 0; i <= x; i++) {
          negativeVectors.push_back(
              std::vector<T>(componentSize, defaultValue));
        }
      }
      auto it = negativeVectors.begin();
      int j;
      for (j = 0; j < accessIndex / componentSize; j++) {
        ++it;
      }
      return (*it)[accessIndex - j * componentSize];
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

  // Const version of operator[]
  const T &operator[](int accessIndex) const {
    // For const access, we can't modify the structure, so we need to handle
    // this carefully
    if (accessIndex >= 0) {
      if (accessIndex >= this->getPositiveSize()) {
        // Return default value for out-of-bounds const access
        static T defaultVal = defaultValue;
        return defaultVal;
      }
      auto it = positiveVectors.begin();
      int j;
      for (j = 0; j < accessIndex / componentSize; j++) {
        ++it;
      }
      return (*it)[accessIndex - j * componentSize];
    } else {
      accessIndex = -1 - accessIndex;
      if (accessIndex >= this->getNegativeSize()) {
        // Return default value for out-of-bounds const access
        static T defaultVal = defaultValue;
        return defaultVal;
      }
      auto it = negativeVectors.begin();
      int j;
      for (j = 0; j < accessIndex / componentSize; j++) {
        ++it;
      }
      return (*it)[accessIndex - j * componentSize];
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
bilvector<T> makeZeroLike(const bilvector<T> &proto,
                          int minExponent,
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
    if (a == T{}) continue;

    for (int j = rhsMin; j <= rhsMax; ++j) {
      T b = rhs[j];
      if (b == T{}) continue;
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
