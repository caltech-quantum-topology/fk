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