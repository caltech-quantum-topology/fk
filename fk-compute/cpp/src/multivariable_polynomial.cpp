#include "fk/multivariable_polynomial.hpp"

MultivariablePolynomial::MultivariablePolynomial(
    int numVariables, int degree, const std::vector<int> &maxDegrees)
    : numXVariables(numVariables), maxDegree(degree) {

  if (maxDegrees.empty()) {
    maxXDegrees = std::vector<int>(numVariables, degree);
  } else {
    if (maxDegrees.size() != numVariables) {
      throw std::invalid_argument(
          "Max degrees vector size must match number of variables");
    }
    maxXDegrees = maxDegrees;
  }

  // Calculate block sizes for indexing
  blockSizes.resize(numVariables);
  blockSizes[0] = 1;
  for (int i = 1; i < numVariables; i++) {
    blockSizes[i] = (maxXDegrees[i - 1] + 1) * blockSizes[i - 1];
  }

  // Calculate total size needed
  int totalSize = 1;
  for (int i = 0; i < numVariables; i++) {
    totalSize *= (maxXDegrees[i] + 1);
  }

  // Initialize coefficient storage
  coefficients.resize(totalSize, bilvector<int>(0, 1, 20, 0));
}

int MultivariablePolynomial::multiIndexToLinear(
    const std::vector<int> &xPowers) const {
  if (xPowers.size() != numXVariables) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  int linearIndex = 0;
  for (int i = 0; i < numXVariables; i++) {
    if (xPowers[i] < 0 || xPowers[i] > maxXDegrees[i]) {
      throw std::out_of_range("X power out of bounds for variable " +
                              std::to_string(i));
    }
    linearIndex += xPowers[i] * blockSizes[i];
  }
  return linearIndex;
}

std::vector<int>
MultivariablePolynomial::linearToMultiIndex(int linearIndex) const {
  std::vector<int> xPowers(numXVariables);
  for (int i = numXVariables - 1; i >= 0; i--) {
    xPowers[i] = (linearIndex / blockSizes[i]) % (maxXDegrees[i] + 1);
  }
  return xPowers;
}

int MultivariablePolynomial::getCoefficient(
    int qPower, const std::vector<int> &xPowers) const {
  int linearIndex = multiIndexToLinear(xPowers);
  return coefficients[linearIndex][qPower];
}

void MultivariablePolynomial::setCoefficient(int qPower,
                                             const std::vector<int> &xPowers,
                                             int coefficient) {
  int linearIndex = multiIndexToLinear(xPowers);
  coefficients[linearIndex][qPower] = coefficient;
}

void MultivariablePolynomial::addToCoefficient(int qPower,
                                               const std::vector<int> &xPowers,
                                               int coefficient) {
  int linearIndex = multiIndexToLinear(xPowers);
  coefficients[linearIndex][qPower] += coefficient;
}

bilvector<int> &
MultivariablePolynomial::getQPolynomial(const std::vector<int> &xPowers) {
  int linearIndex = multiIndexToLinear(xPowers);
  return coefficients[linearIndex];
}

const bilvector<int> &
MultivariablePolynomial::getQPolynomial(const std::vector<int> &xPowers) const {
  int linearIndex = multiIndexToLinear(xPowers);
  return coefficients[linearIndex];
}

std::vector<bilvector<int>> &MultivariablePolynomial::getCoefficients() {
  return coefficients;
}

const std::vector<bilvector<int>> &
MultivariablePolynomial::getCoefficients() const {
  return coefficients;
}

int MultivariablePolynomial::getNumXVariables() const { return numXVariables; }

const std::vector<int> &MultivariablePolynomial::getMaxXDegrees() const {
  return maxXDegrees;
}

const std::vector<int> &MultivariablePolynomial::getBlockSizes() const {
  return blockSizes;
}

void MultivariablePolynomial::clear() {
  for (auto &bilvec : coefficients) {
    bilvec = bilvector<int>(0, 1, 20, 0);
  }
}

bool MultivariablePolynomial::isZero() const {
  for (const auto &bilvec : coefficients) {
    for (int j = bilvec.getMaxNegativeIndex();
         j <= bilvec.getMaxPositiveIndex(); j++) {
      if (bilvec[j] != 0) {
        return false;
      }
    }
  }
  return true;
}

void MultivariablePolynomial::exportToJson(const std::string &fileName) const {
  std::ofstream outputFile;
  outputFile.open(fileName + ".json");
  outputFile << "{\n\t\"coefficient_q_powers\":[\n";

  for (size_t i = 0; i < coefficients.size(); i++) {
    outputFile << "\t\t[";
    bool isFirstWrite = true;

    for (int j = coefficients[i].getMaxNegativeIndex();
         j <= coefficients[i].getMaxPositiveIndex(); j++) {
      if (coefficients[i][j] != 0) {
        if (!isFirstWrite) {
          outputFile << ",[";
        } else {
          outputFile << "[";
          isFirstWrite = false;
        }
        outputFile << j << "," << coefficients[i][j] << "]";
      }
    }

    if (i < coefficients.size() - 1) {
      outputFile << "],\n";
    } else {
      outputFile << "]\n";
    }
  }

  outputFile << "\t],\n";
  outputFile << "\t\"metadata\": {\n";
  outputFile << "\t\t\"num_x_variables\": " << numXVariables << ",\n";
  outputFile << "\t\t\"max_x_degrees\": [";
  for (int i = 0; i < numXVariables; i++) {
    outputFile << maxXDegrees[i];
    if (i < numXVariables - 1)
      outputFile << ",";
  }
  outputFile << "]\n";
  outputFile << "\t}\n}";
  outputFile.close();
}

void MultivariablePolynomial::print(int maxTerms) const {
  std::cout << "Multivariable Polynomial P(q";
  for (int i = 0; i < numXVariables; i++) {
    std::cout << ", x" << (i + 1);
  }
  std::cout << "):\n";

  int termCount = 0;
  for (size_t linearIndex = 0;
       linearIndex < coefficients.size() && termCount < maxTerms;
       linearIndex++) {
    std::vector<int> xPowers = linearToMultiIndex(linearIndex);

    for (int j = coefficients[linearIndex].getMaxNegativeIndex();
         j <= coefficients[linearIndex].getMaxPositiveIndex() &&
         termCount < maxTerms;
         j++) {
      int coeff = coefficients[linearIndex][j];
      if (coeff != 0) {
        if (termCount > 0)
          std::cout << " + ";
        std::cout << coeff << "*q^" << j;
        for (int k = 0; k < numXVariables; k++) {
          if (xPowers[k] > 0) {
            std::cout << "*x" << (k + 1) << "^" << xPowers[k];
          }
        }
        termCount++;
      }
    }
  }
  if (termCount == maxTerms) {
    std::cout << " + ...";
  }
  std::cout << std::endl;
}

// Arithmetic operations implementation

void MultivariablePolynomial::checkCompatibility(
    const MultivariablePolynomial &other) const {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }

  for (int i = 0; i < numXVariables; i++) {
    if (maxXDegrees[i] != other.maxXDegrees[i]) {
      throw std::invalid_argument(
          "Polynomials must have the same maximum degrees for each variable");
    }
  }
}

MultivariablePolynomial &
MultivariablePolynomial::operator+=(const MultivariablePolynomial &other) {
  checkCompatibility(other);

  // Iterate through all coefficient positions
  for (size_t linearIndex = 0; linearIndex < coefficients.size();
       linearIndex++) {
    const auto &otherBilvec = other.coefficients[linearIndex];
    auto &thisBilvec = coefficients[linearIndex];

    // Add coefficients from other polynomial's bilvector
    for (int j = otherBilvec.getMaxNegativeIndex();
         j <= otherBilvec.getMaxPositiveIndex(); j++) {
      int otherCoeff = otherBilvec[j];
      if (otherCoeff != 0) {
        thisBilvec[j] += otherCoeff;
      }
    }
  }

  return *this;
}

MultivariablePolynomial &
MultivariablePolynomial::operator-=(const MultivariablePolynomial &other) {
  checkCompatibility(other);

  // Iterate through all coefficient positions
  for (size_t linearIndex = 0; linearIndex < coefficients.size();
       linearIndex++) {
    const auto &otherBilvec = other.coefficients[linearIndex];
    auto &thisBilvec = coefficients[linearIndex];

    // Subtract coefficients from other polynomial's bilvector
    for (int j = otherBilvec.getMaxNegativeIndex();
         j <= otherBilvec.getMaxPositiveIndex(); j++) {
      int otherCoeff = otherBilvec[j];
      if (otherCoeff != 0) {
        thisBilvec[j] -= otherCoeff;
      }
    }
  }

  return *this;
}

MultivariablePolynomial &
MultivariablePolynomial::operator*=(const MultivariablePolynomial &other) {
  checkCompatibility(other);

  // Create a new polynomial to store the result
  MultivariablePolynomial result(numXVariables, maxDegree, maxXDegrees);

  // Multiply each term in this polynomial with each term in other polynomial
  for (size_t thisLinearIndex = 0; thisLinearIndex < coefficients.size();
       thisLinearIndex++) {
    std::vector<int> thisXPowers = linearToMultiIndex(thisLinearIndex);
    const auto &thisBilvec = coefficients[thisLinearIndex];

    for (size_t otherLinearIndex = 0;
         otherLinearIndex < other.coefficients.size(); otherLinearIndex++) {
      std::vector<int> otherXPowers =
          other.linearToMultiIndex(otherLinearIndex);
      const auto &otherBilvec = other.coefficients[otherLinearIndex];

      // Check if the product x-powers would be within bounds
      std::vector<int> productXPowers(numXVariables);
      bool withinBounds = true;
      for (int i = 0; i < numXVariables; i++) {
        productXPowers[i] = thisXPowers[i] + otherXPowers[i];
        if (productXPowers[i] > maxXDegrees[i]) {
          withinBounds = false;
          break;
        }
      }

      if (withinBounds) {
        // Multiply all q-coefficient combinations
        for (int thisQ = thisBilvec.getMaxNegativeIndex();
             thisQ <= thisBilvec.getMaxPositiveIndex(); thisQ++) {
          int thisCoeff = thisBilvec[thisQ];
          if (thisCoeff != 0) {
            for (int otherQ = otherBilvec.getMaxNegativeIndex();
                 otherQ <= otherBilvec.getMaxPositiveIndex(); otherQ++) {
              int otherCoeff = otherBilvec[otherQ];
              if (otherCoeff != 0) {
                int productQ = thisQ + otherQ;
                int productCoeff = thisCoeff * otherCoeff;
                result.addToCoefficient(productQ, productXPowers, productCoeff);
              }
            }
          }
        }
      }
    }
  }

  // Replace this polynomial's coefficients with the result
  *this = std::move(result);
  return *this;
}

// Friend function implementations
MultivariablePolynomial operator+(const MultivariablePolynomial &lhs,
                                  const MultivariablePolynomial &rhs) {
  MultivariablePolynomial result = lhs;
  result += rhs;
  return result;
}

MultivariablePolynomial operator-(const MultivariablePolynomial &lhs,
                                  const MultivariablePolynomial &rhs) {
  MultivariablePolynomial result = lhs;
  result -= rhs;
  return result;
}

MultivariablePolynomial operator*(const MultivariablePolynomial &lhs,
                                  const MultivariablePolynomial &rhs) {
  MultivariablePolynomial result = lhs;
  result *= rhs;
  return result;
}