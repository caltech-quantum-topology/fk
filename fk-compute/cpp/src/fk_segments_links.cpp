// implement multithreading at top level of recursion

#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "fk/bilvector.hpp"
#include "fk/linalg.hpp"
#include "fk/multivariable_polynomial.hpp"
#include "fk/qalg_links.hpp"
#include "fk/solution_pool_1a_double_links.hpp"
#include "fk/string_to_int.hpp"

class FK {
private:
  std::vector<int> angles;
  std::vector<int> angle_signs;
  std::vector<int> braid;
  std::vector<int> signs;
  std::vector<int> segments;
  std::vector<int> accumulatorBlockSizes = {1};
  std::vector<int> closed_strand_components = {};
  std::vector<std::vector<std::array<int, 2>>> crossingMatrices = {};
  std::vector<int> crossingRelationTypes = {};
  std::vector<std::vector<int>> transitions;
  std::vector<bool> trivial_angles_;
  std::vector<int> nontrivial_map;
  std::vector<int> inversion_data;
  MultivariablePolynomial result{
      1, 1}; // Initialize with dummy values, will be reassigned
  std::vector<std::vector<std::vector<int>>> variableAssignments;
  std::vector<std::vector<int>> numericalAssignments;
  std::vector<int> top_crossing_components = {};
  std::vector<int> bottom_crossing_components = {};
  int components;
  int writhe = 0;
  int prefactors;
  int crossings;
  int degree;
  void computeNumericalAssignment(const std::vector<int> &angles) {
    for (int i = 0; i < crossings + 1; i++) {
      for (int j = 0; j < prefactors + 1; j++) {
        numericalAssignments[i][j] =
            computeDotProduct(variableAssignments[i][j], angles);
      }
    }

    double qPowerAccumulatorDouble = (writhe - prefactors) / 2.0;
    std::vector<double> xPowerAccumulatorDouble(components, 0);
    int initialCoefficient = 1;
    for (int i = 0; i < prefactors; i++) {
      qPowerAccumulatorDouble -= numericalAssignments[0][i + 1];
      xPowerAccumulatorDouble[closed_strand_components[i]] -= 0.5;
    }
    xPowerAccumulatorDouble[0] -= 0.5;

    for (int crossingIndex = 0; crossingIndex < crossings; crossingIndex++) {
      int matrixParam_i =
          numericalAssignments[crossingMatrices[crossingIndex][0][0]]
                              [crossingMatrices[crossingIndex][0][1]];
      int matrixParam_k =
          numericalAssignments[crossingMatrices[crossingIndex][2][0]]
                              [crossingMatrices[crossingIndex][2][1]];
      int matrixParam_j =
          numericalAssignments[crossingMatrices[crossingIndex][1][0]]
                              [crossingMatrices[crossingIndex][1][1]];
      int matrixParam_m =
          numericalAssignments[crossingMatrices[crossingIndex][3][0]]
                              [crossingMatrices[crossingIndex][3][1]];
      int topComponent = top_crossing_components[crossingIndex];
      int bottomComponent = bottom_crossing_components[crossingIndex];
      if (crossingRelationTypes[crossingIndex] == 1 ||
          crossingRelationTypes[crossingIndex] == 2) {
        qPowerAccumulatorDouble += (matrixParam_j + matrixParam_m) / 2.0 +
                                   matrixParam_j * matrixParam_m;
        xPowerAccumulatorDouble[topComponent] +=
            ((matrixParam_j + matrixParam_k + 1) / 4.0);
        xPowerAccumulatorDouble[bottomComponent] +=
            ((3 * matrixParam_m - matrixParam_i + 1) / 4.0);
      } else if (crossingRelationTypes[crossingIndex] == 4) {
        qPowerAccumulatorDouble -= (matrixParam_i + matrixParam_k +
                                    matrixParam_m * (matrixParam_m + 1) -
                                    matrixParam_i * (matrixParam_i + 1)) /
                                       2.0 +
                                   matrixParam_i * matrixParam_k;
        xPowerAccumulatorDouble[topComponent] -=
            ((3 * matrixParam_j - matrixParam_k + 1) / 4.0);
        xPowerAccumulatorDouble[bottomComponent] -=
            ((matrixParam_i + matrixParam_m + 1) / 4.0);
        if ((matrixParam_j - matrixParam_k) % 2 == 0) {
          initialCoefficient *= -1;
        }
      } else if (crossingRelationTypes[crossingIndex] == 3) {
        qPowerAccumulatorDouble -= (matrixParam_i + matrixParam_k +
                                    matrixParam_m * (matrixParam_m + 1) -
                                    matrixParam_i * (matrixParam_i + 1)) /
                                       2.0 +
                                   matrixParam_i * matrixParam_k;
        xPowerAccumulatorDouble[topComponent] -=
            ((3 * matrixParam_j - matrixParam_k + 1) / 4.0);
        xPowerAccumulatorDouble[bottomComponent] -=
            ((matrixParam_i + matrixParam_m + 1) / 4.0);
        if ((matrixParam_k - matrixParam_j) % 2 == 1) {
          initialCoefficient *= -1;
        }
      }
    }
    std::cout << "qPowerAccumulatorDouble : " << qPowerAccumulatorDouble
              << "\n";
    int qPowerAccumulator = static_cast<int>(std::floor(
        qPowerAccumulatorDouble)); // currently, we are losing the actual q
                                   // powers because rounding down; later,
                                   // implement output with the exact q powers
    std::vector<int> xPowerAccumulator(components);
    std::vector<int> maxXDegrees(components);
    std::vector<int> blockSizes(components);
    blockSizes[0] = 1;
    for (int n = 0; n < components; n++) {
      xPowerAccumulator[n] = xPowerAccumulatorDouble[n];
      maxXDegrees[n] = degree - xPowerAccumulator[n];
      if (n != 0) {
        blockSizes[n] = (maxXDegrees[n - 1] + 1) * blockSizes[n - 1];
      }
    }
    int totalProductSize =
        blockSizes[components - 1] * (maxXDegrees[components - 1] + 1);
    // std::cout << "x_acc: " << xPowerAccumulatorDouble[0] << "\n\n";
    std::cout << "x_acc: " << xPowerAccumulatorDouble[0] << " "
              << xPowerAccumulatorDouble[1]
              << "\n\n"; // modifying: x_acc's are coming out to half-integers
                         // in general
    // exit(0);
    // if (xPowerAccumulatorDouble[0] > 15) {
    //     std::cout << "exiting due to x_acc large enough\n";
    //     exit(0);
    // }
    std::vector<bilvector<int>> polynomialTerms(
        totalProductSize,
        bilvector<int>(0, 1, 20, 0)); // error is the degree issue; modifying
    polynomialTerms[0][0] = initialCoefficient;
    for (int crossingIndex = 0; crossingIndex < crossings; crossingIndex++) {
      if (crossingRelationTypes[crossingIndex] == 1) {
        int param_i =
            numericalAssignments[crossingMatrices[crossingIndex][0][0]]
                                [crossingMatrices[crossingIndex][0][1]];
        int param_m =
            numericalAssignments[crossingMatrices[crossingIndex][3][0]]
                                [crossingMatrices[crossingIndex][3][1]];
        if (param_i > 0) {
          computePositiveQBinomial(polynomialTerms, param_i, param_i - param_m,
                                   false);
        } else {
          computeNegativeQBinomial(polynomialTerms, param_i, param_i - param_m,
                                   false);
        }
      } else if (crossingRelationTypes[crossingIndex] == 2) {
        int param_i =
            numericalAssignments[crossingMatrices[crossingIndex][0][0]]
                                [crossingMatrices[crossingIndex][0][1]];
        int param_m =
            numericalAssignments[crossingMatrices[crossingIndex][3][0]]
                                [crossingMatrices[crossingIndex][3][1]];
        computeNegativeQBinomial(polynomialTerms, param_i, param_m, false);
      } else if (crossingRelationTypes[crossingIndex] == 3) {
        int param_j =
            numericalAssignments[crossingMatrices[crossingIndex][1][0]]
                                [crossingMatrices[crossingIndex][1][1]];
        int param_k =
            numericalAssignments[crossingMatrices[crossingIndex][2][0]]
                                [crossingMatrices[crossingIndex][2][1]];
        computeNegativeQBinomial(polynomialTerms, param_j, param_k, true);
      } else {
        int param_j =
            numericalAssignments[crossingMatrices[crossingIndex][1][0]]
                                [crossingMatrices[crossingIndex][1][1]];
        int param_k =
            numericalAssignments[crossingMatrices[crossingIndex][2][0]]
                                [crossingMatrices[crossingIndex][2][1]];
        if (param_j > 0) {
          computePositiveQBinomial(polynomialTerms, param_j, param_j - param_k,
                                   true);
        } else {
          computeNegativeQBinomial(polynomialTerms, param_j, param_j - param_k,
                                   true);
        }
      }
    }
    for (int crossingIndex = 0; crossingIndex < crossings; crossingIndex++) {
      if (crossingRelationTypes[crossingIndex] == 1) {
        int bottomComp = bottom_crossing_components[crossingIndex];
        int param_j =
            numericalAssignments[crossingMatrices[crossingIndex][1][0]]
                                [crossingMatrices[crossingIndex][1][1]];
        int param_k =
            numericalAssignments[crossingMatrices[crossingIndex][2][0]]
                                [crossingMatrices[crossingIndex][2][1]];
        computeXQPochhammer(polynomialTerms, param_k, param_j + 1, bottomComp,
                            components, maxXDegrees, blockSizes);
      } else if (crossingRelationTypes[crossingIndex] == 2) {
        int bottomComp = bottom_crossing_components[crossingIndex];
        int param_j =
            numericalAssignments[crossingMatrices[crossingIndex][1][0]]
                                [crossingMatrices[crossingIndex][1][1]];
        int param_k =
            numericalAssignments[crossingMatrices[crossingIndex][2][0]]
                                [crossingMatrices[crossingIndex][2][1]];
        computeXQInversePochhammer(polynomialTerms, param_j, param_k + 1,
                                   bottomComp, components, maxXDegrees,
                                   blockSizes);
      } else if (crossingRelationTypes[crossingIndex] == 3) {
        int topComp = top_crossing_components[crossingIndex];
        int param_i =
            numericalAssignments[crossingMatrices[crossingIndex][0][0]]
                                [crossingMatrices[crossingIndex][0][1]];
        int param_m =
            numericalAssignments[crossingMatrices[crossingIndex][3][0]]
                                [crossingMatrices[crossingIndex][3][1]];
        computeXQInversePochhammer(polynomialTerms, param_i, param_m + 1,
                                   topComp, components, maxXDegrees,
                                   blockSizes);
      } else {
        int topComp = top_crossing_components[crossingIndex];
        int param_i =
            numericalAssignments[crossingMatrices[crossingIndex][0][0]]
                                [crossingMatrices[crossingIndex][0][1]];
        int param_m =
            numericalAssignments[crossingMatrices[crossingIndex][3][0]]
                                [crossingMatrices[crossingIndex][3][1]];
        computeXQPochhammer(polynomialTerms, param_m, param_i + 1, topComp,
                            components, maxXDegrees, blockSizes);
      }
    }
    auto& resultCoeffs = result.getCoefficients();
    performOffsetAddition(resultCoeffs, polynomialTerms,
                          xPowerAccumulator, qPowerAccumulator, components,
                          maxXDegrees, 1, accumulatorBlockSizes, blockSizes);
    result.syncFromDenseVector(resultCoeffs);
  }
  void writeResultsToJson(std::string fileName) {
    result.exportToJson(fileName);
  }

public:
  std::vector<std::vector<double>> inequalities;
  std::vector<std::vector<double>> criteria;
  std::vector<std::vector<int>> extensions;
  std::string metadata;
  FK(std::string infile_, std::string outfile_) {

    std::ifstream infile;
    infile.open(infile_ + ".csv");
    if (infile.is_open()) {
      std::string line;

      std::getline(infile, line, '\n');
      int index = line.find(",");
      degree = parseStringToInteger(line.substr(0, index));

      std::getline(infile, line, '\n');
      index = line.find(",");
      components = parseStringToInteger(line.substr(0, index));

      std::getline(infile, line, '\n');
      index = line.find(",");
      writhe = parseStringToInteger(line.substr(0, index));

      // Initialize the polynomial after we know the degree and components
      result = MultivariablePolynomial(components, degree);

      std::getline(infile, line, '\n');
      int height = 0;
      index = line.find(",");
      bool done = false;
      while (true) {
        if (index == -1) {
          break;
        }

        int c = parseStringToInteger(line.substr(0, index));
        crossingMatrices.push_back({{height, c - 1},
                                    {height, c},
                                    {height + 1, c - 1},
                                    {height + 1, c}});
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");

        crossingRelationTypes.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);

        height++;
        index = line.find(",");
      }
      std::getline(infile, line, '\n');
      index = line.find(",");
      done = false;
      while (true) {
        if (index == -1) {
          break;
        }
        closed_strand_components.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");
      }
      prefactors = closed_strand_components.size();
      crossings = crossingRelationTypes.size();
      variableAssignments.resize(crossings + 1,
                                 std::vector<std::vector<int>>(prefactors + 1));
      numericalAssignments.resize(crossings + 1,
                                  std::vector<int>(prefactors + 1));
      std::getline(infile, line, '\n');
      index = line.find(",");
      done = false;
      while (true) {
        if (index == -1) {
          break;
        }
        top_crossing_components.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");
        bottom_crossing_components.push_back(
            parseStringToInteger(line.substr(0, index)));
        line = line.substr(index + 1, line.size() - index - 1);
        index = line.find(",");
      }
      int stage = 0;
      int criteria_index = 0;
      int inequality_index = 0;
      int extension_index = 0;
      while (std::getline(infile, line, '\n')) {
        if (line[0] == '/') {
          stage++;
        } else if (stage == 0) {
          criteria.push_back({});
          done = false;
          index = line.find(",");
          while (true) {
            if (index == -1) {
              break;
            }
            criteria[criteria_index].push_back(
                parseStringToDouble(line.substr(0, index)));
            line = line.substr(index + 1, line.size() - index - 1);
            index = line.find(",");
          }
          criteria_index++;
        } else if (stage == 1) {
          inequalities.push_back({});
          done = false;
          index = line.find(",");
          while (true) {
            if (index == -1) {
              break;
            }
            inequalities[inequality_index].push_back(
                parseStringToInteger(line.substr(0, index)));
            line = line.substr(index + 1, line.size() - index - 1);
            index = line.find(",");
          }
          inequality_index++;
        } else if (stage == 2) {
          done = false;
          index = line.find(",");
          while (true) {
            if (index == -1) {
              break;
            }
            variableAssignments[extension_index % (crossings + 1)]
                               [extension_index / (crossings + 1)]
                                   .push_back(parseStringToInteger(
                                       line.substr(0, index)));
            line = line.substr(index + 1, line.size() - index - 1);
            index = line.find(",");
          }
          extension_index++;
        }
      }
    } else {
      std::cout << "ERROR: Unable to open file '" + infile_ + ".csv'!";
      exit(0);
    }
    for (int i = 1; i < components; i++) {
      accumulatorBlockSizes.push_back(accumulatorBlockSizes[i - 1] *
                                      (degree + 1));
    }
    std::function<void(const std::vector<int> &)> f_wrapper =
        [this](const std::vector<int> &v) { computeNumericalAssignment(v); };
    pooling(criteria, inequalities, f_wrapper);
    std::vector<int> increment_offset(components);
    increment_offset[0] = 1;
    std::vector<int> maxima(components, degree - 1);

    // for (int w = 0; w < degree + 1; w++) {
    //     for (int a = 0; a < degree + 1; a++) {
    //         for (int j = acc[w * (degree + 1) + a].get_max_nindex(); j <=
    //         acc[w * (degree + 1) + a].get_max_pindex(); j++) {
    //             std::cout << acc[w * (degree + 1) + a][j] << " ";
    //         }
    //         std::cout << "\t";
    //     }
    //     std::cout << "\n";
    // }
    // std::cout << std::endl;
    // std::cout << accumulatorBlockSizes.size() << " " <<
    // increment_offset.size() << " " << components << " " << maxima.size() <<
    // "\n";
    auto& resultCoeffs1 = result.getCoefficients();
    auto& resultCoeffs2 = result.getCoefficients();
    performOffsetAddition(resultCoeffs1, resultCoeffs2,
                          increment_offset, 0, components, maxima, -1,
                          accumulatorBlockSizes, accumulatorBlockSizes);
    result.syncFromDenseVector(resultCoeffs1);
    writeResultsToJson(outfile_);

    // for (int w = 0; w < degree + 1; w++) {
    //     for (int a = 0; a < degree + 1; a++) {
    //         for (int j = acc[w * (degree + 1) + a].get_max_nindex(); j <=
    //         acc[w * (degree + 1) + a].get_max_pindex(); j++) {
    //             std::cout << acc[w * (degree + 1) + a][j] << " ";
    //         }
    //         std::cout << "\t";
    //     }
    //     std::cout << "\n";
    // }
    std::cout << std::endl;
  }
};

int main(int argc, char *argv[]) {
  std::cout << argc << std::endl;
  /*
  if (argc < 4){
    std::cerr << "Usage error: input and output file expected in argv" <<
  std::endl; return 1;
  }
  */

  const std::string in = argv[1];
  const std::string out = argv[2];
  std::cout << in << std::endl;
  std::cout << out << std::endl;
  try {
    FK(in, out);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

// implement multithreading at top level of recursion
