#ifndef FK_FK_HPP
#define FK_FK_HPP

#include <array>
#include <string>
#include <vector>

#include "bilvector.hpp"
#include "multivariable_polynomial.hpp"

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

  void computeNumericalAssignment(const std::vector<int> &angles);
  void writeResultsToJson(std::string fileName);

public:
  std::vector<std::vector<double>> inequalities;
  std::vector<std::vector<double>> criteria;
  std::vector<std::vector<int>> extensions;
  std::string metadata;

  FK(std::string infile_, std::string outfile_);
};

#endif // FK_FK_HPP