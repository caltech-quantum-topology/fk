#include "fk/bilvector.hpp"
#include "fk/qalg_links.hpp"
#include <fstream>
#include <iostream>
#include <vector>

void writeBilvectorToJson(std::vector<bilvector<int>> &bilvectorTerms,
                          const std::string &outputFileName,
                          const std::string &testName) {
  std::ofstream outputFile;
  outputFile.open(outputFileName);
  outputFile << "{\n\t\"test_name\": \"" << testName << "\",\n";
  outputFile << "\t\"coefficient_q_powers\":[\n";

  for (size_t i = 0; i < bilvectorTerms.size(); i++) {
    outputFile << "\t\t[";
    bool isFirstWrite = true;
    for (int j = bilvectorTerms[i].getMaxNegativeIndex();
         j <= bilvectorTerms[i].getMaxPositiveIndex(); j++) {
      if (bilvectorTerms[i][j] != 0) {
        if (!isFirstWrite) {
          outputFile << ",[";
        } else {
          outputFile << "[";
          isFirstWrite = false;
        }
        outputFile << j << "," << bilvectorTerms[i][j] << "]";
      }
    }
    if (i < bilvectorTerms.size() - 1) {
      outputFile << "],\n";
    } else {
      outputFile << "]\n";
    }
  }
  outputFile << "\t]\n}";
  outputFile.close();
}

int main() {
  std::cout << "Testing computeXQPochhammer function\n";
  std::cout << "================================\n\n";

  // Test case 1: Simple case
  {
    std::cout << "Test 1: up=3, low=1, component=0, components=2\n";
    int components = 1;
    std::vector<int> lengths = {5}; // Component sizes
    std::vector<int> blocks = {1};  // Block structure
    int prod = 20;                  // Large enough size for the term vector

    // Initialize term with bilvectors - need enough elements for indexing
    std::vector<bilvector<int>> term(prod, bilvector<int>(0, 1, 20, 0));

    // Initialize first element as in the original code
    term[0][0] = 1;

    std::cout << "Running computeXQPochhammer(term, 3, 1, 0, components, "
                 "lengths, blocks)...\n";
    computeXQPochhammer(term, 3, 1, 0, components, lengths, blocks);

    writeBilvectorToJson(term, "test1_output.json",
                         "computeXQPochhammer(3,1,0)");
  }

  // Test case 2: Different parameters
  {
    std::cout << "Test 2: up=2, low=0, component=1, components=3\n";
    int components = 1;
    std::vector<int> lengths = {5}; // Component sizes
    std::vector<int> blocks = {1};  // Block structure
    int prod = 20;                  // Large enough size for the term vector

    // Initialize term with bilvectors - need enough elements for indexing
    std::vector<bilvector<int>> term(prod, bilvector<int>(0, 1, 20, 0));

    // Initialize first element as in the original code
    term[0][0] = 1;

    std::cout << "Running computeXQPochhammer(term, 2, 0, 1, components, "
                 "lengths, blocks)...\n";
    computeXQPochhammer(term, 12, 1, 0, components, lengths, blocks);

    writeBilvectorToJson(term, "test2_output.json",
                         "computeXQPochhammer(2,0,1)");
  }

  // Test case 3: Edge case with up = low
  {
    std::cout << "Test 3: up=1, low=1, component=0, components=1\n";
    int components = 1;
    std::vector<int> lengths = {5};
    std::vector<int> blocks = {2};
    int prod = 20; // Large enough size for the term vector

    // Initialize term with bilvectors - need enough elements for indexing
    std::vector<bilvector<int>> term(prod, bilvector<int>(0, 1, 20, 0));

    // Initialize first element as in the original code
    term[0][0] = 1;

    std::cout << "Running computeXQPochhammer(term, 1, 1, 0, components, "
                 "lengths, blocks)...\n";
    computeXQPochhammer(term, 3, 1, 0, components, lengths, blocks);

    writeBilvectorToJson(term, "test3_output.json",
                         "computeXQPochhammer(1,1,0)");
  }

  std::cout << "\nCompleted! Results saved to:\n";
  std::cout << "- test1_output.json\n";
  std::cout << "- test2_output.json\n";
  std::cout << "- test3_output.json\n";

  return 0;
}
