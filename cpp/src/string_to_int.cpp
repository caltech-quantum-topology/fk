#include "fk/string_to_int.hpp"
#include <cmath>

int parseStringToInteger(std::string inputString) {
  int accumulator = 0;
  bool isNegative = false;
  int startIndex = 0;
  if (inputString[0] == '-') {
    isNegative = true;
    startIndex = 1;
  }
  for (size_t i = startIndex; i < inputString.size(); i++) {
    accumulator +=
        (inputString[i] - '0') * std::pow(10, inputString.size() - i - 1);
  }
  if (isNegative) {
    return -1 * accumulator;
  }
  return accumulator;
}

double parseStringToDouble(std::string inputString) {
  double accumulator = 0;
  bool isNegative = false;
  int startIndex = 0;
  if (inputString[0] == '-') {
    isNegative = true;
    startIndex = 1;
  }
  int decimalPointIndex = inputString.find('.');
  if (decimalPointIndex == -1) {
    decimalPointIndex = inputString.size();
  }
  for (int i = startIndex; i < decimalPointIndex; i++) {
    accumulator +=
        (inputString[i] - '0') * std::pow(10, decimalPointIndex - i - 1);
  }
  for (int i = decimalPointIndex + 1; i < (int)inputString.size(); i++) {
    accumulator +=
        (inputString[i] - '0') * std::pow(10, -(i - decimalPointIndex));
  }
  if (isNegative) {
    return -1 * accumulator;
  }
  return accumulator;
}