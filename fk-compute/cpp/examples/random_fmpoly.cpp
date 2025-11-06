#include "fk/fmpoly.hpp"
#include <functional>
#include <iostream>
#include <random>

/**
 * Create a random polynomial in q, x₁, …, xₙ.
 *
 * The exponents for each x-variable range from minXPower to maxXPower
 * (inclusive), independently for each variable. For each multi-index in x, we
 * choose roughly `avgQPerXPower` different q-exponents between minQPower and
 * maxQPower, and assign them random integer coefficients.
 *
 * @param numXVariables   Number of x-variables (x₁,…,xₙ)
 * @param minXPower       Minimum exponent for each xᵢ
 * @param maxXPower       Maximum exponent for each xᵢ
 * @param minQPower       Minimum exponent for q
 * @param maxQPower       Maximum exponent for q
 * @param minCoeff        Minimum (integer) coefficient value
 * @param maxCoeff        Maximum (integer) coefficient value
 * @param avgQPerXPower   Desired average number of q-powers per x-multi-index
 * @param seed            Optional RNG seed (for reproducibility)
 *
 * @return                Random FMPoly
 */
FMPoly randomFMPoly(int numXVariables, int minXPower, int maxXPower,
                    int minQPower, int maxQPower, int minCoeff, int maxCoeff,
                    double avgQPerXPower,
                    std::uint64_t seed = std::random_device{}()) {
  if (numXVariables <= 0)
    throw std::invalid_argument("numXVariables must be positive");
  if (minXPower > maxXPower)
    throw std::invalid_argument("minXPower > maxXPower");
  if (minQPower > maxQPower)
    throw std::invalid_argument("minQPower > maxQPower");
  if (minCoeff > maxCoeff)
    throw std::invalid_argument("minCoeff > maxCoeff");

  // Use a degree large enough to comfortably cover the requested x-exponents.
  int maxAbsX = std::max(std::abs(minXPower), std::abs(maxXPower));
  FMPoly poly(numXVariables, maxAbsX);

  // RNG setup
  std::mt19937_64 rng(seed);

  const int numQ = maxQPower - minQPower + 1;
  double p = 0.0;
  if (numQ > 0) {
    p = avgQPerXPower / static_cast<double>(numQ);
    if (p < 0.0)
      p = 0.0;
    if (p > 1.0)
      p = 1.0;
  }

  std::bernoulli_distribution chooseTerm(p);
  std::uniform_int_distribution<int> coeffDist(minCoeff, maxCoeff);

  auto sampleNonZeroCoeff = [&]() -> int {
    if (minCoeff == 0 && maxCoeff == 0) {
      // Degenerate case: only zero available
      return 0;
    }
    int c;
    do {
      c = coeffDist(rng);
    } while (c == 0);
    return c;
  };

  // Enumerate all x-multi-indices within [minXPower, maxXPower] for each
  // variable
  std::vector<int> xPowers(numXVariables, minXPower);

  std::function<void(int)> recurse;
  recurse = [&](int var) {
    if (var == numXVariables) {
      // For this x-multi-index, decide q-powers
      for (int q = minQPower; q <= maxQPower; ++q) {
        if (chooseTerm(rng)) {
          int c = sampleNonZeroCoeff();
          if (c != 0) {
            poly.addToCoefficient(q, xPowers, c);
          }
        }
      }
      return;
    }

    for (int e = minXPower; e <= maxXPower; ++e) {
      xPowers[var] = e;
      recurse(var + 1);
    }
  };

  recurse(0);
  return poly;
}

int main() {
  try {
    // Example parameters
    int numX = 2; // two x variables: x1, x2
    int minX = -1;
    int maxX = 2;
    int minQ = 0;
    int maxQ = 3;
    int minCoeff = -5;
    int maxCoeff = 5;
    double avgQPerX = 1.5;

    // Create random polynomial
    FMPoly poly(2);
    poly.setQPolynomial({2, 3}, {1, 2, 3, 4, 5}, -3);

    // Print result
    std::cout << "Random polynomial:\n";
    poly.print(50); // or: std::cout << poly;
    std::cout << std::endl;
    auto qpoly = poly.getQPolynomial({2, 3});
    for (auto i : qpoly) {
      std::cout << i << std::endl;
    }
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
}
