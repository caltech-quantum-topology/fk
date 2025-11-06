#include "fk/fmpoly.hpp"
#include "fk/multivariable_polynomial.hpp"
#include <chrono>
#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <algorithm>

class PerformanceTimer {
private:
    std::chrono::high_resolution_clock::time_point start_time;

public:
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    double elapsed_ms() const {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        return duration.count() / 1000.0;  // Convert to milliseconds
    }
};

// Create a random polynomial for testing
void populateRandomPolynomial(FMPoly& poly, int numTerms, int maxQPower, int maxXPower, std::mt19937& gen) {
    std::uniform_int_distribution<int> qDist(-maxQPower, maxQPower);
    std::uniform_int_distribution<int> xDist(0, maxXPower);
    std::uniform_int_distribution<int> coeffDist(1, 100);

    int numVars = poly.getNumXVariables();

    for (int i = 0; i < numTerms; ++i) {
        int qPower = qDist(gen);
        std::vector<int> xPowers(numVars);
        for (int j = 0; j < numVars; ++j) {
            xPowers[j] = xDist(gen);
        }
        int coeff = coeffDist(gen);

        poly.setCoefficient(qPower, xPowers, coeff);
    }
}

void populateRandomPolynomial(MultivariablePolynomial& poly, int numTerms, int maxQPower, int maxXPower, std::mt19937& gen) {
    std::uniform_int_distribution<int> qDist(-maxQPower, maxQPower);
    std::uniform_int_distribution<int> xDist(0, maxXPower);
    std::uniform_int_distribution<int> coeffDist(1, 100);

    int numVars = poly.getNumXVariables();

    for (int i = 0; i < numTerms; ++i) {
        int qPower = qDist(gen);
        std::vector<int> xPowers(numVars);
        for (int j = 0; j < numVars; ++j) {
            xPowers[j] = xDist(gen);
        }
        int coeff = coeffDist(gen);

        poly.setCoefficient(qPower, xPowers, coeff);
    }
}

struct TestParams {
    int numVariables;
    int numTerms;
    int maxQPower;
    int maxXPower;
    int numIterations;
    std::string description;
};

void runPerformanceTest(const TestParams& params) {
    std::cout << "\n=== " << params.description << " ===\n";
    std::cout << "Variables: " << params.numVariables << ", Terms: " << params.numTerms
              << ", Max Q Power: " << params.maxQPower << ", Max X Power: " << params.maxXPower
              << ", Iterations: " << params.numIterations << "\n\n";

    // Setup random number generator
    std::mt19937 gen(42);  // Fixed seed for reproducibility

    PerformanceTimer timer;

    // Test FMPoly
    std::cout << "Testing FMPoly (FLINT-based):\n";
    std::vector<double> fmpolyTimes;

    for (int iter = 0; iter < params.numIterations; ++iter) {
        // Create polynomials
        FMPoly poly1(params.numVariables, params.maxXPower);
        FMPoly poly2(params.numVariables, params.maxXPower);

        // Populate with random data
        populateRandomPolynomial(poly1, params.numTerms, params.maxQPower, params.maxXPower, gen);
        populateRandomPolynomial(poly2, params.numTerms, params.maxQPower, params.maxXPower, gen);

        // Time the multiplication
        timer.start();
        FMPoly result = poly1 * poly2;
        double elapsed = timer.elapsed_ms();
        fmpolyTimes.push_back(elapsed);

        if (iter == 0) {
            std::cout << "  Sample result - is zero: " << (result.isZero() ? "yes" : "no") << "\n";
        }
    }

    // Test MultivariablePolynomial
    std::cout << "Testing MultivariablePolynomial (sparse hash-based):\n";
    std::vector<double> multipolyTimes;

    // Reset generator for fair comparison
    gen.seed(42);

    for (int iter = 0; iter < params.numIterations; ++iter) {
        // Create polynomials
        MultivariablePolynomial poly1(params.numVariables, params.maxXPower);
        MultivariablePolynomial poly2(params.numVariables, params.maxXPower);

        // Populate with random data
        populateRandomPolynomial(poly1, params.numTerms, params.maxQPower, params.maxXPower, gen);
        populateRandomPolynomial(poly2, params.numTerms, params.maxQPower, params.maxXPower, gen);

        // Time the multiplication
        timer.start();
        MultivariablePolynomial result = poly1 * poly2;
        double elapsed = timer.elapsed_ms();
        multipolyTimes.push_back(elapsed);

        if (iter == 0) {
            std::cout << "  Sample result - is zero: " << (result.isZero() ? "yes" : "no") << "\n";
        }
    }

    // Calculate statistics
    auto calcStats = [](const std::vector<double>& times) {
        double sum = 0, min_time = times[0], max_time = times[0];
        for (double t : times) {
            sum += t;
            min_time = std::min(min_time, t);
            max_time = std::max(max_time, t);
        }
        double avg = sum / times.size();

        // Calculate median
        std::vector<double> sorted_times = times;
        std::sort(sorted_times.begin(), sorted_times.end());
        double median = sorted_times[sorted_times.size() / 2];

        return std::make_tuple(avg, min_time, max_time, median);
    };

    auto [fmpoly_avg, fmpoly_min, fmpoly_max, fmpoly_median] = calcStats(fmpolyTimes);
    auto [multipoly_avg, multipoly_min, multipoly_max, multipoly_median] = calcStats(multipolyTimes);

    // Print results
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nResults (times in milliseconds):\n";
    std::cout << std::setw(25) << "Implementation" << std::setw(12) << "Average" << std::setw(12) << "Median"
              << std::setw(12) << "Min" << std::setw(12) << "Max" << "\n";
    std::cout << std::string(73, '-') << "\n";
    std::cout << std::setw(25) << "FMPoly (FLINT)" << std::setw(12) << fmpoly_avg << std::setw(12) << fmpoly_median
              << std::setw(12) << fmpoly_min << std::setw(12) << fmpoly_max << "\n";
    std::cout << std::setw(25) << "MultivariablePolynomial" << std::setw(12) << multipoly_avg << std::setw(12) << multipoly_median
              << std::setw(12) << multipoly_min << std::setw(12) << multipoly_max << "\n";

    // Calculate speedup
    double speedup = multipoly_avg / fmpoly_avg;
    std::cout << "\nSpeedup (MultivariablePoly / FMPoly): " << std::setprecision(2) << speedup << "x\n";

    if (speedup > 1.0) {
        std::cout << "✓ FMPoly is " << speedup << "x faster than MultivariablePolynomial\n";
    } else {
        std::cout << "⚠ MultivariablePolynomial is " << (1.0/speedup) << "x faster than FMPoly\n";
    }
}

int main() {
    std::cout << "=== Performance Comparison: FMPoly vs MultivariablePolynomial ===\n";
    std::cout << "Testing multiplication performance with various polynomial sizes\n";

    // Define test cases
    std::vector<TestParams> testCases = {
        {2, 10, 5, 3, 100, "Small polynomials (2 vars, 10 terms)"},
        {2, 50, 10, 5, 50, "Medium polynomials (2 vars, 50 terms)"},
        {3, 25, 8, 4, 50, "3D medium polynomials (3 vars, 25 terms)"},
        {2, 100, 15, 8, 20, "Large polynomials (2 vars, 100 terms)"},
        {4, 30, 6, 3, 30, "4D polynomials (4 vars, 30 terms)"},
        {2, 200, 20, 10, 10, "Very large polynomials (2 vars, 200 terms)"}
    };

    try {
        for (const auto& testCase : testCases) {
            runPerformanceTest(testCase);
        }

        std::cout << "\n=== Summary ===\n";
        std::cout << "Performance comparison completed successfully!\n";
        std::cout << "Results show relative performance of FLINT-based FMPoly\n";
        std::cout << "versus hash-based MultivariablePolynomial for multiplication.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error during performance testing: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}