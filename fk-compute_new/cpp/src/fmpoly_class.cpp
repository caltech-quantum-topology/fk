#include "fk/fmpoly_class.hpp"
#include <algorithm>
#include <atomic>
#include <flint/fmpz.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace {

// FLINT mpoly contexts are immutable after initialization, so one context per
// variable count is shared by all FMPoly instances and threads. Contexts are
// created lazily and intentionally never freed (they live for the whole run).
constexpr int MAX_FAST_CTX_VARS = 64;
std::atomic<const fmpz_mpoly_ctx_struct *> fast_ctx_cache[MAX_FAST_CTX_VARS + 1];
std::mutex ctx_mutex;
std::unordered_map<int, const fmpz_mpoly_ctx_struct *> slow_ctx_cache;

const fmpz_mpoly_ctx_struct *sharedContext(int nvars) {
  if (nvars <= MAX_FAST_CTX_VARS) {
    const fmpz_mpoly_ctx_struct *c =
        fast_ctx_cache[nvars].load(std::memory_order_acquire);
    if (c) {
      return c;
    }
  }
  std::lock_guard<std::mutex> lock(ctx_mutex);
  auto it = slow_ctx_cache.find(nvars);
  if (it == slow_ctx_cache.end()) {
    auto *c = new fmpz_mpoly_ctx_struct[1];
    fmpz_mpoly_ctx_init(c, nvars, ORD_LEX);
    it = slow_ctx_cache.emplace(nvars, c).first;
    if (nvars <= MAX_FAST_CTX_VARS) {
      fast_ctx_cache[nvars].store(c, std::memory_order_release);
    }
  }
  return it->second;
}

// Build `out = in * monomial(shift)`, re-expressing `in` after its ground
// powers were lowered by `shift` (all entries nonnegative).
void buildShiftedCopy(fmpz_mpoly_t out, const fmpz_mpoly_t in,
                      const std::vector<ulong> &shift,
                      const fmpz_mpoly_ctx_struct *ctx) {
  fmpz_mpoly_t mono;
  fmpz_mpoly_init(mono, ctx);
  fmpz_mpoly_set_coeff_si_ui(mono, 1, shift.data(), ctx);
  fmpz_mpoly_mul(out, in, mono, ctx);
  fmpz_mpoly_clear(mono, ctx);
}

} // namespace


FMPoly::FMPoly(int numVariables, int degree, const std::vector<int> &maxDegrees)
    : numXVariables(numVariables), allGroundPowers(numVariables + 1, 0) {

  if (maxDegrees.empty()) {
    maxXDegrees = std::vector<int>(numVariables, degree);
  } else {
    if (maxDegrees.size() != static_cast<size_t>(numVariables)) {
      throw std::invalid_argument(
          "Max degrees vector size must match number of variables");
    }
    maxXDegrees = maxDegrees;
  }

  setupContext();
}

FMPoly::FMPoly(const FMPoly &source, int newNumVariables,
               int targetVariableIndex, int degree,
               const std::vector<int> &maxDegrees)
    : numXVariables(newNumVariables), allGroundPowers(newNumVariables + 1, 0) {

  // Copy q ground power from source
  allGroundPowers[0] = source.allGroundPowers[0];

  if (newNumVariables < source.numXVariables) {
    throw std::invalid_argument(
        "New number of variables must be >= source's number of variables");
  }

  if (targetVariableIndex < 0 || targetVariableIndex >= newNumVariables) {
    throw std::invalid_argument(
        "Target variable index must be in range [0, newNumVariables)");
  }

  if (source.numXVariables != 1) {
    throw std::invalid_argument(
        "Source polynomial must have exactly 1 x variable");
  }

  // Set up max degrees
  if (maxDegrees.empty()) {
    maxXDegrees = std::vector<int>(newNumVariables, degree);
  } else {
    if (maxDegrees.size() != static_cast<size_t>(newNumVariables)) {
      throw std::invalid_argument(
          "Max degrees vector size must match new number of variables");
    }
    maxXDegrees = maxDegrees;
  }

  setupContext();

  // Copy terms from source, mapping the single x-variable to
  // targetVariableIndex. Source terms are in canonical LEX order over
  // (q, x0); mapping x0 to a fixed destination slot preserves that order,
  // so the terms can be pushed directly.
  slong numTerms = fmpz_mpoly_length(source.poly, source.ctx);
  int dst_nvars = numXVariables + 1;

  ulong src_exp[2];
  std::vector<ulong> dst_exp(dst_nvars, 0);
  fmpz_t coeff;
  fmpz_init(coeff);

  for (slong t = 0; t < numTerms; t++) {
    fmpz_mpoly_get_term_coeff_fmpz(coeff, source.poly, t, source.ctx);
    fmpz_mpoly_get_term_exp_ui(src_exp, source.poly, t, source.ctx);

    dst_exp[0] = src_exp[0];                      // q exponent
    dst_exp[targetVariableIndex + 1] = src_exp[1]; // x exponent

    fmpz_mpoly_push_term_fmpz_ui(poly, coeff, dst_exp.data(), ctx);
  }

  fmpz_clear(coeff);
}

FMPoly::FMPoly(const FMPoly &other)
    : numXVariables(other.numXVariables),
      allGroundPowers(other.allGroundPowers), maxXDegrees(other.maxXDegrees),
      blockSizes(other.blockSizes) {

  setupContext();
  fmpz_mpoly_set(poly, other.poly, ctx);
}

FMPoly::FMPoly(FMPoly &&other) noexcept
    : numXVariables(other.numXVariables),
      allGroundPowers(std::move(other.allGroundPowers)),
      ctx(other.ctx), maxXDegrees(std::move(other.maxXDegrees)),
      blockSizes(std::move(other.blockSizes)) {
  fmpz_mpoly_init(poly, ctx);
  fmpz_mpoly_swap(poly, other.poly, ctx);
}

FMPoly &FMPoly::operator=(FMPoly &&other) noexcept {
  if (this != &other) {
    fmpz_mpoly_clear(poly, ctx);
    numXVariables = other.numXVariables;
    allGroundPowers = std::move(other.allGroundPowers);
    ctx = other.ctx;
    maxXDegrees = std::move(other.maxXDegrees);
    blockSizes = std::move(other.blockSizes);
    fmpz_mpoly_init(poly, ctx);
    fmpz_mpoly_swap(poly, other.poly, ctx);
  }
  return *this;
}

FMPoly &FMPoly::operator=(const FMPoly &other) {
  if (this != &other) {
    if (numXVariables != other.numXVariables) {
      fmpz_mpoly_clear(poly, ctx);
      numXVariables = other.numXVariables;
      ctx = sharedContext(numXVariables + 1);
      fmpz_mpoly_init(poly, ctx);
    }
    allGroundPowers = other.allGroundPowers;
    maxXDegrees = other.maxXDegrees;
    fmpz_mpoly_set(poly, other.poly, ctx);
  }
  return *this;
}

FMPoly::~FMPoly() { fmpz_mpoly_clear(poly, ctx); }

void FMPoly::setupContext() {
  // q, x1, ..., xn share one immutable FLINT context per variable count.
  ctx = sharedContext(numXVariables + 1);
  fmpz_mpoly_init(poly, ctx);
}

void FMPoly::storedExponents(int qPower, const std::vector<int> &xPowers,
                             ulong *exp) const {
  exp[0] = static_cast<ulong>(qPower - allGroundPowers[0]);
  for (int i = 0; i < numXVariables; i++) {
    exp[i + 1] = static_cast<ulong>(xPowers[i] - allGroundPowers[i + 1]);
  }
}

bool FMPoly::needsGroundAdjustment(int qPower,
                                   const std::vector<int> &xPowers) const {
  if (qPower < allGroundPowers[0]) {
    return true;
  }
  for (int i = 0; i < numXVariables; i++) {
    if (xPowers[i] < allGroundPowers[i + 1]) {
      return true;
    }
  }
  return false;
}

bool FMPoly::getExponentsFromMonomial(const fmpz *exps, int &qPower,
                                      std::vector<int> &xPowers) const {
  xPowers.resize(numXVariables);

  // Extract q power (handle offset)
  qPower = fmpz_get_si(&exps[0]) + allGroundPowers[0];

  // Extract x powers (reverse ground powers offset)
  for (int i = 0; i < numXVariables; i++) {
    xPowers[i] = fmpz_get_si(&exps[i + 1]) + allGroundPowers[i + 1];
  }

  return true;
}

void FMPoly::adjustGroundPowersIfNeeded(int qPower,
                                        const std::vector<int> &xPowers) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  // Check if we need to adjust any ground powers
  bool needToAdjust = false;
  std::vector<int> newGroundPowers = allGroundPowers;

  // Check q power
  if (qPower < allGroundPowers[0]) {
    newGroundPowers[0] = qPower;
    needToAdjust = true;
  }

  // Check x powers
  for (int i = 0; i < numXVariables; i++) {
    if (xPowers[i] < allGroundPowers[i + 1]) {
      newGroundPowers[i + 1] = xPowers[i];
      needToAdjust = true;
    }
  }

  if (!needToAdjust) {
    return;
  }

  // Lowering the ground powers raises every stored exponent by the same
  // shift, which is a single monomial multiplication.
  std::vector<ulong> shift(numXVariables + 1);
  for (int j = 0; j <= numXVariables; j++) {
    shift[j] = static_cast<ulong>(allGroundPowers[j] - newGroundPowers[j]);
  }

  fmpz_mpoly_t newPoly;
  fmpz_mpoly_init(newPoly, ctx);
  buildShiftedCopy(newPoly, poly, shift, ctx);
  fmpz_mpoly_swap(poly, newPoly, ctx);
  fmpz_mpoly_clear(newPoly, ctx);
  allGroundPowers = newGroundPowers;
}

void FMPoly::setCoefficient(int qPower, const std::vector<int> &xPowers,
                            int coefficient) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  // Adjust ground powers if needed to handle negative exponents
  adjustGroundPowersIfNeeded(qPower, xPowers);

  std::vector<ulong> exp(numXVariables + 1);
  storedExponents(qPower, xPowers, exp.data());

  fmpz_t coeff;
  fmpz_init(coeff);
  fmpz_set_si(coeff, coefficient);
  fmpz_mpoly_set_coeff_fmpz_ui(poly, coeff, exp.data(), ctx);
  fmpz_clear(coeff);
}

void FMPoly::addToCoefficient(int qPower, const std::vector<int> &xPowers,
                              int coefficient) {
  if (coefficient == 0) {
    return;
  }

  fmpz_t coeff;
  fmpz_init(coeff);
  fmpz_set_si(coeff, coefficient);
  addToCoefficientFmpz(qPower, xPowers, coeff);
  fmpz_clear(coeff);
}

void FMPoly::addToCoefficientFmpz(int qPower, const std::vector<int> &xPowers,
                                  const fmpz_t coefficient) {
  if (fmpz_is_zero(coefficient)) return;

  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  if (needsGroundAdjustment(qPower, xPowers)) {
    adjustGroundPowersIfNeeded(qPower, xPowers);
  }

  std::vector<ulong> exp(numXVariables + 1);
  storedExponents(qPower, xPowers, exp.data());

  fmpz_t current;
  fmpz_init(current);
  fmpz_mpoly_get_coeff_fmpz_ui(current, poly, exp.data(), ctx);
  fmpz_add(current, current, coefficient);
  fmpz_mpoly_set_coeff_fmpz_ui(poly, current, exp.data(), ctx);
  fmpz_clear(current);
}

void FMPoly::addQPolynomial(const std::vector<int> &xPowers,
                            const QPolynomial &qPoly) {
  if (xPowers.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "X powers vector size must match number of variables");
  }

  if (qPoly.isZero())
    return;

  // One ground adjustment up front covers every q power in the range
  if (needsGroundAdjustment(qPoly.getMinPower(), xPowers)) {
    adjustGroundPowersIfNeeded(qPoly.getMinPower(), xPowers);
  }

  std::vector<ulong> exp(numXVariables + 1);
  fmpz_t coeff, current;
  fmpz_init(coeff);
  fmpz_init(current);

  for (int qPower = qPoly.getMinPower(); qPower <= qPoly.getMaxPower();
       qPower++) {
    qPoly.getCoefficientFmpz(coeff, qPower);
    if (!fmpz_is_zero(coeff)) {
      storedExponents(qPower, xPowers, exp.data());
      fmpz_mpoly_get_coeff_fmpz_ui(current, poly, exp.data(), ctx);
      fmpz_add(current, current, coeff);
      fmpz_mpoly_set_coeff_fmpz_ui(poly, current, exp.data(), ctx);
    }
  }

  fmpz_clear(coeff);
  fmpz_clear(current);
}

FMPoly FMPoly::truncate(const std::vector<int> &maxXdegrees) const {
  if (maxXdegrees.size() != static_cast<size_t>(numXVariables)) {
    throw std::invalid_argument(
        "Max degrees vector size must match number of variables");
  }

  // Result: same x-structure & ground powers, but coefficients truncated
  FMPoly result(numXVariables, 0, maxXdegrees);
  result.allGroundPowers = allGroundPowers;

  slong numTerms = fmpz_mpoly_length(poly, ctx);
  if (numTerms == 0) {
    return result;
  }

  fmpz_t coeff;
  fmpz_init(coeff);
  std::vector<ulong> exp(numXVariables + 1);

  // stored_exp[j+1] = x_j - allGroundPowers[j+1], so the stored cutoff is
  // maxXdegrees[j] - allGroundPowers[j+1]
  std::vector<long> max_stored_x(numXVariables);
  for (int j = 0; j < numXVariables; ++j) {
    max_stored_x[j] = maxXdegrees[j] - allGroundPowers[j + 1];
  }

  // Keeping a subsequence of canonically ordered terms preserves the
  // order, so terms can be pushed directly without sorted insertion.
  for (slong i = 0; i < numTerms; ++i) {
    fmpz_mpoly_get_term_exp_ui(exp.data(), poly, i, ctx);

    bool include = true;
    for (int j = 0; j < numXVariables; ++j) {
      if (static_cast<long>(exp[j + 1]) > max_stored_x[j]) {
        include = false;
        break;
      }
    }
    if (!include) {
      continue;
    }

    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);
    fmpz_mpoly_push_term_fmpz_ui(result.poly, coeff, exp.data(), ctx);
  }

  fmpz_clear(coeff);
  return result;
}

FMPoly FMPoly::truncate(int maxDegree) const {
  std::vector<int> maxXdegrees(numXVariables, maxDegree);
  return truncate(maxXdegrees);
}

int FMPoly::getNumXVariables() const { return numXVariables; }

void FMPoly::clear() { fmpz_mpoly_zero(poly, ctx); }

bool FMPoly::isZero() const { return fmpz_mpoly_is_zero(poly, ctx); }

void FMPoly::exportToJson(const std::string &fileName) const {
  std::ofstream outputFile(fileName + ".json");
  outputFile << "{\n\t\"terms\":[\n";

  // Get number of terms in the polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);

  // Collect all terms grouped by x-powers
  std::map<std::vector<int>, std::vector<std::pair<int, std::string>>> xToQTerms;

  for (slong i = 0; i < numTerms; ++i) {
    // Get coefficient of this term
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // Skip zero coefficients (shouldn't happen in FLINT, but be defensive)
    if (fmpz_is_zero(coeff)) {
      fmpz_clear(coeff);
      continue;
    }

    // Get exponent vector for this term
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; ++j) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Extract q-power and x-powers from exponent vector
    int qPower;
    std::vector<int> xPowers;
    getExponentsFromMonomial(exps, qPower, xPowers);

    // Convert coefficient to string to preserve full precision
    char *coeffStr = fmpz_get_str(NULL, 10, coeff);
    xToQTerms[xPowers].emplace_back(qPower, std::string(coeffStr));
    flint_free(coeffStr);

    // Clean up
    fmpz_clear(coeff);
    for (int j = 0; j <= numXVariables; ++j) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  // Output terms in sorted order
  bool firstTerm = true;
  for (const auto &entry : xToQTerms) {
    const std::vector<int> &xPowers = entry.first;
    const auto &qTerms = entry.second;

    if (!qTerms.empty()) {
      if (!firstTerm) {
        outputFile << ",\n";
      }
      firstTerm = false;

      outputFile << "\t\t{\"x\": [";
      for (size_t k = 0; k < xPowers.size(); k++) {
        outputFile << xPowers[k];
        if (k < xPowers.size() - 1)
          outputFile << ",";
      }
      outputFile << "], \"q_terms\": [";

      for (size_t i = 0; i < qTerms.size(); i++) {
        // Output coefficient as string to preserve precision
        outputFile << "{\"q\": " << qTerms[i].first
                   << ", \"c\": \"" << qTerms[i].second << "\"}";
        if (i < qTerms.size() - 1)
          outputFile << ", ";
      }
      outputFile << "]}";
    }
  }

  outputFile << "\n\t],\n";
  outputFile << "\t\"metadata\": {\n";
  outputFile << "\t\t\"num_x_variables\": " << numXVariables << ",\n";
  outputFile << "\t\t\"max_x_degrees\": [";
  for (int i = 0; i < numXVariables; i++) {
    outputFile << maxXDegrees[i];
    if (i < numXVariables - 1)
      outputFile << ",";
  }
  outputFile << "],\n";
  outputFile << "\t\t\"storage_type\": \"flint\"\n";
  outputFile << "\t}\n}";
  outputFile.close();
}

void FMPoly::exportToJson(const std::string &fileName,
                          const std::vector<double> &overall_x_powers,
                          double overall_q_power) const {
  std::ofstream outputFile(fileName + ".json");
  outputFile << "{\n\t\"terms\":[\n";

  slong numTerms = fmpz_mpoly_length(poly, ctx);

  std::map<std::vector<int>, std::vector<std::pair<int, std::string>>> xToQTerms;

  for (slong i = 0; i < numTerms; ++i) {
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    if (fmpz_is_zero(coeff)) {
      fmpz_clear(coeff);
      continue;
    }

    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; ++j) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    int qPower;
    std::vector<int> xPowers;
    getExponentsFromMonomial(exps, qPower, xPowers);

    char *coeffStr = fmpz_get_str(NULL, 10, coeff);
    xToQTerms[xPowers].emplace_back(qPower, std::string(coeffStr));
    flint_free(coeffStr);

    fmpz_clear(coeff);
    for (int j = 0; j <= numXVariables; ++j) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  bool firstTerm = true;
  for (const auto &entry : xToQTerms) {
    const std::vector<int> &xPowers = entry.first;
    const auto &qTerms = entry.second;

    if (!qTerms.empty()) {
      if (!firstTerm)
        outputFile << ",\n";
      firstTerm = false;

      outputFile << "\t\t{\"x\": [";
      for (size_t k = 0; k < xPowers.size(); k++) {
        outputFile << xPowers[k];
        if (k < xPowers.size() - 1)
          outputFile << ",";
      }
      outputFile << "], \"q_terms\": [";

      for (size_t i = 0; i < qTerms.size(); i++) {
        outputFile << "{\"q\": " << qTerms[i].first
                   << ", \"c\": \"" << qTerms[i].second << "\"}";
        if (i < qTerms.size() - 1)
          outputFile << ", ";
      }
      outputFile << "]}";
    }
  }

  outputFile << "\n\t],\n";
  outputFile << "\t\"metadata\": {\n";
  outputFile << "\t\t\"num_x_variables\": " << numXVariables << ",\n";
  outputFile << "\t\t\"max_x_degrees\": [";
  for (int i = 0; i < numXVariables; i++) {
    outputFile << maxXDegrees[i];
    if (i < numXVariables - 1)
      outputFile << ",";
  }
  outputFile << "],\n";
  outputFile << "\t\t\"overall_x_powers\": [";
  for (size_t i = 0; i < overall_x_powers.size(); i++) {
    outputFile << overall_x_powers[i];
    if (i < overall_x_powers.size() - 1)
      outputFile << ",";
  }
  outputFile << "],\n";
  outputFile << "\t\t\"overall_q_power\": " << overall_q_power << ",\n";
  outputFile << "\t\t\"storage_type\": \"flint\"\n";
  outputFile << "\t}\n}";
  outputFile.close();
}

void FMPoly::print(int maxTerms) const {
  std::cout << "FMPoly P(q";
  for (int i = 0; i < numXVariables; i++) {
    std::cout << ", x" << (i + 1);
  }
  std::cout << "):\n";

  if (isZero()) {
    std::cout << "0\n";
    return;
  }

  // Get number of terms in the polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);
  if (numTerms == 0) {
    std::cout << "0\n";
    return;
  }

  // Limit the number of terms to display
  slong termsToShow =
      (maxTerms > 0 && maxTerms < numTerms) ? maxTerms : numTerms;

  bool first = true;

  // Iterate through terms in the FLINT polynomial
  for (slong i = 0; i < termsToShow; i++) {
    // Get coefficient of this term
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // Get exponent vector for this term
    fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
    fmpz **exp_ptrs =
        (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_init(&exps[j]);
      exp_ptrs[j] = &exps[j];
    }

    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Extract q-power and x-powers from exponent vector
    int qPower;
    std::vector<int> xPowers;
    getExponentsFromMonomial(exps, qPower, xPowers);

    if (!fmpz_is_zero(coeff)) {
      // Determine sign and get absolute value
      int sgn = fmpz_sgn(coeff);
      fmpz_t absCoeff;
      fmpz_init(absCoeff);
      fmpz_abs(absCoeff, coeff);

      // Print sign
      if (!first) {
        std::cout << (sgn > 0 ? " + " : " - ");
      } else if (sgn < 0) {
        std::cout << "-";
      }
      first = false;

      // Check if this is a constant term
      bool isConstantTerm = (qPower == 0);
      for (int j = 0; j < numXVariables; j++) {
        if (xPowers[j] != 0) {
          isConstantTerm = false;
          break;
        }
      }

      // Print coefficient if not 1 or if it's a constant term
      bool coeffIsOne = fmpz_is_one(absCoeff);
      if (!coeffIsOne || isConstantTerm) {
        // Convert to string to handle big coefficients
        char *coeffStr = fmpz_get_str(NULL, 10, absCoeff);
        std::cout << coeffStr;
        flint_free(coeffStr);
      }

      // Print q term
      if (qPower != 0) {
        if (!coeffIsOne || isConstantTerm)
          std::cout << "*";
        std::cout << "q";
        if (qPower != 1) {
          std::cout << "^" << qPower;
        }
      }

      // Print x terms
      for (int j = 0; j < numXVariables; j++) {
        if (xPowers[j] != 0) {
          if (!coeffIsOne || isConstantTerm || qPower != 0)
            std::cout << "*";
          std::cout << "x" << (j + 1);
          if (xPowers[j] != 1) {
            std::cout << "^" << xPowers[j];
          }
        }
      }

      fmpz_clear(absCoeff);
    }

    // Cleanup
    fmpz_clear(coeff);
    for (int j = 0; j <= numXVariables; j++) {
      fmpz_clear(&exps[j]);
    }
    flint_free(exp_ptrs);
    flint_free(exps);
  }

  if (termsToShow < numTerms) {
    std::cout << " + ... (" << (numTerms - termsToShow) << " more terms)";
  }

  std::cout << "\n";
}

std::vector<std::pair<std::vector<int>, QPolynomial>>
FMPoly::getCoefficients() const {
  // Map x-powers -> QPolynomial in q
  std::map<std::vector<int>, QPolynomial> xToQPoly;

  // Number of terms in the FLINT multivariate polynomial
  slong numTerms = fmpz_mpoly_length(poly, ctx);

  // Allocate buffers ONCE outside the loop
  fmpz_t coeff;
  fmpz_init(coeff);
  fmpz *exps = (fmpz *)flint_malloc((numXVariables + 1) * sizeof(fmpz));
  fmpz **exp_ptrs = (fmpz **)flint_malloc((numXVariables + 1) * sizeof(fmpz *));
  for (int j = 0; j <= numXVariables; ++j) {
    fmpz_init(&exps[j]);
    exp_ptrs[j] = &exps[j];
  }

  // Loop through terms, reusing buffers
  for (slong i = 0; i < numTerms; ++i) {
    // Get coefficient of this term
    fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

    // FLINT shouldn't store zero terms, but be defensive
    if (fmpz_is_zero(coeff)) {
      continue;
    }

    // Get exponent vector for this term: [q, x1, ..., xn]
    fmpz_mpoly_get_term_exp_fmpz(exp_ptrs, poly, i, ctx);

    // Convert FLINT exponents (with ground offsets) back to logical powers
    int qPower;
    std::vector<int> xPowers;
    getExponentsFromMonomial(exps, qPower, xPowers);

    // Accumulate into the QPolynomial for this xPowers using fmpz directly
    QPolynomial &qp = xToQPoly[xPowers]; // default-constructs if not present
    qp.addToCoefficientFmpz(qPower, coeff);
  }

  // Clean up buffers ONCE after the loop
  fmpz_clear(coeff);
  for (int j = 0; j <= numXVariables; ++j) {
    fmpz_clear(&exps[j]);
  }
  flint_free(exp_ptrs);
  flint_free(exps);

  // Convert the map into the sparse vector of (xPowers, QPolynomial),
  // skipping any zero polynomials (to emulate BMPoly's behavior).
  std::vector<std::pair<std::vector<int>, QPolynomial>> result;
  result.reserve(xToQPoly.size());

  for (auto &entry : xToQPoly) {
    const auto &xPowers = entry.first;
    const auto &qp = entry.second;

    if (!qp.isZero()) {
      result.emplace_back(xPowers, qp);
    }
  }

  return result;
}

void FMPoly::checkCompatibility(const FMPoly &other) const {
  if (numXVariables != other.numXVariables) {
    throw std::invalid_argument(
        "Polynomials must have the same number of x variables");
  }
}

FMPoly &FMPoly::operator+=(const FMPoly &other) {
  checkCompatibility(other);

  if (other.isZero()) {
    return *this;
  }

  // Reconcile ground powers to the componentwise minimum, re-expressing
  // stored exponents where needed, then use FLINT's native merge-add.
  const int n = numXVariables + 1;
  std::vector<int> target(n);
  for (int i = 0; i < n; i++) {
    target[i] = std::min(allGroundPowers[i], other.allGroundPowers[i]);
  }

  if (target != allGroundPowers) {
    std::vector<ulong> shift(n);
    for (int i = 0; i < n; i++) {
      shift[i] = static_cast<ulong>(allGroundPowers[i] - target[i]);
    }
    fmpz_mpoly_t shifted;
    fmpz_mpoly_init(shifted, ctx);
    buildShiftedCopy(shifted, poly, shift, ctx);
    fmpz_mpoly_swap(poly, shifted, ctx);
    fmpz_mpoly_clear(shifted, ctx);
  }

  if (target == other.allGroundPowers) {
    fmpz_mpoly_add(poly, poly, other.poly, ctx);
  } else {
    std::vector<ulong> shift(n);
    for (int i = 0; i < n; i++) {
      shift[i] = static_cast<ulong>(other.allGroundPowers[i] - target[i]);
    }
    fmpz_mpoly_t shifted;
    fmpz_mpoly_init(shifted, ctx);
    buildShiftedCopy(shifted, other.poly, shift, ctx);
    fmpz_mpoly_add(poly, poly, shifted, ctx);
    fmpz_mpoly_clear(shifted, ctx);
  }

  allGroundPowers = std::move(target);
  return *this;
}

FMPoly &FMPoly::operator*=(const FMPoly &other) {
  checkCompatibility(other);

  // Quick zero checks: 0 * anything = 0
  if (isZero() || other.isZero()) {
    clear();
    // For a zero polynomial the ground powers don't really matter;
    // leaving them as-is is fine.
    return *this;
  }

  // Multiply underlying FLINT polynomials using the existing context.
  // We use a temporary to avoid aliasing issues (e.g. self-multiplication).
  fmpz_mpoly_t prod;
  fmpz_mpoly_init(prod, ctx);
  fmpz_mpoly_mul(prod, poly, other.poly, ctx);

  // Swap result into this->poly and clean up
  fmpz_mpoly_swap(poly, prod, ctx);
  fmpz_mpoly_clear(prod, ctx);

  // Update ground powers: product exponents are sums of absolute exponents,
  // and stored exponents are sums of stored exponents, so the new ground
  // vector is just the componentwise sum.
  for (int i = 0; i <= numXVariables; ++i) {
    allGroundPowers[i] += other.allGroundPowers[i];
  }

  // maxXDegrees / blockSizes are just "metadata" compatibility; we can
  // leave them unchanged, since FLINT doesn't use them.
  return *this;
}

FMPoly FMPoly::operator*(const QPolynomial &qPoly) const {
  // Result has the same x-structure metadata as *this
  FMPoly result(numXVariables, 0, maxXDegrees);

  // Zero short-circuits
  if (isZero() || qPoly.isZero()) {
    return result; // already zero poly with same shape
  }

  // Build an FMPoly that represents "qPoly in the first variable, no
  // x-dependence"
  FMPoly factor(numXVariables, 0, maxXDegrees);
  std::vector<int> zeroXPowers(numXVariables, 0);
  factor.addQPolynomial(zeroXPowers, qPoly);

  // Now use the fast FLINT-based FMPoly multiplication
  result = *this;
  result *= factor; // this is your fmpz_mpoly_mul-based operator*=

  return result;
}
