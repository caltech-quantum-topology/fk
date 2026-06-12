#pragma once

/**
 * Polynomial Configuration Header
 *
 * The FK computation is implemented against the FLINT-backed FMPoly /
 * QPolynomial pair. Historical alternative backends (MultivariablePolynomial,
 * BMPoly, HMPoly, ZMPoly) were never compatible with the current fk_main
 * pipeline (they lack the exportToJson overloads and fmpz coefficient
 * methods) and live in _attic/cpp-backup-pre-refactor for reference.
 */

#define POLYNOMIAL_TYPE 1

#include "fk/fmpoly.hpp"
using PolynomialType = FMPoly;
using QPolynomialType = QPolynomial;
#define POLYNOMIAL_CLASS_NAME "FMPoly"
