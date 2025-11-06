#pragma once


/**
 * Polynomial Configuration Header
 *
 * This header allows easy switching between different polynomial implementations.
 * Simply change the USE_FMPOLY macro to switch between:
 * - MultivariablePolynomial (sparse implementation with unordered_map)
 * - FMPoly (FLINT-based dense implementation for performance)
 *
 * Both classes have identical public interfaces, making them interchangeable.
 */

// Configuration: Set to 1 to use FMPoly, 0 to use MultivariablePolynomial
#ifndef USE_FMPOLY
#define USE_FMPOLY 0
#endif

#if USE_FMPOLY
    #include "fk/fmpoly.hpp"
    using PolynomialType = FMPoly;
    #define POLYNOMIAL_CLASS_NAME "FMPoly"
#else
    #include "fk/multivariable_polynomial.hpp"
    using PolynomialType = MultivariablePolynomial;
    #define POLYNOMIAL_CLASS_NAME "MultivariablePolynomial"
#endif
