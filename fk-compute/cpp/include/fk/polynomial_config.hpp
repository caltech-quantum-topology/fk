#pragma once


/**
 * Polynomial Configuration Header
 *
 * This header allows easy switching between different polynomial implementations.
 * Set the POLYNOMIAL_TYPE macro to choose between:
 * - 0: MultivariablePolynomial (sparse implementation with unordered_map)
 * - 1: FMPoly (FLINT-based dense implementation for performance)
 * - 2: BMPoly (Basic vector-based implementation with negative exponent support)
 *
 * All classes have identical public interfaces, making them interchangeable.
 */

// Configuration: Set to 0, 1, or 2 to choose polynomial implementation
#ifndef POLYNOMIAL_TYPE
#define POLYNOMIAL_TYPE 1  // Default to FMPoly
#endif

#if POLYNOMIAL_TYPE == 0
    #include "fk/multivariable_polynomial.hpp"
    using PolynomialType = MultivariablePolynomial;
    #define POLYNOMIAL_CLASS_NAME "MultivariablePolynomial"
#elif POLYNOMIAL_TYPE == 1
    #include "fk/fmpoly.hpp"
    using PolynomialType = FMPoly;
    #define POLYNOMIAL_CLASS_NAME "FMPoly"
#elif POLYNOMIAL_TYPE == 2
    #include "fk/bmpoly.hpp"
    using PolynomialType = BMPoly;
    #define POLYNOMIAL_CLASS_NAME "BMPoly"
#else
    #error "Invalid POLYNOMIAL_TYPE: must be 0 (MultivariablePolynomial), 1 (FMPoly), or 2 (BMPoly)"
#endif

// Backward compatibility: Support old USE_FMPOLY macro
#ifdef USE_FMPOLY
    #if USE_FMPOLY
        #undef POLYNOMIAL_TYPE
        #define POLYNOMIAL_TYPE 1
    #else
        #undef POLYNOMIAL_TYPE
        #define POLYNOMIAL_TYPE 0
    #endif
#endif
