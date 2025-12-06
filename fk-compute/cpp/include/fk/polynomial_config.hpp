#pragma once


/**
 * Polynomial Configuration Header
 *
 * This header allows easy switching between different polynomial implementations.
 * Set the POLYNOMIAL_TYPE macro to choose between:
 * - 0: MultivariablePolynomial (sparse implementation with unordered_map)
 * - 1: FMPoly (FLINT-based dense implementation for performance)
 * - 2: BMPoly (Basic vector-based implementation with negative exponent support)
 * - 3: HMPoly (Hash map polynomial implementation)
 * - 4: ZMPoly (Sparse implementation with arbitrary precision integer coefficients)
 *
 * All classes have identical public interfaces, making them interchangeable.
 */

// Configuration: Set to 0, 1, 2, 3, or 4 to choose polynomial implementation
#ifndef POLYNOMIAL_TYPE
#define POLYNOMIAL_TYPE 0  // Default to MultivariablePolynomial
#endif

#if POLYNOMIAL_TYPE == 0
    #include "fk/multivariable_polynomial.hpp"
    using PolynomialType = MultivariablePolynomial;
    using QPolynomialType = bilvector<int>;
    #define POLYNOMIAL_CLASS_NAME "MultivariablePolynomial"
#elif POLYNOMIAL_TYPE == 1
    #include "fk/fmpoly.hpp"
    using PolynomialType = FMPoly;
    using QPolynomialType = QPolynomial;
    #define POLYNOMIAL_CLASS_NAME "FMPoly"
#elif POLYNOMIAL_TYPE == 2
    #include "fk/bmpoly.hpp"
    using PolynomialType = BMPoly;
    using QPolynomialType = bilvector<int>;
    #define POLYNOMIAL_CLASS_NAME "BMPoly"
#elif POLYNOMIAL_TYPE == 3
    #include "fk/hmpoly.hpp"
    using PolynomialType = HMPoly;
    using QPolynomialType = bilvector<int>;
    #define POLYNOMIAL_CLASS_NAME "HMPoly"
#elif POLYNOMIAL_TYPE == 4
    #include "fk/zmpoly.hpp"
    #include "fk/qpolynomial.hpp"
    using PolynomialType = ZMPoly;
    using QPolynomialType = QPolynomial;
    #define POLYNOMIAL_CLASS_NAME "ZMPoly"
#else
    #error "Invalid POLYNOMIAL_TYPE: must be 0 (MultivariablePolynomial), 1 (FMPoly), 2 (BMPoly), 3 (HMPoly), or 4 (ZMPoly)"
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
