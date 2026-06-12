#pragma once

#include "fk/polynomial_config.hpp"

// Cached q-algebra building blocks used by the FK crossing factors.

QPolynomialType QBinomialPositive(int upperLimit, int lowerLimit);
QPolynomialType QBinomialNegative(int upperLimit, int lowerLimit);
QPolynomialType QBinomial(int upperLimit, int lowerLimit);

PolynomialType qpochhammer_xq_q(int n, int qpow);
PolynomialType inverse_qpochhammer_xq_q(int n, int qpow, int xMax);
