=== Performance Comparison: FMPoly vs MultivariablePolynomial vs BMPoly ===
Testing multiplication performance with various polynomial sizes

=== Small polynomials (2 vars, 10 terms) ===
Variables: 2, Terms: 10, Max Q Power: 5, Max X Power: 3, Iterations: 100

Testing FMPoly (FLINT-based):
  Sample result - is zero: no
Testing MultivariablePolynomial (sparse hash-based):
  Sample result - is zero: no
Testing BMPoly (vector-based):
  Sample result - is zero: no

Results (times in milliseconds):
           Implementation     Average      Median         Min         Max
-------------------------------------------------------------------------
           FMPoly (FLINT)       0.003       0.003       0.002       0.015
  MultivariablePolynomial       0.006       0.006       0.004       0.017
          BMPoly (vector)       0.008       0.008       0.006       0.014

Relative Performance:
                              Comparison        Speedup
-------------------------------------------------------
             MultivariablePoly / FMPoly:           1.94x
                        BMPoly / FMPoly:           2.86x
             MultivariablePoly / BMPoly:           0.68x

Ranking (fastest to slowest):
  1. FMPoly (fastest)
  2. MultivariablePolynomial (0.51x slower)
  3. BMPoly (0.35x slower)

=== Medium polynomials (2 vars, 50 terms) ===
Variables: 2, Terms: 50, Max Q Power: 10, Max X Power: 5, Iterations: 50

Testing FMPoly (FLINT-based):
  Sample result - is zero: no
Testing MultivariablePolynomial (sparse hash-based):
  Sample result - is zero: no
Testing BMPoly (vector-based):
  Sample result - is zero: no

Results (times in milliseconds):
           Implementation     Average      Median         Min         Max
-------------------------------------------------------------------------
           FMPoly (FLINT)       0.049       0.061       0.022       0.076
  MultivariablePolynomial       0.096       0.097       0.086       0.109
          BMPoly (vector)       0.090       0.088       0.077       0.112

Relative Performance:
                              Comparison        Speedup
-------------------------------------------------------
             MultivariablePoly / FMPoly:           1.97x
                        BMPoly / FMPoly:           1.86x
             MultivariablePoly / BMPoly:           1.06x

Ranking (fastest to slowest):
  1. FMPoly (fastest)
  2. BMPoly (0.54x slower)
  3. MultivariablePolynomial (0.51x slower)

=== 3D medium polynomials (3 vars, 25 terms) ===
Variables: 3, Terms: 25, Max Q Power: 8, Max X Power: 4, Iterations: 50

Testing FMPoly (FLINT-based):
  Sample result - is zero: no
Testing MultivariablePolynomial (sparse hash-based):
  Sample result - is zero: no
Testing BMPoly (vector-based):
  Sample result - is zero: no

Results (times in milliseconds):
           Implementation     Average      Median         Min         Max
-------------------------------------------------------------------------
           FMPoly (FLINT)       0.016       0.016       0.014       0.023
  MultivariablePolynomial       0.068       0.068       0.061       0.077
          BMPoly (vector)       0.079       0.077       0.071       0.107

Relative Performance:
                              Comparison        Speedup
-------------------------------------------------------
             MultivariablePoly / FMPoly:           4.23x
                        BMPoly / FMPoly:           4.92x
             MultivariablePoly / BMPoly:           0.86x

Ranking (fastest to slowest):
  1. FMPoly (fastest)
  2. MultivariablePolynomial (0.24x slower)
  3. BMPoly (0.20x slower)

=== Large polynomials (2 vars, 100 terms) ===
Variables: 2, Terms: 100, Max Q Power: 15, Max X Power: 8, Iterations: 20

Testing FMPoly (FLINT-based):
  Sample result - is zero: no
Testing MultivariablePolynomial (sparse hash-based):
  Sample result - is zero: no
Testing BMPoly (vector-based):
  Sample result - is zero: no

Results (times in milliseconds):
           Implementation     Average      Median         Min         Max
-------------------------------------------------------------------------
           FMPoly (FLINT)       0.074       0.073       0.072       0.086
  MultivariablePolynomial       0.429       0.428       0.413       0.450
          BMPoly (vector)       0.417       0.416       0.388       0.449

Relative Performance:
                              Comparison        Speedup
-------------------------------------------------------
             MultivariablePoly / FMPoly:           5.77x
                        BMPoly / FMPoly:           5.62x
             MultivariablePoly / BMPoly:           1.03x

Ranking (fastest to slowest):
  1. FMPoly (fastest)
  2. BMPoly (0.18x slower)
  3. MultivariablePolynomial (0.17x slower)

=== 4D polynomials (4 vars, 30 terms) ===
Variables: 4, Terms: 30, Max Q Power: 6, Max X Power: 3, Iterations: 30

Testing FMPoly (FLINT-based):
  Sample result - is zero: no
Testing MultivariablePolynomial (sparse hash-based):
  Sample result - is zero: no
Testing BMPoly (vector-based):
  Sample result - is zero: no

Results (times in milliseconds):
           Implementation     Average      Median         Min         Max
-------------------------------------------------------------------------
           FMPoly (FLINT)       0.024       0.024       0.022       0.038
  MultivariablePolynomial       0.113       0.114       0.094       0.120
          BMPoly (vector)       0.227       0.246       0.166       0.295

Relative Performance:
                              Comparison        Speedup
-------------------------------------------------------
             MultivariablePoly / FMPoly:           4.70x
                        BMPoly / FMPoly:           9.49x
             MultivariablePoly / BMPoly:           0.50x

Ranking (fastest to slowest):
  1. FMPoly (fastest)
  2. MultivariablePolynomial (0.21x slower)
  3. BMPoly (0.11x slower)

=== Very large polynomials (2 vars, 200 terms) ===
Variables: 2, Terms: 200, Max Q Power: 20, Max X Power: 10, Iterations: 10

Testing FMPoly (FLINT-based):
  Sample result - is zero: no
Testing MultivariablePolynomial (sparse hash-based):
  Sample result - is zero: no
Testing BMPoly (vector-based):
  Sample result - is zero: no

Results (times in milliseconds):
           Implementation     Average      Median         Min         Max
-------------------------------------------------------------------------
           FMPoly (FLINT)       0.250       0.260       0.178       0.266
  MultivariablePolynomial       1.643       1.652       1.557       1.758
          BMPoly (vector)       1.552       1.566       1.454       1.660

Relative Performance:
                              Comparison        Speedup
-------------------------------------------------------
             MultivariablePoly / FMPoly:           6.57x
                        BMPoly / FMPoly:           6.21x
             MultivariablePoly / BMPoly:           1.06x

Ranking (fastest to slowest):
  1. FMPoly (fastest)
  2. BMPoly (0.16x slower)
  3. MultivariablePolynomial (0.15x slower)

=== Summary ===
Performance comparison completed successfully!
Results show relative performance of three polynomial implementations:
  - FMPoly (FLINT-based)
  - MultivariablePolynomial (sparse hash-based)
  - BMPoly (dense vector-based)
