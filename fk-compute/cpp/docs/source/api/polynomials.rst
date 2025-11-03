Multivariable Polynomials API
============================

.. contents:: Table of Contents
   :local:
   :depth: 3

Overview
--------

The :class:`MultivariablePolynomial` class provides a memory-efficient sparse representation
of polynomials in variables q, x₁, x₂, ..., xₙ with support for negative exponents.

.. class:: feature-check

   Thread-safe for read operations, requires external synchronization for writes

.. class:: feature-check

   Memory efficient: 90%+ reduction vs dense representation

.. class:: feature-check

   Supports arbitrary positive/negative q powers

Class Reference
---------------

MultivariablePolynomial
~~~~~~~~~~~~~~~~~~~~~~~

.. cpp:class:: MultivariablePolynomial

   Represents polynomials of the form:

   .. math::

      P(q, x_1, x_2, \ldots, x_n) = \sum_{i,j} c_{i,j} \times q^j \times x_1^{a_{1,i}} \times x_2^{a_{2,i}} \times \cdots \times x_n^{a_{n,i}}

   **Storage**: Uses sparse representation with :cpp:type:`std::unordered_map` for coefficients.

   **Thread Safety**: :class:`thread-safe` for read operations, requires synchronization for writes.

Constructors
~~~~~~~~~~~~

.. cpp:function:: MultivariablePolynomial::MultivariablePolynomial(int numVariables, int degree = 10, const std::vector<int>& maxDegrees = {})

   Creates a new multivariable polynomial.

   :param numVariables: Number of x variables (x₁, x₂, ..., xₙ)
   :param degree: Default maximum degree hint (advisory only)
   :param maxDegrees: Per-variable degree hints (advisory only)

   **Example**:

   .. code-block:: cpp

      // Polynomial in 3 variables with default settings
      MultivariablePolynomial poly1(3);

      // Polynomial with degree hints
      MultivariablePolynomial poly2(2, 5);

      // Polynomial with per-variable degree hints
      std::vector<int> degrees = {3, 5, 7};
      MultivariablePolynomial poly3(3, 0, degrees);

   **Complexity**: :class:`complexity-o` O(1)

Arithmetic Operations
~~~~~~~~~~~~~~~~~~~~

.. cpp:function:: MultivariablePolynomial MultivariablePolynomial::operator+(const MultivariablePolynomial& other) const

   Polynomial addition.

   :param other: Polynomial to add
   :returns: New polynomial representing the sum
   :throws std::invalid_argument: If polynomials have incompatible dimensions

   **Example**:

   .. code-block:: cpp

      MultivariablePolynomial p1(2), p2(2);
      p1.getCoefficient({1, 0})[0] = 3;  // 3*x₁
      p2.getCoefficient({0, 1})[0] = 2;  // 2*x₂

      auto sum = p1 + p2;  // 3*x₁ + 2*x₂

   **Complexity**: :class:`complexity-o` O(|terms₁| + |terms₂|)

.. cpp:function:: MultivariablePolynomial MultivariablePolynomial::operator-(const MultivariablePolynomial& other) const

   Polynomial subtraction.

   :param other: Polynomial to subtract
   :returns: New polynomial representing the difference

   **Complexity**: :class:`complexity-o` O(|terms₁| + |terms₂|)

.. cpp:function:: MultivariablePolynomial MultivariablePolynomial::operator*(const MultivariablePolynomial& other) const

   Polynomial multiplication.

   :param other: Polynomial to multiply
   :returns: New polynomial representing the product

   **Example**:

   .. code-block:: cpp

      MultivariablePolynomial p1(2), p2(2);
      p1.getCoefficient({1, 0})[0] = 2;  // 2*x₁
      p2.getCoefficient({0, 1})[1] = 3;  // 3*q*x₂

      auto product = p1 * p2;  // 6*q*x₁*x₂
      // Result: coefficient 6 at exponents {1,1} for q^1

   **Complexity**: :class:`complexity-o` O(|terms₁| × |terms₂|)

In-Place Operations
~~~~~~~~~~~~~~~~~~

.. cpp:function:: MultivariablePolynomial& MultivariablePolynomial::operator+=(const MultivariablePolynomial& other)

   In-place addition (more memory efficient).

   :param other: Polynomial to add
   :returns: Reference to this polynomial

   **Performance Note**: Preferred over ``operator+`` for performance-critical code.

.. cpp:function:: MultivariablePolynomial& MultivariablePolynomial::operator-=(const MultivariablePolynomial& other)

   In-place subtraction.

.. cpp:function:: MultivariablePolynomial& MultivariablePolynomial::operator*=(const MultivariablePolynomial& other)

   In-place multiplication.

Coefficient Access
~~~~~~~~~~~~~~~~~

.. cpp:function:: bilvector<int>& MultivariablePolynomial::getCoefficient(const std::vector<int>& xExponents)

   Access coefficient bilvector for given x-variable exponents.

   :param xExponents: Vector of exponents for x₁, x₂, ..., xₙ
   :returns: Reference to bilvector containing q-polynomial coefficients
   :throws std::invalid_argument: If exponent vector size doesn't match numVariables

   **Example**:

   .. code-block:: cpp

      MultivariablePolynomial poly(3);

      // Set coefficient of q^2 * x₁^1 * x₂^0 * x₃^2
      poly.getCoefficient({1, 0, 2})[2] = 5;

      // Set coefficient of q^(-1) * x₁^0 * x₂^1 * x₃^0
      poly.getCoefficient({0, 1, 0})[-1] = 3;

   **Thread Safety**: :class:`thread-unsafe` - requires external synchronization

.. cpp:function:: const bilvector<int>& MultivariablePolynomial::getCoefficient(const std::vector<int>& xExponents) const

   Const access to coefficient bilvector.

   **Thread Safety**: :class:`thread-safe`

Backward Compatibility
~~~~~~~~~~~~~~~~~~~~~

.. cpp:function:: std::vector<std::vector<bilvector<int>>> MultivariablePolynomial::getCoefficients() const

   Convert to dense format for backward compatibility.

   :returns: Dense coefficient array indexed by [x₁][x₂]...[xₙ]

   **Performance Warning**: This operation is expensive for sparse polynomials.
   Only use when interfacing with legacy code.

   **Example**:

   .. code-block:: cpp

      auto dense = poly.getCoefficients();
      // Access coefficient of x₁^i * x₂^j * ... * xₙ^k:
      // dense[i][j]...[k] gives bilvector for q-polynomial

.. cpp:function:: void MultivariablePolynomial::syncFromDenseVector(const std::vector<std::vector<bilvector<int>>>& dense)

   Load coefficients from dense format.

   :param dense: Dense coefficient array
   :throws std::invalid_argument: If dimensions don't match

Utility Methods
~~~~~~~~~~~~~~

.. cpp:function:: bool MultivariablePolynomial::isEmpty() const

   Check if polynomial has no non-zero terms.

   **Complexity**: :class:`complexity-o` O(1)

.. cpp:function:: size_t MultivariablePolynomial::getTermCount() const

   Get number of non-zero terms.

   **Complexity**: :class:`complexity-o` O(1)

.. cpp:function:: void MultivariablePolynomial::clear()

   Remove all terms, making polynomial zero.

   **Complexity**: :class:`complexity-o` O(|terms|)

.. cpp:function:: void MultivariablePolynomial::pruneZeros()

   Remove terms with zero coefficients to save memory.

   **Example**:

   .. code-block:: cpp

      // After many operations, some terms might become zero
      poly += other_poly;
      poly -= other_poly;  // May leave zero terms

      poly.pruneZeros();   // Clean up memory

   **Complexity**: :class:`complexity-o` O(|terms|)

Export Functions
~~~~~~~~~~~~~~~

.. cpp:function:: void MultivariablePolynomial::exportToJson(const std::string& filename) const

   Export polynomial to JSON format.

   :param filename: Output filename (without .json extension)

   **JSON Format**:

   .. code-block:: json

      {
        "metadata": {
          "numVariables": 3,
          "maxDegrees": [10, 10, 10]
        },
        "terms": [
          {
            "xExponents": [1, 0, 2],
            "qCoefficients": {
              "2": 5,
              "-1": 3
            }
          }
        ]
      }

Helper Classes
--------------

VectorHash
~~~~~~~~~~

.. cpp:struct:: VectorHash

   Custom hash function for :cpp:type:`std::vector<int>` keys in sparse storage.

   .. cpp:function:: std::size_t VectorHash::operator()(const std::vector<int>& v) const

      Compute hash of integer vector.

      **Hash Quality**: Designed for mathematical vectors with small integers.
      Handles negative values correctly.

Usage Patterns
--------------

Pattern 1: Building Polynomials Term by Term
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   MultivariablePolynomial buildPolynomial() {
       MultivariablePolynomial poly(3);

       // Add terms: 2*q^1*x₁ + 3*q^(-1)*x₂ + q^0*x₃^2
       poly.getCoefficient({1, 0, 0})[1] = 2;   // 2*q*x₁
       poly.getCoefficient({0, 1, 0})[-1] = 3;  // 3*q⁻¹*x₂
       poly.getCoefficient({0, 0, 2})[0] = 1;   // x₃²

       return poly;
   }

Pattern 2: Polynomial Arithmetic Chains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   MultivariablePolynomial compute(const std::vector<MultivariablePolynomial>& polys) {
       MultivariablePolynomial result(polys[0].getNumVariables());

       // Efficient in-place operations
       for (const auto& p : polys) {
           result += p;  // In-place addition
       }

       result.pruneZeros();  // Clean up
       return result;
   }

Pattern 3: Thread-Safe Reading
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   class PolynomialAnalyzer {
       const MultivariablePolynomial& poly_;
   public:
       PolynomialAnalyzer(const MultivariablePolynomial& p) : poly_(p) {}

       // Thread-safe operations
       size_t countTerms() const { return poly_.getTermCount(); }
       bool hasTermAt(const std::vector<int>& exp) const {
           return !poly_.getCoefficient(exp).isEmpty();
       }
   };

Performance Considerations
-------------------------

**Memory Usage**

- Sparse storage uses ~8-24 bytes per non-zero term
- Dense storage would use ``degree^numVariables * sizeof(bilvector<int>)``
- For polynomials with < 10% non-zero terms, sparse is more efficient

**Operation Costs**

.. table:: Operation Complexity
   :class: performance-table

   +------------------+------------------+-------------------+
   | Operation        | Time Complexity  | Memory Usage      |
   +==================+==================+===================+
   | Addition         | O(|t₁| + |t₂|)   | O(|t₁| + |t₂|)   |
   +------------------+------------------+-------------------+
   | Multiplication   | O(|t₁| × |t₂|)   | O(|t₁| × |t₂|)   |
   +------------------+------------------+-------------------+
   | Coefficient      | O(1) average     | O(1)              |
   | Access           | O(log |terms|)   |                   |
   |                  | worst case       |                   |
   +------------------+------------------+-------------------+

**Best Practices**

1. **Use in-place operations** (``+=``, ``-=``, ``*=``) when possible
2. **Call pruneZeros()** periodically after many operations
3. **Pre-allocate expected capacity** for predictable workloads
4. **Use const references** to enable thread-safe sharing

Thread Safety Details
---------------------

.. class:: warning-box

   **Thread Safety Model**: The polynomial class follows a multiple-reader/single-writer pattern.
   Multiple threads can safely call const methods simultaneously, but write operations
   require external synchronization.

**Safe Operations** (multiple threads):

- ``getCoefficient()`` const version
- ``isEmpty()``, ``getTermCount()``
- Copy constructor, ``operator+``, ``operator*``

**Unsafe Operations** (require synchronization):

- ``getCoefficient()`` non-const version
- ``operator+=``, ``operator-=``, ``operator*=``
- ``clear()``, ``pruneZeros()``

**Example Thread-Safe Usage**:

.. code-block:: cpp

   class ThreadSafePolynomialProcessor {
       mutable std::shared_mutex mutex_;
       MultivariablePolynomial poly_;

   public:
       // Multiple readers can call this simultaneously
       size_t getTermCount() const {
           std::shared_lock lock(mutex_);
           return poly_.getTermCount();
       }

       // Only one writer at a time
       void addTerm(const std::vector<int>& exp, int qPower, int coeff) {
           std::unique_lock lock(mutex_);
           poly_.getCoefficient(exp)[qPower] += coeff;
       }
   };

Error Handling
--------------

The polynomial class uses exceptions for error reporting:

.. cpp:exception:: std::invalid_argument

   Thrown when:
   - Exponent vector size doesn't match number of variables
   - Attempting arithmetic on incompatible polynomials
   - Invalid parameters to constructors

.. cpp:exception:: std::bad_alloc

   Thrown when:
   - System runs out of memory during large operations
   - Attempting to create polynomials too large for available memory

**Example Error Handling**:

.. code-block:: cpp

   try {
       MultivariablePolynomial poly(3);
       poly.getCoefficient({1, 2})[0] = 5;  // Error: wrong dimensions
   } catch (const std::invalid_argument& e) {
       std::cerr << "Dimension error: " << e.what() << std::endl;
   } catch (const std::bad_alloc& e) {
       std::cerr << "Memory allocation failed: " << e.what() << std::endl;
   }