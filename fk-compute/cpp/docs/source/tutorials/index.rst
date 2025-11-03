Tutorials
=========

.. contents:: Table of Contents
   :local:
   :depth: 2

This section provides step-by-step tutorials for learning the FK Computation library.

Tutorial Series
---------------

Tutorial 1: Your First Polynomial
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Goal**: Learn the basics of creating and manipulating polynomials.

**Prerequisites**: Basic C++ knowledge, library installed.

**Time**: 15 minutes

.. toctree::
   :maxdepth: 1

   tutorial_01_first_polynomial

**What you'll learn:**

- Creating multivariable polynomials
- Setting coefficients for different terms
- Basic arithmetic operations
- Understanding sparse vs dense storage

Tutorial 2: Constraint Systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Goal**: Solve your first constraint satisfaction problem.

**Prerequisites**: Tutorial 1 completed.

**Time**: 20 minutes

.. toctree::
   :maxdepth: 1

   tutorial_02_constraints

**What you'll learn:**

- Setting up linear inequality constraints
- Understanding solution callbacks
- Interpreting results
- Debugging constraint systems

Tutorial 3: Parallel Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Goal**: Speed up computations with parallel processing.

**Prerequisites**: Tutorial 2 completed.

**Time**: 25 minutes

.. toctree::
   :maxdepth: 1

   tutorial_03_parallel

**What you'll learn:**

- When to use parallel processing
- Thread configuration and optimization
- Monitoring parallel performance
- Troubleshooting parallel issues

Tutorial 4: Fault Tolerance
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Goal**: Build robust long-running computations.

**Prerequisites**: Tutorial 3 completed.

**Time**: 30 minutes

.. toctree::
   :maxdepth: 1

   tutorial_04_checkpointing

**What you'll learn:**

- Setting up automatic checkpointing
- Handling interruptions gracefully
- Resuming computations
- Best practices for reliability

Tutorial 5: Advanced Data Structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Goal**: Master the library's specialized data structures.

**Prerequisites**: All previous tutorials.

**Time**: 35 minutes

.. toctree::
   :maxdepth: 1

   tutorial_05_data_structures

**What you'll learn:**

- Using bilvectors for bidirectional indexing
- Binary trees with substring matching
- Custom hash functions
- Memory optimization techniques

Tutorial 6: Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Goal**: Optimize your code for maximum performance.

**Prerequisites**: All previous tutorials.

**Time**: 45 minutes

.. toctree::
   :maxdepth: 1

   tutorial_06_optimization

**What you'll learn:**

- Profiling and benchmarking
- Memory usage optimization
- CPU optimization techniques
- Scaling to large problems

Mini-Tutorials
--------------

Quick Start Guides (5-10 minutes each)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   mini/setup_development_environment
   mini/first_parallel_computation
   mini/export_results_json
   mini/debug_constraint_system
   mini/monitor_memory_usage
   mini/integrate_with_cmake

Specialized Topics (15-20 minutes each)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   specialized/numerical_stability
   specialized/custom_solution_filters
   specialized/batch_processing
   specialized/distributed_computing
   specialized/web_api_integration
   specialized/mathematical_validation

Problem-Solving Workshops
-------------------------

Workshop 1: Small Research Problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Scenario**: Enumerate lattice points in a polytope.

**Skills Applied**: Basic constraint solving, result analysis.

**Duration**: 1 hour

**Structure**:

1. Problem definition and mathematical background
2. Constraint formulation
3. Implementation and testing
4. Result interpretation
5. Extensions and variations

Workshop 2: Medium Computational Challenge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Scenario**: Parameter sweep for mathematical model.

**Skills Applied**: Parallel processing, data management.

**Duration**: 2 hours

**Structure**:

1. Problem scaling analysis
2. Parallel implementation design
3. Performance optimization
4. Result aggregation and analysis
5. Visualization and reporting

Workshop 3: Large-Scale Application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Scenario**: High-performance scientific computation.

**Skills Applied**: All library features, optimization.

**Duration**: 3 hours

**Structure**:

1. Architecture design for scalability
2. Implementation with checkpointing
3. Performance tuning and profiling
4. Reliability testing
5. Production deployment considerations

Learning Paths
--------------

Beginner Path: Basic Usage
~~~~~~~~~~~~~~~~~~~~~~~~~

**Total Time**: 2-3 hours

**Sequence**:

1. Tutorial 1: Your First Polynomial (15 min)
2. Tutorial 2: Constraint Systems (20 min)
3. Mini: Setup Development Environment (10 min)
4. Mini: Debug Constraint System (10 min)
5. Workshop 1: Small Research Problem (1 hour)

**Learning Outcomes**:

- Understand basic library concepts
- Solve simple mathematical problems
- Debug common issues
- Apply knowledge to real problems

Intermediate Path: Performance and Reliability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Total Time**: 4-5 hours

**Prerequisites**: Beginner path completed

**Sequence**:

1. Tutorial 3: Parallel Processing (25 min)
2. Tutorial 4: Fault Tolerance (30 min)
3. Mini: Monitor Memory Usage (10 min)
4. Specialized: Numerical Stability (20 min)
5. Workshop 2: Medium Computational Challenge (2 hours)

**Learning Outcomes**:

- Master parallel processing techniques
- Build fault-tolerant applications
- Optimize performance effectively
- Handle medium-scale problems

Advanced Path: Expert Usage
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Total Time**: 6-8 hours

**Prerequisites**: Intermediate path completed

**Sequence**:

1. Tutorial 5: Advanced Data Structures (35 min)
2. Tutorial 6: Performance Optimization (45 min)
3. Specialized: Distributed Computing (20 min)
4. Specialized: Web API Integration (20 min)
5. Workshop 3: Large-Scale Application (3 hours)

**Learning Outcomes**:

- Expert-level library usage
- Advanced optimization techniques
- Production system development
- Large-scale problem solving

Tutorial Support Materials
--------------------------

Code Templates
~~~~~~~~~~~~~

Each tutorial includes downloadable code templates:

.. code-block:: bash

   # Download tutorial templates
   git clone https://github.com/fk-computation/tutorial-code
   cd tutorial-code

   # Each tutorial has its own directory
   ls tutorials/
   # tutorial_01/  tutorial_02/  tutorial_03/  ...

   # Starter code and solutions provided
   ls tutorials/tutorial_01/
   # starter.cpp  solution.cpp  Makefile  README.md

Data Sets
~~~~~~~~

Sample data sets for tutorials:

- **Small datasets**: Tutorial problems with known solutions
- **Medium datasets**: Realistic problems for learning
- **Large datasets**: Performance testing and optimization
- **Benchmark datasets**: Standardized problems for comparison

Reference Solutions
~~~~~~~~~~~~~~~~~~

Complete reference implementations:

- **Starter code**: Skeleton code to begin each tutorial
- **Intermediate checkpoints**: Partial solutions for complex tutorials
- **Complete solutions**: Full implementations with explanations
- **Alternative approaches**: Different solution strategies

Interactive Tools
~~~~~~~~~~~~~~~~

**Jupyter Notebooks** (via Python bindings):

.. code-block:: bash

   # Install Python bindings
   pip install fk-computation-python

   # Launch interactive tutorials
   jupyter notebook tutorials/interactive/

**Web-based Tutorials**:

- Interactive code editing and execution
- Real-time result visualization
- Progress tracking and hints
- Community discussion forums

Assessment and Certification
---------------------------

Knowledge Checks
~~~~~~~~~~~~~~~

Each tutorial includes:

- **Pre-assessment**: Check prerequisite knowledge
- **Progress checks**: Verify understanding during tutorial
- **Post-assessment**: Confirm learning objectives met

**Example Knowledge Check**:

.. code-block:: text

   Question: When should you use parallel processing?

   A) Always, for any problem size
   B) Only for problems requiring > 100ms sequential time
   C) Never, sequential is always faster
   D) Only when you have 16+ CPU cores

   Correct Answer: B
   Explanation: Parallel processing has overhead costs that only
   pay off for sufficiently large problems.

Practical Exercises
~~~~~~~~~~~~~~~~~~

Hands-on coding exercises with automated checking:

.. code-block:: cpp

   // Exercise: Implement solution filtering
   class SolutionFilter {
   public:
       // TODO: Implement this method
       bool shouldAccept(const std::vector<int>& solution) {
           // Your code here
           return false;
       }
   };

   // Automated test
   TEST_CASE("Solution filtering works correctly") {
       SolutionFilter filter;
       REQUIRE(filter.shouldAccept({1, 2, 3}) == true);
       REQUIRE(filter.shouldAccept({0, 0, 0}) == false);
   }

Certification Levels
~~~~~~~~~~~~~~~~~~~

**Bronze Certification**: Basic Usage

- Complete Beginner path
- Pass knowledge assessments (80%+)
- Submit one working project

**Silver Certification**: Performance & Reliability

- Complete Intermediate path
- Demonstrate performance optimization
- Build fault-tolerant application

**Gold Certification**: Expert Usage

- Complete Advanced path
- Contribute to library documentation or examples
- Solve challenging computational problem

Community Features
-----------------

Discussion Forums
~~~~~~~~~~~~~~~

Tutorial-specific discussion areas:

- **Question & Answer**: Get help with specific problems
- **Code Reviews**: Share code for feedback
- **Project Showcase**: Display completed tutorial projects
- **Study Groups**: Form learning groups with other students

Contribution Opportunities
~~~~~~~~~~~~~~~~~~~~~~~~~

Ways to contribute while learning:

- **Improve tutorials**: Suggest clarifications and improvements
- **Add examples**: Contribute new example problems
- **Test on platforms**: Verify tutorials work on different systems
- **Translate content**: Help make tutorials available in other languages

Mentorship Program
~~~~~~~~~~~~~~~~

Connect with experienced users:

- **Beginner mentoring**: Get paired with intermediate users
- **Expert guidance**: Advanced users can connect with library maintainers
- **Project collaboration**: Work on real problems with experienced developers

Getting Help
-----------

Tutorial Support
~~~~~~~~~~~~~~~

When you need help with tutorials:

1. **Check the FAQ**: Common issues and solutions
2. **Search forums**: Look for similar questions
3. **Ask specific questions**: Provide code, error messages, and context
4. **Join office hours**: Weekly live Q&A sessions

Technical Support
~~~~~~~~~~~~~~~~

For technical issues:

1. **Documentation**: Refer to :doc:`../api/index` for detailed API docs
2. **Examples**: Check :doc:`../examples/index` for working code
3. **Performance**: See :doc:`../performance/index` for optimization help
4. **Bug reports**: Use GitHub issues for suspected bugs

Learning Resources
~~~~~~~~~~~~~~~~~

Additional learning materials:

- **Video tutorials**: Step-by-step video walkthroughs
- **Webinars**: Live presentations on advanced topics
- **Conference talks**: Recordings from mathematical computing conferences
- **Research papers**: Academic papers using the library

Prerequisites and Setup
-----------------------

System Requirements
~~~~~~~~~~~~~~~~~~

**Minimum Requirements**:

- C++17 compatible compiler
- 4GB RAM
- 2 CPU cores

**Recommended for Tutorials**:

- GCC 9+ or Clang 10+
- 8GB RAM
- 4+ CPU cores
- SSD storage

Development Environment
~~~~~~~~~~~~~~~~~~~~~~

**Required Tools**:

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt install build-essential cmake git

   # macOS
   xcode-select --install
   brew install cmake git

   # Verify installation
   g++ --version
   cmake --version

**Optional Tools** (recommended):

.. code-block:: bash

   # Code formatting
   sudo apt install clang-format

   # Debugging and profiling
   sudo apt install gdb valgrind

   # Documentation building
   pip install sphinx sphinx-rtd-theme

Environment Verification
~~~~~~~~~~~~~~~~~~~~~~~

Test your setup before starting tutorials:

.. code-block:: bash

   # Clone tutorial repository
   git clone https://github.com/fk-computation/fk-compute
   cd fk-compute/cpp

   # Build library
   make all

   # Run verification test
   make run-test

   # Expected output
   echo "✓ All tests passed - ready for tutorials!"

Next Steps
----------

Ready to start learning? Choose your path:

**New to the Library?**
   Start with :doc:`tutorial_01_first_polynomial`

**Want Quick Results?**
   Try :doc:`mini/first_parallel_computation`

**Specific Problem to Solve?**
   Browse :doc:`../examples/index` for similar use cases

**Performance Focus?**
   Jump to :doc:`tutorial_06_optimization`

**Need Help Choosing?**
   Take the `learning path assessment <#learning-paths>`_ above.