Parallel Processing API
=======================

.. contents:: Table of Contents
   :local:
   :depth: 3

Overview
--------

The parallel processing framework provides high-performance multi-threaded computation
with work stealing, automatic load balancing, and thread-safe data structures.

.. class:: feature-check

   Automatic CPU core detection and optimal thread allocation

.. class:: feature-check

   Work-stealing queues for dynamic load balancing

.. class:: feature-check

   Thread-safe solution collection with atomic operations

.. class:: feature-check

   Progress monitoring and performance statistics

High-Level Functions
-------------------

parallelPooling
~~~~~~~~~~~~~~

.. cpp:function:: void parallelPooling(std::vector<std::vector<double>> main_inequalities, std::vector<std::vector<double>> supporting_inequalities, const std::function<void(const std::vector<int>&)>& callback, int num_threads = 0)

   Execute solution pool enumeration using parallel processing.

   :param main_inequalities: Primary constraint matrix
   :param supporting_inequalities: Additional constraint matrix
   :param callback: Function called for each solution found
   :param num_threads: Number of threads (0 = auto-detect)

   **Auto-Detection**: When ``num_threads = 0``, uses :cpp:func:`std::thread::hardware_concurrency()`

   **Example**:

   .. code-block:: cpp

      auto callback = [](const std::vector<int>& solution) {
          std::cout << "Solution: ";
          for (int x : solution) std::cout << x << " ";
          std::cout << std::endl;
      };

      std::vector<std::vector<double>> main_ineq = {
          {1, 1, -5},  // x + y ≤ 5
          {1, 0, 0},   // x ≥ 0
          {0, 1, 0}    // y ≥ 0
      };

      // Use all available cores
      parallelPooling(main_ineq, {}, callback);

      // Use specific thread count
      parallelPooling(main_ineq, {}, callback, 8);

   **Thread Safety**: :class:`thread-safe` - callback may be called concurrently

   **Performance**: Linear speedup for problems with >1000 solutions

parallelPoolingWithCheckpoints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. cpp:function:: void parallelPoolingWithCheckpoints(std::vector<std::vector<double>> main_inequalities, std::vector<std::vector<double>> supporting_inequalities, const std::function<void(const std::vector<int>&)>& callback, const std::string& checkpoint_file, bool resume_from_checkpoint = false, size_t checkpoint_interval = 1000, int num_threads = 0)

   Execute parallel computation with automatic checkpointing for fault tolerance.

   :param main_inequalities: Primary constraint matrix
   :param supporting_inequalities: Additional constraint matrix
   :param callback: Function called for each solution found
   :param checkpoint_file: Path to checkpoint file
   :param resume_from_checkpoint: Resume from existing checkpoint if available
   :param checkpoint_interval: Save checkpoint every N solutions
   :param num_threads: Number of threads (0 = auto-detect)

   **Example**:

   .. code-block:: cpp

      #include <csignal>

      // Setup graceful interruption
      std::signal(SIGINT, [](int) {
          std::cout << "Interrupting computation..." << std::endl;
          // Checkpoint will be saved automatically
          std::exit(0);
      });

      auto callback = [](const std::vector<int>& solution) {
          // Process solution
      };

      try {
          parallelPoolingWithCheckpoints(
              main_inequalities,
              supporting_inequalities,
              callback,
              "computation.ckpt",  // checkpoint file
              true,                // resume if exists
              500,                 // checkpoint every 500 solutions
              16                   // use 16 threads
          );
      } catch (const std::exception& e) {
          std::cout << "Computation saved to checkpoint" << std::endl;
      }

   **Fault Tolerance**: Computation can be resumed from any checkpoint

   **File Format**: Binary format optimized for size and speed

Core Infrastructure Classes
---------------------------

WorkQueue
~~~~~~~~~

.. cpp:class:: WorkQueue

   Thread-safe work distribution system with work stealing capabilities.

   **Design**: Lock-free queue using atomic operations for high performance.

   .. cpp:function:: void WorkQueue::addWork(const WorkItem& item)

      Add work item to the queue.

      **Thread Safety**: :class:`thread-safe`

   .. cpp:function:: bool WorkQueue::getWork(WorkItem& item)

      Attempt to retrieve work item from queue.

      :returns: true if work was available, false if queue is empty

      **Thread Safety**: :class:`thread-safe`

   .. cpp:function:: void WorkQueue::markFinished()

      Signal that no more work will be added.

   .. cpp:function:: bool WorkQueue::isFinished() const

      Check if queue is finished and empty.

   .. cpp:function:: size_t WorkQueue::getQueueSize() const

      Get current number of queued items.

      **Performance**: :class:`complexity-o` O(1) with atomic operations

   **Example Usage**:

   .. code-block:: cpp

      WorkQueue queue;

      // Producer thread
      std::thread producer([&queue]() {
          for (int i = 0; i < 1000; ++i) {
              WorkItem item = createWorkItem(i);
              queue.addWork(item);
          }
          queue.markFinished();
      });

      // Worker thread
      std::thread worker([&queue]() {
          WorkItem item;
          while (queue.getWork(item) || !queue.isFinished()) {
              processWorkItem(item);
          }
      });

ThreadSafeSolutionCollector
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. cpp:class:: ThreadSafeSolutionCollector

   Collects solutions from multiple threads with thread-safe operations.

   .. cpp:function:: ThreadSafeSolutionCollector::ThreadSafeSolutionCollector(const std::function<void(const std::vector<int>&)>& func)

      Constructor with solution processing function.

      :param func: Callback function for each solution

   .. cpp:function:: void ThreadSafeSolutionCollector::addSolution(const std::vector<int>& solution)

      Add solution from worker thread.

      **Thread Safety**: :class:`thread-safe` with internal synchronization

   .. cpp:function:: size_t ThreadSafeSolutionCollector::getSolutionCount() const

      Get total number of solutions collected.

      **Atomicity**: Uses atomic counter for lock-free access

   **Example**:

   .. code-block:: cpp

      auto processor = [](const std::vector<int>& sol) {
          std::cout << "Solution: ";
          for (int x : sol) std::cout << x << " ";
          std::cout << std::endl;
      };

      ThreadSafeSolutionCollector collector(processor);

      // Multiple worker threads can safely call:
      collector.addSolution({1, 2, 3});
      collector.addSolution({4, 5, 6});

      std::cout << "Total: " << collector.getSolutionCount() << std::endl;

ThreadPool
~~~~~~~~~~

.. cpp:class:: ThreadPool

   Manages worker thread lifecycle and coordination.

   .. cpp:function:: ThreadPool::ThreadPool(int num_threads)

      Create thread pool with specified number of threads.

      :param num_threads: Number of worker threads to create

   .. cpp:function:: ThreadPool::~ThreadPool()

      Destructor - automatically stops all threads and cleans up resources.

   .. cpp:function:: void ThreadPool::start()

      Start all worker threads.

   .. cpp:function:: void ThreadPool::stop()

      Stop all worker threads gracefully.

   .. cpp:function:: bool ThreadPool::isRunning() const

      Check if thread pool is currently running.

   .. cpp:function:: int ThreadPool::getThreadCount() const

      Get number of threads in pool.

   **Resource Management**: RAII-compliant with automatic cleanup

ParallelWorker
~~~~~~~~~~~~~

.. cpp:class:: ParallelWorker

   Individual worker thread implementation.

   **Design**: Each worker runs independently, stealing work when local queue is empty.

   **Work Stealing**: Automatically balances load across threads for optimal performance.

Work Distribution Strategies
---------------------------

Automatic Load Balancing
~~~~~~~~~~~~~~~~~~~~~~~~

The framework uses several strategies to ensure optimal work distribution:

1. **Initial Partitioning**: Problem space divided roughly equally among threads
2. **Work Stealing**: Idle threads steal work from busy threads
3. **Dynamic Subdivision**: Large work items automatically subdivided
4. **Queue Monitoring**: Real-time adjustment based on queue sizes

**Algorithm**:

.. code-block:: text

   for each worker thread:
       while computation not finished:
           if local_queue not empty:
               process local work
           else if global_queue not empty:
               steal work from global queue
           else:
               try to steal from other workers
               if no work found:
                   brief sleep before retry

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

**Thread Count Selection**:

.. code-block:: cpp

   int selectOptimalThreads(size_t expected_solutions) {
       int hw_threads = std::thread::hardware_concurrency();

       if (expected_solutions < 100) {
           return 1;  // Sequential for small problems
       } else if (expected_solutions < 10000) {
           return std::min(4, hw_threads);  // Limited parallelism
       } else {
           return hw_threads;  // Full parallelism
       }
   }

**Memory Considerations**:

.. code-block:: cpp

   // Each thread needs stack space and local data
   size_t memory_per_thread = 8 * 1024 * 1024;  // 8MB estimate
   size_t available_memory = getAvailableMemory();
   int max_threads = available_memory / memory_per_thread;

   int threads = std::min(max_threads, hw_threads);

Progress Monitoring
------------------

Real-Time Statistics
~~~~~~~~~~~~~~~~~~~

.. cpp:class:: ParallelProgressMonitor

   Provides real-time progress information for parallel computations.

   **Metrics Tracked**:
   - Solutions found per second
   - Work queue sizes
   - Thread utilization
   - Memory usage
   - Estimated completion time

   **Example Usage**:

   .. code-block:: cpp

      class ProgressReporter {
          std::atomic<size_t> solutions_found{0};
          std::chrono::steady_clock::time_point start_time;

      public:
          void reportSolution() {
              size_t count = ++solutions_found;
              if (count % 1000 == 0) {
                  auto elapsed = std::chrono::steady_clock::now() - start_time;
                  auto seconds = std::chrono::duration<double>(elapsed).count();
                  double rate = count / seconds;

                  std::cout << "Progress: " << count
                           << " solutions (" << rate << "/sec)" << std::endl;
              }
          }
      };

Performance Benchmarking
~~~~~~~~~~~~~~~~~~~~~~~~

Built-in benchmarking capabilities:

.. code-block:: cpp

   struct BenchmarkResults {
       size_t solutions_found;
       std::chrono::milliseconds sequential_time;
       std::chrono::milliseconds parallel_time;
       double speedup_ratio;
       double efficiency;
       int thread_count;
   };

   BenchmarkResults benchmark(const ProblemDefinition& problem) {
       // Run sequential version
       auto seq_start = std::chrono::steady_clock::now();
       size_t seq_solutions = runSequential(problem);
       auto seq_time = std::chrono::steady_clock::now() - seq_start;

       // Run parallel version
       auto par_start = std::chrono::steady_clock::now();
       size_t par_solutions = runParallel(problem);
       auto par_time = std::chrono::steady_clock::now() - par_start;

       return {
           par_solutions,
           std::chrono::duration_cast<std::chrono::milliseconds>(seq_time),
           std::chrono::duration_cast<std::chrono::milliseconds>(par_time),
           static_cast<double>(seq_time.count()) / par_time.count(),
           /* efficiency calculation */,
           std::thread::hardware_concurrency()
       };
   }

Thread Synchronization
----------------------

Lock-Free Design
~~~~~~~~~~~~~~~~

The framework minimizes locks for maximum performance:

**Atomic Operations**:

.. code-block:: cpp

   class LockFreeCounter {
       std::atomic<size_t> count_{0};
   public:
       size_t increment() { return count_.fetch_add(1); }
       size_t get() const { return count_.load(); }
   };

**Memory Ordering**:

.. code-block:: cpp

   // Relaxed ordering for counters (performance critical)
   solutions_found.fetch_add(1, std::memory_order_relaxed);

   // Acquire-release for work queue synchronization
   bool work_available = queue_empty.load(std::memory_order_acquire);
   queue_empty.store(false, std::memory_order_release);

Condition Variables
~~~~~~~~~~~~~~~~~~

For coordination between threads:

.. code-block:: cpp

   class WorkNotifier {
       mutable std::mutex mutex_;
       std::condition_variable cv_;
       bool work_available_ = false;

   public:
       void notifyWork() {
           std::lock_guard<std::mutex> lock(mutex_);
           work_available_ = true;
           cv_.notify_all();
       }

       void waitForWork() {
           std::unique_lock<std::mutex> lock(mutex_);
           cv_.wait(lock, [this] { return work_available_; });
           work_available_ = false;
       }
   };

Error Handling in Parallel Context
----------------------------------

Exception Propagation
~~~~~~~~~~~~~~~~~~~~

Exceptions in worker threads are captured and re-thrown in main thread:

.. code-block:: cpp

   class ParallelExceptionHandler {
       std::exception_ptr stored_exception_;
       std::mutex exception_mutex_;

   public:
       void captureException(std::exception_ptr e) {
           std::lock_guard<std::mutex> lock(exception_mutex_);
           if (!stored_exception_) {
               stored_exception_ = e;
           }
       }

       void rethrowIfException() {
           std::lock_guard<std::mutex> lock(exception_mutex_);
           if (stored_exception_) {
               std::rethrow_exception(stored_exception_);
           }
       }
   };

Graceful Shutdown
~~~~~~~~~~~~~~~~

.. code-block:: cpp

   class GracefulShutdown {
       std::atomic<bool> shutdown_requested_{false};

   public:
       void requestShutdown() {
           shutdown_requested_.store(true);
       }

       bool shouldContinue() const {
           return !shutdown_requested_.load();
       }
   };

Usage Patterns
--------------

Pattern 1: Simple Parallel Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   void simpleParallelSearch() {
       std::vector<std::vector<int>> solutions;
       std::mutex solutions_mutex;

       auto callback = [&](const std::vector<int>& sol) {
           std::lock_guard<std::mutex> lock(solutions_mutex);
           solutions.push_back(sol);
       };

       parallelPooling(main_inequalities, supporting_inequalities, callback);

       std::cout << "Found " << solutions.size() << " solutions" << std::endl;
   }

Pattern 2: Progress Monitoring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   void monitoredParallelSearch() {
       std::atomic<size_t> progress{0};
       auto start_time = std::chrono::steady_clock::now();

       auto callback = [&](const std::vector<int>& sol) {
           size_t count = ++progress;
           if (count % 1000 == 0) {
               auto elapsed = std::chrono::steady_clock::now() - start_time;
               auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed);
               std::cout << "Progress: " << count << " solutions in "
                        << ms.count() << "ms" << std::endl;
           }
       };

       parallelPooling(main_inequalities, supporting_inequalities, callback);
   }

Pattern 3: Resource-Constrained Execution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   void resourceConstrainedSearch() {
       // Limit threads based on available memory
       size_t available_memory_gb = getAvailableMemoryGB();
       int max_threads = std::min(
           static_cast<int>(available_memory_gb / 2),  // 2GB per thread
           static_cast<int>(std::thread::hardware_concurrency())
       );

       auto callback = [](const std::vector<int>& sol) {
           // Lightweight processing only
       };

       parallelPooling(main_inequalities, supporting_inequalities,
                      callback, max_threads);
   }

Best Practices
--------------

**Performance**

1. **Use appropriate thread counts**: More threads ≠ faster for small problems
2. **Minimize callback overhead**: Keep solution processing lightweight
3. **Enable checkpointing**: For computations > 1 hour
4. **Monitor memory usage**: Large problems can exhaust memory

**Reliability**

1. **Handle interruptions gracefully**: Use signal handlers with checkpointing
2. **Validate results**: Compare solution counts between runs
3. **Test with different thread counts**: Ensure correctness across configurations
4. **Use timeouts**: For computations with unknown bounds

**Debugging**

1. **Start with single thread**: Easier to debug algorithmic issues
2. **Use thread-safe logging**: Avoid interleaved output
3. **Monitor queue sizes**: Detect load balancing issues
4. **Profile with tools**: Use thread profilers for optimization