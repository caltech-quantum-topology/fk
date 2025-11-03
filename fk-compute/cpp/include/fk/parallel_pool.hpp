#pragma once

#include "fk/solution_pool_1a_double_links.hpp"

#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

/**
 * Thread-safe work item for parallel processing
 */
struct WorkItem {
  std::vector<std::vector<double>> criteria;
  std::vector<std::vector<double>> main_inequalities;
  std::vector<std::vector<double>> supporting_inequalities;
  int size;
  size_t work_id;
};

/**
 * Thread-safe solution collector
 */
class ThreadSafeSolutionCollector {
private:
  std::mutex solutions_mutex_;
  std::atomic<size_t> solution_count_;
  std::function<void(const std::vector<int>&)> user_function_;

public:
  ThreadSafeSolutionCollector(const std::function<void(const std::vector<int>&)>& func);

  void addSolution(const std::vector<int>& solution);
  size_t getSolutionCount() const;
};

/**
 * Thread-safe work queue with work stealing
 */
class WorkQueue {
private:
  std::queue<WorkItem> queue_;
  mutable std::mutex queue_mutex_;
  std::condition_variable queue_cv_;
  std::atomic<bool> finished_;
  std::atomic<size_t> active_workers_;
  std::atomic<size_t> work_id_counter_;

public:
  WorkQueue();

  void addWork(const WorkItem& item);
  bool getWork(WorkItem& item, std::chrono::milliseconds timeout = std::chrono::milliseconds(100));
  void setFinished();
  bool isFinished() const;
  void incrementActiveWorkers();
  void decrementActiveWorkers();
  size_t getActiveWorkers() const;
  size_t generateWorkId();
  size_t getQueueSize() const;
};


/**
 * Parallel worker thread class
 */
class ParallelWorker {
private:
  int worker_id_;
  WorkQueue* work_queue_;
  ThreadSafeSolutionCollector* solution_collector_;
  std::atomic<size_t>* operations_processed_;
  std::atomic<size_t>* solutions_found_;
  std::atomic<bool>* should_stop_;

public:
  ParallelWorker(int id, WorkQueue* queue, ThreadSafeSolutionCollector* collector,
                 std::atomic<size_t>* operations_processed, std::atomic<size_t>* solutions_found,
                 std::atomic<bool>* stop_flag);

  void run();

private:
  void processWorkItem(const WorkItem& item);
};

/**
 * Thread pool manager for parallel computation
 */
class ThreadPool {
private:
  std::vector<std::thread> workers_;
  std::vector<std::unique_ptr<ParallelWorker>> worker_objects_;
  WorkQueue work_queue_;
  std::atomic<bool> should_stop_;
  int num_threads_;

public:
  ThreadPool(int num_threads = -1);  // -1 = use hardware concurrency
  ~ThreadPool();

  void start(ThreadSafeSolutionCollector* collector, std::atomic<size_t>* operations_processed,
             std::atomic<size_t>* solutions_found);
  void stop();
  void addWork(const WorkItem& item);
  bool isFinished() const;
  size_t getQueueSize() const;
  size_t getActiveWorkers() const;
  int getNumThreads() const;
};


/**
 * Parallel version without checkpointing for simple cases
 */
void parallelPooling(
    std::vector<std::vector<double>> main_inequalities,
    std::vector<std::vector<double>> supporting_inequalities,
    const std::function<void(const std::vector<int>&)>& function,
    int num_threads = -1  // -1 = use hardware concurrency
);

/**
 * Utility function to get optimal number of threads
 */
int getOptimalThreadCount();