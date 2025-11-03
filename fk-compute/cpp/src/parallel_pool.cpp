#include "fk/parallel_pool.hpp"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>

// ThreadSafeSolutionCollector implementation
ThreadSafeSolutionCollector::ThreadSafeSolutionCollector(
    const std::function<void(const std::vector<int>&)>& func)
    : solution_count_(0), user_function_(func) {}

void ThreadSafeSolutionCollector::addSolution(const std::vector<int>& solution) {
  {
    std::lock_guard<std::mutex> lock(solutions_mutex_);
    user_function_(solution);
  }
  solution_count_.fetch_add(1);
}

size_t ThreadSafeSolutionCollector::getSolutionCount() const {
  return solution_count_.load();
}

// WorkQueue implementation
WorkQueue::WorkQueue() : finished_(false), active_workers_(0), work_id_counter_(0) {}

void WorkQueue::addWork(const WorkItem& item) {
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    queue_.push(item);
  }
  queue_cv_.notify_one();
}

bool WorkQueue::getWork(WorkItem& item, std::chrono::milliseconds timeout) {
  std::unique_lock<std::mutex> lock(queue_mutex_);

  if (queue_cv_.wait_for(lock, timeout, [this] { return !queue_.empty() || finished_.load(); })) {
    if (!queue_.empty()) {
      item = queue_.front();
      queue_.pop();
      return true;
    }
  }
  return false;
}

void WorkQueue::setFinished() {
  finished_.store(true);
  queue_cv_.notify_all();
}

bool WorkQueue::isFinished() const {
  return finished_.load() && queue_.empty();
}

void WorkQueue::incrementActiveWorkers() {
  active_workers_.fetch_add(1);
}

void WorkQueue::decrementActiveWorkers() {
  active_workers_.fetch_sub(1);
}

size_t WorkQueue::getActiveWorkers() const {
  return active_workers_.load();
}

size_t WorkQueue::generateWorkId() {
  return work_id_counter_.fetch_add(1);
}

size_t WorkQueue::getQueueSize() const {
  std::lock_guard<std::mutex> lock(queue_mutex_);
  return queue_.size();
}

// ParallelCheckpointState implementation
ParallelCheckpointState::ParallelCheckpointState()
    : operations_processed_(0), solutions_found_(0) {
  checkpoint_state_.version = 1;
  checkpoint_state_.total_processed = 0;
  checkpoint_state_.solutions_found = 0;
}

void ParallelCheckpointState::updateState(
    const std::queue<std::vector<std::vector<double>>>& processing_queue,
    const std::vector<std::vector<std::vector<double>>>& visited_solutions,
    const std::vector<std::vector<double>>& main_inequalities,
    const std::vector<std::vector<double>>& supporting_inequalities) {

  std::lock_guard<std::mutex> lock(state_mutex_);

  // Update checkpoint state
  checkpoint_state_.main_inequalities = main_inequalities;
  checkpoint_state_.supporting_inequalities = supporting_inequalities;
  checkpoint_state_.visited_solutions = visited_solutions;
  checkpoint_state_.total_processed = operations_processed_.load();
  checkpoint_state_.solutions_found = solutions_found_.load();

  // Clear and rebuild processing queue
  while (!checkpoint_state_.processing_queue.empty()) {
    checkpoint_state_.processing_queue.pop();
  }

  std::queue<std::vector<std::vector<double>>> temp_queue = processing_queue;
  while (!temp_queue.empty()) {
    checkpoint_state_.processing_queue.push(temp_queue.front());
    temp_queue.pop();
  }
}

CheckpointState ParallelCheckpointState::getCheckpointState() const {
  std::lock_guard<std::mutex> lock(state_mutex_);
  return checkpoint_state_;
}

void ParallelCheckpointState::incrementOperations() {
  operations_processed_.fetch_add(1);
}

void ParallelCheckpointState::incrementSolutions() {
  solutions_found_.fetch_add(1);
}

size_t ParallelCheckpointState::getOperationsProcessed() const {
  return operations_processed_.load();
}

size_t ParallelCheckpointState::getSolutionsFound() const {
  return solutions_found_.load();
}

// ParallelWorker implementation
ParallelWorker::ParallelWorker(int id, WorkQueue* queue, ThreadSafeSolutionCollector* collector,
                               ParallelCheckpointState* state, std::atomic<bool>* stop_flag)
    : worker_id_(id), work_queue_(queue), solution_collector_(collector),
      checkpoint_state_(state), should_stop_(stop_flag) {}

void ParallelWorker::run() {
  WorkItem item;

  while (!should_stop_->load()) {
    if (work_queue_->getWork(item, std::chrono::milliseconds(100))) {
      work_queue_->incrementActiveWorkers();
      processWorkItem(item);
      work_queue_->decrementActiveWorkers();
      checkpoint_state_->incrementOperations();
    } else if (work_queue_->isFinished() && work_queue_->getActiveWorkers() == 0) {
      break;
    }
  }
}

void ParallelWorker::processWorkItem(const WorkItem& item) {
  // Create a local solution collector function that forwards to the thread-safe collector
  auto local_solution_handler = [this](const std::vector<int>& solution) {
    solution_collector_->addSolution(solution);
    checkpoint_state_->incrementSolutions();
  };

  // Process the criteria using the existing processCriteria function
  if (processCriteria(item.criteria, item.main_inequalities,
                     item.supporting_inequalities, local_solution_handler, item.size)) {
    // If processCriteria returns true, it means computation is complete
    return;
  }

  // Generate new work items by trying different supporting inequalities
  int mains = item.main_inequalities.size();
  int support = item.supporting_inequalities.size();

  for (int i = 0; i < support && !should_stop_->load(); i++) {
    for (int q = 0; q < mains; q++) {
      for (int j = 1; j < item.size; j++) {
        if (item.criteria[q][j] > 0 && item.supporting_inequalities[i][j] < 0) {
          // Create new criteria by combining with supporting inequality
          std::vector<std::vector<double>> new_criteria = item.criteria;
          for (int k = 0; k < item.size; k++) {
            new_criteria[q][k] = item.criteria[q][k] + item.supporting_inequalities[i][k] / 2.0;
          }

          // Create new work item and add to queue
          WorkItem new_item;
          new_item.criteria = new_criteria;
          new_item.main_inequalities = item.main_inequalities;
          new_item.supporting_inequalities = item.supporting_inequalities;
          new_item.size = item.size;
          new_item.work_id = work_queue_->generateWorkId();

          work_queue_->addWork(new_item);
          break;
        }
      }
    }
  }
}

// ThreadPool implementation
ThreadPool::ThreadPool(int num_threads) : should_stop_(false) {
  if (num_threads <= 0) {
    num_threads_ = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
  } else {
    num_threads_ = num_threads;
  }
}

ThreadPool::~ThreadPool() {
  stop();
}

void ThreadPool::start(ThreadSafeSolutionCollector* collector, ParallelCheckpointState* checkpoint_state) {
  should_stop_.store(false);

  // Create worker objects and threads
  for (int i = 0; i < num_threads_; i++) {
    worker_objects_.emplace_back(
        std::make_unique<ParallelWorker>(i, &work_queue_, collector, checkpoint_state, &should_stop_));

    workers_.emplace_back([this, i]() {
      worker_objects_[i]->run();
    });
  }
}

void ThreadPool::stop() {
  should_stop_.store(true);
  work_queue_.setFinished();

  for (auto& worker : workers_) {
    if (worker.joinable()) {
      worker.join();
    }
  }

  workers_.clear();
  worker_objects_.clear();
}

void ThreadPool::addWork(const WorkItem& item) {
  work_queue_.addWork(item);
}

bool ThreadPool::isFinished() const {
  return work_queue_.isFinished();
}

size_t ThreadPool::getQueueSize() const {
  return work_queue_.getQueueSize();
}

size_t ThreadPool::getActiveWorkers() const {
  return work_queue_.getActiveWorkers();
}

int ThreadPool::getNumThreads() const {
  return num_threads_;
}

// Utility function implementation
int getOptimalThreadCount() {
  int hardware_threads = std::thread::hardware_concurrency();
  return hardware_threads > 0 ? hardware_threads : 4;  // fallback to 4 if detection fails
}

// Simple parallel version without checkpointing
void parallelPooling(
    std::vector<std::vector<double>> main_inequalities,
    std::vector<std::vector<double>> supporting_inequalities,
    const std::function<void(const std::vector<int>&)>& function,
    int num_threads) {

  int mains = main_inequalities.size();
  int size = main_inequalities[0].size();

  // Add main inequalities to supporting inequalities
  for (const auto& x : main_inequalities) {
    supporting_inequalities.push_back(x);
  }

  // Create thread pool and collectors
  ThreadPool thread_pool(num_threads);
  ThreadSafeSolutionCollector solution_collector(function);
  ParallelCheckpointState checkpoint_state;

  std::cout << "Starting parallel computation with " << thread_pool.getNumThreads()
            << " threads..." << std::endl;

  // Start thread pool
  thread_pool.start(&solution_collector, &checkpoint_state);

  // Try to process the initial criteria
  if (processCriteria(main_inequalities, main_inequalities, supporting_inequalities,
                     function, size)) {
    thread_pool.stop();
    return;
  }

  // Create initial work item
  WorkItem initial_work;
  initial_work.criteria = main_inequalities;
  initial_work.main_inequalities = main_inequalities;
  initial_work.supporting_inequalities = supporting_inequalities;
  initial_work.size = size;
  initial_work.work_id = 0;

  thread_pool.addWork(initial_work);

  // Monitor progress
  auto start_time = std::chrono::steady_clock::now();
  size_t last_operations = 0;
  size_t last_solutions = 0;

  while (!thread_pool.isFinished()) {
    std::this_thread::sleep_for(std::chrono::seconds(5));

    size_t current_operations = checkpoint_state.getOperationsProcessed();
    size_t current_solutions = checkpoint_state.getSolutionsFound();

    auto current_time = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();

    std::cout << "Progress: " << current_operations << " operations ("
              << (current_operations - last_operations) << " new), "
              << current_solutions << " solutions ("
              << (current_solutions - last_solutions) << " new), "
              << thread_pool.getQueueSize() << " queued, "
              << thread_pool.getActiveWorkers() << " active workers, "
              << elapsed << "s elapsed" << std::endl;

    last_operations = current_operations;
    last_solutions = current_solutions;
  }

  thread_pool.stop();

  auto end_time = std::chrono::steady_clock::now();
  auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

  std::cout << "Parallel computation completed:" << std::endl;
  std::cout << "  Total operations: " << checkpoint_state.getOperationsProcessed() << std::endl;
  std::cout << "  Total solutions: " << checkpoint_state.getSolutionsFound() << std::endl;
  std::cout << "  Total time: " << total_time << "ms" << std::endl;
  std::cout << "  Threads used: " << thread_pool.getNumThreads() << std::endl;
}

// Parallel version with checkpointing
void parallelPoolingWithCheckpoints(
    std::vector<std::vector<double>> main_inequalities,
    std::vector<std::vector<double>> supporting_inequalities,
    const std::function<void(const std::vector<int>&)>& function,
    const std::string& checkpoint_file,
    bool resume_from_checkpoint,
    size_t checkpoint_interval,
    int num_threads) {

  // For now, delegate to the simpler version and add checkpoint integration later
  // This is a complex integration that would require careful synchronization
  std::cout << "Note: Checkpointing with parallel execution is complex and not yet fully implemented." << std::endl;
  std::cout << "Using parallel execution without checkpointing..." << std::endl;

  parallelPooling(main_inequalities, supporting_inequalities, function, num_threads);
}