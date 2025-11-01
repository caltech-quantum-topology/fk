#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <functional>

/**
 * Checkpoint system for solution pool computation
 * Allows saving and resuming long-running computations
 */

struct CheckpointState {
  // Processing queue state
  std::queue<std::vector<std::vector<double>>> processing_queue;

  // Visited solutions tracking
  std::vector<std::vector<std::vector<double>>> visited_solutions;

  // Configuration
  std::vector<std::vector<double>> main_inequalities;
  std::vector<std::vector<double>> supporting_inequalities;

  // Progress tracking
  size_t total_processed;
  size_t solutions_found;

  // Metadata
  std::string timestamp;
  int version;
};

class CheckpointManager {
private:
  std::string checkpoint_file_;
  size_t save_interval_;
  size_t operations_since_save_;
  bool auto_save_enabled_;

public:
  /**
   * Constructor
   * @param checkpoint_file Path to checkpoint file
   * @param save_interval How often to save (in number of operations)
   */
  CheckpointManager(const std::string& checkpoint_file, size_t save_interval = 1000);

  /**
   * Save current state to checkpoint file
   */
  bool saveCheckpoint(const CheckpointState& state);

  /**
   * Load state from checkpoint file
   */
  bool loadCheckpoint(CheckpointState& state);

  /**
   * Check if checkpoint file exists
   */
  bool checkpointExists() const;

  /**
   * Delete checkpoint file
   */
  bool deleteCheckpoint();

  /**
   * Enable/disable automatic saving
   */
  void setAutoSave(bool enabled);

  /**
   * Increment operation counter and auto-save if needed
   */
  void incrementOperation(const CheckpointState& state);

  /**
   * Force save regardless of interval
   */
  void forceSave(const CheckpointState& state);

  /**
   * Get progress information
   */
  std::string getProgressInfo(const CheckpointState& state) const;
};

/**
 * Checkpointed version of pooling function
 */
void poolingWithCheckpoints(
    std::vector<std::vector<double>> main_inequalities,
    std::vector<std::vector<double>> supporting_inequalities,
    const std::function<void(const std::vector<int>&)>& function,
    const std::string& checkpoint_file = "",
    bool resume_from_checkpoint = false,
    size_t checkpoint_interval = 1000
);