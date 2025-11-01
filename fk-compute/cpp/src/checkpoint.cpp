#include "fk/checkpoint.hpp"
#include "fk/solution_pool_1a_double_links.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

CheckpointManager::CheckpointManager(const std::string& checkpoint_file, size_t save_interval)
    : checkpoint_file_(checkpoint_file), save_interval_(save_interval),
      operations_since_save_(0), auto_save_enabled_(true) {}

bool CheckpointManager::saveCheckpoint(const CheckpointState& state) {
    std::ofstream file(checkpoint_file_, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open checkpoint file for writing: " << checkpoint_file_ << std::endl;
        return false;
    }

    try {
        // Write version
        file.write(reinterpret_cast<const char*>(&state.version), sizeof(state.version));

        // Write timestamp
        size_t timestamp_size = state.timestamp.size();
        file.write(reinterpret_cast<const char*>(&timestamp_size), sizeof(timestamp_size));
        file.write(state.timestamp.c_str(), timestamp_size);

        // Write progress tracking
        file.write(reinterpret_cast<const char*>(&state.total_processed), sizeof(state.total_processed));
        file.write(reinterpret_cast<const char*>(&state.solutions_found), sizeof(state.solutions_found));

        // Write main inequalities
        size_t main_rows = state.main_inequalities.size();
        file.write(reinterpret_cast<const char*>(&main_rows), sizeof(main_rows));
        if (main_rows > 0) {
            size_t main_cols = state.main_inequalities[0].size();
            file.write(reinterpret_cast<const char*>(&main_cols), sizeof(main_cols));
            for (const auto& row : state.main_inequalities) {
                file.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(double));
            }
        }

        // Write supporting inequalities
        size_t support_rows = state.supporting_inequalities.size();
        file.write(reinterpret_cast<const char*>(&support_rows), sizeof(support_rows));
        if (support_rows > 0) {
            size_t support_cols = state.supporting_inequalities[0].size();
            file.write(reinterpret_cast<const char*>(&support_cols), sizeof(support_cols));
            for (const auto& row : state.supporting_inequalities) {
                file.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(double));
            }
        }

        // Write processing queue
        size_t queue_size = state.processing_queue.size();
        file.write(reinterpret_cast<const char*>(&queue_size), sizeof(queue_size));

        // Create a copy of the queue to iterate through it
        std::queue<std::vector<std::vector<double>>> queue_copy = state.processing_queue;
        while (!queue_copy.empty()) {
            const auto& criteria = queue_copy.front();
            size_t criteria_rows = criteria.size();
            file.write(reinterpret_cast<const char*>(&criteria_rows), sizeof(criteria_rows));
            if (criteria_rows > 0) {
                size_t criteria_cols = criteria[0].size();
                file.write(reinterpret_cast<const char*>(&criteria_cols), sizeof(criteria_cols));
                for (const auto& row : criteria) {
                    file.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(double));
                }
            }
            queue_copy.pop();
        }

        // Write visited solutions
        size_t visited_size = state.visited_solutions.size();
        file.write(reinterpret_cast<const char*>(&visited_size), sizeof(visited_size));
        for (const auto& solution_set : state.visited_solutions) {
            size_t set_size = solution_set.size();
            file.write(reinterpret_cast<const char*>(&set_size), sizeof(set_size));
            for (const auto& solution : solution_set) {
                size_t solution_size = solution.size();
                file.write(reinterpret_cast<const char*>(&solution_size), sizeof(solution_size));
                file.write(reinterpret_cast<const char*>(solution.data()), solution.size() * sizeof(double));
            }
        }

        file.close();
        operations_since_save_ = 0;
        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error saving checkpoint: " << e.what() << std::endl;
        return false;
    }
}

bool CheckpointManager::loadCheckpoint(CheckpointState& state) {
    std::ifstream file(checkpoint_file_, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }

    try {
        // Read version
        file.read(reinterpret_cast<char*>(&state.version), sizeof(state.version));

        // Read timestamp
        size_t timestamp_size;
        file.read(reinterpret_cast<char*>(&timestamp_size), sizeof(timestamp_size));
        state.timestamp.resize(timestamp_size);
        file.read(&state.timestamp[0], timestamp_size);

        // Read progress tracking
        file.read(reinterpret_cast<char*>(&state.total_processed), sizeof(state.total_processed));
        file.read(reinterpret_cast<char*>(&state.solutions_found), sizeof(state.solutions_found));

        // Read main inequalities
        size_t main_rows;
        file.read(reinterpret_cast<char*>(&main_rows), sizeof(main_rows));
        state.main_inequalities.resize(main_rows);
        if (main_rows > 0) {
            size_t main_cols;
            file.read(reinterpret_cast<char*>(&main_cols), sizeof(main_cols));
            for (auto& row : state.main_inequalities) {
                row.resize(main_cols);
                file.read(reinterpret_cast<char*>(row.data()), row.size() * sizeof(double));
            }
        }

        // Read supporting inequalities
        size_t support_rows;
        file.read(reinterpret_cast<char*>(&support_rows), sizeof(support_rows));
        state.supporting_inequalities.resize(support_rows);
        if (support_rows > 0) {
            size_t support_cols;
            file.read(reinterpret_cast<char*>(&support_cols), sizeof(support_cols));
            for (auto& row : state.supporting_inequalities) {
                row.resize(support_cols);
                file.read(reinterpret_cast<char*>(row.data()), row.size() * sizeof(double));
            }
        }

        // Read processing queue
        size_t queue_size;
        file.read(reinterpret_cast<char*>(&queue_size), sizeof(queue_size));

        // Clear the queue and rebuild it
        while (!state.processing_queue.empty()) {
            state.processing_queue.pop();
        }

        for (size_t i = 0; i < queue_size; i++) {
            size_t criteria_rows;
            file.read(reinterpret_cast<char*>(&criteria_rows), sizeof(criteria_rows));
            std::vector<std::vector<double>> criteria(criteria_rows);
            if (criteria_rows > 0) {
                size_t criteria_cols;
                file.read(reinterpret_cast<char*>(&criteria_cols), sizeof(criteria_cols));
                for (auto& row : criteria) {
                    row.resize(criteria_cols);
                    file.read(reinterpret_cast<char*>(row.data()), row.size() * sizeof(double));
                }
            }
            state.processing_queue.push(criteria);
        }

        // Read visited solutions
        size_t visited_size;
        file.read(reinterpret_cast<char*>(&visited_size), sizeof(visited_size));
        state.visited_solutions.resize(visited_size);
        for (auto& solution_set : state.visited_solutions) {
            size_t set_size;
            file.read(reinterpret_cast<char*>(&set_size), sizeof(set_size));
            solution_set.resize(set_size);
            for (auto& solution : solution_set) {
                size_t solution_size;
                file.read(reinterpret_cast<char*>(&solution_size), sizeof(solution_size));
                solution.resize(solution_size);
                file.read(reinterpret_cast<char*>(solution.data()), solution.size() * sizeof(double));
            }
        }

        file.close();
        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error loading checkpoint: " << e.what() << std::endl;
        return false;
    }
}

bool CheckpointManager::checkpointExists() const {
    std::ifstream file(checkpoint_file_);
    return file.good();
}

bool CheckpointManager::deleteCheckpoint() {
    return std::remove(checkpoint_file_.c_str()) == 0;
}

void CheckpointManager::setAutoSave(bool enabled) {
    auto_save_enabled_ = enabled;
}

void CheckpointManager::incrementOperation(const CheckpointState& state) {
    operations_since_save_++;
    if (auto_save_enabled_ && operations_since_save_ >= save_interval_) {
        saveCheckpoint(state);
    }
}

void CheckpointManager::forceSave(const CheckpointState& state) {
    saveCheckpoint(state);
}

std::string CheckpointManager::getProgressInfo(const CheckpointState& state) const {
    std::ostringstream oss;
    oss << "Progress: " << state.total_processed << " operations processed, "
        << state.solutions_found << " solutions found, "
        << state.processing_queue.size() << " items in queue";
    if (!state.timestamp.empty()) {
        oss << " (last saved: " << state.timestamp << ")";
    }
    return oss.str();
}

// Helper function to get current timestamp
std::string getCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

void poolingWithCheckpoints(
    std::vector<std::vector<double>> main_inequalities,
    std::vector<std::vector<double>> supporting_inequalities,
    const std::function<void(const std::vector<int>&)>& function,
    const std::string& checkpoint_file,
    bool resume_from_checkpoint,
    size_t checkpoint_interval) {

    CheckpointState state;
    CheckpointManager checkpoint_manager(checkpoint_file, checkpoint_interval);

    // Initialize state
    state.version = 1;
    state.total_processed = 0;
    state.solutions_found = 0;
    state.timestamp = getCurrentTimestamp();

    bool should_resume = false;

    // Try to load from checkpoint if requested
    if (resume_from_checkpoint && !checkpoint_file.empty() && checkpoint_manager.checkpointExists()) {
        std::cout << "Loading checkpoint from: " << checkpoint_file << std::endl;
        if (checkpoint_manager.loadCheckpoint(state)) {
            std::cout << "Checkpoint loaded successfully." << std::endl;
            std::cout << checkpoint_manager.getProgressInfo(state) << std::endl;
            should_resume = true;
        } else {
            std::cout << "Failed to load checkpoint, starting fresh." << std::endl;
        }
    }

    if (!should_resume) {
        // Initialize fresh state
        state.main_inequalities = main_inequalities;
        state.supporting_inequalities = supporting_inequalities;

        // Add main inequalities to supporting inequalities
        for (const auto& x : main_inequalities) {
            state.supporting_inequalities.push_back(x);
        }

        // Initialize visited solutions tracking
        int mains = main_inequalities.size();
        state.visited_solutions.resize(mains);
        for (int i = 0; i < mains; i++) {
            state.visited_solutions[i].push_back(main_inequalities[i]);
        }

        // Initialize processing queue
        state.processing_queue.push(main_inequalities);
    }

    // Create wrapper function that tracks solutions found
    auto tracking_function = [&](const std::vector<int>& point) {
        state.solutions_found++;
        function(point);  // Call original function
    };

    // Main processing loop with checkpointing
    int mains = state.main_inequalities.size();
    int size = state.main_inequalities[0].size();
    int support = state.supporting_inequalities.size();

    // Try to process the initial criteria if not resumed
    if (!should_resume) {
        if (processCriteria(state.main_inequalities, state.main_inequalities,
                           state.supporting_inequalities, tracking_function, size)) {
            // Save final checkpoint
            if (!checkpoint_file.empty()) {
                state.timestamp = getCurrentTimestamp();
                checkpoint_manager.forceSave(state);
                std::cout << "Computation completed. Final checkpoint saved." << std::endl;
            }
            return;
        }
    }

    std::vector<std::vector<double>> criteria(mains, std::vector<double>(size));
    std::vector<std::vector<double>> new_criteria(mains, std::vector<double>(size));

    // Main processing loop
    while (!state.processing_queue.empty()) {
        criteria = state.processing_queue.front();
        state.processing_queue.pop();
        state.total_processed++;

        // Try different supporting inequalities
        for (int i = 0; i < support; i++) {
            bool found_improvement = false;

            // Check each main criterion
            for (int q = 0; q < mains && !found_improvement; q++) {
                // Check each variable coefficient
                for (int j = 1; j < size; j++) {
                    if (criteria[q][j] > 0 && state.supporting_inequalities[i][j] < 0) {
                        // Create new criteria by combining with supporting inequality
                        for (int k = 0; k < size; k++) {
                            new_criteria[q][k] = criteria[q][k] + state.supporting_inequalities[i][k] / 2.0;
                        }

                        // Check if this new criterion set has been visited
                        bool is_new = false;
                        for (int visdex = 0; visdex < mains; visdex++) {
                            bool found_in_visited = false;
                            for (const auto& visited : state.visited_solutions[visdex]) {
                                if (visited == new_criteria[visdex]) {
                                    found_in_visited = true;
                                    break;
                                }
                            }
                            if (!found_in_visited) {
                                is_new = true;
                                break;
                            }
                        }

                        if (is_new) {
                            // Mark as visited
                            for (int s = 0; s < mains; s++) {
                                state.visited_solutions[s].push_back(new_criteria[s]);
                            }

                            // Try to process new criteria
                            if (processCriteria(new_criteria, state.main_inequalities,
                                              state.supporting_inequalities, tracking_function, size)) {
                                // Save final checkpoint
                                if (!checkpoint_file.empty()) {
                                    state.timestamp = getCurrentTimestamp();
                                    checkpoint_manager.forceSave(state);
                                    std::cout << "Computation completed. Final checkpoint saved." << std::endl;
                                }
                                return;
                            }

                            // Add to queue for further processing
                            state.processing_queue.push(new_criteria);
                            found_improvement = true;
                        }
                        break;
                    }
                }
            }
        }

        // Periodic checkpointing
        if (!checkpoint_file.empty()) {
            state.timestamp = getCurrentTimestamp();
            checkpoint_manager.incrementOperation(state);

            // Print progress periodically
            if (state.total_processed % (checkpoint_interval * 10) == 0) {
                std::cout << checkpoint_manager.getProgressInfo(state) << std::endl;
            }
        }
    }

    // Save final checkpoint
    if (!checkpoint_file.empty()) {
        state.timestamp = getCurrentTimestamp();
        checkpoint_manager.forceSave(state);
        std::cout << "Computation completed. Final checkpoint saved." << std::endl;
        std::cout << checkpoint_manager.getProgressInfo(state) << std::endl;
    }
}
