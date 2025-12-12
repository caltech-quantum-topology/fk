#pragma once

#include <queue>
#include <set>
#include <unordered_set>
#include <stack>
#include <vector>
#include <functional>
#include <optional>

namespace graph_search {

template<typename Node>
class BreadthFirstSearch {
public:
    using NodeValidator = std::function<bool(const Node&)>;
    using NeighborGenerator = std::function<std::vector<Node>(const Node&)>;
    using VisitedChecker = std::function<bool(const Node&)>;
    using VisitedMarker = std::function<void(const Node&)>;

    struct Config {
        NodeValidator is_valid_goal;
        NeighborGenerator get_neighbors;
        VisitedChecker has_been_visited;
        VisitedMarker mark_visited;
    };

    explicit BreadthFirstSearch(const Config& config) : config_(config) {}

    std::optional<Node> search(const Node& start_node) {
        std::queue<Node> queue;

        config_.mark_visited(start_node);
        queue.push(start_node);

        while (!queue.empty()) {
            Node current = std::move(queue.front());
            queue.pop();

            if (config_.is_valid_goal(current)) {
                return current;
            }

            for (const auto& neighbor : config_.get_neighbors(current)) {
                if (!config_.has_been_visited(neighbor)) {
                    config_.mark_visited(neighbor);
                    queue.push(neighbor);
                }
            }
        }

        return std::nullopt;
    }

private:
    Config config_;
};

template<typename Node>
class BreadthFirstSearchWithVisited {
public:
    using NodeValidator = std::function<bool(const Node&)>;
    using NeighborGenerator = std::function<std::vector<Node>(const Node&)>;

    struct Config {
        NodeValidator is_valid_goal;
        NeighborGenerator get_neighbors;
    };

    explicit BreadthFirstSearchWithVisited(const Config& config) : config_(config) {}

    std::optional<Node> search(const Node& start_node) {
        std::queue<Node> queue;
        std::set<Node> visited;

        visited.insert(start_node);
        queue.push(start_node);

        while (!queue.empty()) {
            Node current = std::move(queue.front());
            queue.pop();

            if (config_.is_valid_goal(current)) {
                return current;
            }

            for (const auto& neighbor : config_.get_neighbors(current)) {
                if (visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    queue.push(neighbor);
                }
            }
        }

        return std::nullopt;
    }

private:
    Config config_;

};

template<typename Node>
class DepthFirstSearch {
public:
    using NodeValidator = std::function<bool(const Node&)>;
    using NeighborGenerator = std::function<std::vector<Node>(const Node&)>;
    using VisitedChecker = std::function<bool(const Node&)>;
    using VisitedMarker = std::function<void(const Node&)>;
    using NodeProcessor = std::function<void(const Node&)>;

    struct Config {
        NodeValidator is_valid_goal;
        NeighborGenerator get_neighbors;
        VisitedChecker has_been_visited;
        VisitedMarker mark_visited;
        NodeProcessor process_node;
    };

    explicit DepthFirstSearch(const Config& config) : config_(config) {}

    std::vector<Node> search_all(const Node& start_node) {
        std::vector<Node> results;
        std::stack<Node> stack;

        config_.mark_visited(start_node);
        stack.push(start_node);

        while (!stack.empty()) {
            Node current = std::move(stack.top());
            stack.pop();

            if (config_.is_valid_goal(current)) {
                results.push_back(current);
                if (config_.process_node) {
                    config_.process_node(current);
                }
                continue;
            }

            for (const auto& neighbor : config_.get_neighbors(current)) {
                if (!config_.has_been_visited(neighbor)) {
                    config_.mark_visited(neighbor);
                    stack.push(neighbor);
                }
            }
        }

        return results;
    }

    std::optional<Node> search_first(const Node& start_node) {
        std::stack<Node> stack;

        config_.mark_visited(start_node);
        stack.push(start_node);

        while (!stack.empty()) {
            Node current = std::move(stack.top());
            stack.pop();

            if (config_.is_valid_goal(current)) {
                if (config_.process_node) {
                    config_.process_node(current);
                }
                return current;
            }

            for (const auto& neighbor : config_.get_neighbors(current)) {
                if (!config_.has_been_visited(neighbor)) {
                    config_.mark_visited(neighbor);
                    stack.push(neighbor);
                }
            }
        }

        return std::nullopt;
    }

private:
    Config config_;
};

template<typename Node>
class DepthFirstSearchWithVisited {
public:
    using NodeValidator = std::function<bool(const Node&)>;
    using NeighborGenerator = std::function<std::vector<Node>(const Node&)>;
    using NodeProcessor = std::function<void(const Node&)>;

    struct Config {
        NodeValidator is_valid_goal;
        NeighborGenerator get_neighbors;
        NodeProcessor process_node;
    };

    explicit DepthFirstSearchWithVisited(const Config& config) : config_(config) {}

    std::vector<Node> search_all(const Node& start_node) {
        std::vector<Node> results;
        std::stack<Node> stack;
        std::unordered_set<Node> visited;

        visited.insert(start_node);
        stack.push(start_node);

        while (!stack.empty()) {
            Node current = std::move(stack.top());
            stack.pop();

            if (config_.is_valid_goal(current)) {
                results.push_back(current);
                if (config_.process_node) {
                    config_.process_node(current);
                }
                continue;
            }

            for (const auto& neighbor : config_.get_neighbors(current)) {
                if (visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    stack.push(neighbor);
                }
            }
        }

        return results;
    }

    std::optional<Node> search_first(const Node& start_node) {
        std::stack<Node> stack;
        std::unordered_set<Node> visited;

        visited.insert(start_node);
        stack.push(start_node);

        while (!stack.empty()) {
            Node current = std::move(stack.top());
            stack.pop();

            if (config_.is_valid_goal(current)) {
                if (config_.process_node) {
                    config_.process_node(current);
                }
                return current;
            }

            for (const auto& neighbor : config_.get_neighbors(current)) {
                if (visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    stack.push(neighbor);
                }
            }
        }

        return std::nullopt;
    }

private:
    Config config_;
};

} // namespace graph_search