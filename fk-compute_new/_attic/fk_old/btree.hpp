/**
 * @file btree.hpp
 * @brief Specialized B-Tree implementation for vector storage and BFS history tracking
 * 
 * This file implements a reverse-order B-Tree optimized for storing and querying
 * vectors in breadth-first search algorithms. The tree is used as a visited-state
 * history to prevent revisiting the same vectors during BFS traversal in 
 * Gukov-Manolescu invariant computations.
 * 
 * Key Features:
 * - Reverse-order storage: vectors stored from last element to first
 * - Substring matching: shorter vectors match if they are suffixes of stored vectors
 * - Memory-efficient: uses dynamic allocation for tree nodes
 * - Fast lookup: O(vector_length) containment queries
 * 
 * Usage Pattern:
 * 1. Insert vectors using update() during BFS exploration
 * 2. Check if vectors already visited using contains() before processing
 * 3. Supports partial matching for suffix-based queries
 * 
 * WARNING: Current implementation has memory leak - custom destructor needed!
 */

#pragma once

#include <vector>
/**
 * @brief Reverse-order B-Tree for efficient vector storage and BFS history tracking
 * @tparam T Element type of vectors to be stored (typically int for integer vectors)
 * 
 * This tree stores vectors in reverse order (last element to first) to enable
 * efficient suffix matching. Each node represents one element of a vector path,
 * with children representing subsequent elements in the reverse traversal.
 * 
 * Tree Structure Example:
 * For vectors [1,2,3] and [4,2,3]:
 * Root -> 3 -> 2 -> {1, 4}
 * 
 * This enables fast suffix queries: [2,3] matches both stored vectors.
 */
template <typename T>
struct btree {
    private:
        std::vector<T> children_values = {};      // Values at current tree level
        std::vector<btree<T>*> children_pointers = {}; // Pointers to child subtrees
        /**
         * @brief Recursive helper to insert vector into tree (reverse-order traversal)
         * @param v Vector to insert
         * @param size Current position in reverse traversal (v.size() down to 0)
         * 
         * Algorithm:
         * 1. Process vector elements from last to first (reverse order)
         * 2. Search existing children for matching value at current position
         * 3. If match found, recursively insert into existing subtree
         * 4. If no match, create new child node and continue recursion
         * 
         * Base case: size == 0 indicates complete vector path has been stored
         */
        void private_update (std::vector<T> v, int size) {
            if (size != 0) {
                size--;  // Move to next element in reverse order
                
                // Search for existing child with current element value
                for (int i = 0; i < children_values.size(); i++) {
                    if (children_values[i] == v[size]) {
                        // Found existing path, continue recursion in existing subtree
                        (*children_pointers[i]).private_update(v, size);
                        return;
                    }
                }
                
                // No existing child found, create new branch
                // WARNING: Memory leak! Need custom destructor to manage child nodes
                btree<T>* btr_ptr = new btree<T>; 
                children_values.push_back(v[size]);
                children_pointers.push_back(btr_ptr);
                
                // Continue recursive insertion in new child
                (*btr_ptr).private_update(v, size);
            }
            // Base case: size == 0, complete vector path stored
        }
        /**
         * @brief Recursive helper to check if vector exists in tree (reverse-order search)
         * @param v Vector to search for
         * @param size Current position in reverse traversal (v.size() down to 0)
         * @return True if vector (or its suffix) exists in tree, false otherwise
         * 
         * Algorithm:
         * 1. Search vector elements from last to first (matching insertion order)
         * 2. At each level, look for child with matching element value
         * 3. If match found, continue search in corresponding subtree
         * 4. If no match found at any level, vector is not stored
         * 
         * Base case: size == 0 indicates complete path found (vector exists)
         * 
         * Note: This enables suffix matching - shorter vectors will match if they
         * correspond to suffixes of any stored vectors.
         */
        bool private_contains (std::vector<T> v, int size) {
            if (size != 0) {
                size--;  // Move to next element in reverse order
                
                // Search children for matching element value
                for (int i = 0; i < children_values.size(); i++) {
                    if (children_values[i] == v[size]) {
                        // Found matching path, continue search in subtree
                        return (*children_pointers[i]).private_contains(v, size);
                    }
                }
                
                // No matching child found, vector not in tree
                return false;
            }
            
            // Base case: completed full traversal, vector exists
            return true;
        }
    public:
        /**
         * @brief Insert vector into the tree for BFS history tracking
         * @param v Vector to insert (typically representing a state in BFS)
         * 
         * This is the main insertion interface used during BFS traversal.
         * Call this method to mark a vector as "visited" to prevent
         * revisiting the same state in future BFS iterations.
         * 
         * Time Complexity: O(vector_length * children_per_level)
         * Space Complexity: O(vector_length) for new unique paths
         */
        void update (std::vector<T> v) {
            private_update(v, v.size());
        }
        
        /**
         * @brief Check if vector already exists in the tree
         * @param v Vector to search for
         * @return True if vector has been previously inserted, false otherwise
         * 
         * This is the main query interface used during BFS traversal.
         * Call this method before processing a vector to determine if
         * it has already been visited.
         * 
         * Important: Due to reverse-order storage, this also returns true
         * for suffix matches. For example, if [1,2,3] is stored, then
         * contains([2,3]) will return true.
         * 
         * Time Complexity: O(vector_length * children_per_level)
         */
        bool contains (std::vector<T> v) {
            return private_contains(v, v.size());
        }
        
        /**
         * @brief Get values stored at current tree level (for debugging/inspection)
         * @return Copy of children_values vector
         */
        std::vector<T> get_children_values () {
            return children_values;
        }
        
        /**
         * @brief Get pointers to child subtrees (for debugging/inspection)
         * @return Copy of children_pointers vector
         * 
         * WARNING: Be careful with these pointers - they are managed internally
         * and should not be deleted externally.
         */
        std::vector<btree<T>*> get_children_pointers () {
            return children_pointers;
        }
};