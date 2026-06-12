/**
 * @file bilvector.hpp
 * @brief Bi-directional linked-list vector container for efficient sparse array operations
 * 
 * This file implements a specialized container that avoids the O(n²) cost of dynamic
 * vector reallocation by using linked lists of fixed-size vector chunks. The container
 * supports both positive and negative indexing, making it ideal for sparse polynomial
 * coefficient storage in Gukov-Manolescu invariant computations.
 * 
 * Key advantages over std::vector:
 * - Avoids quadratic reallocation costs for append-heavy workloads
 * - Supports negative indexing without offset calculations
 * - Memory-efficient for sparse data with clustered access patterns
 * - Amortized O(1) access for sequential operations
 */

#pragma once

#include <vector>
#include <list>

/**
 * @brief Bi-directional vector container using linked lists of vector chunks
 * @tparam T Element type stored in the container
 * 
 * This container manages two separate linked lists of vectors:
 * - pvectors: for positive indices [0, 1, 2, ...]
 * - nvectors: for negative indices [-1, -2, -3, ...]
 * 
 * Each vector chunk has fixed size `component_size`, avoiding the reallocation
 * costs of dynamically growing arrays. New chunks are allocated on-demand when
 * accessing indices beyond the current capacity.
 */
template <typename T>
struct bilvector {
    private:
        int component_size;    // Fixed size of each vector chunk
        int nnvectors = 0;     // Number of negative vector chunks allocated
        int npvectors = 0;     // Number of positive vector chunks allocated
        int max_nindex = 0;    // Most negative index accessed (for tracking)
        int max_pindex = 0;    // Most positive index accessed (for tracking)
        T default_;            // Default value for new elements
        std::list<std::vector<T>> nvectors = {};  // Chunks for negative indices
        std::list<std::vector<T>> pvectors = {};  // Chunks for positive indices
    public:
        /**
         * @brief Constructor initializing the bilvector with pre-allocated chunks
         * @param initial_nnvectors Number of negative vector chunks to pre-allocate
         * @param initial_npvectors Number of positive vector chunks to pre-allocate  
         * @param component_size_ Fixed size for each vector chunk
         * @param default__ Default value for initializing new elements
         * 
         * Pre-allocating chunks can improve performance for known access patterns,
         * avoiding initial allocation overhead during computation.
         */
        bilvector (int initial_nnvectors, int initial_npvectors, int component_size_, T default__) {
            component_size = component_size_;
            default_ = default__;
            
            // Pre-allocate negative vector chunks
            nnvectors = initial_nnvectors;
            for (int i = 0; i < initial_nnvectors; i++) {
                nvectors.push_back(std::vector<T>(component_size, default_));
            }
            
            // Pre-allocate positive vector chunks
            npvectors = initial_npvectors;
            for (int i = 0; i < initial_npvectors; i++) {
                pvectors.push_back(std::vector<T>(component_size, default_));
            }
        }
        /**
         * @brief Get total capacity for negative indices
         * @return Total number of elements that can be stored with negative indices
         */
        int nsize () {
            return nnvectors * component_size;
        }
        
        /**
         * @brief Get total capacity for positive indices
         * @return Total number of elements that can be stored with positive indices
         */
        int psize () {
            return npvectors * component_size;
        }
        /**
         * @brief Access element at given index with automatic capacity expansion
         * @param index Integer index (can be positive or negative)
         * @return Reference to element at the specified index
         * 
         * This operator supports bi-directional indexing:
         * - Positive indices [0, 1, 2, ...] access pvectors
         * - Negative indices [-1, -2, -3, ...] access nvectors
         * 
         * Capacity expansion occurs on-demand when accessing out-of-bounds indices.
         * The container automatically allocates additional vector chunks as needed,
         * avoiding the O(n) reallocation cost of std::vector::push_back().
         * 
         * Complexity: O(chunk_count) for chunk traversal, amortized O(1) for sequential access
         */
        T& operator[] (int index) {
            // Update access tracking for debugging/optimization
            if (index > max_pindex) {
                max_pindex = index;
            }
            else if (index < max_nindex) {
                max_nindex = index;
            }
            
            if (index >= 0) {
                // Handle positive indices using pvectors
                if (index >= (*this).psize()) {
                    // Calculate number of additional chunks needed
                    int x = (index - (*this).psize()) / component_size;
                    npvectors += x + 1;
                    
                    // Allocate new vector chunks to accommodate the index
                    for (int i = 0; i <= x; i++) {
                        pvectors.push_back(std::vector<T>(component_size, default_));
                    }
                }
                
                // Navigate to the appropriate vector chunk
                auto it = pvectors.begin();
                int j;
                for (j = 0; j < index / component_size; j++) {
                    ++it;
                }
                
                // Return reference to element within the chunk
                return (*it)[index - j * component_size];
            }
            else {
                // Handle negative indices using nvectors
                // Convert negative index to positive offset: -1 → 0, -2 → 1, etc.
                index = -1 - index;
                
                if (index >= (*this).nsize()) {
                    // Calculate number of additional chunks needed
                    int x = (index - (*this).nsize()) / component_size;
                    nnvectors += x + 1;
                    
                    // Allocate new vector chunks to accommodate the index
                    for (int i = 0; i <= x; i++) {
                        nvectors.push_back(std::vector<T>(component_size, default_));
                    }
                }
                
                // Navigate to the appropriate vector chunk
                auto it = nvectors.begin();
                int j;
                for (j = 0; j < index / component_size; j++) {
                    ++it;
                }
                
                // Return reference to element within the chunk
                return (*it)[index - j * component_size];
            }
        }
        /**
         * @brief Get number of allocated negative vector chunks
         * @return Current count of nvectors list elements
         */
        int get_nnvectors() {
            return nnvectors;
        }
        
        /**
         * @brief Get number of allocated positive vector chunks  
         * @return Current count of pvectors list elements
         */
        int get_npvectors() {
            return npvectors;
        }
        
        /**
         * @brief Get most negative index that has been accessed
         * @return Minimum index value accessed via operator[] (for profiling)
         */
        int get_max_nindex() {
            return max_nindex;
        }
        
        /**
         * @brief Get most positive index that has been accessed
         * @return Maximum index value accessed via operator[] (for profiling)
         */
        int get_max_pindex() {
            return max_pindex;
        }
        
        /**
         * @brief Get the fixed size of each vector chunk
         * @return Component size used for all vector chunks in the container
         */
        int get_component_size() {
            return component_size;
        }
};