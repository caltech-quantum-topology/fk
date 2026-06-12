/**
 * @file string_to_numeric.hpp
 * @brief String to numeric conversion utilities for Gukov-Manolescu invariant computation
 * 
 * This file provides efficient string-to-number conversion functions for parsing
 * numerical data from CSV files and input streams. The functions handle both
 * integer and floating-point conversions with proper sign handling.
 * 
 * Note: These are custom implementations that may be more efficient than std::stoi/std::stod
 * for simple numeric strings without scientific notation or advanced formatting.
 */

#pragma once

#include <string>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <limits>

/**
 * @brief Convert string to integer
 * @param s Input string containing integer representation
 * @return Parsed integer value
 * 
 * Supports:
 * - Positive and negative integers
 * - Leading minus sign for negative numbers
 * 
 * Note: No validation - assumes input is well-formed integer string
 */
int string_to_int(const std::string& s) {
    if (s.empty()) return 0;
    
    int acc = 0;
    bool negate = false;
    size_t start = 0;
    
    // Handle negative sign
    if (s[0] == '-') {
        negate = true;
        start = 1;
    }
    
    // Parse digits from left to right
    for (size_t i = start; i < s.size(); i++) {
        acc += (s[i] - '0') * static_cast<int>(std::pow(10, s.size() - i - 1));
    }
    
    return negate ? -acc : acc;
}

/**
 * @brief Convert string to double-precision floating point
 * @param s Input string containing floating-point representation
 * @return Parsed double value
 * 
 * Supports:
 * - Positive and negative numbers
 * - Decimal point notation (e.g., "3.14159")
 * - Integer strings (treated as whole numbers)
 * 
 * Note: No validation - assumes input is well-formed numeric string
 */
double string_to_double(const std::string& s) {
    if (s.empty()) return 0.0;
    
    double acc = 0.0;
    bool negate = false;
    size_t start = 0;
    
    // Handle negative sign
    if (s[0] == '-') {
        negate = true;
        start = 1;
    }
    
    // Find decimal point position
    size_t decimal_pos = s.find('.');
    if (decimal_pos == std::string::npos) {
        decimal_pos = s.size();  // No decimal point found
    }
    
    // Parse integer part
    for (size_t i = start; i < decimal_pos; i++) {
        acc += (s[i] - '0') * std::pow(10, decimal_pos - i - 1);
    }
    
    // Parse fractional part
    for (size_t i = decimal_pos + 1; i < s.size(); i++) {
        acc += (s[i] - '0') * std::pow(10, -(static_cast<int>(i - decimal_pos)));
    }
    
    return negate ? -acc : acc;
}

/**
 * @brief Generic template function for string to numeric conversion
 * @tparam T Target numeric type (int, long, float, double, etc.)
 * @param s Input string containing numeric representation
 * @return Parsed value of type T
 * 
 * Supports all fundamental numeric types including:
 * - Integer types: int, long, long long, unsigned variants
 * - Floating-point types: float, double, long double
 * - Automatic sign detection and proper type conversion
 * 
 * Note: Uses SFINAE to enable appropriate conversion logic based on type
 */
template<typename T>
typename std::enable_if<std::is_integral<T>::value, T>::type
string_to_numeric(const std::string& s) {
    if (s.empty()) return T(0);
    
    T acc = 0;
    bool negate = false;
    size_t start = 0;
    
    // Handle negative sign
    if (s[0] == '-') {
        negate = true;
        start = 1;
    }
    
    // Parse digits from left to right
    for (size_t i = start; i < s.size(); i++) {
        acc += (s[i] - '0') * static_cast<T>(std::pow(10, s.size() - i - 1));
    }
    
    return negate ? -acc : acc;
}

template<typename T>
typename std::enable_if<std::is_floating_point<T>::value, T>::type
string_to_numeric(const std::string& s) {
    if (s.empty()) return T(0.0);
    
    T acc = T(0.0);
    bool negate = false;
    size_t start = 0;
    
    // Handle negative sign
    if (s[0] == '-') {
        negate = true;
        start = 1;
    }
    
    // Find decimal point position
    size_t decimal_pos = s.find('.');
    if (decimal_pos == std::string::npos) {
        decimal_pos = s.size();  // No decimal point found
    }
    
    // Parse integer part
    for (size_t i = start; i < decimal_pos; i++) {
        acc += (s[i] - '0') * std::pow(T(10), decimal_pos - i - 1);
    }
    
    // Parse fractional part
    for (size_t i = decimal_pos + 1; i < s.size(); i++) {
        acc += (s[i] - '0') * std::pow(T(10), -(static_cast<int>(i - decimal_pos)));
    }
    
    return negate ? -acc : acc;
}

/**
 * @brief Convenience wrapper functions for specific numeric types
 */

// Integer type wrappers
inline long string_to_long(const std::string& s) {
    return string_to_numeric<long>(s);
}

inline long long string_to_long_long(const std::string& s) {
    return string_to_numeric<long long>(s);
}

inline unsigned int string_to_unsigned_int(const std::string& s) {
    return string_to_numeric<unsigned int>(s);
}

inline unsigned long string_to_unsigned_long(const std::string& s) {
    return string_to_numeric<unsigned long>(s);
}

// Floating-point type wrappers
inline float string_to_float(const std::string& s) {
    return string_to_numeric<float>(s);
}

inline long double string_to_long_double(const std::string& s) {
    return string_to_numeric<long double>(s);
}
