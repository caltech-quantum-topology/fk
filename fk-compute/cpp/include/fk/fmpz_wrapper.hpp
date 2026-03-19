#pragma once

#include <flint/fmpz.h>

/**
 * Minimal RAII wrapper for fmpz_t to allow use in std::vector and std::map
 * with correct copy/move semantics.
 */
struct fmpz_wrapper {
    fmpz_t val;

    fmpz_wrapper()                      { fmpz_init(val); }
    explicit fmpz_wrapper(slong v)      { fmpz_init_set_si(val, v); }
    fmpz_wrapper(const fmpz_wrapper& o) { fmpz_init_set(val, o.val); }
    fmpz_wrapper(fmpz_wrapper&& o)      { fmpz_init(val); fmpz_swap(val, o.val); }

    fmpz_wrapper& operator=(const fmpz_wrapper& o) { fmpz_set(val, o.val); return *this; }
    fmpz_wrapper& operator=(slong v)               { fmpz_set_si(val, v); return *this; }

    ~fmpz_wrapper() { fmpz_clear(val); }

    bool is_zero() const               { return fmpz_is_zero(val); }
    void add(const fmpz_wrapper& o)    { fmpz_add(val, val, o.val); }
    void sub(const fmpz_wrapper& o)    { fmpz_sub(val, val, o.val); }
    void negate()                      { fmpz_neg(val, val); }
};
