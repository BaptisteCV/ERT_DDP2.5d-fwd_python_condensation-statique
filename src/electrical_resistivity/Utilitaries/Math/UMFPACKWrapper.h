/*
 * UMFPACKWrapper.h
 *
 * C++ wrapper for UMFPACK C API (SuiteSparse 5.10.1)
 * Replaces boost::numeric::bindings::umfpack to remove external dependency
 * Uses native UMFPACK C interface directly
 */

#ifndef UMFPACKWRAPPER_H_
#define UMFPACKWRAPPER_H_

#include <umfpack.h>
#include <cstring>
#include <stdexcept>

namespace umfpack_wrapper {

// ========================================================================
// Opaque types for symbolic and numeric factorizations
// ========================================================================

class SymbolicFactorization {
public:
    SymbolicFactorization() : symbolic_(nullptr) {}
    
    ~SymbolicFactorization() {
        if (symbolic_ != nullptr) {
            umfpack_di_free_symbolic(&symbolic_);
        }
    }
    
    void* get() const { return symbolic_; }
    void** get_ptr() { return &symbolic_; }
    
private:
    void* symbolic_;
    
    // Prevent copying
    SymbolicFactorization(const SymbolicFactorization&) = delete;
    SymbolicFactorization& operator=(const SymbolicFactorization&) = delete;
};

class NumericFactorization {
public:
    NumericFactorization() : numeric_(nullptr) {}
    
    ~NumericFactorization() {
        if (numeric_ != nullptr) {
            umfpack_di_free_numeric(&numeric_);
        }
    }
    
    void* get() const { return numeric_; }
    void** get_ptr() { return &numeric_; }
    
private:
    void* numeric_;
    
    // Prevent copying
    NumericFactorization(const NumericFactorization&) = delete;
    NumericFactorization& operator=(const NumericFactorization&) = delete;
};

// ========================================================================
// Symbolic factorization: compute the column-form matrix directly from ublas
// ========================================================================

template<typename MatrixType>
int symbolic(
    int n_row,
    int n_col,
    const int* Ap,
    const int* Ai,
    const double* Ax,
    SymbolicFactorization& symbolic,
    const double* control = nullptr,
    double* info = nullptr)
{
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    
    // Use default control if not provided
    if (control == nullptr) {
        umfpack_di_defaults(Control);
    } else {
        std::memcpy(Control, control, UMFPACK_CONTROL * sizeof(double));
    }
    
    // Initialize info
    if (info == nullptr) {
        info = Info;
    }
    
    int status = umfpack_di_symbolic(
        n_row, n_col, Ap, Ai, Ax,
        symbolic.get_ptr(), Control, info);
    
    return status;
}

// ========================================================================
// Numeric factorization: factor matrix numerically
// ========================================================================

template<typename MatrixType>
int numeric(
    const int* Ap,
    const int* Ai,
    const double* Ax,
    SymbolicFactorization& symbolic,
    NumericFactorization& numeric,
    const double* control = nullptr,
    double* info = nullptr)
{
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    
    // Use default control if not provided
    if (control == nullptr) {
        umfpack_di_defaults(Control);
    } else {
        std::memcpy(Control, control, UMFPACK_CONTROL * sizeof(double));
    }
    
    // Initialize info
    if (info == nullptr) {
        info = Info;
    }
    
    int status = umfpack_di_numeric(
        Ap, Ai, Ax,
        symbolic.get(), numeric.get_ptr(),
        Control, info);
    
    return status;
}

// ========================================================================
// Solve: solve Ax = b using factorized matrix
// ========================================================================

template<typename VectorIn, typename VectorOut>
int solve(
    const int* Ap,
    const int* Ai,
    const double* Ax,
    double* x,
    const double* b,
    NumericFactorization& numeric,
    int sys = UMFPACK_A,
    const double* control = nullptr,
    double* info = nullptr)
{
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    
    // Use default control if not provided
    if (control == nullptr) {
        umfpack_di_defaults(Control);
    } else {
        std::memcpy(Control, control, UMFPACK_CONTROL * sizeof(double));
    }
    
    // Initialize info
    if (info == nullptr) {
        info = Info;
    }
    
    int status = umfpack_di_solve(
        sys, Ap, Ai, Ax,
        x, b,
        numeric.get(), Control, info);
    
    return status;
}

} // namespace umfpack_wrapper

#endif /* UMFPACKWRAPPER_H_ */
