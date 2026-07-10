#include "LinearSystemEigen.h"
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <iostream>

using namespace std;

namespace eigen_cs {

bool LinearSystemFactorize(const SparseMatrix& matrix, LinearSystemSparseFactorization& factorization) {
    factorization.solver.analyzePattern(matrix);
    factorization.solver.factorize(matrix);
    factorization.ready = (factorization.solver.info() == Eigen::Success);
    if (!factorization.ready) {
        cout << "ERROR: Sparse factorization failed" << endl;
    }
    return factorization.ready;
}

Vector LinearSystemSolveFactorized(const LinearSystemSparseFactorization& factorization, const Vector& rhs) {
    return factorization.solver.solve(rhs);
}

bool LinearSystemFactorizeLLT(const SparseMatrix& matrix, LinearSystemSparseLLTFactorization& factorization) {
    factorization.solver.analyzePattern(matrix);
    factorization.solver.factorize(matrix);
    factorization.ready = (factorization.solver.info() == Eigen::Success);
    if (!factorization.ready) {
        cout << "ERROR: Sparse LLT factorization failed" << endl;
    }
    return factorization.ready;
}

Vector LinearSystemSolveFactorizedLLT(const LinearSystemSparseLLTFactorization& factorization, const Vector& rhs) {
    return factorization.solver.solve(rhs);
}

Vector LinearSystemSolving(const SparseMatrix& matrix, const Vector& rhs) {
    Eigen::SparseLU<SparseMatrix> solver;
    solver.analyzePattern(matrix);
    solver.factorize(matrix);
    if (solver.info() != Eigen::Success) {
        cout << "ERROR: Sparse solve failed" << endl;
        return Vector::Zero(rhs.size());
    }
    return solver.solve(rhs);
}

Vector LinearSystemSolvingDense(DenseMatrix matrix, Vector rhs) {
    // Eigen::PartialPivLU<DenseMatrix> solver(matrix);
    // cout << "Dense LU factorization succeeded" << endl;
    // return solver.solve(rhs);

    Eigen::LLT<DenseMatrix> solver(matrix);
    if (solver.info() != Eigen::Success) {
        cout << "ERROR: Dense LLT factorization failed" << endl;
        return Vector::Zero(rhs.size());
    }

    Vector solution = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
        cout << "ERROR: Dense LLT solve failed" << endl;
        return Vector::Zero(rhs.size());
    }
    // cout << "Dense LLT solve succeeded" << endl;

    return solution;
}

}  // namespace eigen_cs