#ifndef LINEARSYSTEMEIGEN_H_
#define LINEARSYSTEMEIGEN_H_

#include "../Storage/EigenCSStructures.h"
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

namespace eigen_cs {

struct LinearSystemSparseFactorization {
    Eigen::SparseLU<SparseMatrix> solver;
    bool ready;

    LinearSystemSparseFactorization() : ready(false) {}
};

struct LinearSystemSparseLLTFactorization {
    Eigen::SimplicialLLT<SparseMatrix> solver;
    bool ready;

    LinearSystemSparseLLTFactorization() : ready(false) {}
};

bool LinearSystemFactorize(const SparseMatrix& matrix, LinearSystemSparseFactorization& factorization);
Vector LinearSystemSolveFactorized(const LinearSystemSparseFactorization& factorization, const Vector& rhs);
bool LinearSystemFactorizeLLT(const SparseMatrix& matrix, LinearSystemSparseLLTFactorization& factorization);
Vector LinearSystemSolveFactorizedLLT(const LinearSystemSparseLLTFactorization& factorization, const Vector& rhs);
Vector LinearSystemSolving(const SparseMatrix& matrix, const Vector& rhs);
Vector LinearSystemSolvingDense(DenseMatrix matrix, Vector rhs);

}  // namespace eigen_cs

#endif /* LINEARSYSTEMEIGEN_H_ */