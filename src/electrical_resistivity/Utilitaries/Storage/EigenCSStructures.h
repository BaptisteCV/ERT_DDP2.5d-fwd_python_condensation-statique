#ifndef EIGENCSSTRUCTURES_H_
#define EIGENCSSTRUCTURES_H_

#ifdef G
#pragma push_macro("G")
#undef G
#define EIGEN_CS_RESTORE_G
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <map>
#include <utility>

namespace eigen_cs {

using DenseMatrix = Eigen::MatrixXd;
using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
using Vector = Eigen::VectorXd;
using Triplet = Eigen::Triplet<double>;

inline DenseMatrix InitializationDenseMatrix(int rows, int cols, double value = 0.0) {
    return DenseMatrix::Constant(rows, cols, value);
}

inline Vector InitializationVector(int size, double value = 0.0) {
    return Vector::Constant(size, value);
}

inline SparseMatrix InitializationSparseMatrix(int rows, int cols) {
    SparseMatrix matrix(rows, cols);
    matrix.setZero();
    return matrix;
}

inline int return_index(int i, int j, int Ny) {
    return i * Ny + j;
}

inline std::pair<int, int> return_indices(int index, int Ny) {
    return std::make_pair(index / Ny, index % Ny);
}

}  // namespace eigen_cs

#ifdef EIGEN_CS_RESTORE_G
#pragma pop_macro("G")
#undef EIGEN_CS_RESTORE_G
#endif

#endif /* EIGENCSSTRUCTURES_H_ */