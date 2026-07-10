/*
 * MatrixConversion.h
 *
 * Utilities for converting ublas matrices to UMFPACK CSC format
 */

#ifndef MATRIXCONVERSION_H_
#define MATRIXCONVERSION_H_

#include <vector>
#include <algorithm>

namespace matrix_convert {

/*
 * Convert ublas::compressed_matrix to UMFPACK CSC format (Ap, Ai, Ax)
 * 
 * The ublas compressed_matrix is already in CSC format internally,
 * so we just extract the raw pointers
 */
template<typename MatrixType>
void ublas_to_csc(
    const MatrixType& mat,
    std::vector<int>& Ap,
    std::vector<int>& Ai,
    std::vector<double>& Ax)
{
    int n_col = mat.size2();
    int n_row = mat.size1();
    int nnz = mat.nnz();
    
    // Allocate arrays
    Ap.resize(n_col + 1);
    Ai.resize(nnz);
    Ax.resize(nnz);
    
    // Extract CSC format directly from ublas compressed_matrix
    const int* data_ip1 = mat.index1_data().begin();
    const int* data_ip2 = mat.index2_data().begin();
    const double* data_val = mat.value_data().begin();
    
    // Copy column pointers
    std::copy(data_ip1, data_ip1 + (n_col + 1), Ap.begin());
    
    // Copy row indices and values
    std::copy(data_ip2, data_ip2 + nnz, Ai.begin());
    std::copy(data_val, data_val + nnz, Ax.begin());
}

} // namespace matrix_convert

#endif /* MATRIXCONVERSION_H_ */
