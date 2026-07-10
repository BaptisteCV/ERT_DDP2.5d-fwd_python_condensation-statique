/*
 * LinearSystem.h
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet
 */

#ifndef LINEARSYSTEM_H_
#define LINEARSYSTEM_H_

#include "../Storage/UblasStructures.h"
#include "UMFPACKWrapper.h"
#include "MatrixConversion.h"
#include <vector>

// ===================== BEGIN ADDED (CS optimization) =====================
namespace ublas = boost::numeric::ublas;

// Use native UMFPACK C API wrapper instead of boost::numeric::bindings
typedef ublas::compressed_matrix<double, ublas::column_major, 0,
    ublas::unbounded_array<int>, ublas::unbounded_array<double> > ublas_matrix_syst;

struct LinearSystemUMFFactorization {
	ublas_matrix_syst mat;
	umfpack_wrapper::SymbolicFactorization symbolic;
	umfpack_wrapper::NumericFactorization numeric;
	std::vector<int> Ap, Ai;
	std::vector<double> Ax;
	bool ready;

	LinearSystemUMFFactorization() : ready(false) {}
};

ublas_vector LinearSystemSolving(ublas_matrix, ublas_vector);
ublas_vector LinearSystemSolvingDense(ublas_dense_matrix, ublas_vector);
void LinearSystemFactorize(ublas_matrix, LinearSystemUMFFactorization&);
ublas_vector LinearSystemSolveFactorized(const LinearSystemUMFFactorization&, ublas_vector);
// ====================== END ADDED (CS optimization) ======================


#endif /* LINEARSYSTEM_H_ */
