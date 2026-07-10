/*
 * UblasStructures.h
 *
 *  Created on: Apr 23, 2013
 *      Author: delphineroubinet
 */

#ifndef UBLASSTRUCTURES_H_
#define UBLASSTRUCTURES_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <map>

// Removed: boost::numeric::bindings is no longer needed 
// Now using native UMFPACK C API via UMFPACKWrapper.h

typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::column_major, 0,boost::numeric::ublas::unbounded_array<int>, boost::numeric::ublas::unbounded_array<double> > ublas_matrix;
//typedef boost::numeric::ublas::matrix<double> ublas_matrix;
typedef boost::numeric::ublas::matrix<double> ublas_dense_matrix;
typedef boost::numeric::ublas::vector<double> ublas_vector;

void print_matrix(ublas_matrix,int option=0);
void print_vector(ublas_vector);
ublas_vector InitializationVector(int);
ublas_vector InitializationVector(int,double);
ublas_matrix InitializationMatrix(int,int);
ublas_matrix InitializationMatrix(int,int,double);
ublas_matrix RelativeMatrix(ublas_matrix,ublas_matrix);

// in the first input matrix, the indicated line (first input is the line id) is replaced by a matrix block (second input matrix)
void LinesReplacementMatrix1(int,ublas_matrix&,ublas_matrix);
void LinesReplacementMatrix2(int,ublas_matrix&,ublas_matrix);
// in the first input vector, the indicated line (first input is the line id) is replaced by a vector block (second input vector)
void LinesReplacementVector(int,ublas_vector&,ublas_vector);

double min_vector(ublas_vector);
double max_vector(ublas_vector);
double min_matrix(ublas_matrix);
double max_matrix(ublas_matrix);
double ErrorComputation(ublas_matrix,ublas_matrix);
double ErrorComputationL2(ublas_matrix,ublas_matrix);
bool PositiveVector(ublas_vector);
bool PositiveMatrix(ublas_matrix);
ublas_matrix ExtractMatrix(ublas_matrix,int,int,int,int);
double AveragedMatrix(ublas_matrix);
std::map<int,double> MatrixToMap(ublas_matrix,int);
bool ConstantValue(ublas_matrix,double&);

#endif /* UBLASSTRUCTURES_H_ */
