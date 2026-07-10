/*
 * LinearSystem.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet
 */

// ============================================================
// SOLVER SELECTION — uncomment exactly ONE of the following:
// ============================================================
#define USE_UMFPACK   // LU direct (non-symétrique)   — SuiteSparse
// #define USE_CHOLMOD  // Cholesky (matrices SPD only) — SuiteSparse
// #define USE_EIGEN    // solveurs Eigen (conversion ublas <-> Eigen)
// ============================================================

#include "LinearSystem.h"
#include <iostream>
#include <cstring>
#include <boost/numeric/ublas/lu.hpp>
#include <umfpack.h>   // toujours requis pour LinearSystemFactorize/SolveFactorized

#if defined(USE_CHOLMOD)
#  include <cholmod.h>
#elif defined(USE_EIGEN)
#  include <Eigen/Sparse>
#  include <Eigen/SparseLU>
#  include <Eigen/SparseCholesky>
#  include <Eigen/IterativeLinearSolvers>
#endif

namespace ublas = boost::numeric::ublas;
using namespace std;

// ===================== BEGIN ADDED (CS optimization) =====================
// Ces fonctions utilisent toujours UMFPACK (pré-factorisation réutilisable pour le chemin CS).
void LinearSystemFactorize(ublas_matrix mat_syst, LinearSystemUMFFactorization& fac){
	fac.mat = mat_syst;

	int n_col = fac.mat.size2();
	int nnz = fac.mat.nnz();

	fac.Ap.resize(n_col + 1);
	fac.Ai.resize(nnz);
	fac.Ax.resize(nnz);

	const int* data_ip1 = fac.mat.index1_data().begin();
	const int* data_ip2 = fac.mat.index2_data().begin();
	const double* data_val = fac.mat.value_data().begin();

	std::copy(data_ip1, data_ip1 + (n_col + 1), fac.Ap.begin());
	std::copy(data_ip2, data_ip2 + nnz, fac.Ai.begin());
	std::copy(data_val, data_val + nnz, fac.Ax.begin());

	double Control[UMFPACK_CONTROL];
	double Info[UMFPACK_INFO];
	umfpack_di_defaults(Control);

	int status = umfpack_di_symbolic(
		fac.mat.size1(), fac.mat.size2(),
		fac.Ap.data(), fac.Ai.data(), fac.Ax.data(),
		fac.symbolic.get_ptr(), Control, Info);
	if (status < 0) {
		cerr << "ERROR in LinearSystemFactorize: umfpack_di_symbolic failed with status " << status << endl;
		fac.ready = false;
		return;
	}

	status = umfpack_di_numeric(
		fac.Ap.data(), fac.Ai.data(), fac.Ax.data(),
		fac.symbolic.get(), fac.numeric.get_ptr(),
		Control, Info);
	if (status < 0) {
		cerr << "ERROR in LinearSystemFactorize: umfpack_di_numeric failed with status " << status << endl;
		fac.ready = false;
		return;
	}

	fac.ready = true;
}

ublas_vector LinearSystemSolveFactorized(const LinearSystemUMFFactorization& fac, ublas_vector vect_syst){
	int Ny = vect_syst.size();
	ublas::vector<double> X(Ny);

	if (!fac.ready) {
		cout << "WARNING in LinearSystemSolveFactorized (LinearSystem.cpp): factorization not ready" << endl;
		return X;
	}

	if ((unsigned)fac.mat.size2() == vect_syst.size()) {
		double Control[UMFPACK_CONTROL];
		double Info[UMFPACK_INFO];
		umfpack_di_defaults(Control);
		int status = umfpack_di_solve(
			UMFPACK_A,
			fac.Ap.data(), fac.Ai.data(), fac.Ax.data(),
			X.data().begin(), vect_syst.data().begin(),
			fac.numeric.get(), Control, Info);
		if (status < 0) {
			cerr << "ERROR in LinearSystemSolveFactorized: umfpack_di_solve failed with status " << status << endl;
		}
	} else {
		cout << "WARNING in LinearSystemSolveFactorized (LinearSystem.cpp): system not well defined" << endl;
	}

	return X;
}
// ====================== END ADDED (CS optimization) ======================


// Résout Ax = b en une passe (factorisation + solve).
// Le solveur actif est contrôlé par le #define en tête de fichier.
ublas_vector LinearSystemSolving(ublas_matrix mat_syst, ublas_vector vect_syst){

	int n = (int)vect_syst.size();
	ublas::vector<double> X(n, 0.0);

	if ((unsigned)mat_syst.size2() != vect_syst.size()) {
		cout << "WARNING in LinearSystemSolving: system not well defined" << endl;
		return X;
	}

	// Conversion CSC commune à tous les solveurs creux
	int nnz = (int)mat_syst.nnz();
	std::vector<int>    Ap(n + 1), Ai(nnz);
	std::vector<double> Ax(nnz);
	std::copy(mat_syst.index1_data().begin(), mat_syst.index1_data().begin() + (n + 1), Ap.begin());
	std::copy(mat_syst.index2_data().begin(), mat_syst.index2_data().begin() + nnz,     Ai.begin());
	std::copy(mat_syst.value_data().begin(),  mat_syst.value_data().begin()  + nnz,     Ax.begin());

#ifdef USE_UMFPACK
	// --- UMFPACK : LU direct non-symétrique ---
	void* symbolic = nullptr;
	void* numeric  = nullptr;
	double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
	umfpack_di_defaults(Control);
	umfpack_di_symbolic(n, n, Ap.data(), Ai.data(), Ax.data(), &symbolic, Control, Info);
	umfpack_di_numeric (Ap.data(), Ai.data(), Ax.data(), symbolic, &numeric, Control, Info);
	umfpack_di_solve(UMFPACK_A, Ap.data(), Ai.data(), Ax.data(),
	                 X.data().begin(), vect_syst.data().begin(), numeric, Control, Info);
	umfpack_di_free_symbolic(&symbolic);
	umfpack_di_free_numeric (&numeric);

#elif defined(USE_CHOLMOD)
	// --- Cholmod : Cholesky (matrice SPD obligatoire) ---
	cholmod_common c;
	cholmod_start(&c);
	cholmod_sparse* A = cholmod_allocate_sparse(n, n, nnz, true, true, -1, CHOLMOD_REAL, &c);
	std::copy(Ap.begin(), Ap.end(), (int*)A->p);
	std::copy(Ai.begin(), Ai.end(), (int*)A->i);
	std::copy(Ax.begin(), Ax.end(), (double*)A->x);
	cholmod_dense* b_chol = cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, &c);
	std::copy(vect_syst.begin(), vect_syst.end(), (double*)b_chol->x);
	cholmod_factor* L = cholmod_analyze(A, &c);
	cholmod_factorize(A, L, &c);
	cholmod_dense* x_chol = cholmod_solve(CHOLMOD_A, L, b_chol, &c);
	std::copy((double*)x_chol->x, (double*)x_chol->x + n, X.begin());
	cholmod_free_sparse(&A, &c);
	cholmod_free_dense(&b_chol, &c);
	cholmod_free_dense(&x_chol, &c);
	cholmod_free_factor(&L, &c);
	cholmod_finish(&c);

#elif defined(USE_EIGEN)
	// --- Eigen : conversion ublas CSC -> Eigen -> résolution -> ublas ---
	//
	// Conversion ublas -> Eigen SparseMatrix (format CSC déjà disponible)
	Eigen::SparseMatrix<double> A_eigen(n, n);
	A_eigen = Eigen::Map<const Eigen::SparseMatrix<double>>(
	    n, n, nnz, Ap.data(), Ai.data(), Ax.data());
	A_eigen.makeCompressed();

	// Conversion ublas_vector -> Eigen VectorXd
	Eigen::VectorXd b_eigen = Eigen::Map<const Eigen::VectorXd>(
	    vect_syst.data().begin(), n);

	// Résolution — choisir UNE ligne de déclaration du solveur :
	// Eigen::SparseLU<Eigen::SparseMatrix<double>>             solver; // LU direct       (non-symétrique)
	// Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>     solver; // Cholesky LLT    (SPD)
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>    solver; // Cholesky LDLT   (SPD, plus stable que LLT)
	// Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver; // CG itératif     (SPD, grandes matrices)
	// Eigen::BiCGSTAB<Eigen::SparseMatrix<double>>          solver; // BiCGSTAB itér.  (non-symétrique)

	// compute() = analyzePattern + factorize pour les solveurs directs,
	//            initialisation du préconditionneur pour les itératifs.
	clock_t start = clock();
	solver.compute(A_eigen);
	clock_t end = clock() - start;
	cout << "Eigen compute time: " << (double)end / CLOCKS_PER_SEC << " seconds" << endl;
	if (solver.info() != Eigen::Success) {
		cout << "ERROR in LinearSystemSolving: Eigen compute failed" << endl;
		return X;
	}
	Eigen::VectorXd x_eigen = solver.solve(b_eigen);
	if (solver.info() != Eigen::Success) {
		cout << "ERROR in LinearSystemSolving: Eigen solve failed" << endl;
		return X;
	}

	// Conversion Eigen VectorXd -> ublas_vector
	std::copy(x_eigen.data(), x_eigen.data() + n, X.begin());

#else
	#error "LinearSystem.cpp: no solver selected — define USE_UMFPACK, USE_CHOLMOD or USE_EIGEN"
#endif

	return X;
}


ublas_vector LinearSystemSolvingDense(ublas_dense_matrix mat_syst, ublas_vector vect_syst){
	int n = vect_syst.size();
	ublas::vector<double> X(n);

	if ((unsigned)mat_syst.size1() == vect_syst.size() && mat_syst.size1() == mat_syst.size2()){
		ublas::permutation_matrix<std::size_t> pm(mat_syst.size1());
		int res = boost::numeric::ublas::lu_factorize(mat_syst, pm);
		if (res != 0) {
			cout << "WARNING in LinearSystemSolvingDense (LinearSystem.cpp): factorization failed" << endl;
			return X;
		}
		X = vect_syst;
		boost::numeric::ublas::lu_substitute(mat_syst, pm, X);
	} else {
		cout << "WARNING in LinearSystemSolvingDense (LinearSystem.cpp): system not well defined" << endl;
	}

	return X;
}
