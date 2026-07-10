#ifndef EPMDEFINITIONCSEIGEN_H_
#define EPMDEFINITIONCSEIGEN_H_

#include "../Utilitaries/Storage/EigenCSStructures.h"
#include "../Utilitaries/Math/LinearSystemEigen.h"
#include "EPMComputationCSEigen.h"
#include "../DomainDefinition/Domain.h"
#include "../Utilitaries/Storage/UtilitaryStorage.h"
#include <map>
#include <vector>
#include <string>

namespace eigen_cs {

class EPMSystemCSEigen {
public:
    int Nx;
    int Ny;
    double delta_x;
    double delta_y;
    double w;

    DenseMatrix porous_cond;
    BoundaryConditionsCSEigen bc_def_epm;
    SourceTermsCSEigen source_terms_epm;
    DenseMatrix potential;

    SparseMatrix block_Aii;
    SparseMatrix block_Abb;
    SparseMatrix block_Aib;
    SparseMatrix block_Abi;

    DenseMatrix matrix_Y;
    DenseMatrix matrix_S;
    DenseMatrix matrix_S_precomputed;
    Vector vector_delta_diag;

    Vector vector_z;
    Vector vector_rhs;
    Vector vector_phi_i;
    Vector vector_phi_b;
    Vector vector_bi;
    Vector vector_bb;

    LinearSystemSparseFactorization factor_Aii;
    bool factor_Aii_ready;

    LinearSystemSparseLLTFactorization factor_Aii_llt;
    bool factor_Aii_llt_ready;

    std::vector<int> interior_global_indices;
    std::vector<int> border_global_indices;
    std::map<int, int> global_to_interior;
    std::map<int, int> global_to_border;
    int Ni;
    int Nb;

public:
    EPMSystemCSEigen();
    EPMSystemCSEigen(int Nx_, int Ny_, double Lx, double Ly, const DenseMatrix& conductivity);
    EPMSystemCSEigen(const EPMSystemCSEigen& other);
    EPMSystemCSEigen& operator=(const EPMSystemCSEigen& other);

    void SetConductivity(const DenseMatrix& conductivity);
    void InitializeSCIndexing();
    void DetermineBasicLinearSystemCS();
    void UpdateSCBlocksFourier();
    void UpdateSCBoundaryConditionsAndSources(bool fourier);
    void PrecomputeSCForCurrentFrequency();
    double ComputeCondensedSystemResidualNorm() const;
    void SolveCondensedSystem();
    void ReconstructInteriorFromBorder();
    Vector BuildGlobalSolutionVectorFromSC() const;
    DenseMatrix PotentialComputationSC();
    // ExportFirstRowOfAiiForFrequency is commented out (A_ii export optional)
    // void ExportFirstRowOfAiiForFrequency(const std::string& output_dir, double current_w);
    void ExportFirstRowsToFile(const std::string& output_dir, int elec_nb, double current_w);
};

}  // namespace eigen_cs

#endif /* EPMDEFINITIONCSEIGEN_H_ */