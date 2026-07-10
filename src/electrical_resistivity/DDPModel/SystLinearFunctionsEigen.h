#ifndef SYSTLINEARFUNCTIONSEIGEN_H_
#define SYSTLINEARFUNCTIONSEIGEN_H_

#include "../Utilitaries/Storage/EigenCSStructures.h"
#include "../Utilitaries/Storage/UblasStructures.h"
#include "../InputOutput/BoundaryConditions.h"
#include "../DFNModel/DFNComputation.h"
#include "DDPModel.h"

// ── Sélection du solveur pour DDP_EIGEN ──────────────────────────────────────
// Décommenter UNE ligne :
#define DDP_EIGEN_SOLVER_LLT       // SimplicialLLT  — Cholesky simplicial (SPD requis)
// #define DDP_EIGEN_SOLVER_UMFPACK   // UmfPackLU      — LU multifrontal (UMFPACK, rapide)
// ─────────────────────────────────────────────────────────────────────────────

#if defined(DDP_EIGEN_SOLVER_LLT)
#  include <Eigen/SparseCholesky>
   using DDPEigenSolver = Eigen::SimplicialLLT<eigen_cs::SparseMatrix>;
#elif defined(DDP_EIGEN_SOLVER_UMFPACK)
#  include <Eigen/UmfPackSupport>
   using DDPEigenSolver = Eigen::UmfPackLU<eigen_cs::SparseMatrix>;
#else
#  error "SystLinearFunctionsEigen.h : aucun solveur sélectionné — définir DDP_EIGEN_SOLVER_LLT ou DDP_EIGEN_SOLVER_UMFPACK"
#endif

// ─────────────────────────────────────────────────────────────────────────────
// DDPEigenBase : données du système de base (indépendantes de w et de la source)
// ─────────────────────────────────────────────────────────────────────────────
struct DDPEigenBase {
    eigen_cs::SparseMatrix matrix_base; // stencil poreux + pattern DFN (valeurs DFN = 0)
    int n;
    int nb_nodes_fract, Nx, Ny;
    double delta_x, delta_y;
    Eigen::MatrixXd porous_cond;           // copie dense (O(1) accès)
    SubNetworkMap* subnetwork_map_ptr;     // pointeur — pas de copie du SubNetworkMap
};

// ─────────────────────────────────────────────────────────────────────────────
// Construit le système de base Eigen. Prend DDPModel par référence non-const
// pour pouvoir stocker un pointeur non-const vers son SubNetworkMap.
// ─────────────────────────────────────────────────────────────────────────────
DDPEigenBase DetermineBasicLinearSystemDDPEigen(DDPModel& m);

// ─────────────────────────────────────────────────────────────────────────────
// Remplit les termes w-dépendants. matrix_syst doit être une copie fraîche de
// base.matrix_base. vector_syst doit être initialisé à zéro.
// ─────────────────────────────────────────────────────────────────────────────
void UpdateLinearSystemFourierEigen(
    DDPEigenBase& base,
    double current_w,
    const BoundaryConditionsEPM& bc_def_epm,
    const BoundaryConditionsDFN& bc_map,
    const SourceTermsEPM& source_terms_epm,
    bool singularity,
    eigen_cs::SparseMatrix& matrix_syst,
    eigen_cs::Vector& vector_syst);

// ─────────────────────────────────────────────────────────────────────────────
// Résolution en réutilisant le pattern symbolique pré-calculé.
// analyzePattern doit avoir été appelé une fois sur le solver avant le premier appel.
// Le type du solver dépend du #define DDP_EIGEN_SOLVER_* actif.
// Sorties 100% Eigen : DFNPotential (vecteur dense), PorousPotential (matrice dense Nx×Ny).
// ─────────────────────────────────────────────────────────────────────────────
void EvaluatePotentialEigen(
    int nb_nodes_fract, int Nx, int Ny,
    DDPEigenSolver& solver,
    eigen_cs::SparseMatrix& matrix_syst,
    const eigen_cs::Vector& vector_syst,
    eigen_cs::Vector& DFNPotential,
    eigen_cs::DenseMatrix& PorousPotential);

#endif /* SYSTLINEARFUNCTIONSEIGEN_H_ */
