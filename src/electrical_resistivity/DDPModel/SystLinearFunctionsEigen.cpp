#include "SystLinearFunctionsEigen.h"
#include "DDPDFNFunctions.h"
#include "ExchangeFunctions.h"
#include "../DFNModel/DFNComputation.h"
#include "../Utilitaries/Math/Operations.h"
#include "../EPMModel/EPMComputation.h"
#include <boost/math/special_functions/bessel.hpp>
#include <iostream>
#include <omp.h>

using namespace std;

using EigenSparse  = eigen_cs::SparseMatrix;
using EigenVec     = eigen_cs::Vector;
using EigenTriplet = eigen_cs::Triplet;

// Wrappers locaux sans ambiguïté (évite le conflit return_index global / eigen_cs)
static inline int ri(int i, int j, int Ny)  { return i * Ny + j; }
static inline pair<int,int> rij(int idx, int Ny) { return {idx / Ny, idx % Ny}; }

// ─────────────────────────────────────────────────────────────────────────────
// Construction du système de base Eigen
// ─────────────────────────────────────────────────────────────────────────────

DDPEigenBase DetermineBasicLinearSystemDDPEigen(DDPModel& m) {
    //double t_eigen_total = omp_get_wtime();
    DDPEigenBase base;
    base.nb_nodes_fract     = m.nb_nodes_fract;
    base.Nx                 = m.Nx;
    base.Ny                 = m.Ny;
    base.delta_x            = m.delta_x;
    base.delta_y            = m.delta_y;
    base.subnetwork_map_ptr = &m.subnetwork_map;   // pointeur — zéro copie
    base.n                  = m.nb_nodes_fract + m.Nx * m.Ny;

    // Conversion ublas → Eigen dense pour accès O(1)
    base.porous_cond.resize(m.Nx, m.Ny);
    for (int i = 0; i < m.Nx; i++)
        for (int j = 0; j < m.Ny; j++)
            base.porous_cond(i, j) = m.porous_cond(i, j);

    const int Nx = m.Nx, Ny = m.Ny, nf = m.nb_nodes_fract;
    const double dx = m.delta_x, dy = m.delta_y;

    // ── 1. Triplets : stencil poreux (valeurs réelles) + pattern DFN (valeurs 0) ──
    vector<EigenTriplet> triplets;
    triplets.reserve(base.n * 6);

    //double t_stencil = omp_get_wtime();
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = nf + ri(i, j, Ny);
            double diag = 0.0;
            if (i > 0) {
                double c = Harmonic_Mean(base.porous_cond(i,j), base.porous_cond(i-1,j)) * dy / dx;
                triplets.emplace_back(idx, nf + ri(i-1, j, Ny), -c);
                diag += c;
            }
            if (i < Nx-1) {
                double c = Harmonic_Mean(base.porous_cond(i,j), base.porous_cond(i+1,j)) * dy / dx;
                triplets.emplace_back(idx, nf + ri(i+1, j, Ny), -c);
                diag += c;
            }
            if (j > 0) {
                double c = Harmonic_Mean(base.porous_cond(i,j), base.porous_cond(i,j-1)) * dx / dy;
                triplets.emplace_back(idx, nf + ri(i, j-1, Ny), -c);
                diag += c;
            }
            if (j < Ny-1) {
                double c = Harmonic_Mean(base.porous_cond(i,j), base.porous_cond(i,j+1)) * dx / dy;
                triplets.emplace_back(idx, nf + ri(i, j+1, Ny), -c);
                diag += c;
            }
            triplets.emplace_back(idx, idx, diag);
        }
    }

    //cout << "[PROFILE][DetermineBasicLinearSystemDDPEigen] stencil poreux (triplets): "
    //     << omp_get_wtime() - t_stencil << " sec" << endl;

    //double t_dfn_pattern = omp_get_wtime();
    for (SubNetworkMap::const_iterator it = m.subnetwork_map.begin();
         it != m.subnetwork_map.end(); ++it) {
        int mat_idx = nf + it->first;
        for (const FractureMesh& frac : it->second.meshes) {
            int n1 = frac.p_ori.index, n2 = frac.p_tar.index;
            triplets.emplace_back(n1, n1, 0.0);
            triplets.emplace_back(n1, n2, 0.0);
            triplets.emplace_back(n2, n1, 0.0);
            triplets.emplace_back(n2, n2, 0.0);
            triplets.emplace_back(n1, mat_idx, 0.0);
            triplets.emplace_back(n2, mat_idx, 0.0);
            triplets.emplace_back(mat_idx, n1, 0.0);
            triplets.emplace_back(mat_idx, n2, 0.0);
        }
    }

    //cout << "[PROFILE][DetermineBasicLinearSystemDDPEigen] pattern DFN (triplets zero): "
    //     << omp_get_wtime() - t_dfn_pattern << " sec" << endl;

    //double t_compress = omp_get_wtime();
    base.matrix_base.resize(base.n, base.n);
    base.matrix_base.setFromTriplets(triplets.begin(), triplets.end());
    base.matrix_base.makeCompressed();
    //cout << "[PROFILE][DetermineBasicLinearSystemDDPEigen] setFromTriplets + makeCompressed: "
    //     << omp_get_wtime() - t_compress << " sec" << endl;
    //cout << "[PROFILE][DetermineBasicLinearSystemDDPEigen] total: "
    //     << omp_get_wtime() - t_eigen_total << " sec"
    //     << "  (n=" << base.n << " nnz=" << base.matrix_base.nonZeros() << ")" << endl;

    return base;
}

// ─────────────────────────────────────────────────────────────────────────────
// Sous-fonctions d'assemblage
// ─────────────────────────────────────────────────────────────────────────────

static void DFNDDPLinearSystemFourierEigen(
    const vector<FractureMesh>& meshes, int block_idx, int nf,
    double alpha, double w, EigenSparse& mat)
{
    double a1, a2, a3;
    for (const FractureMesh& frac : meshes) {
        double coeff   = frac.aperture * frac.conductivity;
        double L       = frac.ReturnLength();
        double k_fract = alpha / (frac.conductivity * frac.aperture);
        int n1 = frac.p_ori.index, n2 = frac.p_tar.index, nm = nf + block_idx;

        ReturnDerivativeCoeffFullFourier(0, L, k_fract, w, a1, a2, a3);
        mat.coeffRef(n1, n1) += coeff * frac.aperture * w * w - coeff * a1;
        mat.coeffRef(n1, n2) -= coeff * a2;
        mat.coeffRef(n1, nm) -= coeff * a3;

        ReturnDerivativeCoeffFullFourier(L, L, k_fract, w, a1, a2, a3);
        mat.coeffRef(n2, n1) += coeff * a1;
        mat.coeffRef(n2, n2) += coeff * a2;
        mat.coeffRef(n2, nm) += coeff * a3;
    }
}

static void PorousDDPLinearSystemUpdateFourierEigen(
    const vector<FractureMesh>& meshes, int block_idx, int nf,
    double alpha, int Ny, double w, double dx, double dy,
    const Eigen::MatrixXd& porous_cond, EigenSparse& mat)
{
    int idx = nf + block_idx;
    auto ij = rij(block_idx, Ny);
    mat.coeffRef(idx, idx) += porous_cond(ij.first, ij.second) * w * w * dx * dy;

    double c1, c2, c3;
    for (const FractureMesh& frac : meshes) {
        CumulatePotFractCoeffFourier(frac.ReturnLength(),
            alpha / (frac.conductivity * frac.aperture), w, c1, c2, c3);
        mat.coeffRef(idx, frac.p_ori.index) -= alpha * c1;
        mat.coeffRef(idx, frac.p_tar.index) -= alpha * c2;
        mat.coeffRef(idx, idx)              += alpha * (frac.ReturnLength() - c3);
    }
}

static void BCDDPLinearSystemFourierEigen(
    int nf, double w, double dx, double dy,
    const SourceTermsEPM& src,
    SubNetworkMap& snmap,
    const BoundaryConditionsDFN& bc_map,
    EigenSparse& mat, EigenVec& vec)
{
    for (auto it = bc_map.begin(); it != bc_map.end(); ++it) {
        int node = it->first;
        if (it->second.bc_type == DIRICHLET) {
            for (EigenSparse::InnerIterator jt(mat, node); jt; ++jt)
                jt.valueRef() = 0.0;
            mat.coeffRef(node, node) = 1.0;
            vec(node) = it->second.bc_value;
        } else if (it->second.bc_type == NEUMANN) {
            vec(node) = it->second.bc_value;
        } else if (it->second.bc_type == MIXED_BC) {
            NetworkMeshes& sn = snmap.at(it->second.block_index);
            CgalPoint2D npos  = sn.ReturnNodeCoord(node);
            CgalPoint2D sterm = CenterPosition(src.begin()->first, {dx, dy});
            double r = distance_2D(sterm, npos);
            string border;
            if (sn.border_map.find(node) != sn.border_map.end())
                border = sn.border_map.at(node);
            double alpha = w * boost::math::cyl_bessel_k(1, w*r)
                             / boost::math::cyl_bessel_k(0, w*r);
            double proj = 1.0;
            if      (border == LEFT_BORDER  || border == RIGHT_BORDER)
                proj = abs(npos.x() - sterm.x()) / r;
            else if (border == BOTTOM_BORDER || border == TOP_BORDER)
                proj = abs(npos.y() - sterm.y()) / r;
            mat.coeffRef(node, node) += 2.0 * sn.return_aperture_from_node(node)
                                            * sn.return_conductivity_from_node(node)
                                            * alpha * proj / (2.0 + alpha * dx);
        }
    }
}

static void UpdateLinearSystemBC_EPMEigen(
    int nf, int Nx, int Ny,
    const BoundaryConditionsEPM& bc,
    const Eigen::MatrixXd& pc,
    double dx, double dy,
    EigenSparse& mat, EigenVec& vec)
{
    double c1, c3;
    map<string, BoundaryConditions> cell;
    auto apply = [&](int i, int j, const string& border) {
        int idx = nf + ri(i, j, Ny);
        cell = bc.at({i, j});
        BoundCondCoeffs(border, cell.at(border), pc(i,j), dx, dy, c1, c3);
        mat.coeffRef(idx, idx) += c1;
        vec(idx) += c3;
    };
    for (int j = 0; j < Ny; j++) apply(0,    j, LEFT_BORDER);
    for (int j = 0; j < Ny; j++) apply(Nx-1, j, RIGHT_BORDER);
    for (int i = 0; i < Nx; i++) apply(i,    0, BOTTOM_BORDER);
    for (int i = 0; i < Nx; i++) apply(i, Ny-1, TOP_BORDER);
}

static void SourceTermFourierEigen(
    const SourceTermsEPM& src, int nf, int Ny, bool singularity, EigenVec& vec)
{
    if (!singularity)
        for (auto it = src.begin(); it != src.end(); ++it)
            vec(ri(it->first.first, it->first.second, Ny) + nf) += it->second * 0.5;
}

// ─────────────────────────────────────────────────────────────────────────────
// UpdateLinearSystemFourierEigen : point d'entrée public
// ─────────────────────────────────────────────────────────────────────────────

void UpdateLinearSystemFourierEigen(
    DDPEigenBase& base, double w,
    const BoundaryConditionsEPM& bc_def_epm,
    const BoundaryConditionsDFN& bc_map,
    const SourceTermsEPM& src,
    bool singularity,
    EigenSparse& mat, EigenVec& vec)
{
    SubNetworkMap& snmap = *base.subnetwork_map_ptr;
    for (auto it = snmap.begin(); it != snmap.end(); ++it) {
        auto ij    = rij(it->first, base.Ny);
        double sig = base.porous_cond(ij.first, ij.second);
        double alp = ReturnAlphaBlock(sig, base.delta_x, base.delta_y, it->second);
        DFNDDPLinearSystemFourierEigen(it->second.meshes, it->first,
                                       base.nb_nodes_fract, alp, w, mat);
        PorousDDPLinearSystemUpdateFourierEigen(it->second.meshes, it->first,
                                                base.nb_nodes_fract, alp,
                                                base.Ny, w, base.delta_x, base.delta_y,
                                                base.porous_cond, mat);
    }
    BCDDPLinearSystemFourierEigen(base.nb_nodes_fract, w, base.delta_x, base.delta_y,
                                   src, snmap, bc_map, mat, vec);
    UpdateLinearSystemBC_EPMEigen(base.nb_nodes_fract, base.Nx, base.Ny,
                                   bc_def_epm, base.porous_cond,
                                   base.delta_x, base.delta_y, mat, vec);
    SourceTermFourierEigen(src, base.nb_nodes_fract, base.Ny, singularity, vec);
}

// ─────────────────────────────────────────────────────────────────────────────
// EvaluatePotentialEigen : factorisation numérique + résolution
// (analyzePattern doit avoir été appelé une fois avant le premier appel)
// Le solveur actif est contrôlé par le #define DDP_EIGEN_SOLVER_* dans le .h
// ─────────────────────────────────────────────────────────────────────────────

void EvaluatePotentialEigen(
    int nf, int Nx, int Ny,
    DDPEigenSolver& solver,
    EigenSparse& mat, const EigenVec& vec,
    EigenVec& DFNPotential, Eigen::MatrixXd& PorousPotential)
{
    solver.factorize(mat);
    if (solver.info() != Eigen::Success) {
#if defined(DDP_EIGEN_SOLVER_LLT)
        cout << "[ERROR][DDP_EIGEN] SimplicialLLT factorization failed (matrice non SPD ?)" << endl;
#elif defined(DDP_EIGEN_SOLVER_UMFPACK)
        cout << "[ERROR][DDP_EIGEN] UmfPackLU factorization failed" << endl;
#endif
        return;
    }
    EigenVec phi = solver.solve(vec);

    // Recopie sans boucle ublas : affectation directe dans des conteneurs Eigen denses
    DFNPotential = phi.head(nf);
    // phi[nf + i*Ny + j] → PorousPotential(i,j) : mapping row-major via Eigen::Map
    PorousPotential = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        phi.data() + nf, Nx, Ny);
}
