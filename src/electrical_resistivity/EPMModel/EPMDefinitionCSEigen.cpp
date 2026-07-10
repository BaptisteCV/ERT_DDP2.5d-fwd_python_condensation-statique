#include "EPMDefinitionCSEigen.h"
#include "../Utilitaries/Math/LinearSystemEigen.h"
#include "../Utilitaries/Math/Operations.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <iomanip>

// using namespace std; removed — use std:: qualifiers

namespace eigen_cs {

namespace {

inline void AddBlockCoefficient(std::vector<Triplet>& triplets, int row, int col, double value) {
    if (value != 0.0) {
        triplets.emplace_back(row, col, value);
    }
}

}  // namespace

EPMSystemCSEigen::EPMSystemCSEigen()
        : Nx(0), Ny(0), delta_x(0.0), delta_y(0.0), w(0.0), factor_Aii_ready(false), factor_Aii_llt_ready(false), Ni(0), Nb(0) {
}

EPMSystemCSEigen::EPMSystemCSEigen(int Nx_, int Ny_, double Lx, double Ly, const DenseMatrix& conductivity)
    : Nx(Nx_), Ny(Ny_), delta_x(Lx / Nx_), delta_y(Ly / Ny_), w(0.0),
            porous_cond(conductivity), factor_Aii_ready(false), factor_Aii_llt_ready(false), Ni(0), Nb(0) {
    potential = DenseMatrix::Zero(Nx, Ny);
}

EPMSystemCSEigen::EPMSystemCSEigen(const EPMSystemCSEigen& other)
    : Nx(other.Nx), Ny(other.Ny), delta_x(other.delta_x), delta_y(other.delta_y), w(other.w),
      porous_cond(other.porous_cond), bc_def_epm(other.bc_def_epm), source_terms_epm(other.source_terms_epm),
      potential(other.potential), block_Aii(other.block_Aii), block_Abb(other.block_Abb),
      block_Aib(other.block_Aib), block_Abi(other.block_Abi), matrix_Y(other.matrix_Y), matrix_S(other.matrix_S),
      matrix_S_precomputed(other.matrix_S_precomputed), vector_z(other.vector_z), vector_rhs(other.vector_rhs),
      vector_phi_i(other.vector_phi_i), vector_phi_b(other.vector_phi_b), vector_bi(other.vector_bi),
            vector_delta_diag(other.vector_delta_diag),
        vector_bb(other.vector_bb), factor_Aii_ready(false), factor_Aii_llt_ready(false),
      interior_global_indices(other.interior_global_indices), border_global_indices(other.border_global_indices),
      global_to_interior(other.global_to_interior), global_to_border(other.global_to_border),
      Ni(other.Ni), Nb(other.Nb) {
}

EPMSystemCSEigen& EPMSystemCSEigen::operator=(const EPMSystemCSEigen& other) {
    if (this != &other) {
        Nx = other.Nx;
        Ny = other.Ny;
        delta_x = other.delta_x;
        delta_y = other.delta_y;
        w = other.w;
        porous_cond = other.porous_cond;
        bc_def_epm = other.bc_def_epm;
        source_terms_epm = other.source_terms_epm;
        potential = other.potential;
        block_Aii = other.block_Aii;
        block_Abb = other.block_Abb;
        block_Aib = other.block_Aib;
        block_Abi = other.block_Abi;
        matrix_Y = other.matrix_Y;
        matrix_S = other.matrix_S;
        matrix_S_precomputed = other.matrix_S_precomputed;
        vector_delta_diag = other.vector_delta_diag;
        vector_z = other.vector_z;
        vector_rhs = other.vector_rhs;
        vector_phi_i = other.vector_phi_i;
        vector_phi_b = other.vector_phi_b;
        vector_bi = other.vector_bi;
        vector_bb = other.vector_bb;
        interior_global_indices = other.interior_global_indices;
        border_global_indices = other.border_global_indices;
        global_to_interior = other.global_to_interior;
        global_to_border = other.global_to_border;
        Ni = other.Ni;
        Nb = other.Nb;
        factor_Aii_ready = false;
        factor_Aii_llt_ready = false;
    }
    return *this;
}

void EPMSystemCSEigen::SetConductivity(const DenseMatrix& conductivity) {
    porous_cond = conductivity;
}

void EPMSystemCSEigen::InitializeSCIndexing() {
    interior_global_indices.clear();
    border_global_indices.clear();
    global_to_interior.clear();
    global_to_border.clear();
    Ni = 0;
    Nb = 0;

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            int global_idx = return_index(i, j, Ny);
            bool is_border = (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1);
            if (is_border) {
                global_to_border[global_idx] = Nb;
                border_global_indices.push_back(global_idx);
                ++Nb;
            } else {
                global_to_interior[global_idx] = Ni;
                interior_global_indices.push_back(global_idx);
                ++Ni;
            }
        }
    }
}

void EPMSystemCSEigen::DetermineBasicLinearSystemCS() {
    std::vector<Triplet> triplets_Aii;
    std::vector<Triplet> triplets_Aib;
    std::vector<Triplet> triplets_Abi;
    std::vector<Triplet> triplets_Abb;
    triplets_Aii.reserve(4 * Nx * Ny);
    triplets_Aib.reserve(4 * Nx * Ny);
    triplets_Abi.reserve(4 * Nx * Ny);
    triplets_Abb.reserve(4 * Nx * Ny);

    vector_bi = Vector::Zero(Ni);
    vector_bb = Vector::Zero(Nb);

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            int index_center = return_index(i, j, Ny);
            bool row_interior = (global_to_interior.count(index_center) > 0);
            int local_row = row_interior ? global_to_interior[index_center] : global_to_border[index_center];
            double coeff;

            auto handle_neighbor = [&](int i_neigh, int j_neigh, double factor, double sign) {
                int index_neigh = return_index(i_neigh, j_neigh, Ny);
                bool col_interior = (global_to_interior.count(index_neigh) > 0);
                int local_col = col_interior ? global_to_interior[index_neigh] : global_to_border[index_neigh];
                if (row_interior && col_interior) {
                    AddBlockCoefficient(triplets_Aii, local_row, local_row, factor);
                    AddBlockCoefficient(triplets_Aii, local_row, local_col, -factor);
                } else if (row_interior && !col_interior) {
                    AddBlockCoefficient(triplets_Aii, local_row, local_row, factor);
                    AddBlockCoefficient(triplets_Aib, local_row, local_col, -factor);
                } else if (!row_interior && col_interior) {
                    AddBlockCoefficient(triplets_Abb, local_row, local_row, factor);
                    AddBlockCoefficient(triplets_Abi, local_row, local_col, -factor);
                } else {
                    AddBlockCoefficient(triplets_Abb, local_row, local_row, factor);
                    AddBlockCoefficient(triplets_Abb, local_row, local_col, -factor);
                }
            };

            if (i != 0) {
                coeff = Harmonic_Mean(porous_cond(i, j), porous_cond(i - 1, j)) * delta_y / delta_x;
                handle_neighbor(i - 1, j, coeff, -1.0);
            }
            if (i != Nx - 1) {
                coeff = Harmonic_Mean(porous_cond(i, j), porous_cond(i + 1, j)) * delta_y / delta_x;
                handle_neighbor(i + 1, j, coeff, -1.0);
            }
            if (j != 0) {
                coeff = Harmonic_Mean(porous_cond(i, j), porous_cond(i, j - 1)) * delta_x / delta_y;
                handle_neighbor(i, j - 1, coeff, -1.0);
            }
            if (j != Ny - 1) {
                coeff = Harmonic_Mean(porous_cond(i, j), porous_cond(i, j + 1)) * delta_x / delta_y;
                handle_neighbor(i, j + 1, coeff, -1.0);
            }
        }
    }

    block_Aii.resize(Ni, Ni);
    block_Aib.resize(Ni, Nb);
    block_Abi.resize(Nb, Ni);
    block_Abb.resize(Nb, Nb);
    block_Aii.setFromTriplets(triplets_Aii.begin(), triplets_Aii.end());
    block_Aib.setFromTriplets(triplets_Aib.begin(), triplets_Aib.end());
    block_Abi.setFromTriplets(triplets_Abi.begin(), triplets_Abi.end());
    block_Abb.setFromTriplets(triplets_Abb.begin(), triplets_Abb.end());
    block_Aii.makeCompressed();
    block_Aib.makeCompressed();
    block_Abi.makeCompressed();
    block_Abb.makeCompressed();
}

void EPMSystemCSEigen::UpdateSCBlocksFourier() {
    factor_Aii_ready = false;
    factor_Aii_llt_ready = false;
    
    // Must uncompress matrices before using coeffRef() to modify them
    // (they were makeCompressed() in DetermineBasicLinearSystemCS)
    // block_Aii.uncompress();
    // block_Abb.uncompress();
    
    for (int li = 0; li < Ni; ++li) {
        int g = interior_global_indices[li];
        std::pair<int, int> ij = return_indices(g, Ny);
        block_Aii.coeffRef(li, li) += porous_cond(ij.first, ij.second) * w * w * delta_x * delta_y;
    }
    for (int lb = 0; lb < Nb; ++lb) {
        int g = border_global_indices[lb];
        std::pair<int, int> ij = return_indices(g, Ny);
        block_Abb.coeffRef(lb, lb) += porous_cond(ij.first, ij.second) * w * w * delta_x * delta_y;
    }
    
    // Recompress matrices for efficient storage and future factorizations
    // block_Aii.makeCompressed();
    // block_Abb.makeCompressed();
}

void EPMSystemCSEigen::UpdateSCBoundaryConditionsAndSources(bool fourier) {
    std::map<std::string, BoundaryConditions> bc_cell;
    double coeff1 = 0.0;
    double coeff3 = 0.0;
    Vector delta_diag = Vector::Zero(Nb);

    for (int lb = 0; lb < Nb; ++lb) {
        int g = border_global_indices[lb];
        std::pair<int, int> ij = return_indices(g, Ny);
        int i = ij.first;
        int j = ij.second;
        bc_cell = bc_def_epm[std::make_pair(i, j)];

        if (i == 0) {
            BoundCondCoeffsCSEigen(LEFT_BORDER, bc_cell[LEFT_BORDER], porous_cond(i, j), delta_x, delta_y, coeff1, coeff3);
            vector_bb(lb) += coeff3;
            delta_diag(lb) += coeff1;
        }
        if (i == Nx - 1) {
            BoundCondCoeffsCSEigen(RIGHT_BORDER, bc_cell[RIGHT_BORDER], porous_cond(i, j), delta_x, delta_y, coeff1, coeff3);
            vector_bb(lb) += coeff3;
            delta_diag(lb) += coeff1;
        }
        if (j == 0) {
            BoundCondCoeffsCSEigen(BOTTOM_BORDER, bc_cell[BOTTOM_BORDER], porous_cond(i, j), delta_x, delta_y, coeff1, coeff3);
            vector_bb(lb) += coeff3;
            delta_diag(lb) += coeff1;
        }
        if (j == Ny - 1) {
            BoundCondCoeffsCSEigen(TOP_BORDER, bc_cell[TOP_BORDER], porous_cond(i, j), delta_x, delta_y, coeff1, coeff3);
            vector_bb(lb) += coeff3;
            delta_diag(lb) += coeff1;
        }
    }

    bool bi_updated = false;
    for (SourceTermsCSEigen::iterator it = source_terms_epm.begin(); it != source_terms_epm.end(); ++it) {
        int i = it->first.first;
        int j = it->first.second;
        int g = return_index(i, j, Ny);
        double val = it->second;
        if (fourier) {
            val *= 0.5;
        }
        if (global_to_interior.count(g) > 0) {
            vector_bi(global_to_interior[g]) += val;
            bi_updated = true;
        } else if (global_to_border.count(g) > 0) {
            vector_bb(global_to_border[g]) += val;
        }
    }
    matrix_S = matrix_S_precomputed;
    vector_delta_diag = delta_diag;
    for (int r = 0; r < Nb; ++r) {
        matrix_S(r, r) += delta_diag(r);
    }

    if (bi_updated) {
        if (Ni > 0) {
            // Ancien chemin conservé pour référence : factorisation LU creuse.
            // if (!factor_Aii_ready) {
            //     LinearSystemFactorize(block_Aii, factor_Aii);
            //     factor_Aii_ready = true;
            // }
            // vector_z = LinearSystemSolveFactorized(factor_Aii, vector_bi);

            if (!factor_Aii_llt_ready) {
                LinearSystemFactorizeLLT(block_Aii, factor_Aii_llt);
                factor_Aii_llt_ready = true;
            }
            vector_z = LinearSystemSolveFactorizedLLT(factor_Aii_llt, vector_bi);
        } else {
            vector_z = Vector::Zero(0);
        }
    }

    vector_rhs = vector_bb - block_Abi * vector_z;
}

void EPMSystemCSEigen::PrecomputeSCForCurrentFrequency() {
    clock_t t_total = clock();
    
    // Step 1: Factorize Aii
    clock_t t_factorize = clock();
    // Ancien chemin conservé pour référence : factorisation LU creuse.
    // if (Ni > 0) {
    //     LinearSystemFactorize(block_Aii, factor_Aii);
    //     factor_Aii_ready = true;
    // }
    if (Ni > 0) {
        LinearSystemFactorizeLLT(block_Aii, factor_Aii_llt);
        factor_Aii_llt_ready = true;
    }
    double t_factorize_total = (double)(clock() - t_factorize) / CLOCKS_PER_SEC;
    std::cout << "[PROFILE][EigenCS]     Step 1 - Factorize Aii: " << t_factorize_total << " sec" << std::endl;

    // Step 2: Compute matrix_Y (Aii^{-1} * Aib)
    clock_t t_matrix_y = clock();
    matrix_Y = DenseMatrix::Zero(Ni, Nb);
    for (int c = 0; c < Nb; ++c) {
        Vector rhs_col = Vector::Zero(Ni);
        for (SparseMatrix::InnerIterator it(block_Aib, c); it; ++it) {
            rhs_col(it.row()) = it.value();
        }
        Vector x_col = (Ni > 0) ? LinearSystemSolveFactorizedLLT(factor_Aii_llt, rhs_col) : Vector::Zero(0);
        for (int r = 0; r < Ni; ++r) {
            matrix_Y(r, c) = x_col(r);
        }
    }
    double t_matrix_y_total = (double)(clock() - t_matrix_y) / CLOCKS_PER_SEC;
    std::cout << "[PROFILE][EigenCS]     Step 2 - Compute matrix_Y: " << t_matrix_y_total << " sec" << std::endl;

    // Step 3: Compute condensed system matrix S = Abb - Abi*Y
    clock_t t_matrix_s = clock();
    matrix_S = DenseMatrix::Zero(Nb, Nb);
    for (int i = 0; i < Nb; ++i) {
        for (int j = 0; j < Nb; ++j) {
            matrix_S(i, j) = block_Abb.coeff(i, j);
        }
    }
    matrix_S.noalias() -= (block_Abi * matrix_Y);
    matrix_S_precomputed = matrix_S;
    double t_matrix_s_total = (double)(clock() - t_matrix_s) / CLOCKS_PER_SEC;
    std::cout << "[PROFILE][EigenCS]     Step 3 - Compute matrix_S: " << t_matrix_s_total << " sec" << std::endl;

    // Step 4: Compute vector_z (Aii^{-1} * bi)
    clock_t t_vector_z = clock();
    vector_z = (Ni > 0) ? LinearSystemSolveFactorizedLLT(factor_Aii_llt, vector_bi) : Vector::Zero(0);
    double t_vector_z_total = (double)(clock() - t_vector_z) / CLOCKS_PER_SEC;
    std::cout << "[PROFILE][EigenCS]     Step 4 - Compute vector_z: " << t_vector_z_total << " sec" << std::endl;

    // Step 5: Compute vector_rhs
    clock_t t_vector_rhs = clock();
    vector_rhs = vector_bb - block_Abi * vector_z;
    double t_vector_rhs_total = (double)(clock() - t_vector_rhs) / CLOCKS_PER_SEC;
    std::cout << "[PROFILE][EigenCS]     Step 5 - Compute vector_rhs: " << t_vector_rhs_total << " sec" << std::endl;

    double t_total_time = (double)(clock() - t_total) / CLOCKS_PER_SEC;
    std::cout << "[PROFILE][EigenCS]   Total PrecomputeSCForCurrentFrequency: " << t_total_time << " sec" << std::endl;
}

double EPMSystemCSEigen::ComputeCondensedSystemResidualNorm() const {
    if (Ni == 0 && Nb == 0) {
        return 0.0;
    }

    if (Nb > 0 && vector_phi_b.size() != Nb) {
        return 0.0;
    }

    Vector phi_i = Vector::Zero(Ni);
    if (Ni > 0) {
        if (factor_Aii_llt_ready) {
            Vector correction = block_Aib * vector_phi_b;
            phi_i = vector_z - LinearSystemSolveFactorizedLLT(factor_Aii_llt, correction);
        } else if (vector_phi_i.size() == Ni) {
            phi_i = vector_phi_i;
        } else {
            return 0.0;
        }
    }

    Vector residual_i = Vector::Zero(Ni);
    if (Ni > 0) {
        residual_i = block_Aii * phi_i + block_Aib * vector_phi_b - vector_bi;
    }

    Vector residual_b = Vector::Zero(Nb);
    if (Nb > 0) {
        residual_b = block_Abi * phi_i + block_Abb * vector_phi_b - vector_bb;
    }

    Vector residual_b_with_bc = residual_b;
    if (Nb > 0 && matrix_S.rows() == Nb && matrix_S_precomputed.rows() == Nb && vector_phi_b.size() == Nb) {
        // BC diagonal correction is injected into matrix_S (not block_Abb).
        // Reinject it in the border residual to compare apples-to-apples.
        for (int i = 0; i < Nb; ++i) {
            double delta_diag = matrix_S(i, i) - matrix_S_precomputed(i, i);
            residual_b_with_bc(i) += delta_diag * vector_phi_b(i);
        }
    }

    Vector residual_condensed = Vector::Zero(Nb);
    if (Nb > 0 && matrix_S.rows() == Nb && matrix_S.cols() == Nb && vector_rhs.size() == Nb && vector_phi_b.size() == Nb) {
        residual_condensed = matrix_S * vector_phi_b - vector_rhs;
    }


    double norm_aii = (Ni > 0) ? residual_i.norm() : 0.0;
    double denom_aii = 0.0;
    if (Ni > 0) {
        Vector rhs_aii = vector_bi - block_Aib * vector_phi_b;
        denom_aii = rhs_aii.norm();
    }
    double relative_aii = (denom_aii > 0.0) ? (norm_aii / denom_aii) : norm_aii;

    double norm_condensed = (Nb > 0) ? residual_condensed.norm() : 0.0;
    double denom_condensed = (Nb > 0) ? vector_rhs.norm() : 0.0;
    double relative_condensed = (denom_condensed > 0.0) ? (norm_condensed / denom_condensed) : norm_condensed;

    // Debug logging per source loop: residuals of A_ii and S inversions.
    std::cout << "[PROFILE][EigenCS]   Residual of A_ii inversion ||A_ii*phi_i-rhs_i||_2 = " << norm_aii
              << " | Relative ||A_ii*phi_i-rhs_i||_2/||rhs_i||_2 = " << relative_aii
              << std::endl;

    std::cout << "[PROFILE][EigenCS]   Residual of S inversion ||S*phi_b-rhs||_2 = " << norm_condensed
              << " | Relative ||S*phi_b-rhs||_2/||rhs||_2 = " << relative_condensed
              << std::endl;

    Vector residual(Ni + Nb);
    if (Ni > 0) {
        residual.head(Ni) = residual_i;
    }
    if (Nb > 0) {
        residual.tail(Nb) = residual_b;
    }
    return relative_condensed;
}

void EPMSystemCSEigen::SolveCondensedSystem() {
    if (Nb > 0) {
        vector_phi_b = LinearSystemSolvingDense(matrix_S, vector_rhs);
    } else {
        vector_phi_b = Vector::Zero(0);
    }
}

void EPMSystemCSEigen::ReconstructInteriorFromBorder() {
    vector_phi_i = Vector::Zero(Ni);
    if (Ni == 0) {
        return;
    }

    Vector vector_X = block_Aib * vector_phi_b;
    if (!factor_Aii_llt_ready) {
        LinearSystemFactorizeLLT(block_Aii, factor_Aii_llt);
        factor_Aii_llt_ready = true;
    }
    Vector solve_Aii_X = LinearSystemSolveFactorizedLLT(factor_Aii_llt, vector_X);
    for (int r = 0; r < Ni; ++r) {
        vector_phi_i(r) = vector_z(r) - solve_Aii_X(r);
    }
}

Vector EPMSystemCSEigen::BuildGlobalSolutionVectorFromSC() const {
    Vector global_phi = Vector::Zero(Nx * Ny);
    for (int i = 0; i < Ni; ++i) {
        int global_idx = interior_global_indices[i];
        global_phi(global_idx) = vector_phi_i(i);
    }
    for (int i = 0; i < Nb; ++i) {
        int global_idx = border_global_indices[i];
        global_phi(global_idx) = vector_phi_b(i);
    }
    return global_phi;
}

DenseMatrix EPMSystemCSEigen::PotentialComputationSC() {
    Vector phi_global = BuildGlobalSolutionVectorFromSC();
    DenseMatrix potential_sc(Nx, Ny);
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            potential_sc(i, j) = phi_global(return_index(i, j, Ny));
        }
    }
    potential = potential_sc;
    return potential_sc;
}

/*
void EPMSystemCSEigen::ExportFirstRowOfAiiForFrequency(const std::string& output_dir, double current_w) {
    // Create matrices subdirectory if it doesn't exist
    std::string matrices_dir = output_dir + "/matrices";
    std::filesystem::create_directories(matrices_dir);
    
    // Create filename based on frequency
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(10) << current_w;
    std::string freq_filename = matrices_dir + "/freq_" + oss.str() + ".txt";
    
    // Open file in write mode (overwrite - called once per frequency)
    std::ofstream out_file(freq_filename.c_str(), std::ios::out);
    if (!out_file.is_open()) {
        std::cerr << "[ERROR][EigenCS] Cannot open output file: " << freq_filename << std::endl;
        return;
    }
    
    // Write first row of block_Aii (same for all sources)
    out_file << "=== Frequency " << std::scientific << std::setprecision(10) << current_w << " ===" << std::endl;
    out_file << "A_ii_row_0: ";
    for (int col = 0; col < block_Aii.cols(); ++col) {
        out_file << std::scientific << std::setprecision(12) << block_Aii.coeff(0, col) << " ";
    }
    out_file << std::endl << std::endl;
    
    out_file.close();
}
*/

void EPMSystemCSEigen::ExportFirstRowsToFile(const std::string& output_dir, int elec_nb, double current_w) {
    // Create matrices subdirectory if it doesn't exist
    std::string matrices_dir = output_dir + "/matrices";
    std::filesystem::create_directories(matrices_dir);
    
    // Create filename based on frequency
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(10) << current_w;
    std::string freq_filename = matrices_dir + "/freq_" + oss.str() + ".txt";
    
    // Open file in append mode (sources are appended to the file)
    std::ofstream out_file(freq_filename.c_str(), std::ios::app);
    if (!out_file.is_open()) {
        std::cerr << "[ERROR][EigenCS] Cannot open output file: " << freq_filename << std::endl;
        return;
    }
    
    // Write source/electrode number as header
    out_file << "=== Source " << elec_nb << " ===" << std::endl;

    out_file << "delta_diag (all coefficients with indices)" << std::endl;
    for (int r = 0; r < (int)vector_delta_diag.size(); ++r) {
        // Get the global index and convert to (i,j)
        int global_idx = border_global_indices[r];
        std::pair<int, int> ij = return_indices(global_idx, Ny);
        int i = ij.first;
        int j = ij.second;
        
        out_file << "delta_diag_" << r 
                 << " (i=" << i << ",j=" << j << ")"
                 << ": "
                 << std::scientific << std::setprecision(12) << vector_delta_diag(r)
                 << std::endl;
    }
    out_file << std::endl;
    
    out_file.close();
}

}  // namespace eigen_cs