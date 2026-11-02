#include "FourierModelFunctionsCSEigen.h"
#include "../EPMModel/EPMDefinitionCSEigen.h"
#include "../EPMModel/EPMComputation.h"
#include "../EPMModel/EPMComputationCSEigen.h"
#include "../DDPModel/DDPModelCS.h"
#include "../DDPModel/DDPDFNFunctions.h"
#include "../DFNModel/DFNComputation.h"
#include "../DFNModel/NetworkMeshes.h"
#include "../Utilitaries/Storage/UblasStructures.h"
#include "../Utilitaries/Storage/UtilitaryStorage.h"
#include <sstream>
#include <iostream>
#include <ctime>
#include <fstream>
#include <filesystem>
#include <omp.h>

// using namespace std; removed — use std:: qualifiers

namespace eigen_cs {

namespace {

SparseMatrix DenseMatrixToSparse(const DenseMatrix& dense) {
    std::vector<Triplet> triplets;
    triplets.reserve(static_cast<std::size_t>(dense.rows() * dense.cols()));
    for (int i = 0; i < dense.rows(); ++i) {
        for (int j = 0; j < dense.cols(); ++j) {
            if (dense(i, j) != 0.0) {
                triplets.emplace_back(i, j, dense(i, j));
            }
        }
    }
    SparseMatrix sparse(dense.rows(), dense.cols());
    sparse.setFromTriplets(triplets.begin(), triplets.end());
    sparse.makeCompressed();
    return sparse;
}

void WriteSparseMatrixResult(const std::string& output_file, const SparseMatrix& matrix) {
    std::ofstream output(output_file.c_str(), std::ofstream::out);
    if (!output.is_open()) {
        std::cerr << "[ERROR][EigenCS] Cannot open output file: " << output_file << std::endl;
        return;
    }
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            output << matrix.coeff(i, j) << " ";
        }
        output << std::endl;
    }
}

void WriteDenseMatrixResult(const std::string& output_file, const DenseMatrix& matrix) {
    std::ofstream output(output_file.c_str(), std::ofstream::out);
    if (!output.is_open()) {
        std::cerr << "[ERROR][EigenCS] Cannot open output file: " << output_file << std::endl;
        return;
    }
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            output << matrix(i, j) << " ";
        }
        output << std::endl;
    }
}

}  // namespace

// CS-typed version of ReturnBoundCondDFNFourier: uses BoundaryConditionsCSEigen
// (keyed by pair<int,int>) instead of BoundaryConditionsEPM (keyed by pair<double,double>).
static BoundaryConditionsDFN ReturnBoundCondDFNFourierCS(
    Domain domain,
    BoundaryConditionsDef bc_map_classic,
    BordersMap border_map,
    NodesMap nodes_map,
    int Nx, int Ny,
    double delta_x, double delta_y,
    BoundaryConditionsCSEigen bc_def_epm_cs)
{
    BoundaryConditionsDFN bc_map;
    for (BordersMap::const_iterator it = border_map.begin(); it != border_map.end(); ++it) {
        const int index = it->first;
        const std::string& bc_type = bc_map_classic.at(it->second).bc_type;
        if (bc_type == DIRICHLET || bc_type == NEUMANN) {
            bc_map[index] = bc_map_classic.at(it->second);
        } else if (bc_type == DIR_GRADIENT) {
            double border_dist = 0.0, position = 0.0;
            std::string border1, border2;
            position = domain.ReturnBorderCoordinate(nodes_map.at(index), border_dist, border1, border2);
            if (bc_map_classic.at(border1).bc_type == DIRICHLET && bc_map_classic.at(border2).bc_type == DIRICHLET) {
                double bc_value = (bc_map_classic.at(border2).bc_value - bc_map_classic.at(border1).bc_value)
                                  / border_dist * position + bc_map_classic.at(border1).bc_value;
                bc_map[index] = BoundaryConditions(DIRICHLET, bc_value);
            }
        } else if (bc_type == MIXED_BC) {
            CgalPoint2D node_pos = nodes_map.at(index);
            double x = node_pos.x(), y = node_pos.y();
            if (x == domain.max_pt.i) x -= EPSILON;
            if (y == domain.max_pt.j) y -= EPSILON;
            std::pair<int,int> ij = ::return_indices(x, y, delta_x, delta_y, Nx, Ny);
            BoundaryConditionsCSEigen::const_iterator it_bc = bc_def_epm_cs.find(ij);
            if (it_bc != bc_def_epm_cs.end()) {
                const std::map<std::string, BoundaryConditions>& cell_bc = it_bc->second;
                std::map<std::string, BoundaryConditions>::const_iterator it_cell = cell_bc.find(it->second);
                if (it_cell != cell_bc.end()) {
                    bc_map[index] = BoundaryConditions(MIXED_BC, it_cell->second.bc_value);
                    bc_map[index].block_index = return_index(ij.first, ij.second, Ny);
                }
            }
        }
    }
    return bc_map;
}

std::map<int, SparseMatrix> FourierPotentialEPMSourceTermCSEigen(
    Parameters param,
    Domain domain,
    bool singularity,
    bool bc_on_border,
    std::map<int, std::pair<double, double> > source_term_positions,
    bool build_result_map) {
    (void)singularity;
    (void)bc_on_border;

    std::map<int, SparseMatrix> res_matrices;
    const double source_value = 1.0;
    const bool fourier = true;
    const BoundaryConditionsDef bc_map = param.simu_param.bc_map;
    const std::vector<double> w_vect = param.simu_param.EPM_param.w_vect;

    clock_t t_global_start = clock();
    clock_t t_step_start = clock();

    t_step_start = clock();
    DenseMatrix conductivity(param.simu_param.EPM_param.Nx, param.simu_param.EPM_param.Ny);
    if (param.domain_param.mat_cond != -1) {
        conductivity.setConstant(param.domain_param.mat_cond);
    } else {
        conductivity.setZero();
        for (int i = 0; i < param.simu_param.EPM_param.Nx; ++i) {
            for (int j = 0; j < param.simu_param.EPM_param.Ny; ++j) {
                conductivity(i, j) = param.domain_param.mat_cond_dist(i, j);
            }
        }
    }
        std::cout << "[PROFILE][EigenCS] Build conductivity: "
            << (double)(clock() - t_step_start) / CLOCKS_PER_SEC << " sec" << std::endl;

        t_step_start = clock();
    std::map<int, std::map<double, Vector> > potential_by_elec;
    for (std::map<int, std::pair<double, double> >::const_iterator it = source_term_positions.begin(); it != source_term_positions.end(); ++it) {
        potential_by_elec[it->first] = std::map<double, Vector>();
    }
        std::cout << "[PROFILE][EigenCS] Init potential_by_elec: "
            << (double)(clock() - t_step_start) / CLOCKS_PER_SEC << " sec" << std::endl;

        t_step_start = clock();
    EPMSystemCSEigen epd_model_base(static_cast<int>(conductivity.rows()), static_cast<int>(conductivity.cols()), domain.domain_size_x(), domain.domain_size_y(), conductivity);
        std::cout << "[PROFILE][EigenCS] Construct EPMSystemCSEigen base: "
            << (double)(clock() - t_step_start) / CLOCKS_PER_SEC << " sec" << std::endl;

        t_step_start = clock();
    epd_model_base.SetTopBorderCondensation(bc_map);
    epd_model_base.InitializeSCIndexing();
        std::cout << "[PROFILE][EigenCS] InitializeSCIndexing: "
            << (double)(clock() - t_step_start) / CLOCKS_PER_SEC << " sec" << std::endl;

        t_step_start = clock();
    epd_model_base.DetermineBasicLinearSystemCS();
        std::cout << "[PROFILE][EigenCS] DetermineBasicLinearSystemCS: "
            << (double)(clock() - t_step_start) / CLOCKS_PER_SEC << " sec" << std::endl;

    // Accumulateurs sur toutes les (fréquences × électrodes)
    double t_acc_copy      = 0.0;
    double t_acc_fourier   = 0.0;
    double t_acc_precompute= 0.0;
    double t_acc_bcdef     = 0.0;
    double t_acc_update    = 0.0;
    double t_acc_solve     = 0.0;
    double t_acc_recon     = 0.0;
    double t_acc_build     = 0.0;
    int nb_solves_cs = 0;

    for (std::vector<double>::const_iterator it_w = w_vect.begin(); it_w != w_vect.end(); ++it_w) {
        double current_w = *it_w;

        t_step_start = clock();
        EPMSystemCSEigen epd_model_w = epd_model_base;
        epd_model_w.w = current_w;
        t_acc_copy += (double)(clock() - t_step_start) / CLOCKS_PER_SEC;

        t_step_start = clock();
        epd_model_w.UpdateSCBlocksFourier();
        t_acc_fourier += (double)(clock() - t_step_start) / CLOCKS_PER_SEC;

        t_step_start = clock();
        epd_model_w.PrecomputeSCForCurrentFrequency();
        t_acc_precompute += (double)(clock() - t_step_start) / CLOCKS_PER_SEC;

        Vector vector_bi_base = epd_model_w.vector_bi;
        Vector vector_bb_base = epd_model_w.vector_bb;

        for (std::map<int, std::pair<double, double> >::const_iterator it_e = source_term_positions.begin(); it_e != source_term_positions.end(); ++it_e) {
            int elec_nb = it_e->first;
            double source_position_x = it_e->second.first;
            double source_position_y = it_e->second.second;

            epd_model_w.vector_bi = vector_bi_base;
            epd_model_w.vector_bb = vector_bb_base;

            t_step_start = clock();
            SourceTermsDef source_terms_local;
            source_terms_local.push_back(SourceTerm(source_position_x, source_position_y, source_value));
            SourceTermsCSEigen source_terms_epm_local;
            epd_model_w.bc_def_epm = BoundaryConditionsDefinitionCSEigen(
                domain, epd_model_w.Nx, epd_model_w.Ny, bc_map,
                source_terms_local, source_terms_epm_local, fourier, epd_model_w.w);
            epd_model_w.source_terms_epm = source_terms_epm_local;
            t_acc_bcdef += (double)(clock() - t_step_start) / CLOCKS_PER_SEC;

            t_step_start = clock();
            epd_model_w.UpdateSCBoundaryConditionsAndSources(fourier);
            t_acc_update += (double)(clock() - t_step_start) / CLOCKS_PER_SEC;

            t_step_start = clock();
            epd_model_w.SolveCondensedSystem();
            t_acc_solve += (double)(clock() - t_step_start) / CLOCKS_PER_SEC;

            t_step_start = clock();
            epd_model_w.ReconstructInteriorFromBorder();
            t_acc_recon += (double)(clock() - t_step_start) / CLOCKS_PER_SEC;

            t_step_start = clock();
            Vector phi_global = epd_model_w.BuildGlobalSolutionVectorFromSC();
            t_acc_build += (double)(clock() - t_step_start) / CLOCKS_PER_SEC;

            potential_by_elec[elec_nb][current_w] = phi_global;
            ++nb_solves_cs;
        }
    }

    t_step_start = clock();
    // Écriture disque uniquement en mode exécutable (MAIN.cpp) : en mode
    // librairie (PERFORM.cpp, build_result_map==true), le potentiel est
    // transmis à Python via res_matrices, l'écriture serait un coût pur.
    if (!build_result_map) {
        std::filesystem::create_directories(param.files_param.output_path);
        for (std::map<int, std::pair<double, double> >::const_iterator it_e = source_term_positions.begin(); it_e != source_term_positions.end(); ++it_e) {
            int elec_nb = it_e->first;
            const std::map<double, Vector>& potential_map = potential_by_elec[elec_nb];
            if (potential_map.empty()) {
                continue;
            }

            std::stringstream ss;
            ss << elec_nb;
            const std::string output_file = param.files_param.output_path + std::string("Results") + ss.str() + ".txt";
            std::ofstream output(output_file.c_str(), std::ofstream::out);
            if (!output.is_open()) {
                std::cerr << "[ERROR][EigenCS] Cannot open output file: " << output_file << std::endl;
                continue;
            }

            for (std::map<double, Vector>::const_iterator it_vec = potential_map.begin(); it_vec != potential_map.end(); ++it_vec) {
                const Vector& current_vector = it_vec->second;
                for (int j = 0; j < epd_model_base.Ny; ++j) {
                    for (int i = 0; i < epd_model_base.Nx; ++i) {
                        output << current_vector(return_index(i, (epd_model_base.Ny - 1) - j, epd_model_base.Ny)) << " ";
                    }
                    output << "\n";
                }
            }
            output.close();
        }
    }
    double t_write_cs = (double)(clock() - t_step_start) / CLOCKS_PER_SEC;

    // Empilement des potentiels (toutes fréquences) en une matrice creuse Eigen
    // par électrode, uniquement nécessaire quand le résultat est renvoyé à l'appelant
    // (PERFORM.cpp / pybind) : MAIN.cpp n'en a pas besoin (résultats déjà écrits sur disque).
    if (build_result_map) {
        const int n_freq_res = (int)w_vect.size();
        const int Nx_res = epd_model_base.Nx;
        const int Ny_res = epd_model_base.Ny;
        for (std::map<int, std::pair<double, double> >::const_iterator it_e = source_term_positions.begin();
             it_e != source_term_positions.end(); ++it_e) {
            int elec_nb = it_e->first;
            const std::map<double, Vector>& potential_map = potential_by_elec[elec_nb];
            if (potential_map.empty()) continue;

            DenseMatrix A_dense(n_freq_res * Nx_res, Ny_res);
            A_dense.setZero();
            int cpt_w = 0;
            for (std::map<double, Vector>::const_iterator it_vec = potential_map.begin();
                 it_vec != potential_map.end(); ++it_vec) {
                const Vector& phi = it_vec->second;
                for (int i = 0; i < Nx_res; ++i)
                    for (int j = 0; j < Ny_res; ++j)
                        A_dense(j + Ny_res * cpt_w, i) = phi(return_index(i, (Ny_res - 1) - j, Ny_res));
                cpt_w++;
            }
            res_matrices[elec_nb] = A_dense.sparseView();
        }
    }

    const int n_freq_cs        = (int)w_vect.size();
    const double t_boucle_freq_cs   = t_acc_copy + t_acc_fourier + t_acc_precompute;
    const double t_boucle_source_cs = t_acc_bcdef + t_acc_update + t_acc_solve + t_acc_recon + t_acc_build;

    std::cout << "\n[PROFILE] ============ FourierPotentialEPMSourceTermCSEigen ============" << std::endl;
    std::cout << "[PROFILE] Hors boucles :" << std::endl;
    std::cout << "[PROFILE]   Build conductivity / EPMSystemCSEigen init : (inclus dans GRAND TOTAL)" << std::endl;

    std::cout << "[PROFILE] Boucle frequences (" << n_freq_cs << " freq) :" << std::endl;
    std::cout << "[PROFILE]   Copie base + set w          : " << t_acc_copy       << " sec  (moy: " << t_acc_copy      /n_freq_cs << " sec/freq)" << std::endl;
    std::cout << "[PROFILE]   UpdateSCBlocksFourier       : " << t_acc_fourier    << " sec  (moy: " << t_acc_fourier   /n_freq_cs << " sec/freq)" << std::endl;
    std::cout << "[PROFILE]   PrecomputeSCForFrequency    : " << t_acc_precompute << " sec  (moy: " << t_acc_precompute/n_freq_cs << " sec/freq)" << std::endl;
    std::cout << "[PROFILE]   Sous-total boucle freq      : " << t_boucle_freq_cs << " sec" << std::endl;

    std::cout << "[PROFILE] Boucle sources (" << nb_solves_cs << " inversions) :" << std::endl;
    std::cout << "[PROFILE]   BCDef (BoundaryConditions)  : " << t_acc_bcdef  << " sec  (moy: " << t_acc_bcdef /nb_solves_cs << " sec/inv)" << std::endl;
    std::cout << "[PROFILE]   UpdateSCBoundaryConditions  : " << t_acc_update << " sec  (moy: " << t_acc_update/nb_solves_cs << " sec/inv)" << std::endl;
    std::cout << "[PROFILE]   SolveCondensedSystem        : " << t_acc_solve  << " sec  (moy: " << t_acc_solve /nb_solves_cs << " sec/inv)" << std::endl;
    std::cout << "[PROFILE]   ReconstructInterior         : " << t_acc_recon  << " sec  (moy: " << t_acc_recon /nb_solves_cs << " sec/inv)" << std::endl;
    std::cout << "[PROFILE]   BuildGlobalSolutionVector   : " << t_acc_build  << " sec  (moy: " << t_acc_build /nb_solves_cs << " sec/inv)" << std::endl;
    std::cout << "[PROFILE]   Sous-total boucle sources   : " << t_boucle_source_cs << " sec" << std::endl;

    std::cout << "[PROFILE] Ecriture resultats            : " << t_write_cs << " sec" << std::endl;
    std::cout << "[PROFILE] GRAND TOTAL                   : " << (double)(clock() - t_global_start) / CLOCKS_PER_SEC << " sec" << std::endl;
    std::cout << "[PROFILE] ====================================================================\n" << std::endl;

    return res_matrices;
}

// ─────────────────────────────────────────────────────────────────────────────
// DDP_CS : parallélisation sur les sources/électrodes (boucle interne).
// La boucle colonnes de matrix_Y est aussi parallélisée dans PrecomputeDDPSCForCurrentFrequency.
//   → Commenter la ligne suivante pour revenir en mode séquentiel.
#define DDP_CS_PARALLEL_SOURCES
static const int DDP_CS_NB_THREADS = 14;
// ─────────────────────────────────────────────────────────────────────────────
std::map<int, SparseMatrix> FourierPotentialDDPSourceTermCSEigen(
    Parameters param,
    Domain domain,
    bool singularity,
    bool bc_on_border,
    std::map<int, std::pair<double, double> > source_term_positions,
    bool build_result_map)
{
    (void)bc_on_border;

    std::map<int, SparseMatrix> res_matrices;
    const double source_value = 1.0;
    const bool fourier = true;
    const BoundaryConditionsDef bc_map = param.simu_param.bc_map;
    const std::vector<double> w_vect = param.simu_param.DDP_param.w_vect;
    const int Nx = param.simu_param.DDP_param.Nx;
    const int Ny = param.simu_param.DDP_param.Ny;
    const int n_freq = (int)w_vect.size();

    double t_global_start = omp_get_wtime();
    double t0;

    // 1. Build conductivity matrix
    t0 = omp_get_wtime();
    DenseMatrix conductivity(Nx, Ny);
    if (param.domain_param.mat_cond != -1) {
        conductivity.setConstant(param.domain_param.mat_cond);
    } else {
        conductivity.setZero();
        for (int i = 0; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
                conductivity(i, j) = param.domain_param.mat_cond_dist(i, j);
    }
    // std::cout << "[PROFILE][DDPCS] Build conductivity: "
    //           << omp_get_wtime() - t0 << " sec" << std::endl;

    // 2. Generate fracture network
    t0 = omp_get_wtime();
    NetworkMeshes net_mesh(domain, param, param.simu_param.simu_option);
    // std::cout << "[PROFILE][DDPCS] Generate fracture network: "
    //           << omp_get_wtime() - t0 << " sec" << std::endl;

    // 3. Construct base DDP model and initialise indexing + basic blocks
    t0 = omp_get_wtime();
    DDPModelCS ddp_base(param, domain, net_mesh, fourier, singularity);
    ddp_base.InitializeDFNSCIndexing();
    ddp_base.DetermineBasicLinearSystemDDPCS();
    // std::cout << "[PROFILE][DDPCS] Construct base DDPModelCS + indexing + basic system: "
    //           << omp_get_wtime() - t0 << " sec" << std::endl;

    // Prepare potential storage (pre-insert outer keys).
    // Thread-safe dans la boucle électrodes : chaque thread écrit dans un elec_nb distinct.
    std::map<int, std::map<double, Vector> > potential_by_elec;
    for (std::map<int, std::pair<double,double> >::const_iterator it = source_term_positions.begin();
         it != source_term_positions.end(); ++it) {
        potential_by_elec[it->first] = std::map<double, Vector>();
    }

    // Convert source_term_positions map → vector for indexed parallel for
    std::vector<std::pair<int, std::pair<double,double> > > elec_vec(
        source_term_positions.begin(), source_term_positions.end());
    const int nb_elec_total = (int)elec_vec.size();

#ifdef DDP_CS_PARALLEL_SOURCES
    const int nb_threads_max = DDP_CS_NB_THREADS;
#else
    const int nb_threads_max = 1;
#endif

    // Accumulator timers (wall-time, summed over threads via reduction)
    double t_acc_ddp_copy       = 0.0;
    double t_acc_ddp_fourier    = 0.0;
    double t_acc_ddp_precompute = 0.0;
    double t_acc_ddp_bc_CS      = 0.0;
    double t_acc_ddp_update     = 0.0;
    double t_acc_ddp_solve      = 0.0;
    double t_acc_ddp_recon      = 0.0;
    double t_acc_ddp_build      = 0.0;
    int    nb_solves_ddp = 0;

    double t_wall_loops = omp_get_wtime();

    // 4. Boucle fréquences (séquentielle) — précomputation Y parallélisée en interne.
    for (int wi = 0; wi < n_freq; ++wi) {
        const double current_w = w_vect[wi];
        // std::cout << "\n[PROFILE][DDPCS] === Frequency w = " << current_w << " ===" << std::endl;

        // 4a. Copie du modèle de base et affectation de la fréquence
        t0 = omp_get_wtime();
        DDPModelCS ddp_w = ddp_base;
        ddp_w.current_w = current_w;
        ddp_w.w = current_w;
        t_acc_ddp_copy += omp_get_wtime() - t0;

        // 4b. Mise à jour des blocs dépendant de w
        t0 = omp_get_wtime();
        ddp_w.UpdateDDPSCBlocksFourier();
        t_acc_ddp_fourier += omp_get_wtime() - t0;

        // 4c. Précomputation de la condensation statique.
        // La boucle colonnes de Y est parallélisée en interne (DDPModelCS.cpp).
        // matrix_Y_DDP (~266 MB) est libérée à la fin → copies ddp_per_thread légères.
        t0 = omp_get_wtime();
        ddp_w.PrecomputeDDPSCForCurrentFrequency();
        t_acc_ddp_precompute += omp_get_wtime() - t0;

        // vector_zi_DDP : seul champ cumulatif modifié par électrode → à restaurer.
        // factor_Ai1plusi2 (shared_ptr) : solve() est const → lecture concurrente sûre.
        const Vector vector_zi_DDP_base = ddp_w.vector_zi_DDP;

        // Un DDPModelCS par thread, initialisé depuis l'état gelé après précomputation.
        std::vector<DDPModelCS> ddp_per_thread(nb_threads_max, ddp_w);

        // Timers par-fréquence (réduits sur les threads via reduction)
        double t_bc_par=0, t_upd_par=0, t_sol_par=0, t_rec_par=0, t_bld_par=0;
        int    nb_sol_par = 0;
        double t_wall_elec = omp_get_wtime();

        // 4d. Boucle électrodes — parallélisée si DDP_CS_PARALLEL_SOURCES.
        // Pas de race condition : chaque thread écrit dans potential_by_elec[elec_nb]
        // avec un elec_nb distinct (un seul thread par électrode).
#ifdef DDP_CS_PARALLEL_SOURCES
#       pragma omp parallel for num_threads(DDP_CS_NB_THREADS) \
            reduction(+:t_bc_par) reduction(+:t_upd_par) reduction(+:t_sol_par) \
            reduction(+:t_rec_par) reduction(+:t_bld_par) reduction(+:nb_sol_par)
#endif
        for (int ei = 0; ei < nb_elec_total; ++ei) {
            const int    elec_nb = elec_vec[ei].first;
            const double sx      = elec_vec[ei].second.first;
            const double sy      = elec_vec[ei].second.second;

            DDPModelCS& ddp_local = ddp_per_thread[omp_get_thread_num()];

            // Restauration du seul champ cumulatif
            ddp_local.vector_zi_DDP = vector_zi_DDP_base;

            double t_loc;

            t_loc = omp_get_wtime();
            SourceTermsDef source_terms_local;
            source_terms_local.push_back(SourceTerm(sx, sy, source_value));
            SourceTermsCSEigen source_terms_epm_local;
            ddp_local.bc_def_epm = BoundaryConditionsDefinitionCSEigen(
                domain, Nx, Ny, bc_map,
                source_terms_local, source_terms_epm_local,
                fourier, current_w);
            ddp_local.source_terms_epm = source_terms_epm_local;
            t_bc_par += omp_get_wtime() - t_loc;

            t_loc = omp_get_wtime();
            ddp_local.UpdateDDPSCBoundaryConditionsAndSources(fourier);
            t_upd_par += omp_get_wtime() - t_loc;

            t_loc = omp_get_wtime();
            ddp_local.SolveCondensedDDPSystem();
            t_sol_par += omp_get_wtime() - t_loc;

            t_loc = omp_get_wtime();
            ddp_local.ReconstructDDPInteriorFromBorder();
            t_rec_par += omp_get_wtime() - t_loc;

            t_loc = omp_get_wtime();
            potential_by_elec.at(elec_nb)[current_w] = ddp_local.BuildGlobalSolutionVectorFromSC();
            t_bld_par += omp_get_wtime() - t_loc;

            ++nb_sol_par;
        }

        double t_wall_elec_elapsed = omp_get_wtime() - t_wall_elec;
        double t_sum_elec = t_bc_par + t_upd_par + t_sol_par + t_rec_par + t_bld_par;

        t_acc_ddp_bc_CS  += t_bc_par;
        t_acc_ddp_update += t_upd_par;
        t_acc_ddp_solve  += t_sol_par;
        t_acc_ddp_recon  += t_rec_par;
        t_acc_ddp_build  += t_bld_par;
        nb_solves_ddp    += nb_sol_par;

        // std::cout << "[PROFILE][DDPCS]   Boucle elec w=" << current_w
        //           << " : mur=" << t_wall_elec_elapsed << " s"
        //           << "  somme threads=" << t_sum_elec << " s" << std::endl;
    }

    double t_wall_loops_elapsed = omp_get_wtime() - t_wall_loops;

    // 5. Write results to files — uniquement en mode exécutable (MAIN.cpp) :
    // en mode librairie (PERFORM.cpp, build_result_map==true), le potentiel
    // est transmis à Python via res_matrices, l'écriture serait un coût pur.
    t0 = omp_get_wtime();
    if (!build_result_map) {
    std::filesystem::create_directories(param.files_param.output_path);
#ifdef DDP_CS_PARALLEL_SOURCES
#   pragma omp parallel for num_threads(DDP_CS_NB_THREADS)
#endif
    for (int ei = 0; ei < nb_elec_total; ++ei) {
        const int elec_nb = elec_vec[ei].first;
        const std::map<double, Vector>& potential_map = potential_by_elec.at(elec_nb);
        if (potential_map.empty()) continue;

        std::stringstream ss;
        ss << elec_nb;
        const std::string output_file = param.files_param.output_path + std::string("Matrix_Potential") + ss.str() + ".txt";
        std::ofstream output(output_file.c_str(), std::ofstream::out);
        if (!output.is_open()) {
            std::cerr << "[ERROR][DDPCS] Cannot open output file: " << output_file << std::endl;
            continue;
        }
        for (std::map<double, Vector>::const_iterator it_vec = potential_map.begin();
             it_vec != potential_map.end(); ++it_vec) {
            const Vector& phi = it_vec->second;
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i)
                    output << phi(return_index(i, (Ny - 1) - j, Ny)) << " ";
                output << "\n";
            }
        }
        output.close();
    }
    }
    double t_ddp_write = omp_get_wtime() - t0;

    // 5b. Empilement des potentiels (toutes fréquences) en une matrice creuse Eigen
    // par électrode, uniquement nécessaire quand le résultat est renvoyé à l'appelant
    // (PERFORM.cpp / pybind) : MAIN.cpp n'en a pas besoin (résultats déjà écrits sur disque).
    if (build_result_map) {
#ifdef DDP_CS_PARALLEL_SOURCES
#       pragma omp parallel for num_threads(DDP_CS_NB_THREADS)
#endif
        for (int ei = 0; ei < nb_elec_total; ++ei) {
            const int elec_nb = elec_vec[ei].first;
            const std::map<double, Vector>& potential_map = potential_by_elec.at(elec_nb);
            if (potential_map.empty()) continue;

            DenseMatrix A_dense(n_freq * Nx, Ny);
            A_dense.setZero();
            int cpt_w = 0;
            for (std::map<double, Vector>::const_iterator it_vec = potential_map.begin();
                 it_vec != potential_map.end(); ++it_vec) {
                const Vector& phi = it_vec->second;
                for (int i = 0; i < Nx; ++i)
                    for (int j = 0; j < Ny; ++j)
                        A_dense(j + Ny * cpt_w, i) = phi(return_index(i, (Ny - 1) - j, Ny));
                cpt_w++;
            }
#           pragma omp critical
            { res_matrices[elec_nb] = A_dense.sparseView(); }
        }
    }

    // ── Profiling final ──────────────────────────────────────────────────────
    const double t_boucle_freq   = t_acc_ddp_copy + t_acc_ddp_fourier + t_acc_ddp_precompute;
    const double t_boucle_source = t_acc_ddp_bc_CS + t_acc_ddp_update + t_acc_ddp_solve
                                 + t_acc_ddp_recon + t_acc_ddp_build;
    const double t_grand_total   = omp_get_wtime() - t_global_start;

#ifdef DDP_CS_PARALLEL_SOURCES
    const double efficiency = (t_wall_loops_elapsed > 0 && nb_threads_max > 0)
        ? t_boucle_source / (t_wall_loops_elapsed * nb_threads_max) * 100.0 : 0.0;
#endif

    std::cout << "\n[PROFILE] ============ FourierPotentialDDPSourceTermCSEigen ============" << std::endl;
// #ifdef DDP_CS_PARALLEL_SOURCES
    // std::cout << "[PROFILE] Mode                              : Parallele sources (" << nb_threads_max << " threads)" << std::endl;
// #else
//     std::cout << "[PROFILE] Mode                              : Sequentiel" << std::endl;
// #endif
    std::cout << "[PROFILE] Discretisation :" << std::endl;
    std::cout << "[PROFILE]   Noeuds matrice        : " << Nx << " x " << Ny << " = " << Nx * Ny << std::endl;
    std::cout << "[PROFILE]   Noeuds fractures       : " << ddp_base.nb_nodes_fract
              << "  (int: " << ddp_base.Nfi << ", bord: " << ddp_base.Nfb << ")" << std::endl;
    std::cout << "[PROFILE]   Noeuds interieurs tot  : " << ddp_base.Nfi + ddp_base.Ni
              << "  (fract: " << ddp_base.Nfi << ", matrice: " << ddp_base.Ni << ")" << std::endl;
    std::cout << "[PROFILE]   Noeuds de bord tot     : " << ddp_base.Nfb + ddp_base.Nb
              << "  (fract: " << ddp_base.Nfb << ", matrice: " << ddp_base.Nb << ")" << std::endl;

//     std::cout << "[PROFILE] Boucle frequences (" << n_freq << " freq) :" << std::endl;
//     std::cout << "[PROFILE]   Copie base + set w              : " << t_acc_ddp_copy       << " s  (moy: " << t_acc_ddp_copy      /n_freq << " s/freq)" << std::endl;
//     std::cout << "[PROFILE]   UpdateDDPSCBlocksFourier        : " << t_acc_ddp_fourier    << " s  (moy: " << t_acc_ddp_fourier   /n_freq << " s/freq)" << std::endl;
//     std::cout << "[PROFILE]   PrecomputeDDPSCForFrequency     : " << t_acc_ddp_precompute << " s  (moy: " << t_acc_ddp_precompute/n_freq << " s/freq)" << std::endl;
//     std::cout << "[PROFILE]   Sous-total boucle freq          : " << t_boucle_freq << " s" << std::endl;

//     std::cout << "[PROFILE] Boucle sources (" << nb_solves_ddp << " inversions) :" << std::endl;
//     std::cout << "[PROFILE]   [Les temps ci-dessous sont la SOMME sur tous les threads]" << std::endl;
//     std::cout << "[PROFILE]   BCDef (BoundaryConditions)      : " << t_acc_ddp_bc_CS  << " s  (moy: " << t_acc_ddp_bc_CS /nb_solves_ddp << " s/inv)" << std::endl;
//     std::cout << "[PROFILE]   UpdateDDPSCBoundaryConditions   : " << t_acc_ddp_update << " s  (moy: " << t_acc_ddp_update/nb_solves_ddp << " s/inv)" << std::endl;
//     std::cout << "[PROFILE]   SolveCondensedDDPSystem         : " << t_acc_ddp_solve  << " s  (moy: " << t_acc_ddp_solve /nb_solves_ddp << " s/inv)" << std::endl;
//     std::cout << "[PROFILE]   ReconstructDDPInterior          : " << t_acc_ddp_recon  << " s  (moy: " << t_acc_ddp_recon /nb_solves_ddp << " s/inv)" << std::endl;
//     std::cout << "[PROFILE]   BuildGlobalSolutionVector       : " << t_acc_ddp_build  << " s  (moy: " << t_acc_ddp_build /nb_solves_ddp << " s/inv)" << std::endl;
//     std::cout << "[PROFILE]   Somme temps threads             : " << t_boucle_source << " s" << std::endl;
//     std::cout << "[PROFILE]   Temps mur boucles (freq+elec)   : " << t_wall_loops_elapsed << " s" << std::endl;
// #ifdef DDP_CS_PARALLEL_SOURCES
//     std::cout << "[PROFILE]   Efficacite parallelisation      : " << efficiency << " %" << std::endl;
// #endif

//     std::cout << "[PROFILE] Ecriture resultats                : " << t_ddp_write << " s" << std::endl;
//     std::cout << "[PROFILE] GRAND TOTAL (mur)                 : " << t_grand_total << " s" << std::endl;
//     std::cout << "[PROFILE] ======================================================================\n" << std::endl;

    return res_matrices;
}

}  // namespace eigen_cs