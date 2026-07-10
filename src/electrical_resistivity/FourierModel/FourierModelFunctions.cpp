/*
 * FourierModelFunctions.cpp
 *
 *  Created on: Apr 25, 2013
 *      Author: delphineroubinet
 */

#include "FourierModelFunctions.h"
#include "../EPMModel/EPMDefinition.h"
#include "../InputOutput/Results.h"
#include "../DDPModel/DDPModel.h"
#include "../DDPModel/SystLinearFunctions.h"
#include "../DDPModel/SystLinearFunctionsEigen.h"
#include "../DFNModel/DFNGeneration.h"
#include <omp.h>

using namespace std;

// Building and returning the stacked matrix map is expensive and not used by MAIN.cpp.
// Keep it disabled by default to restore post-processing performance.
static const bool kBuildFourierMatrixMap = false;

std::map<int, ublas_matrix > FourierPotentialEPMSourceTerm(Parameters param, Domain domain, bool singularity, bool bc_on_border,
                                   map<int, pair<double, double> > source_term_positions, bool build_result_map) {
    // Parameters
    std::map<int, ublas_matrix > res_matrices;

    vector<double> w_vect = param.simu_param.EPM_param.w_vect;
    bool fourier = true;
    int elec_nb;
    double source_position_x, source_position_y, source_value = 1.0;
    // --- Timer hors boucles ---
    clock_t t0 = clock();
    EPMSystem epd_model_init(param), epd_model;
    epd_model_init.DetermineBasicLinearSystem();
    double t_model_init = (double)(clock() - t0) / CLOCKS_PER_SEC;

    // --- Timers accumulés sur les boucles ---
    double t_total_copy   = 0.0;
    double t_total_fourier = 0.0;
    double t_total_bcdef  = 0.0;
    double t_total_update = 0.0;
    double t_total_solve  = 0.0;
    double t_total_write  = 0.0;
    int nb_solves = 0;

    for (map<int, pair<double, double> >::iterator it1 = source_term_positions.begin();
         it1 != source_term_positions.end(); it1++) {
        elec_nb = it1->first;
        source_position_x = it1->second.first, source_position_y = it1->second.second;
        if (param.simu_param.source_terms.size() == 0) {
            param.simu_param.source_terms.push_back(SourceTerm(source_position_x, source_position_y, source_value));
        } else { param.simu_param.source_terms[0] = SourceTerm(source_position_x, source_position_y, source_value); }

        map<double, ublas_vector> Potential_Vector;
        for (vector<double>::iterator it2 = w_vect.begin(); it2 != w_vect.end(); it2++) {
            param.simu_param.EPM_param.w_current = *it2;

            t0 = clock();
            epd_model = epd_model_init;
            epd_model.w = param.simu_param.EPM_param.w_current;
            t_total_copy += (double)(clock() - t0) / CLOCKS_PER_SEC;

            t0 = clock();
            epd_model.UpdateLinearSystemFourier();
            t_total_fourier += (double)(clock() - t0) / CLOCKS_PER_SEC;

            t0 = clock();
            SourceTermsEPM source_terms_epm_;
            epd_model.bc_def_epm = BoundaryConditionsDefinition(domain, epd_model.Nx, epd_model.Ny,
                                                                param.simu_param.bc_map, param.simu_param.source_terms,
                                                                source_terms_epm_, fourier, epd_model.w);
            epd_model.source_terms_epm = source_terms_epm_;
            t_total_bcdef += (double)(clock() - t0) / CLOCKS_PER_SEC;

            t0 = clock();
            epd_model.UpdateLinearSystemBoundaryConditionsAndSourceTerms(fourier);
            t_total_update += (double)(clock() - t0) / CLOCKS_PER_SEC;

            t0 = clock();
            Potential_Vector[*it2] = LinearSystemSolving(epd_model.matrix_syst, epd_model.vector_syst);
            t_total_solve += (double)(clock() - t0) / CLOCKS_PER_SEC;
            ++nb_solves;
        }

        t0 = clock();
        stringstream ss;
        string str_cpt;
        ss << elec_nb;
        str_cpt = ss.str();
        string output_file = "Results" + str_cpt + ".txt";
        Results results(param.files_param.output_path + output_file, Potential_Vector, epd_model.Nx, epd_model.Ny);
        t_total_write += (double)(clock() - t0) / CLOCKS_PER_SEC;

        if (kBuildFourierMatrixMap || build_result_map) {
            // Optional compatibility path for external post-processing consumers.
            ublas::matrix<double> A_dense(w_vect.size() * epd_model.Nx, epd_model.Ny, 0.0);
            ublas_vector current_vector;
            int cpt_w = 0;
            for (map<double, ublas_vector>::iterator it = Potential_Vector.begin(); it != Potential_Vector.end(); ++it) {
                current_vector = it->second;
                for (int i = 0; i < epd_model.Nx; i++) {
                    for (int j = 0; j < epd_model.Ny; j++) {
                        A_dense(j + epd_model.Ny * cpt_w, i) =
                            current_vector(return_index(i, (epd_model.Ny - 1) - j, epd_model.Ny));
                    }
                }
                cpt_w++;
            }

            ublas_matrix A(w_vect.size() * epd_model.Nx, epd_model.Ny);
            for (size_t i = 0; i < A_dense.size1(); ++i) {
                for (size_t j = 0; j < A_dense.size2(); ++j) {
                    if (A_dense(i, j) != 0.0) {
                        A.insert_element(i, j, A_dense(i, j));
                    }
                }
            }
            res_matrices[elec_nb] = A;
        }
    }
    double t_loops = t_total_copy + t_total_fourier + t_total_bcdef + t_total_update + t_total_solve + t_total_write;
    cout << "\n[PROFILE] ============ FourierPotentialEPMSourceTerm ============" << endl;
    cout << "[PROFILE] Hors boucles :" << endl;
    cout << "[PROFILE]   EPMSystem init + BasicSystem : " << t_model_init << " sec" << endl;
    cout << "[PROFILE] Dans les boucles (cumul sur " << nb_solves << " inversions) :" << endl;
    cout << "[PROFILE]   Copie modele                 : " << t_total_copy    << " sec" << endl;
    cout << "[PROFILE]   UpdateLinearSystemFourier    : " << t_total_fourier << " sec" << endl;
    cout << "[PROFILE]   BoundaryConditionsDefinition : " << t_total_bcdef   << " sec" << endl;
    cout << "[PROFILE]   UpdateLinearSystem BCs+Src   : " << t_total_update  << " sec" << endl;
    cout << "[PROFILE]   LinearSystemSolving (solve)  : " << t_total_solve   << " sec  (moy: " << t_total_solve/nb_solves << " sec)" << endl;
    cout << "[PROFILE]   Ecriture resultats           : " << t_total_write   << " sec" << endl;
    cout << "[PROFILE]   Sous-total boucles           : " << t_loops << " sec" << endl;
    cout << "[PROFILE]   TOTAL (hors+boucles)         : " << t_model_init + t_loops << " sec" << endl;
    cout << "[PROFILE] ========================================================\n" << endl;

    if (param.simu_param.debug_mode){
        cout << "DEBUG: Fourier matrix-map export enabled: " << (kBuildFourierMatrixMap ? "yes" : "no") << endl;
        if (!res_matrices.empty()) {
            map<int, ublas_matrix>::iterator it_first = res_matrices.begin();
            cout << "DEBUG: first matrix dimensions: " << it_first->second.size1() << "x" << it_first->second.size2() << endl;
        }
    }
    return res_matrices;
}

void FourierPotentialEPMSourceTerm(Parameters param, Domain domain, bool singularity, bool bc_on_border) {

    vector<double> w_vect = param.simu_param.EPM_param.w_vect;
    if (param.simu_param.source_terms.size() != 1) {
        cout << "WARNING in main (PERFORM.cpp): configuration not implemented" << endl;
    }

    double source_position_x = param.simu_param.source_terms[0].x_position;
    double source_position_y = param.simu_param.source_terms[0].y_position;
    double source_value = param.simu_param.source_terms[0].source_value;
    double Lsimu = domain.domain_size_x();
    double source_spacing = param.simu_param.source_terms[0].source_spacing;
    int cpt = 0;
    bool continue_loop = true, fourier = true;

    EPMSystem epd_model_init(param), epd_model;
    epd_model_init.DetermineBasicLinearSystem();

    while (continue_loop) {
        cpt++;
        //cout << "Experiment " << cpt << endl;

        param.simu_param.source_terms[0] = SourceTerm(source_position_x, source_position_y, source_value);

        map<double, ublas_matrix> Potential_Matrices;
        for (vector<double>::iterator it1 = w_vect.begin(); it1 != w_vect.end(); it1++) {
            //cout << "w=" << *it1 << endl;
            param.simu_param.EPM_param.w_current = *it1;

            // update the EPM model and linear system
            epd_model = epd_model_init;
            epd_model.w = param.simu_param.EPM_param.w_current;
            // Fourier terms
            epd_model.UpdateLinearSystemFourier();
            // Define the boundary conditions and source
            SourceTermsEPM source_terms_epm_;
            epd_model.bc_def_epm = BoundaryConditionsDefinition(domain, epd_model.Nx, epd_model.Ny,
                                                                param.simu_param.bc_map, param.simu_param.source_terms,
                                                                source_terms_epm_, fourier, epd_model.w);
            epd_model.source_terms_epm = source_terms_epm_;
            // Update the linear system with the boundary conditions and source terms
            epd_model.UpdateLinearSystemBoundaryConditionsAndSourceTerms(fourier);

            Potential_Matrices[*it1] = epd_model.PotentialComputation();
        }

        stringstream ss;
        string str_cpt;
        ss << cpt;
        str_cpt = ss.str();

        string output_file = "EPM_Xu" + str_cpt;
        output_file = output_file + ".txt";

        Results results(param.files_param.output_path + output_file, Potential_Matrices);

        source_position_x += source_spacing;
        if (source_position_x < Lsimu && source_spacing != -1) { continue_loop = true; }
        else { continue_loop = false; }

    }
}

void FourierPotentialEPM(Parameters param, Domain domain, bool singularity) {

    vector<double> w_vect = param.simu_param.EPM_param.w_vect;
    map<double, ublas_matrix> Potential_Matrices;
    for (vector<double>::iterator it1 = w_vect.begin(); it1 != w_vect.end(); it1++) {

        param.simu_param.EPM_param.w_current = *it1;
        EPMSystem epd_model(param, domain, true, singularity);
        Potential_Matrices[*it1] = epd_model.PotentialComputation();
    }

    Results results(param.files_param.output_path + param.files_param.output_file, Potential_Matrices);
}

void FourierPotentialDDP(Parameters param, Domain domain, bool singularity, bool bc_on_border) {
    // 1.1. Variables
    vector<double> w_vect = param.simu_param.DDP_param.w_vect;
    map<double, ublas_matrix> Potential_Matrices;
    map<double, ublas_vector> Potential_Fractures;
    ublas_vector DFN_potential;
    ublas_matrix matrix_potential;
    // 1.2. Fracture network generation
    NetworkMeshes net_mesh(domain, param, param.simu_param.simu_option);
    bool fourier = true;
    clock_t t1 = clock();
    DDPModel ddp_model(param, domain, net_mesh, fourier, singularity);
    clock_t t2 = clock();
    cout << "time for generating the fracture network = t2-t1" << endl;
    // 2. Simulation for each wave number
    for (vector<double>::iterator it1 = w_vect.begin(); it1 != w_vect.end(); it1++) {
        // 2.1. Update the current wave number
        param.simu_param.DDP_param.w_current = *it1;
        // 2.2. Define the boundary conditions for ddp model (depend on w)
        ddp_model.DDPModelFourierBC(param, domain, net_mesh.domain, net_mesh.border_map, net_mesh.nodes_map);
        // 2.3. Compute electrical potential
        cout << "w=" << ddp_model.current_w << endl;
        ddp_model.ComputePotentialFourier(DFN_potential, matrix_potential, bc_on_border);
        // 2.4. Store the computed values
        Potential_Matrices[*it1] = matrix_potential;
        Potential_Fractures[*it1] = DFN_potential;
    }
    Results results(param.files_param.output_path + param.files_param.output_file, Potential_Matrices);
    Results results1(param.files_param.output_path + param.files_param.output_file1, Potential_Fractures);
}

std::map<int, ublas_matrix > FourierPotentialDDPSourceTerm(Parameters param, Domain domain, bool singularity,
                                   map<int, pair<double, double> > source_term_positions, bool build_result_map) {
    // 1. Global variables and fracture network generation
    // 1.1. Variables
    std::map<int, ublas_matrix > res_matrices;
    vector<double> w_vect = param.simu_param.DDP_param.w_vect;
    // --- Timers hors boucles ---
    clock_t t1 = clock();
    NetworkMeshes net_mesh(domain, param, param.simu_param.simu_option);
    double t_net_mesh = (double)(clock() - t1) / CLOCKS_PER_SEC;

    bool fourier = true;
    double source_value = 1.0, source_position_x, source_position_y;
    int nb_elec;

    t1 = clock();
    DDPModel ddp_model_init(param, domain, net_mesh, fourier, singularity);
    ddp_model_init.DetermineBasicLinearSystemDDP();
    double t_model_init = (double)(clock() - t1) / CLOCKS_PER_SEC;

    // --- Timers accumulés sur les boucles ---
    double t_total_copy   = 0.0;  // copie matrix_syst/vector_syst
    double t_total_bcdef  = 0.0;  // DDPModelFourierBC
    double t_total_update = 0.0;  // UpdateLinearSystemFourier
    double t_total_solve  = 0.0;  // EvaluatePotential
    double t_total_write  = 0.0;  // écriture résultats
    int nb_solves = 0;

    map<int, pair<double, double> >::iterator it1;
    vector<double>::iterator it2;
    map<double, ublas_matrix> Potential_Matrices;
    map<double, ublas_vector> Potential_Fractures;
    ublas_vector DFN_potential;
    ublas_matrix matrix_potential;
    string str_cpt, output_file_matrix;
    Results results;
    BoundaryConditionsEPM bc_def_epm = ddp_model_init.bc_def_epm, bc_def_epm_sing = ddp_model_init.bc_def_epm_sing;
    BoundaryConditionsDFN bc_map = ddp_model_init.bc_map;
    int Nx = ddp_model_init.Nx, Ny = ddp_model_init.Ny, nb_nodes_fract = ddp_model_init.nb_nodes_fract;
    double delta_x = ddp_model_init.delta_x, delta_y = ddp_model_init.delta_y;
    ublas_matrix matrix_syst;
    ublas_vector vector_syst;
    ublas_matrix porous_cond = ddp_model_init.porous_cond;
    SourceTermsEPM source_terms_epm = ddp_model_init.source_terms_epm, source_terms_epm_sing = ddp_model_init.source_terms_epm_sing;
    SubNetworkMap subnetwork_map = ddp_model_init.subnetwork_map;

#pragma omp parallel
    {
#pragma omp single nowait
// #pragma omp single
       {
            for (it1 = source_term_positions.begin(); it1 != source_term_positions.end(); it1++) {
// #pragma omp task firstprivate(res_matrices, it1, it2, nb_elec, source_position_x, source_position_y, param, Potential_Matrices, Potential_Fractures, DFN_potential, matrix_potential, str_cpt, output_file_matrix, results, matrix_syst, vector_syst, bc_def_epm, bc_map, bc_def_epm_sing, source_terms_epm, source_terms_epm_sing, subnetwork_map)
#pragma omp task firstprivate(it1, it2, nb_elec, source_position_x, source_position_y, param, Potential_Matrices, Potential_Fractures, DFN_potential, matrix_potential, str_cpt, output_file_matrix, results, matrix_syst, vector_syst, bc_def_epm, bc_map, bc_def_epm_sing, source_terms_epm, source_terms_epm_sing, subnetwork_map)
                {
                    nb_elec = it1->first;
                    source_position_x = it1->second.first;
                    source_position_y = it1->second.second;
                    if (param.simu_param.source_terms.size() == 0) {
                        param.simu_param.source_terms.push_back(
                                SourceTerm(source_position_x, source_position_y, source_value));
                    } else {
                        param.simu_param.source_terms[0] = SourceTerm(source_position_x, source_position_y,
                                                                      source_value);
                    }

                    source_terms_epm.clear();
                    source_terms_epm_sing.clear();
                    map<double, ublas_matrix> Potential_Matrices;
                    map<double, ublas_vector> Potential_Fractures;
                    ublas_vector DFN_potential;
                    ublas_matrix matrix_potential;

                    double loc_copy = 0.0, loc_bcdef = 0.0, loc_update = 0.0, loc_solve = 0.0;
                    int loc_solves = 0;
                    for (it2 = w_vect.begin(); it2 != w_vect.end(); it2++) {
                        param.simu_param.DDP_param.w_current = *it2;

                        double t0 = omp_get_wtime();
                        matrix_syst = ddp_model_init.matrix_syst;
                        vector_syst = ddp_model_init.vector_syst;
                        loc_copy += omp_get_wtime() - t0;

                        t0 = omp_get_wtime();
                        DDPModelFourierBC(singularity, param, domain, net_mesh.domain, net_mesh.border_map,
                                          net_mesh.nodes_map, source_terms_epm_sing, source_terms_epm, bc_map,
                                          bc_def_epm_sing, bc_def_epm);
                        loc_bcdef += omp_get_wtime() - t0;

                        t0 = omp_get_wtime();
                        UpdateLinearSystemFourier(Nx, Ny, delta_x, delta_y, nb_nodes_fract, *it2, bc_def_epm, bc_map,
                                                  source_terms_epm, porous_cond, subnetwork_map, singularity,
                                                  matrix_syst, vector_syst);
                        loc_update += omp_get_wtime() - t0;

                        t0 = omp_get_wtime();
                        EvaluatePotential(nb_nodes_fract, Nx, Ny, matrix_syst, vector_syst, DFN_potential,
                                          matrix_potential);
                        loc_solve += omp_get_wtime() - t0;
                        ++loc_solves;

                        Potential_Matrices[*it2] = matrix_potential;
                        Potential_Fractures[*it2] = DFN_potential;
                    }

                    double t0 = omp_get_wtime();
                    str_cpt = to_string(nb_elec);
                    output_file_matrix = "Matrix_Potential" + str_cpt + ".txt";
                    results = Results(param.files_param.output_path + output_file_matrix, Potential_Matrices);
                    double loc_write = omp_get_wtime() - t0;

#pragma omp atomic
                    t_total_copy   += loc_copy;
#pragma omp atomic
                    t_total_bcdef  += loc_bcdef;
#pragma omp atomic
                    t_total_update += loc_update;
#pragma omp atomic
                    t_total_solve  += loc_solve;
#pragma omp atomic
                    t_total_write  += loc_write;
#pragma omp atomic
                    nb_solves      += loc_solves;

                    DFN_potential.clear();

                    if (kBuildFourierMatrixMap || build_result_map) {
                        ublas::matrix<double> A_dense(w_vect.size() * Nx, Ny, 0.0);
                        ublas_matrix current_matrix;
                        int cpt_w = 0;
                        for (map<double, ublas_matrix>::iterator it = Potential_Matrices.begin(); it != Potential_Matrices.end(); ++it) {
                            current_matrix = it->second;
                            for (int i = 0; i < Nx; i++)
                                for (int j = 0; j < Ny; j++)
                                    A_dense(j + Ny * cpt_w, i) = current_matrix(i, (Ny - 1) - j);
                            cpt_w++;
                        }
                        ublas_matrix A(w_vect.size() * Nx, Ny);
                        for (size_t i = 0; i < A_dense.size1(); ++i)
                            for (size_t j = 0; j < A_dense.size2(); ++j)
                                if (A_dense(i, j) != 0.0)
                                    A.insert_element(i, j, A_dense(i, j));
#pragma omp critical
                        { res_matrices[nb_elec] = A; }
                    }
                    if (param.simu_param.debug_mode)
                        cout << "DEBUG: res_matrices calculated for simulation number " << nb_elec-1 << endl;
                }
            }
        }
    }
    PrintUpdateLinearSystemFourierStats();
    double t_loops = t_total_copy + t_total_bcdef + t_total_update + t_total_solve + t_total_write;
    cout << "\n[PROFILE] ============ FourierPotentialDDPSourceTerm ============" << endl;
    cout << "[PROFILE] Hors boucles :" << endl;
    cout << "[PROFILE]   NetworkMeshes              : " << t_net_mesh   << " sec" << endl;
    cout << "[PROFILE]   DDPModel init + BasicSystem : " << t_model_init << " sec" << endl;
    cout << "[PROFILE] Dans les boucles (cumul sur " << nb_solves << " inversions) :" << endl;
    cout << "[PROFILE]   Copie matrix/vector_syst    : " << t_total_copy   << " sec" << endl;
    cout << "[PROFILE]   DDPModelFourierBC           : " << t_total_bcdef  << " sec" << endl;
    cout << "[PROFILE]   UpdateLinearSystemFourier   : " << t_total_update << " sec" << endl;
    cout << "[PROFILE]   EvaluatePotential (solve)   : " << t_total_solve  << " sec  (moy: " << t_total_solve/nb_solves << " sec)" << endl;
    cout << "[PROFILE]   Ecriture resultats          : " << t_total_write  << " sec" << endl;
    cout << "[PROFILE]   Sous-total boucles          : " << t_loops << " sec" << endl;
    cout << "[PROFILE]   TOTAL (hors+boucles)        : " << t_net_mesh + t_model_init + t_loops << " sec" << endl;
    cout << "[PROFILE] ========================================================\n" << endl;

    if (param.simu_param.debug_mode){
        cout << "DEBUG: Fourier matrix-map export enabled: " << (kBuildFourierMatrixMap ? "yes" : "no") << endl;
        if (!res_matrices.empty()) {
            map<int, ublas_matrix>::iterator it_first = res_matrices.begin();
            cout << "DEBUG: first matrix dimensions: " << it_first->second.size1() << "x" << it_first->second.size2() << endl;
        }
    }
    return res_matrices;
}

void FourierPotentialDDPSourceTerm(Parameters param, Domain domain, bool singularity) {
    // 1. Global variables and fracture network generation
    // 1.1. Variables
    vector<double> w_vect = param.simu_param.DDP_param.w_vect;
    if (param.simu_param.source_terms.size() != 1) {
        cout << "WARNING in main (PERFORM.cpp): configuration not implemented" << endl;
    }
    // 1.2 Fracture network generation
    NetworkMeshes net_mesh(domain, param, param.simu_param.simu_option);
    /*cout << "NETWORK RESTRICTRED TO ITS BACKBONE AND UPPER VALUES" << endl;
    // Compute the backbone
    net_mesh=net_mesh.return_backbone(param,EPSILON);
    // Cut the domain
    cout << "NETWORK RESTRICTRED TO THE UPPER FRACTURES" << endl;
    DFNCutting(net_mesh,25);*/
    bool fourier = true;
    // 2. Simulation for each source position
    // 2.1. Parameters for the simulations
    // standard parameters and variables
    double source_value = param.simu_param.source_terms[0].source_value;
    double source_position_y = param.simu_param.source_terms[0].y_position;
    double Lsimu = domain.domain_size_x();
    int cpt = 0;
    // definition of the x-coordinate of the source term
    double source_spacing = param.simu_param.source_terms[0].source_spacing;
    vector<double> vect_source_pos_x, vect_values_I;
    ublas_matrix elec_positions;
    // several positions for the source terms defined in a text file (electrical resistivity survey)
    if (source_spacing == -1) { vect_source_pos_x = ReadSourceTermsCoordX(param); }
        // one position for the source terms
    else if (source_spacing == -2) {
        double source_position_x = param.simu_param.source_terms[0].x_position;
        vect_source_pos_x.push_back(source_position_x);
    } else {
        double source_position_x = param.simu_param.source_terms[0].x_position;
        vect_source_pos_x.push_back(source_position_x);
        //while (source_position_x<=Lsimu-4*source_spacing){
        while (source_position_x < Lsimu - source_spacing) {
            source_position_x += source_spacing;
            vect_source_pos_x.push_back(source_position_x);
        }
    }
    // Basic linear system (defined the coefficients depending on the porous domain properties and which do not depend on the Fourier variable, boundary conditions, and source terms)
    DDPModel ddp_model_init(param, domain, net_mesh, fourier, singularity), ddp_model;
    ddp_model_init.DetermineBasicLinearSystemDDP();
    ublas_vector potential;
    // 2.2. Loop for the simulations
    for (vector<double>::iterator it1 = vect_source_pos_x.begin(); it1 != vect_source_pos_x.end(); it1++) {
        cpt++;
        //cout << "Experiment " << cpt << " and x_coord=" << *it1 << endl;
        //2.1. Variables
        param.simu_param.source_terms[0] = SourceTerm(*it1, source_position_y, source_value);
        map<double, ublas_matrix> Potential_Matrices;
        map<double, ublas_vector> Potential_Fractures;
        ublas_vector DFN_potential;
        ublas_matrix matrix_potential;
        // 2.2. Simulation for each wave number
        for (vector<double>::iterator it1 = w_vect.begin(); it1 != w_vect.end(); it1++) {
            ddp_model = ddp_model_init;
            // 2.1. Update the current wave number
            param.simu_param.DDP_param.w_current = *it1;
            ddp_model.current_w = *it1;
            //cout << "w=" << ddp_model.current_w << endl;
            // 2.2. Define the boundary conditions for ddp model (depend on w)
            ddp_model.DDPModelFourierBC(param, domain, net_mesh.domain, net_mesh.border_map, net_mesh.nodes_map);
            // 2.3. Update the linear system with the coefficients which depend on the Fourier variables, boundary conditions, and source terms
            ddp_model.UpdateLinearSystemFourier();
            // 2.4. Compute electrical potential
            ddp_model.EvaluatePotential(DFN_potential, matrix_potential);
            // 2.4. Store the computed values
            Potential_Matrices[*it1] = matrix_potential;
            Potential_Fractures[*it1] = DFN_potential;
        }
        stringstream ss;
        string str_cpt;
        ss << cpt;
        str_cpt = ss.str();
        string output_file_matrix = "Matrix_Potential" + str_cpt;
        output_file_matrix = output_file_matrix + ".txt";
        string output_file_fract = "Fracture_Potential" + str_cpt;
        output_file_fract = output_file_fract + ".txt";
        Results results(param.files_param.output_path + output_file_matrix, Potential_Matrices);
        Results results1(param.files_param.output_path + output_file_fract, Potential_Fractures);
    }
}

// Function to read the x-coordinates of the source term
vector<double> ReadSourceTermsCoordX(Parameters param) {
    vector<double> vect_source_pos_x;
    //string file_name_coord=param.files_param.input_path+param.files_param.file_name_source_xcoord;
    string file_name_coord = param.files_param.input_path + param.files_param.file_name_source_positions;
    ifstream file(file_name_coord.c_str());
    if (file.is_open()) {
        int Nb_coord;
        double x_coord;
        file >> Nb_coord;
        for (int i = 0; i < Nb_coord; i++) {
            file >> x_coord;
            vect_source_pos_x.push_back(x_coord);
        }
    } else {
        cout << "WARNING in ReadSourceTermsCoordX (FourierModelFunctions.cpp): file not found" << endl;
    }
    return vect_source_pos_x;
}

// ─────────────────────────────────────────────────────────────────────────────
// DDP_EIGEN : même algorithme que DDP classique, containers Eigen natifs
// Parallélisation sur les électrodes (sources) :
//   → Décommenter la ligne suivante pour activer, recommenter pour désactiver.
#define DDP_EIGEN_PARALLEL_SOURCES
// ─────────────────────────────────────────────────────────────────────────────
std::map<int, eigen_cs::SparseMatrix> FourierPotentialDDPSourceTermEigen(
    Parameters param, Domain domain, bool singularity,
    map<int, pair<double, double>> source_term_positions,
    bool build_result_map)
{
    std::map<int, eigen_cs::SparseMatrix> res_matrices;
    vector<double> w_vect = param.simu_param.DDP_param.w_vect;
    bool fourier = true;
    double source_value = 1.0;

    // ── 1. Réseau de fractures ──────────────────────────────────────────────
    double t0 = omp_get_wtime();
    NetworkMeshes net_mesh(domain, param, param.simu_param.simu_option);
    double t_net = omp_get_wtime() - t0;

    // ── 2. Système de base Eigen (stencil poreux + pattern DFN) ────────────
    t0 = omp_get_wtime();
    DDPModel ddp_model_init(param, domain, net_mesh, fourier, singularity);
    //double t_ctor = omp_get_wtime() - t0;

    //t0 = omp_get_wtime();
    ddp_model_init.DetermineBasicLinearSystemDDP();
    //double t_basic_ublas = omp_get_wtime() - t0;

    //t0 = omp_get_wtime();
    DDPEigenBase base = DetermineBasicLinearSystemDDPEigen(ddp_model_init);
    //double t_basic_eigen = omp_get_wtime() - t0;

    // ── 3. Matrice complète pour analyzePattern (1ère élec + 1ère fréq) ────
    //t0 = omp_get_wtime();
    // analyzePattern sur une matrice avec BCs donne un ordering AMD optimal :
    // les nœuds Dirichlet sont déjà isolés dans le graphe de dépendance.
    eigen_cs::SparseMatrix mat_pattern = base.matrix_base;
    if (!source_term_positions.empty() && !w_vect.empty()) {
        auto   it0  = source_term_positions.begin();
        double sx0  = it0->second.first, sy0 = it0->second.second;
        Parameters param0 = param;
        if (param0.simu_param.source_terms.empty())
            param0.simu_param.source_terms.push_back(SourceTerm(sx0, sy0, source_value));
        else
            param0.simu_param.source_terms[0] = SourceTerm(sx0, sy0, source_value);
        param0.simu_param.DDP_param.w_current = w_vect.front();

        SourceTermsEPM     src0, src0_sing;
        BoundaryConditionsEPM bc0, bc0_sing;
        BoundaryConditionsDFN bcdfn0;
        DDPModelFourierBC(singularity, param0, domain,
                          net_mesh.domain, net_mesh.border_map, net_mesh.nodes_map,
                          src0_sing, src0, bcdfn0, bc0_sing, bc0);
        eigen_cs::Vector vec0 = eigen_cs::Vector::Zero(base.n);
        UpdateLinearSystemFourierEigen(base, w_vect.front(), bc0, bcdfn0, src0,
                                       singularity, mat_pattern, vec0);
    }

    // ── 4. Un solveur par thread, pattern symbolique pré-calculé ────────────
    // En mode séquentiel : 1 solveur. En mode parallèle : autant que de threads.
#ifdef DDP_EIGEN_PARALLEL_SOURCES
    int nb_threads_max = omp_get_max_threads();
#else
    int nb_threads_max = 1;
#endif
    vector<DDPEigenSolver*> solvers(nb_threads_max);
    for (int t = 0; t < nb_threads_max; t++) {
        solvers[t] = new DDPEigenSolver();
        solvers[t]->analyzePattern(mat_pattern);
    }
    //double t_analyze = omp_get_wtime() - t0;
    double t_init = omp_get_wtime() - t0;

    // ── 5. Conversion map → vecteur pour boucle for indexée (parallel for) ─
    vector<pair<int, pair<double,double>>> elec_vec(
        source_term_positions.begin(), source_term_positions.end());
    int nb_elec_total = (int)elec_vec.size();

    // ── 6. Timers — chaque thread accumule sa propre mesure.
    //    La reduction(+:...) somme les contributions de TOUS les threads :
    //    somme = nb_threads × temps_moyen_par_thread ≈ nb_threads × (temps_mur / nb_threads)
    //    → la somme est invariante par rapport au nb de threads (= travail total réalisé).
    double t_copy=0, t_bc=0, t_assem=0, t_solve=0, t_write=0;
    int nb_solves = 0;
    double t_wall = omp_get_wtime();   // temps mur de la section parallèle

    // ── 7. Boucle électrodes ─────────────────────────────────────────────────
    // Décommenter le #pragma suivant pour activer la parallélisation :
#ifdef DDP_EIGEN_PARALLEL_SOURCES
#   pragma omp parallel for\
        reduction(+:t_copy) reduction(+:t_bc) reduction(+:t_assem) \
        reduction(+:t_solve) reduction(+:t_write) reduction(+:nb_solves)
#endif
    for (int ei = 0; ei < nb_elec_total; ei++) {

        int    nb_elec = elec_vec[ei].first;
        double sx      = elec_vec[ei].second.first;
        double sy      = elec_vec[ei].second.second;

        // Solveur propre au thread courant (thread 0 en mode séquentiel)
        DDPEigenSolver& local_solver = *solvers[omp_get_thread_num()];

        // Copie locale de param : chaque thread a sa propre source
        Parameters local_param = param;
        if (local_param.simu_param.source_terms.empty())
            local_param.simu_param.source_terms.push_back(SourceTerm(sx, sy, source_value));
        else
            local_param.simu_param.source_terms[0] = SourceTerm(sx, sy, source_value);

        // Variables privées par électrode (déclarées dans le corps → toujours privées)
        SourceTermsEPM     src, src_sing;
        BoundaryConditionsEPM bc_epm, bc_epm_sing;
        BoundaryConditionsDFN bc_dfn;
        map<double, eigen_cs::DenseMatrix> Potential_Matrices;
        map<double, eigen_cs::Vector>      Potential_Fractures;
        double t_loc;   // chrono local à cette itération

        // Boucle fréquences (séquentielle à l'intérieur de chaque thread)
        for (int wi = 0; wi < (int)w_vect.size(); wi++) {
            double w = w_vect[wi];
            local_param.simu_param.DDP_param.w_current = w;

            // a. Copie de la matrice de base
            t_loc = omp_get_wtime();
            eigen_cs::SparseMatrix mat = base.matrix_base;
            eigen_cs::Vector       vec = eigen_cs::Vector::Zero(base.n);
            t_copy += omp_get_wtime() - t_loc;

            // b. Calcul des conditions aux limites
            t_loc = omp_get_wtime();
            src.clear(); src_sing.clear();
            DDPModelFourierBC(singularity, local_param, domain,
                              net_mesh.domain, net_mesh.border_map, net_mesh.nodes_map,
                              src_sing, src, bc_dfn, bc_epm_sing, bc_epm);
            t_bc += omp_get_wtime() - t_loc;

            // c. Assemblage des termes w-dépendants
            t_loc = omp_get_wtime();
            UpdateLinearSystemFourierEigen(base, w, bc_epm, bc_dfn, src,
                                           singularity, mat, vec);
            t_assem += omp_get_wtime() - t_loc;

            // d. Résolution (factorisation numérique, pattern symbolique réutilisé)
            t_loc = omp_get_wtime();
            eigen_cs::Vector      phi_dfn;
            eigen_cs::DenseMatrix phi_mat;
            EvaluatePotentialEigen(base.nb_nodes_fract, base.Nx, base.Ny,
                                   local_solver, mat, vec, phi_dfn, phi_mat);
            t_solve += omp_get_wtime() - t_loc;
            ++nb_solves;

            Potential_Matrices[w] = phi_mat;
            Potential_Fractures[w] = phi_dfn;
        }

        // e. Écriture résultats (fichier propre à chaque électrode → pas de contention)
        // Uniquement en mode exécutable (MAIN.cpp) : en mode librairie
        // (PERFORM.cpp, build_result_map==true), le potentiel est transmis
        // à Python via res_matrices, l'écriture serait un coût pur.
        if (!build_result_map) {
            t_loc = omp_get_wtime();
            Results(param.files_param.output_path
                    + "Matrix_Potential" + to_string(nb_elec) + ".txt",
                    Potential_Matrices);
            t_write += omp_get_wtime() - t_loc;
        }

        // f. Empilement des potentiels (toutes fréquences) en une matrice creuse Eigen
        // par électrode, uniquement nécessaire quand le résultat est renvoyé à l'appelant
        // (PERFORM.cpp / pybind) : MAIN.cpp n'en a pas besoin (résultats écrits sur disque).
        if (build_result_map) {
            Eigen::MatrixXd A_dense((int)w_vect.size() * base.Nx, base.Ny);
            A_dense.setZero();
            int cpt_w = 0;
            for (map<double, eigen_cs::DenseMatrix>::iterator it = Potential_Matrices.begin();
                 it != Potential_Matrices.end(); ++it) {
                const eigen_cs::DenseMatrix& current_matrix = it->second;
                for (int i = 0; i < base.Nx; i++)
                    for (int j = 0; j < base.Ny; j++)
                        A_dense(j + base.Ny * cpt_w, i) = current_matrix(i, (base.Ny - 1) - j);
                cpt_w++;
            }
#pragma omp critical
            { res_matrices[nb_elec] = A_dense.sparseView(); }
        }
    }

    double t_wall_elapsed = omp_get_wtime() - t_wall;

    // Nettoyage solveurs
    for (int t = 0; t < nb_threads_max; t++) { delete solvers[t]; }

    // ── Profiling ────────────────────────────────────────────────────────────
    double t_sum = t_copy + t_bc + t_assem + t_solve + t_write;
#ifdef DDP_EIGEN_PARALLEL_SOURCES
    double efficiency = (t_wall_elapsed > 0)
                        ? t_sum / (t_wall_elapsed * nb_threads_max) * 100.0 : 0.0;
#endif
    cout << "\n[PROFILE] ============ FourierPotentialDDPSourceTermEigen ============" << endl;
    std::cout << "[PROFILE] Discretisation :" << std::endl;
    std::cout << "[PROFILE]   Noeuds matrice        : " << base.Nx << " x " << base.Ny << " = " << base.Nx * base.Ny << std::endl;
    std::cout << "[PROFILE]   Noeuds fractures       : " << base.nb_nodes_fract << std::endl;
// #ifdef DDP_EIGEN_PARALLEL_SOURCES
//     cout << "[PROFILE] Mode : PARALLELE (" << nb_threads_max << " threads max, OMP_NUM_THREADS)" << endl;
// #else
//     cout << "[PROFILE] Mode : SEQUENTIEL" << endl;
// #endif
    // cout << "[PROFILE] Hors boucles :" << endl;
    // cout << "[PROFILE]   NetworkMeshes                      : " << t_net  << " sec" << endl;
    // cout << "[PROFILE]   DDPModel + BasicSystem + solveurs  : " << t_init << " sec" << endl;
    // //cout << "[PROFILE]     dont DDPModel ctor (DDPNetworkMeshes+ReturnNbNodes) : " << t_ctor        << " sec" << endl;
    //cout << "[PROFILE]     dont DetermineBasicLinearSystemDDP (ublas)          : " << t_basic_ublas << " sec" << endl;
    //cout << "[PROFILE]     dont DetermineBasicLinearSystemDDPEigen (Eigen)     : " << t_basic_eigen << " sec" << endl;
    //cout << "[PROFILE]     dont analyzePattern (AMD ordering)                  : " << t_analyze     << " sec" << endl;
//     cout << "[PROFILE] Dans les boucles (" << nb_solves << " inversions) :" << endl;
//     cout << "[PROFILE]   [Les temps ci-dessous sont la SOMME sur tous les threads]" << endl;
//     cout << "[PROFILE]   [Diviser par nb_threads pour estimer le temps par thread]" << endl;
//     cout << "[PROFILE]   Copie matrix (Eigen)               : " << t_copy  << " sec" << endl;
//     cout << "[PROFILE]   DDPModelFourierBC                  : " << t_bc    << " sec" << endl;
//     cout << "[PROFILE]   UpdateLinearSystemFourierEigen     : " << t_assem << " sec" << endl;
//     cout << "[PROFILE]   EvaluatePotential (solve)          : " << t_solve
//          << " sec  (moy/inversion: " << t_solve / nb_solves << " sec)" << endl;
//     cout << "[PROFILE]   Ecriture resultats                 : " << t_write << " sec" << endl;
//     cout << "[PROFILE]   -------------------------------------------------------" << endl;
//     cout << "[PROFILE]   Somme threads (travail total)      : " << t_sum          << " sec" << endl;
//     cout << "[PROFILE]   Temps mur (boucles)                : " << t_wall_elapsed << " sec" << endl;
// #ifdef DDP_EIGEN_PARALLEL_SOURCES
//     cout << "[PROFILE]   Acceleration (somme / mur)         : x" << t_sum / t_wall_elapsed << endl;
//     cout << "[PROFILE]   Efficacite parallele               : " << efficiency << " %" << endl;
// #endif
//     cout << "[PROFILE]   TOTAL (hors-boucles + temps mur)   : " << t_net + t_init + t_wall_elapsed << " sec" << endl;
//     cout << "[PROFILE] ====================================================================\n" << endl;

    return res_matrices;
}

// Function to read the source term positions
map<int, pair<double, double> > ReadSourceTermsCoordinates(Parameters param) {
    map<int, pair<double, double> > source_positions;
    string file_name_coord = param.files_param.input_path + param.files_param.file_name_source_positions;
    ifstream file(file_name_coord.c_str());
    int elec_nb;
    double coord_x, coord_y;
    if (file.is_open()) {
        while (!file.eof()) {
            file >> elec_nb >> coord_x >> coord_y;
            source_positions[elec_nb] = make_pair(coord_x, coord_y);
        }
    } else {
        cout << "WARNING in ReadSourceTermsCoordX (FourierModelFunctions.cpp): file not found" << endl;
    }
    return source_positions;
}
