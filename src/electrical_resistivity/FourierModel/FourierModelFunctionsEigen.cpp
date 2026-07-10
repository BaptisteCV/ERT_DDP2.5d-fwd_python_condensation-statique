#include "FourierModelFunctionsEigen.h"
#include "../EPMModel/EPMDefinitionEigen.h"
#include "../InputOutput/Results.h"
#include "FourierModelFunctions.h"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

namespace eigen_epm {

namespace {

eigen_cs::SparseMatrix DenseMatrixToSparse(const eigen_cs::DenseMatrix& dense) {
    std::vector<eigen_cs::Triplet> triplets;
    triplets.reserve(static_cast<size_t>(dense.rows() * dense.cols()));
    for (int i = 0; i < dense.rows(); ++i) {
        for (int j = 0; j < dense.cols(); ++j) {
            if (dense(i, j) != 0.0) {
                triplets.push_back(eigen_cs::Triplet(i, j, dense(i, j)));
            }
        }
    }

    eigen_cs::SparseMatrix sparse(dense.rows(), dense.cols());
    sparse.setFromTriplets(triplets.begin(), triplets.end());
    sparse.makeCompressed();
    return sparse;
}

}  // namespace

std::map<int, eigen_cs::SparseMatrix> FourierPotentialEPMSourceTermEigen(
    Parameters param,
    Domain domain,
    bool singularity,
    bool bc_on_border,
    std::map<int, std::pair<double, double> > source_term_positions) {
    (void)bc_on_border;

    std::map<int, eigen_cs::SparseMatrix> res_matrices;
    const bool fourier = true;
    const double source_value = 1.0;

    std::vector<double> w_vect = param.simu_param.EPM_param.w_vect;

    eigen_epm::EPMSystemEigen epd_model_init(param);
    epd_model_init.DetermineBasicLinearSystem();

    for (std::map<int, std::pair<double, double> >::iterator it1 = source_term_positions.begin();
         it1 != source_term_positions.end();
         ++it1) {
        int elec_nb = it1->first;
        double source_position_x = it1->second.first;
        double source_position_y = it1->second.second;

        if (param.simu_param.source_terms.size() == 0) {
            param.simu_param.source_terms.push_back(SourceTerm(source_position_x, source_position_y, source_value));
        } else {
            param.simu_param.source_terms[0] = SourceTerm(source_position_x, source_position_y, source_value);
        }

        std::map<double, eigen_cs::Vector> potential_vector;

        for (std::vector<double>::iterator it_w = w_vect.begin(); it_w != w_vect.end(); ++it_w) {
            param.simu_param.EPM_param.w_current = *it_w;

            eigen_epm::EPMSystemEigen epd_model = epd_model_init;
            epd_model.w = param.simu_param.EPM_param.w_current;
            epd_model.UpdateLinearSystemFourier();

            SourceTermsEPM source_terms_epm_local;
            epd_model.bc_def_epm = BoundaryConditionsDefinition(
                domain,
                epd_model.Nx,
                epd_model.Ny,
                param.simu_param.bc_map,
                param.simu_param.source_terms,
                source_terms_epm_local,
                fourier,
                epd_model.w);
            epd_model.source_terms_epm = source_terms_epm_local;

            epd_model.UpdateLinearSystemBoundaryConditionsAndSourceTerms(fourier);
            potential_vector[*it_w] = eigen_cs::LinearSystemSolving(epd_model.matrix_syst, epd_model.vector_syst);

            if (singularity) {
                cout << "WARNING in FourierPotentialEPMSourceTermEigen: singularity mode is not implemented in EPM_Eigen" << endl;
            }
        }

        // Reuse existing results exporter by converting Eigen vectors to ublas vectors.
        std::map<double, ublas_vector> potential_vector_ublas;
        for (std::map<double, eigen_cs::Vector>::iterator it_vec = potential_vector.begin(); it_vec != potential_vector.end(); ++it_vec) {
            const eigen_cs::Vector& current = it_vec->second;
            ublas_vector u(current.size());
            for (int k = 0; k < current.size(); ++k) {
                u(k) = current(k);
            }
            potential_vector_ublas[it_vec->first] = u;
        }

        std::stringstream ss;
        ss << elec_nb;
        std::string output_file = "Results" + ss.str() + ".txt";
        Results results(param.files_param.output_path + output_file, potential_vector_ublas, epd_model_init.Nx, epd_model_init.Ny);

        // Optional compatibility map build, same shape/orientation as standard path.
        int n_rows = static_cast<int>(w_vect.size()) * epd_model_init.Ny;
        int n_cols = epd_model_init.Nx;
        eigen_cs::DenseMatrix dense_result = eigen_cs::DenseMatrix::Zero(n_rows, n_cols);
        int cpt_w = 0;
        for (std::map<double, eigen_cs::Vector>::iterator it_vec = potential_vector.begin(); it_vec != potential_vector.end(); ++it_vec) {
            const eigen_cs::Vector& current = it_vec->second;
            for (int i = 0; i < epd_model_init.Nx; ++i) {
                for (int j = 0; j < epd_model_init.Ny; ++j) {
                    dense_result(j + epd_model_init.Ny * cpt_w, i) = current(return_index(i, (epd_model_init.Ny - 1) - j, epd_model_init.Ny));
                }
            }
            ++cpt_w;
        }

        res_matrices[elec_nb] = DenseMatrixToSparse(dense_result);

        // Keep return-map compatibility for Python consumers using sparse outputs.
    }

    return res_matrices;
}

void FourierPotentialEPMSourceTermEigen(Parameters param, Domain domain, bool singularity, bool bc_on_border) {
    std::map<int, std::pair<double, double> > source_term_positions = ReadSourceTermsCoordinates(param);
    FourierPotentialEPMSourceTermEigen(param, domain, singularity, bc_on_border, source_term_positions);
}

void FourierPotentialEPMEigen(Parameters param, Domain domain, bool singularity) {
    std::vector<double> w_vect = param.simu_param.EPM_param.w_vect;
    std::map<double, eigen_cs::DenseMatrix> potential_matrices;

    for (std::vector<double>::iterator it_w = w_vect.begin(); it_w != w_vect.end(); ++it_w) {
        param.simu_param.EPM_param.w_current = *it_w;

        eigen_epm::EPMSystemEigen epd_model(param, domain, true, true, singularity);
        potential_matrices[*it_w] = epd_model.PotentialComputation();
    }

    // Export format compatibility with existing Results class.
    std::map<double, ublas_matrix> potential_matrices_ublas;
    for (std::map<double, eigen_cs::DenseMatrix>::iterator it = potential_matrices.begin(); it != potential_matrices.end(); ++it) {
        const eigen_cs::DenseMatrix& m = it->second;
        ublas_matrix u(m.rows(), m.cols());
        for (int i = 0; i < m.rows(); ++i) {
            for (int j = 0; j < m.cols(); ++j) {
                u(i, j) = m(i, j);
            }
        }
        potential_matrices_ublas[it->first] = u;
    }

    Results results(param.files_param.output_path + param.files_param.output_file, potential_matrices_ublas);
}

}  // namespace eigen_epm
