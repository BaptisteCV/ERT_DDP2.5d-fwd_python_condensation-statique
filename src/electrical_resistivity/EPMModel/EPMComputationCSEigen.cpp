#include "EPMComputationCSEigen.h"
#include "../Utilitaries/Math/Operations.h"
#include <cmath>
#include <iostream>

using namespace std;

namespace eigen_cs {

namespace {

BoundaryConditions ReturnBoundaryConditionsCSEigen(Domain domain,
                                                   std::string border,
                                                   int i,
                                                   int j,
                                                   int Nx,
                                                   int Ny,
                                                   BoundaryConditionsDef bc_map,
                                                   SourceTermsDef source_terms,
                                                   bool fourier,
                                                   double w) {
    if (bc_map.count(border) == 0) {
        cout << "WARNING in ReturnBoundaryConditionsCSEigen: boundary conditions not defined" << endl;
        return BoundaryConditions();
    }

    if (bc_map[border].bc_type == DIRICHLET || bc_map[border].bc_type == NEUMANN) {
        return bc_map[border];
    }

    if (bc_map[border].bc_type == MIXED_BC) {
        double bc_value = 0.0;
        if (!fourier) {
            bc_value = ReturnMixedBoundCondValueCSEigen(source_terms, i, j, Nx, Ny, domain, border, w);
        } else {
            bc_value = ReturnMixedBoundCondValueCSEigen(source_terms, i, j, Nx, Ny, domain, border, w);
        }
        if (std::isnan(bc_value)) {
            return BoundaryConditions(DIRICHLET, 0.0);
        }
        return BoundaryConditions(MIXED_BC, bc_value);
    }

    cout << "WARNING in ReturnBoundaryConditionsCSEigen: boundary condition not implemented" << endl;
    return BoundaryConditions();
}

}  // namespace

double ReturnFVMFactorCSEigen(std::string border, double delta_x, double delta_y) {
    if (border == LEFT_BORDER || border == RIGHT_BORDER) {
        return delta_y / delta_x;
    } else if (border == BOTTOM_BORDER || border == TOP_BORDER) {
        return delta_x / delta_y;
    }
    cout << "WARNING in ReturnFVMFactorCSEigen: border not defined" << endl;
    return -1.0;
}

void BoundCondCoeffsCSEigen(std::string border,
                            BoundaryConditions bound_cond,
                            double current_cond,
                            double delta_x,
                            double delta_y,
                            double& coeff1,
                            double& coeff3) {
    coeff1 = 0.0;
    coeff3 = 0.0;

    if (bound_cond.bc_type == DIRICHLET) {
        double value = 2 * ReturnFVMFactorCSEigen(border, delta_x, delta_y);
        coeff1 = current_cond * value;
        coeff3 = coeff1 * bound_cond.bc_value;
    } else if (bound_cond.bc_type == NEUMANN) {
        coeff1 = 0.0;
        coeff3 = ReturnFVMFactorCSEigen(border, delta_x, delta_y) * current_cond * bound_cond.bc_value;
    } else if (bound_cond.bc_type == MIXED_BC) {
        coeff1 = 2 * current_cond * bound_cond.bc_value;
        if (border == LEFT_BORDER || border == RIGHT_BORDER) {
            coeff1 = coeff1 * delta_y / (2 + bound_cond.bc_value * delta_x);
        } else if (border == BOTTOM_BORDER || border == TOP_BORDER) {
            coeff1 = coeff1 * delta_x / (2 + bound_cond.bc_value * delta_y);
        } else {
            cout << "WARNING in BoundCondCoeffsCSEigen: border not defined" << endl;
        }
        coeff3 = 0.0;
    } else {
        cout << "WARNING in BoundCondCoeffsCSEigen: boundary condition not implemented" << endl;
    }
}

double ReturnMixedBoundCondValueCSEigen(SourceTermsDef source_terms,
                                        int i_cell, int j_cell,
                                        int Nx, int Ny,
                                        Domain domain,
                                        std::string border,
                                        double w) {
    double delta_x = domain.ReturnDeltaX(Nx);
    double delta_y = domain.ReturnDeltaY(Ny);
    CgalPoint2D cell_center = CenterPosition(make_pair(i_cell, j_cell), make_pair(delta_x, delta_y));

    if (source_terms.size() == 0) {
        return 0.0;
    }
    if (source_terms.size() != 1) {
        cout << "WARNING in ReturnMixedBoundCondValueCSEigen: case not implemented (multi-source)" << endl;
        return -1.0;
    }

    double x_source = source_terms.begin()->x_position;
    double y_source = source_terms.begin()->y_position;
    CgalPoint2D pt_source(x_source, y_source);

    double r, alpha;
    if (border == LEFT_BORDER) {
        cell_center = CgalPoint2D(cell_center.x() - 0.5 * delta_x, cell_center.y());
        r = distance_2D(pt_source, cell_center);
        alpha = w * boost::math::cyl_bessel_k(1, w * r) / boost::math::cyl_bessel_k(0, w * r);
        return alpha * fabs(cell_center.x() - x_source) / r;
    } else if (border == RIGHT_BORDER) {
        cell_center = CgalPoint2D(cell_center.x() + 0.5 * delta_x, cell_center.y());
        r = distance_2D(pt_source, cell_center);
        alpha = w * boost::math::cyl_bessel_k(1, w * r) / boost::math::cyl_bessel_k(0, w * r);
        return alpha * fabs(cell_center.x() - x_source) / r;
    } else if (border == BOTTOM_BORDER) {
        cell_center = CgalPoint2D(cell_center.x(), cell_center.y() + 0.5 * delta_y);
        r = distance_2D(pt_source, cell_center);
        alpha = w * boost::math::cyl_bessel_k(1, w * r) / boost::math::cyl_bessel_k(0, w * r);
        return alpha * fabs(cell_center.y() - y_source) / r;
    } else if (border == TOP_BORDER) {
        cell_center = CgalPoint2D(cell_center.x(), cell_center.y() - 0.5 * delta_y);
        r = distance_2D(pt_source, cell_center);
        alpha = w * boost::math::cyl_bessel_k(1, w * r) / boost::math::cyl_bessel_k(0, w * r);
        return alpha * fabs(cell_center.y() - y_source) / r;
    }

    cout << "WARNING in ReturnMixedBoundCondValueCSEigen: border not defined" << endl;
    return -1.0;
}

BoundaryConditionsCSEigen BoundaryConditionsDefinitionCSEigen(Domain domain,
                                                              int Nx,
                                                              int Ny,
                                                              BoundaryConditionsDef bc_map,
                                                              SourceTermsDef source_terms,
                                                              SourceTermsCSEigen& source_terms_epm,
                                                              bool fourier,
                                                              double w) {
    BoundaryConditionsCSEigen bc_def_epm;

    for (int j = 0; j < Ny; j++) {
        bc_def_epm[make_pair(0, j)].insert(make_pair(LEFT_BORDER, ReturnBoundaryConditionsCSEigen(domain, LEFT_BORDER, 0, j, Nx, Ny, bc_map, source_terms, fourier, w)));
        bc_def_epm[make_pair(Nx - 1, j)].insert(make_pair(RIGHT_BORDER, ReturnBoundaryConditionsCSEigen(domain, RIGHT_BORDER, Nx - 1, j, Nx, Ny, bc_map, source_terms, fourier, w)));
    }
    for (int i = 0; i < Nx; i++) {
        bc_def_epm[make_pair(i, 0)].insert(make_pair(BOTTOM_BORDER, ReturnBoundaryConditionsCSEigen(domain, BOTTOM_BORDER, i, 0, Nx, Ny, bc_map, source_terms, fourier, w)));
        bc_def_epm[make_pair(i, Ny - 1)].insert(make_pair(TOP_BORDER, ReturnBoundaryConditionsCSEigen(domain, TOP_BORDER, i, Ny - 1, Nx, Ny, bc_map, source_terms, fourier, w)));
    }

    for (SourceTermsDef::iterator it1 = source_terms.begin(); it1 != source_terms.end(); ++it1) {
        std::pair<int, int> indices = ::return_indices(it1->x_position, it1->y_position, domain.ReturnDeltaX(Nx), domain.ReturnDeltaY(Ny), Nx, Ny);
        source_terms_epm[indices] = it1->source_value;
    }

    return bc_def_epm;
}

}  // namespace eigen_cs