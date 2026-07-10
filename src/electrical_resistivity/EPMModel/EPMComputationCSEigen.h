#ifndef EPMCOMPUTATIONCSEIGEN_H_
#define EPMCOMPUTATIONCSEIGEN_H_

#include "../Utilitaries/Storage/EigenCSStructures.h"
#include "../Utilitaries/Constantes.h"
#include "../Utilitaries/Storage/UtilitaryStorage.h"
#include "../DomainDefinition/Domain.h"
#include "../InputOutput/BoundaryConditions.h"
#include <boost/math/special_functions/bessel.hpp>

namespace eigen_cs {

typedef std::map<std::pair<int, int>, std::map<std::string, BoundaryConditions> > BoundaryConditionsCSEigen;
typedef std::map<std::pair<int, int>, double> SourceTermsCSEigen;

double ReturnFVMFactorCSEigen(std::string border, double delta_x, double delta_y);
void BoundCondCoeffsCSEigen(std::string border,
                            BoundaryConditions bound_cond,
                            double current_cond,
                            double delta_x,
                            double delta_y,
                            double& coeff1,
                            double& coeff3);
double ReturnMixedBoundCondValueCSEigen(SourceTermsDef source_terms,
                                        int i_cell, int j_cell,
                                        int Nx, int Ny,
                                        Domain domain,
                                        std::string border,
                                        double w);
BoundaryConditionsCSEigen BoundaryConditionsDefinitionCSEigen(Domain domain,
                                                              int Nx,
                                                              int Ny,
                                                              BoundaryConditionsDef bc_map,
                                                              SourceTermsDef source_terms,
                                                              SourceTermsCSEigen& source_terms_epm,
                                                              bool fourier = false,
                                                              double w = 0.0);

}  // namespace eigen_cs

#endif /* EPMCOMPUTATIONCSEIGEN_H_ */