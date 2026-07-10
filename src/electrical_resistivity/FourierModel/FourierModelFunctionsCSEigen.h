#ifndef FOURIERMODELFUNCTIONSCSEIGEN_H_
#define FOURIERMODELFUNCTIONSCSEIGEN_H_

#include "../InputOutput/Parameters.h"
#include "../DomainDefinition/Domain.h"
#include "../InputOutput/BoundaryConditions.h"
#include "../Utilitaries/Storage/EigenCSStructures.h"
#include <map>
#include <vector>

namespace eigen_cs {

std::map<int, SparseMatrix> FourierPotentialEPMSourceTermCSEigen(
    Parameters param,
    Domain domain,
    bool singularity,
    bool bc_on_border,
    std::map<int, std::pair<double, double> > source_term_positions,
    bool build_result_map = false);

std::map<int, SparseMatrix> FourierPotentialDDPSourceTermCSEigen(
    Parameters param,
    Domain domain,
    bool singularity,
    bool bc_on_border,
    std::map<int, std::pair<double, double> > source_term_positions,
    bool build_result_map = false);

}  // namespace eigen_cs

#endif /* FOURIERMODELFUNCTIONSCSEIGEN_H_ */