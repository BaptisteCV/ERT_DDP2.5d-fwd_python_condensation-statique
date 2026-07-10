#ifndef FOURIERMODELFUNCTIONSEIGEN_H_
#define FOURIERMODELFUNCTIONSEIGEN_H_

#include "../InputOutput/Parameters.h"
#include "../DomainDefinition/Domain.h"
#include "../Utilitaries/Storage/EigenCSStructures.h"
#include <map>

namespace eigen_epm {

void FourierPotentialEPMEigen(Parameters param, Domain domain, bool singularity);
std::map<int, eigen_cs::SparseMatrix> FourierPotentialEPMSourceTermEigen(
    Parameters param,
    Domain domain,
    bool singularity,
    bool bc_on_border,
    std::map<int, std::pair<double, double> > source_term_positions);
void FourierPotentialEPMSourceTermEigen(Parameters param, Domain domain, bool singularity, bool bc_on_border);

}  // namespace eigen_epm

#endif /* FOURIERMODELFUNCTIONSEIGEN_H_ */
