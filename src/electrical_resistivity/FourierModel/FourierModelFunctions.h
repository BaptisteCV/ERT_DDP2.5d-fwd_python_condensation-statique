/*
 * FourierModelFunctions.h
 *
 *  Created on: Apr 25, 2013
 *      Author: delphineroubinet
 */

#ifndef FOURIERMODELFUNCTIONS_H_
#define FOURIERMODELFUNCTIONS_H_

#include "../InputOutput/Parameters.h"
#include "../DomainDefinition/Domain.h"
#include "../Utilitaries/Storage/UblasStructures.h"
#include "../Utilitaries/Storage/EigenCSStructures.h"

void FourierPotentialEPM(Parameters,Domain,bool);
//void FourierPotentialEPMSourceTerm(Parameters,Domain,bool,bool,std::map<int,std::pair<double,double> >);
std::map<int, ublas_matrix > FourierPotentialEPMSourceTerm(Parameters,Domain,bool,bool,std::map<int,std::pair<double,double> >,bool build_result_map=false);
void FourierPotentialEPMSourceTerm(Parameters,Domain,bool,bool);
std::map<int, ublas_matrix > FourierPotentialDDPSourceTerm(Parameters,Domain,bool,std::map<int,std::pair<double,double> >,bool build_result_map=false);
void FourierPotentialDDPSourceTerm(Parameters,Domain,bool);
void FourierPotentialDDP(Parameters,Domain,bool,bool bc_on_border=false);

std::map<int, eigen_cs::SparseMatrix > FourierPotentialDDPSourceTermEigen(Parameters,Domain,bool,std::map<int,std::pair<double,double> >,bool build_result_map=false);

std::vector<double> ReadSourceTermsCoordX(Parameters);
std::map<int,std::pair<double,double> > ReadSourceTermsCoordinates(Parameters);


#endif /* FOURIERMODELFUNCTIONS_H_ */
