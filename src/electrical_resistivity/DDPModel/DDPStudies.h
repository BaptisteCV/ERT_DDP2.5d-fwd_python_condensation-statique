/*
 * DDPStudies.h
 *
 *  Created on: Apr 29, 2013
 *      Author: delphineroubinet
 */

#ifndef DDPSTUDIES_H_
#define DDPSTUDIES_H_

#include "../Utilitaries/Visualization/VisuStructures.h"
#include "../EPMModel/EPMDefinition.h"


void PotentialComputationDDP(Parameters,Domain,double,double);
void ConductivityPlotDDP(Parameters,Domain);
//void ExchangeCoeffStudy(Parameters,Domain,EPMVisu&,double&);
double ComputeFMFlowMesh(EPMSystem,ublas_matrix,int,int,int,int,double);
void ExchangeCoeffBlockScale(EPMSystem,ublas_matrix,double,double&,double&);
void ExchangeStudy(Parameters,Domain,int,int);
//void ExchangeStudy2(Parameters,Domain);
void ComputationSumFMExchange(std::map<int,Domain>,SubNetworkMap,std::map<int,EPMSystem>,double,std::map<int,double>&,std::map<int,double>&);
void ConductivityPlotEPMDFN(Parameters,Domain);
void PotentialStudy(Parameters,Domain,int,int);

#endif /* DDPSTUDIES_H_ */
