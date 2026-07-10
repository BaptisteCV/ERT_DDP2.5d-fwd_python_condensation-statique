/*
 * DFNComputation.h
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet
 */

#ifndef DFNCOMPUTATION_H_
#define DFNCOMPUTATION_H_

#include "../DFNModel/NetworkMeshes.h"
#include "../InputOutput/Parameters.h"
#include "../Utilitaries/Math/LinearSystem.h"
#include "../Utilitaries/Visualization/VisuStructures.h"
#include "../EPMModel/EPMComputation.h"

typedef std::map<int,double> SourceTermDFN;	//<node index,source term>


//DFNVisu PotentialComputationDFN(Parameters,Domain);
BoundaryConditionsDFN ReturnBoundCondDFN(Domain,BoundaryConditionsDef,BordersMap,NodesMap);
BoundaryConditionsDFN ReturnBoundCondDFNFourier(Domain,BoundaryConditionsDef,BordersMap,NodesMap,BoundaryConditionsEPM,double,double);
ublas_vector DFNComputationClassic(NetworkMeshes,BoundaryConditionsDFN);
ublas_vector DFNComputationClassic(NetworkMeshes,BoundaryConditionsDFN,SourceTermDFN);
ublas_vector DFNComputationClassic(NetworkMeshes,BoundaryConditionsDFN,SourceTermDFN,ublas_vector);
ublas_vector DFNComputationRadialDiffusion(NetworkMeshes,BoundaryConditionsDFN,double,double);
double AnalyticalSolutionFractureClassic(double,double,double,double,double);
double AnalyticalSolutionFractureMatrix(double,double,double,double,double,double);
double AnalyticalSolutionFractureRadialDiffusion(double,double,double,double,double,double,double,double);
double MetricElecFlow(double,double,double,double);
double TotalElectricCurrentDFN(Parameters,NetworkMeshes,ublas_vector,BoundaryConditionsDFN,double&);

ublas_vector DFNMatrixComputation(SubNetworkMap,BoundaryConditionsDFN,ublas_matrix,double,double);
void ReturnDerivativeCoefficients(double,double,double,double,double,double,double,double&,double&,double&);

ublas_vector ComputeFlowVelocities(NetworkMeshes&,Parameters,Domain,std::string option="gradient");
#endif /* DFNCOMPUTATION_H_ */
