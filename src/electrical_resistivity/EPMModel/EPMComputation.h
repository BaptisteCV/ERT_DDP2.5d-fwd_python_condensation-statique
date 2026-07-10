/*
 * EPMDefinition.h
 *
 *  Created on: Apr 25, 2013
 *      Author: delphineroubinet
 */

#ifndef EPMCOMPUTATION_H_
#define EPMCOMPUTATION_H_

#include "../Utilitaries/Math/LinearSystem.h"
#include "../InputOutput/Parameters.h"
#include "../DFNModel/NetworkMeshes.h"
#include "../Utilitaries/Visualization/VisuStructures.h"

// boundary conditions storage for EPM domain
typedef std::map<std::pair<double,double>,std::map<std::string,BoundaryConditions> > BoundaryConditionsEPM;	// map<cell coordinate,map<border,boundary conditions> >
typedef std::map<std::pair<int,int>,double> SourceTermsEPM;	// map<cell indices,source term value>

void CoefficientsLinearSystem(int,int,int,int,std::string,BoundaryConditionsEPM,ublas_matrix,double,double,double&,double&,double&,int&);
void CoefficientsLinearSystemBCOnBorder(int,int,int,int,std::string,BoundaryConditionsEPM,ublas_matrix,double,double,double&,double&,int&);
double ReturnFVMFactor(std::string,double,double);
double ReturnFVMRegularCoeff(std::string,double,double,double,double);
void BoundCondCoeffsBCOnBorder(std::string,BoundaryConditions,double,double,double&,double&,double&);
void BoundCondCoeffs(std::string,BoundaryConditions,double,double,double,double&,double&);
BoundaryConditionsEPM BoundaryConditionsDefinition(Domain,int,int,BoundaryConditionsDef,SourceTermsDef,SourceTermsEPM&,bool fourier=false,double w=0);
BoundaryConditionsEPM BoundaryConditionsDefinition(Domain,int,int,BoundaryConditionsDef,SourceTermsDef);
BoundaryConditions ReturnBoundaryConditions(Domain,std::string,int,int,int,int,BoundaryConditionsDef,SourceTermsDef,bool fourier=false,double w=0);
double ReturnMixedBoundCondValue1(SourceTermsDef,int,int,int,int,Domain,std::string);
double ReturnMixedBoundCondValue2(SourceTermsDef,int,int,int,int,Domain,std::string,double);

void EPMMatrixVectSystem(Parameters,ublas_matrix&,ublas_vector&);
void EPMMatrixVectSystem(Parameters,ublas_matrix&,ublas_vector&);
ublas_matrix EPMDefinition(Parameters,ublas_matrix,ublas_vector,ublas_matrix,ublas_matrix);
ublas_matrix EPMDefinition(Parameters,ublas_matrix,ublas_vector,ublas_matrix);
ublas_matrix EPMDefinition(Parameters,ublas_matrix,ublas_vector);
ublas_matrix EPMDefinition(Parameters,ublas_matrix);
ublas_matrix EPMDefinition(Parameters);
void ElectricCurrentComputation(Parameters,ublas_matrix,ublas_matrix,ublas_matrix&,ublas_matrix&);
double TotalElectricCurrentEPM(Parameters,ublas_matrix,ublas_matrix,std::string);
double EquivalentConductivityEPM(Parameters,ublas_matrix,ublas_matrix);

void print(BoundaryConditionsEPM bc_epm);



#endif /* EPMCOMPUTATION_H_ */
