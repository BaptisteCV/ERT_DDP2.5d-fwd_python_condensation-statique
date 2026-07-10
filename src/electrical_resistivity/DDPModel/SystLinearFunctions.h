/*
 * SystLinearFunctions.h
 *
 *  Created on: Apr 29, 2013
 *      Author: delphineroubinet
 */

#ifndef SYSTLINEARFUNCTIONS_H_
#define SYSTLINEARFUNCTIONS_H_

#include "../Utilitaries/Storage/UblasStructures.h"
#include "../InputOutput/BoundaryConditions.h"
#include "../DFNModel/DFNComputation.h"

void EvaluatePotential(int,int,int,ublas_matrix,ublas_vector,ublas_vector&,ublas_matrix&);
void PrintUpdateLinearSystemFourierStats();
void ResetUpdateLinearSystemFourierStats();
void UpdateLinearSystemFourier(int,int,double,double,int,double,BoundaryConditionsEPM,BoundaryConditionsDFN,SourceTermsEPM,ublas_matrix,SubNetworkMap,bool,ublas_matrix&,ublas_vector&);
void DFNDDPLinearSystemFourier(std::vector<FractureMesh>,int,int,double,double,ublas_matrix&,ublas_vector&);
void PorousDDPLinearSystemUpdateFourier(const std::vector<FractureMesh>&,int,int,double,int,double,double,double,const ublas_matrix&,ublas_matrix&,ublas_vector&);
void BCDDPLinearSystemFourier(int,double,double,double,SourceTermsEPM,SubNetworkMap,BoundaryConditionsDFN,ublas_matrix&,ublas_vector&);
void SourceTermFourier(SourceTermsEPM,int,int,bool,ublas_vector&);
void DDPModelFourierBC(bool,Parameters,Domain,Domain,BordersMap,NodesMap,SourceTermsEPM&,SourceTermsEPM&,BoundaryConditionsDFN&,BoundaryConditionsEPM&,BoundaryConditionsEPM&);
void UpdateLinearSystemBC_EPM(int,int,int,BoundaryConditionsEPM,ublas_matrix,double,double,ublas_matrix&,ublas_vector&);
double ReturnAlphaBCMixed(int,int,double,double,double,SourceTermsEPM,SubNetworkMap);
void PotentialFractureMatrixDistribution(int,int,ublas_vector,ublas_vector&,ublas_matrix&);

/*#include "../Utilitaries/Storage/UblasStructures.h"
#include "../DFNModel/DFNComputation.h"
#include "../InputOutput/BoundaryConditions.h"
#include "../EPMModel/EPMDefinition.h"


class DDPModel{
public:
	// Porous domain definition
	int Nx;	// number of cells (or blocks) discretizing the porous domain
	int Ny;
	double delta_x;	// size of the cells (or blocks) discretizing the porous domain
	double delta_y;
	ublas_matrix porous_cond;	// conductivity of the porous domain
	BoundaryConditionsEPM bc_def_epm;	// boundary conditions for the porous domain
	BoundaryConditionsEPM bc_def_epm_sing;
	SourceTermsEPM source_terms_epm;
	SourceTermsEPM source_terms_epm_sing;
	// DFN domain definition
	SubNetworkMap subnetwork_map;	// division of the fracture network in subnetwork at the block-scale
	BoundaryConditionsDFN bc_map;	// boundary conditions definition for the fracture network
	int nb_nodes_fract;
	// Definition for the fracture-matrix interaction
	ublas_matrix AlphaBlock;
	// Linear system definition
	ublas_matrix matrix_syst;	// matrix of the linear system
	ublas_vector vector_syst;	//  second member of the linear system
	// Fourier variables
	double current_w;
	bool singularity;
public:
	DDPModel(){};
	DDPModel(Parameters,Domain,NetworkMeshes,bool fourier=false,bool singularity=false);
	virtual ~DDPModel(){porous_cond.clear();AlphaBlock.clear();matrix_syst.clear();vector_syst.clear();};

	// New functions to optimize the definition of the linear system (over several experiments and Fourier variables)
	void DetermineBasicLinearSystemDDP();
	void UpdateLinearSystemFourier();
	void PorousDDPLinearSystemUpdateFourier(std::vector<FractureMesh>,int,double);
	void UpdateLinearSystemBC_EPM(int);
	void EvaluatePotential(ublas_vector&,ublas_matrix&);

	void ComputePotential(ublas_vector&,ublas_matrix&);
	void DFNDDPLinearSystem(std::vector<FractureMesh>,int,double);
	void DefinitionLinearSystem(int,std::pair<int,int>,std::string);
	void DefinitionLinearSystemBCOnBorder(int,std::pair<int,int>,std::string);
	void PorousDDPLinearSystem(std::vector<FractureMesh>,int,double);
	void SingularityMethod();
	void BCDDPLinearSystem();
	void BCDDPLinearSystemFourier();
	void PotentialFractureMatrixDistribution(ublas_vector,ublas_vector&,ublas_matrix&);
	void PotentialFractureMatrixDistributionFourier(ublas_vector,ublas_vector&,ublas_matrix&);
	double phi_function(double,int,int);
	void DDPModelFourierBC(Parameters,Domain,Domain,BordersMap,NodesMap);
	// Functions for Fourier model
	void ComputePotentialFourier(ublas_vector&,ublas_matrix&,bool bc_on_border=false);
	void SourceTermFourier();
	void DefineLinearSystemFourier(ublas_matrix&,bool);
	void DFNDDPLinearSystemFourier(std::vector<FractureMesh>,int,double);
	void PorousDDPLinearSystemFourier(std::vector<FractureMesh>,int,double);
	void DetermineLinearSystemFourier(int);
	double ReturnAlphaBCMixed(int,int);
};

void print_subnetwork(SubNetworkMap);
*/

#endif /* SYSTLINEARFUNCTIONS_H_ */
