/*
 * EPMDefinition.h
 *
 *  Created on: Apr 25, 2013
 *      Author: delphineroubinet
 */

#ifndef EPMDEFINITION_H_
#define EPMDEFINITION_H_

#include "../Utilitaries/Storage/UblasStructures.h"
#include "EPMComputation.h"


// class to evaluate potential in Equivalent Porous Medium
class EPMSystem{
public:
	// matrix and vector of the linear system
	ublas_matrix matrix_syst;
	ublas_vector vector_syst;
	// spatial discretization
	int Nx;
	int Ny;
	double delta_x;
	double delta_y;
	double w;//wave number for Fourier transform
	// domain properties
	ublas_matrix porous_cond; // conductivity of the porous domain
	// boundary conditions of the domain
	BoundaryConditionsEPM bc_def_epm;
	BoundaryConditionsEPM bc_def_epm_sing;
	SourceTermsEPM source_terms_epm;
	SourceTermsEPM source_terms_epm_sing;
	std::vector<std::pair<int,int> > fract_indices;	// indices of the meshes corresponding to a projected fracture
	ublas_matrix potential;
public:
	EPMSystem(){};
	EPMSystem(Parameters);
	EPMSystem(Parameters,NetworkMeshes);
	EPMSystem(Parameters,Domain,NetworkMeshes,std::string option=INTERSECTION);
	EPMSystem(Parameters,Domain,bool fourier_syst=false,bool fourier_bc=false,bool singularity=false,bool bc_on_border=false,bool fourier_update=false, bool source_term=true);
	virtual ~EPMSystem(){};
	// New functions to progressively define the linear system in the Fourier domain (implemented on 2018/12/05)
	void DetermineBasicLinearSystem();
	void UpdateLinearSystemFourier();
	void UpdateLinearSystemBoundaryConditionsAndSourceTerms(bool);


	void DetermineLinearSystem(bool,bool,bool,bool fourier_update=false,bool source_term=true);
	void SingularityMethod();
	void SingularityMethodBCOnBorder();
	void DefinitionLinearSystem(int,int,std::string,bool bc_on_border=false);
	void DefinitionLinearSystemStandard(int,int,std::string);
	void DefinitionLinearSystemBCOnBorder(int,int,std::string);
	ublas_matrix PotentialComputation();
	void print_bc_def();
	void SourceTermsImplementation(bool fourier=false);
	double phi_function(double,int,int);
	std::map<int,EPMSystem> DefineSubEPMFract(std::map<int,Domain>,SubNetworkMap,double,std::string option=INTERSECTION);
	std::map<int,EPMSystem> DefineSubEPM(std::map<int,Domain>,double);
	void EPMConductivityDefinition(double,pointcpp<double> min_domain=pointcpp<double>(0,0),NetworkMeshes net_mesh=NetworkMeshes(),std::string option=INTERSECTION);
};

std::map<int,ublas_matrix> ExtractMatrixMap(ublas_matrix,std::map<int,Domain>,std::map<int,EPMSystem>);
//EPMVisu PotentialComputationEPM(Parameters,Domain,double,double);
//EPMVisu PotentialComputationEPMDFN(Parameters,Domain,NetworkMeshes,double,double);
std::map<int,double> AveragedMatrices(std::map<int,EPMSystem>);
ublas_matrix ReturnSigma0(ublas_matrix,double&);

#endif /* EPMDEFINITION_H_ */
