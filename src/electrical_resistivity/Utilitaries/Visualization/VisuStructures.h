

/* * VisuStructures.h
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet*/

/*
#ifndef VISUSTRUCTURES_H_
#define VISUSTRUCTURES_H_

#include "../../DFNModel/NetworkMeshes.h"
#include "../../DomainDefinition/Domain.h"
#include "../../Utilitaries/Math/LinearSystem.h"
#include <boost/numeric/ublas/vector.hpp>

class DFNVisu{
public:
	double L;
	Domain domain;
	NetworkMeshes net_mesh;
	ublas_vector DFN_potential;
	double cond_matrix; // matrix conductivity
	double cond_fract;
	double fract_spacing;
	std::string option;
public:
	DFNVisu();
	~DFNVisu(){};
	DFNVisu(Domain,NetworkMeshes,ublas_vector,double cond_fract=1,double cond_matrix=0,double fract_spacing=0);
	DFNVisu(double,NetworkMeshes,std::string);
	DFNVisu(double,NetworkMeshes,ublas_vector,std::string);
	void display_DFN();
	void display_Potential_Classic();
	//void display_Potential_Radial();
};

class EPMVisu{
public:
	double Lx;
	double Ly;
	ublas_matrix EPM_potential;
	ublas_matrix Current_x;
	ublas_matrix Current_y;
	double min_visu;
	double max_visu;
public:
	EPMVisu();
	~EPMVisu(){};
	EPMVisu(double,double,ublas_matrix,double min_visu=-1,double max_visu=-1);
	EPMVisu(double,double,ublas_matrix,ublas_matrix,ublas_matrix);
	void display_EPM_Potential();
	void display_EPM_Current_x();
	void display_EPM_Current_y();
};

class EPMDFNVisu{
public:
	Parameters param;
	double L;
	double cond_matrix;
	NetworkMeshes net_mesh;
	SubNetworkMap subnetwork_map;
	ublas_vector DFN_potential;
	ublas_matrix EPM_potential;
	ublas_matrix AlphaBlock;
	double min_visu;
	double max_visu;
public:
	EPMDFNVisu();
	~EPMDFNVisu(){};
	EPMDFNVisu(double,NetworkMeshes,ublas_vector,ublas_matrix);
	EPMDFNVisu(double,SubNetworkMap,ublas_vector,ublas_matrix,double,ublas_matrix,double min_visu_=-1,double max_visu_=-1);
	EPMDFNVisu(Parameters,SubNetworkMap,ublas_vector,ublas_matrix);
	void display_EPMDFN_Potential();
	void display_EPMDFN_Potential_Matrix();
	void display_EPMDFN_Current();
};


#endif // VISUSTRUCTURES_H_
*/
