/*
 * Parameters.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

enum{INFINITE_MATRIX,FINITE_MATRIX};

#include "BoundaryConditions.h"
#include "../Utilitaries/Storage/UtilitaryStorage.h"
#include "../Utilitaries/Storage/UblasStructures.h"

#include <string>
#include <fstream>
#include <map>


/*******************************/
/* Functions for DDPParameters */
/*******************************/
class DiscretizationParameters{
public:
	int Nx;
	int Ny;
	int nb_div; //alternative definition
	int nb_seed; //alternative definition
	std::vector<double> w_vect; // vector containing wave values
	double w_current;
public:
	DiscretizationParameters(){nb_div=0;w_current=0.0;nb_seed=0;Nx=0;Ny=0;};
	~DiscretizationParameters(){};
	DiscretizationParameters(int,int);
	double ReturnDeltax(double);
	double ReturnDeltay(double);
};

/*******************************/
/* Functions for DFNParameters */
/*******************************/
class DFNParameters{
public:
	BoundaryConditionsDFN bc_map; // map<node index, boundary conditions>
	double power_law_exp;
	double percolation_param;
	double fract_aperture;
	double fract_conductivity;
	double fract_length;
	int seed;
	int nb_fract;
	std::string generation_option;
public:
	DFNParameters(){fract_length=0.0;fract_conductivity=0.0;seed=0;fract_aperture=0.0;nb_fract=0;percolation_param=0.0;power_law_exp=0.0;};
	~DFNParameters(){};
	DFNParameters(BoundaryConditionsDFN bc_map_){bc_map=bc_map_;fract_length=0.0;fract_conductivity=0.0;seed=0;fract_aperture=0.0;nb_fract=0;percolation_param=0.0;power_law_exp=0.0;};
};

/********************************/
/* Functions for SimuParameters */
/********************************/
class SimuParameters{
public:
	std::string simu_option;
	std::string model_option;
	DFNParameters DFN_param;
	DiscretizationParameters DDP_param;
	DiscretizationParameters EPM_param;
	// standard boundary conditions and source terms
	BoundaryConditionsDef bc_map;
	SourceTermsDef source_terms;
	// boundary conditions and source terms for the secondary potential
	BoundaryConditionsDef bc_map_sing;
	SourceTermsDef source_terms_sing;
	int nb_angle_step;	// number of angle step for ellipse permeability computation
	bool singularity;	// = true: singularity removal method
	bool debug_mode;
public:
	SimuParameters();
	SimuParameters(std::string,std::string,std::string file_name_bc_sing="",std::string file_name_waves="",bool singularity_=false);
	~SimuParameters(){};
};

/**********************************/
/* Functions for DomainParameters */
/**********************************/
class DomainParameters{
public:
	double Lx;	// domain size in the longitudinal direction
	double Ly;	// domain size in the transversal direction
	double mat_cond; // matrix electrical conductivity
	ublas_matrix mat_cond_dist;
public:
	DomainParameters();
	DomainParameters(std::string,std::string file_="",int Nx_=0,int Ny_=0);
	~DomainParameters(){};
};

/*********************************/
/* Functions for FilesParameters */
/*********************************/
class FilesParameters{
public:
	std::string code_path;
	std::string TypeSimu;
	std::string input_path;
	std::string output_path;
	std::string file_name_simu;
	std::string file_name_domain;
	std::string file_name_DFN;
	std::string file_name_bc;
	std::string file_name_bc_sing;
	std::string file_name_waves;
	std::string file_name_source_positions;
	std::string output_file;
	std::string output_file1;
public:
	FilesParameters();
	~FilesParameters(){};
};

/****************************/
/* Functions for Parameters */
/****************************/
class Parameters{
public:
	FilesParameters files_param;
	DomainParameters domain_param;
	SimuParameters simu_param;
public:
	Parameters(bool singularity_=false);
	~Parameters(){};
};



#endif /* PARAMETERS_H_ */
