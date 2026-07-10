/*
 * Results.h
 *
 *  Created on: Jun 29, 2012
 *      Author: delphine
 */

#ifndef RESULTS_H_
#define RESULTS_H_

#include "../InputOutput/Parameters.h"
#include "../DFNModel/NetworkMeshes.h"
#include "../Utilitaries/Storage/UblasStructures.h"
#include "../Utilitaries/Storage/EigenCSStructures.h"
#include <map>
#include <string>


class Results{
public:
	// Output file
	std::string output_file;	// file name to write results
	std::string output_file1;	// file name to write results
	// DFN Generation results
	NetworkMeshes net_meshes;
	// Polar permeability
	std::map<double,double> map_perm;
	// Flow with discretization
	std::map<int,double> map_res1;
	std::map<int,double> map_res2;
	std::map<int,double> map_res3;
	std::map<double,ublas_matrix> map_mat_res;
	std::map<double,ublas_vector> map_fract_res;
	ublas_matrix res1;
public:
	Results(){};
	virtual ~Results(){};
	////////////////////////////////
	// Results for DFN Generation //
	////////////////////////////////
	Results(std::string,NetworkMeshes);
	Results(std::string,std::map<double,double>);
	Results(std::string,std::map<int,double>);
	Results(std::string,std::map<double,ublas_matrix>);
	Results(std::string,std::map<double,eigen_cs::DenseMatrix>);
	Results(std::string,std::map<double,ublas_vector>);
	Results(std::string,std::map<double,ublas_vector>,int,int);
	Results(std::string,ublas_matrix);
	Results(std::string,std::string,std::map<int,double>,std::map<int,double>,std::map<int,double>);
	void WriteDFNDefinition();
};


#endif /* RESULTS_H_ */
