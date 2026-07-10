/*
 * Parameters.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#include "Parameters.h"
#include "../Utilitaries/Constantes.h"
#include <iostream>
//#include <filesystem>

using namespace std;


/*******************************/
/* Functions for DDPParameters */
/*******************************/
DiscretizationParameters::DiscretizationParameters(int Nx_,int Ny_){Nx=Nx_;Ny=Ny_;nb_seed=0;nb_div=0;w_current=0.0;}	// for steady-state simulations
double DiscretizationParameters::ReturnDeltax(double Lx){
	return Lx/Nx;
}
double DiscretizationParameters::ReturnDeltay(double Ly){
	return Ly/Ny;
}

/********************************/
/* Functions for SimuParameters */
/********************************/
SimuParameters::SimuParameters(){singularity=false;nb_angle_step=0;}
SimuParameters::SimuParameters(string file_name,string file_name_bc,string file_name_bc_sing,string file_name_waves,bool singularity_){
	nb_angle_step=0;
	singularity=singularity_;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		fichier >> this->simu_option;
		if (simu_option==FOURIER_POTENTIAL_COMPUTATION||simu_option==POTENTIAL_COMPUTATION||simu_option==CONDUCTIVITY_PLOT){
			fichier >> model_option;
			if (model_option==DDP||model_option==DDP_CS||model_option==DDP_EIGEN){
				int Nx,Ny;
				fichier >> Nx;
				fichier >> Ny;
				fichier >> debug_mode;
				DDP_param=DiscretizationParameters(Nx,Ny);
				EPM_param=DiscretizationParameters(Nx,Ny);
			}
			else if (model_option==EPM||model_option==EPM_Eigen||model_option==EPM_CS||model_option==EPMDFN){
				int Nx,Ny;
				fichier >> Nx;
				fichier >> Ny;
				fichier >> debug_mode;
				EPM_param=DiscretizationParameters(Nx,Ny);
			}
			SourceTermsDef source_terms_def_;
			bc_map=ReadBoundaryConditions(file_name_bc,source_terms_def_);
			source_terms=source_terms_def_;
			if (singularity){
				SourceTermsDef source_terms_def_sing_;
				bc_map_sing=ReadBoundaryConditions(file_name_bc_sing,source_terms_def_sing_);
				source_terms_sing=source_terms_def_sing_;
			}
			if (simu_option==FOURIER_POTENTIAL_COMPUTATION){
				ifstream fichier_waves(file_name_waves.c_str());
				if (fichier_waves.is_open()){
					int Nw;double value;vector<double> vect_waves;
					fichier_waves >> Nw;
					for (int i=0;i<Nw;i++){
						fichier_waves >> value;
						vect_waves.push_back(value);
					}
					fichier_waves.close();
					if (model_option==EPM||model_option==EPM_Eigen||model_option==EPM_CS){
						EPM_param.w_vect=vect_waves;
					}
					else if (model_option==DDP||model_option==DDP_CS||model_option==DDP_EIGEN){
						DDP_param.w_vect=vect_waves;
					}
					else{
						cout << "WARNING in SimuParameters (Parameters.cpp): simu_option not defined" << endl;
					}
				}
			}
		}
		else if (simu_option==POTENTIAL_STUDY){
			fichier >> model_option;
			if (model_option==EPMDFN){
				int Nx,Ny,nb_div,nb_seed;
				fichier >> Nx;
				fichier >> Ny;
				fichier >> nb_div;
				fichier >> nb_seed;
				EPM_param=DiscretizationParameters(Nx,Ny);
				DDP_param=DiscretizationParameters(Nx,Ny);
				DDP_param.nb_div=nb_div;
				DDP_param.nb_seed=nb_seed;
			}
			else{cout << "WARNING in SimuParameters (Parameters.cpp): simu_option not defined" << endl;}
			SourceTermsDef source_terms_def_;
			bc_map=ReadBoundaryConditions(file_name_bc,source_terms_def_);
			source_terms=source_terms_def_;
		}
		else if (simu_option==PERMEABILITY_ELLIPSE){
			fichier >> model_option;
			fichier >> nb_angle_step;
			if (model_option==DDP){
				int Nx,Ny;
				fichier >> Nx;
				fichier >> Ny;
				DDP_param=DiscretizationParameters(Nx,Ny);
			}
			else if (model_option==EPM){
				int Nx,Ny;
				fichier >> Nx;
				fichier >> Ny;
				EPM_param=DiscretizationParameters(Nx,Ny);
			}
		}
		else if (simu_option==EXCHANGE_STUDY){
			int Nx,Ny;
			fichier >> model_option;
			if (model_option==EPMDFN){
				fichier >> Nx;
				fichier >> Ny;
				EPM_param=DiscretizationParameters(Nx,Ny);
			}
			else if (model_option==DDP){
				int nb_div,nb_seed;
				fichier >> Nx;
				fichier >> Ny;
				fichier >> nb_seed;
				DDP_param=DiscretizationParameters(Nx,Ny);
				DDP_param.nb_seed=nb_seed;
			}
			else{cout << "WARNING in SimuParameters (Parameters.cpp): simu_option not defined" << endl;}
			SourceTermsDef source_terms_def_;
			bc_map=ReadBoundaryConditions(file_name_bc,source_terms_def_);
			source_terms=source_terms_def_;
		}
		else if (simu_option!=DFN_GENERATION){
			cout << "WARNING in SimuParameters (Parameters.cpp): simu_option not defined" << endl;
		}
	}
}

/**********************************/
/* Functions for DomainParameters */
/**********************************/
DomainParameters::DomainParameters(){Lx=0;Ly=0;mat_cond=0;}
DomainParameters::DomainParameters(string file_name1,string file_name2,int Nx,int Ny){
	Lx=0;Ly=0;mat_cond=0;
	ifstream fichier(file_name1.c_str());
	if (fichier.is_open()){
		fichier >> this->Lx;
		fichier >> this->Ly;
		fichier >> this->mat_cond;
		fichier.close();
	}
	else{cout << "WARNING in DomainParameters (Parameters.cpp): file1 could not be opened" << endl;}
	// heterogeneous matrix conductivity
	if (mat_cond==-1){
		ifstream fichier1(file_name2.c_str());
		if (fichier1.is_open()){
			ublas_matrix current_mat(Nx,Ny);
			double value;
			for (int j=Ny-1;j>=0;j--){
				for (int i=0;i<Nx;i++){
					fichier1 >> value;
					current_mat(i,j)=value;
				}
			}
			fichier1.close();
			mat_cond_dist=current_mat;
		}
		else{cout << "WARNING in DomainParameters (Parameters.cpp): file2 could not be opened" << endl;}
	}
}

/*********************************/
/* Functions for FilesParameters */
/*********************************/
FilesParameters::FilesParameters(){
	/*cout << "Enter the path for Input and Output folders" << endl;
	cin >> code_path;*/
	// 1. Read the type of simulations and define input path
	//code_path="./../.."; // uncomment if executing ElectricalResistivty from Release folder
	code_path=".";
	string simu_file_name=code_path+"/input/SimuType.txt";
	//cout << simu_file_name << endl;
	ifstream fichier1(simu_file_name.c_str());
	string simu_name;
	if (fichier1.is_open()){fichier1 >> simu_name;fichier1.close();}
	else{cout << "Parameters files not found ERROR1" << endl;}
	TypeSimu = simu_name;
	input_path=code_path+"/input/"+simu_name+"/";
	output_path=code_path+"/output/"+simu_name+"/";

	string file_name=input_path+"FileNames.txt";
	ifstream fichier2(file_name.c_str());
	if (fichier2.is_open()){
		fichier2 >> file_name_simu;
		fichier2 >> file_name_domain;
		fichier2 >> file_name_DFN;
		fichier2 >> file_name_bc;
		fichier2 >> file_name_bc_sing;
		fichier2 >> file_name_waves;
		fichier2 >> file_name_source_positions;
		fichier2 >> output_file;
		fichier2 >> output_file1;
		fichier2.close();
	}
	else{cout << "Parameters files not found" << endl;}
}

/****************************/
/* Functions for Parameters */
/****************************/
Parameters::Parameters(bool singularity_){
	// 1. Simulation parameters
	string file_name1=files_param.input_path+files_param.file_name_simu;
	string file_name_bc=files_param.input_path+files_param.file_name_bc;
	string file_name_bc_sing=files_param.input_path+files_param.file_name_bc_sing;
	string file_name_waves=files_param.input_path+files_param.file_name_waves;
	simu_param=SimuParameters(file_name1,file_name_bc,file_name_bc_sing,file_name_waves,singularity_);
	// 2. Domain description
	file_name1=files_param.input_path+files_param.file_name_domain;
	string file_name2=files_param.input_path+"porous_cond.txt";
	domain_param=DomainParameters(file_name1,file_name2,simu_param.EPM_param.Nx,simu_param.EPM_param.Ny);

}
