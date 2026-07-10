/*
 * FourierModelFunctions.cpp
 *
 *  Created on: Apr 25, 2013
 *      Author: delphineroubinet
 */

#include "FourierModelFunctions.h"
#include "../EPMModel/EPMDefinition.h"
#include "../InputOutput/Results.h"
#include "../DDPModel/DDPModel.h"
#include "../DFNModel/DFNGeneration.h"
#include <omp.h>
using namespace std;

void FourierPotentialEPMSourceTerm(Parameters param,Domain domain,bool singularity,bool bc_on_border,map<int,pair<double,double> > source_term_positions){
	// Parameters
	vector<double> w_vect=param.simu_param.EPM_param.w_vect;
	bool fourier=true;
	int elec_nb;double source_position_x,source_position_y,source_value=1.0;
	// Initial domain and system
	EPMSystem epd_model_init(param),epd_model;
	epd_model_init.DetermineBasicLinearSystem();
	// Loop on the different electrode positions <-> different experiments
	for (map<int,pair<double,double> >::iterator it1=source_term_positions.begin();it1!=source_term_positions.end();it1++){
		// Definition of the considered electrode position
		elec_nb=it1->first;
		cout << "Experiment " << elec_nb << endl;
		source_position_x=it1->second.first,source_position_y=it1->second.second;
		if (param.simu_param.source_terms.size()==0){
			param.simu_param.source_terms.push_back(SourceTerm(source_position_x,source_position_y,source_value));
		}
		else{param.simu_param.source_terms[0]=SourceTerm(source_position_x,source_position_y,source_value);};
		// Each experiment is run for several wave values
		map<double,ublas_matrix> Potential_Matrices;
		for (vector<double>::iterator it1=w_vect.begin();it1!=w_vect.end();it1++){
			//param.simu_param.EPM_param.w_current=*it1;
			// update the EPM model and linear system
			epd_model=epd_model_init;
			//epd_model.w=param.simu_param.EPM_param.w_current;
			// Fourier terms
			epd_model.UpdateLinearSystemFourier(*it1);
			// Define the boundary conditions and source
			SourceTermsEPM source_terms_epm_;
			epd_model.bc_def_epm=BoundaryConditionsDefinition(domain,epd_model.Nx,epd_model.Ny,param.simu_param.bc_map,param.simu_param.source_terms,source_terms_epm_,fourier,*it1);
			epd_model.source_terms_epm=source_terms_epm_;
			// Update the linear system with the boundary conditions and source terms
			epd_model.UpdateLinearSystemBoundaryConditionsAndSourceTerms(fourier);
			Potential_Matrices[*it1]=epd_model.PotentialComputation();
		}
		// Write the results
		stringstream ss;string str_cpt;ss<<elec_nb;str_cpt=ss.str();
		string output_file="Results"+str_cpt;
		output_file=output_file+".txt";
		Results results(param.files_param.output_path+output_file,Potential_Matrices);
	}
}

void FourierPotentialEPMSourceTerm(Parameters param,Domain domain,bool singularity,bool bc_on_border){

	vector<double> w_vect=param.simu_param.EPM_param.w_vect;
	if (param.simu_param.source_terms.size()!=1){
		cout << "WARNING in main (PERFORM.cpp): configuration not implemented" << endl;
	}

	double source_position_x=param.simu_param.source_terms[0].x_position;
	double source_position_y=param.simu_param.source_terms[0].y_position;
	double source_value=param.simu_param.source_terms[0].source_value;
	double Lsimu=domain.domain_size_x();
	double source_spacing=param.simu_param.source_terms[0].source_spacing;
	int cpt=0;bool continue_loop=true,fourier=true;

	EPMSystem epd_model_init(param),epd_model;
	epd_model_init.DetermineBasicLinearSystem();

	while (continue_loop){
		cpt++;
		cout << "Experiment " << cpt << endl;

		param.simu_param.source_terms[0]=SourceTerm(source_position_x,source_position_y,source_value);

		map<double,ublas_matrix> Potential_Matrices;
		for (vector<double>::iterator it1=w_vect.begin();it1!=w_vect.end();it1++){
			//cout << "w=" << *it1 << endl;
			//param.simu_param.EPM_param.w_current=*it1;
			// update the EPM model and linear system
			epd_model=epd_model_init;
			//epd_model.w=param.simu_param.EPM_param.w_current;
			// Fourier terms
			epd_model.UpdateLinearSystemFourier(*it1);
			// Define the boundary conditions and source
			SourceTermsEPM source_terms_epm_;
			epd_model.bc_def_epm=BoundaryConditionsDefinition(domain,epd_model.Nx,epd_model.Ny,param.simu_param.bc_map,param.simu_param.source_terms,source_terms_epm_,fourier,*it1);
			epd_model.source_terms_epm=source_terms_epm_;
			// Update the linear system with the boundary conditions and source terms
			epd_model.UpdateLinearSystemBoundaryConditionsAndSourceTerms(fourier);

			Potential_Matrices[*it1]=epd_model.PotentialComputation();
		}

		stringstream ss;string str_cpt;
		ss<<cpt;
		str_cpt=ss.str();

		string output_file="EPM_Xu"+str_cpt;
		output_file=output_file+".txt";

		Results results(param.files_param.output_path+output_file,Potential_Matrices);

		source_position_x+=source_spacing;
		if (source_position_x<Lsimu&&source_spacing!=-1){continue_loop=true;}
		else{continue_loop=false;}

	}
}

void FourierPotentialEPM(Parameters param,Domain domain,bool singularity){
	vector<double> w_vect=param.simu_param.EPM_param.w_vect;
	map<double,ublas_matrix> Potential_Matrices;
	for (vector<double>::iterator it1=w_vect.begin();it1!=w_vect.end();it1++){
		//param.simu_param.EPM_param.w_current=*it1;
		EPMSystem epd_model(param,domain,*it1,true,singularity);
		Potential_Matrices[*it1]=epd_model.PotentialComputation();
	}

	Results results(param.files_param.output_path+param.files_param.output_file,Potential_Matrices);
}

void FourierPotentialDDP(Parameters param,Domain domain,bool singularity,bool bc_on_border){
	// 1.1. Variables
	vector<double> w_vect=param.simu_param.DDP_param.w_vect;
	map<double,ublas_matrix> Potential_Matrices;map<double,ublas_vector> Potential_Fractures;
	ublas_vector DFN_potential;ublas_matrix matrix_potential;
	// 1.2. Fracture network generation
	NetworkMeshes net_mesh(domain,param,param.simu_param.simu_option);
	bool fourier=true;
	clock_t t1=clock();
	DDPModel ddp_model(param,domain,net_mesh,fourier,singularity);
	clock_t t2=clock();
	cout << "time for generating the fracture network = t2-t1" << endl;
	// 2. Simulation for each wave number
	for (vector<double>::iterator it1=w_vect.begin();it1!=w_vect.end();it1++){
		// 2.1. Update the current wave number
		//param.simu_param.DDP_param.w_current=*it1;
		// 2.2. Define the boundary conditions for ddp model (depend on w)
		ddp_model.DDPModelFourierBC(param,domain,net_mesh.domain,param.simu_param.source_terms,*it1,net_mesh.border_map,net_mesh.nodes_map);
		// 2.3. Compute electrical potential
		//cout << "w=" << ddp_model.current_w << endl;
		ddp_model.ComputePotentialFourier(DFN_potential,matrix_potential,bc_on_border);
		// 2.4. Store the computed values
		Potential_Matrices[*it1]=matrix_potential;
		Potential_Fractures[*it1]=DFN_potential;
	}
	Results results(param.files_param.output_path+param.files_param.output_file,Potential_Matrices);
	Results results1(param.files_param.output_path+param.files_param.output_file1,Potential_Fractures);
}

void FourierPotentialDDPSourceTerm(Parameters param,Domain domain,bool singularity,map<int,pair<double,double> > source_term_positions){
	// 1. Global variables and fracture network generation
	// 1.1. Variables
	//vector<double> w_vect=param.simu_param.DDP_param.w_vect;
	// 1.2 Fracture network generation
	clock_t t1=clock();
	NetworkMeshes net_mesh(domain,param,param.simu_param.simu_option);
	clock_t t2=clock();
        cout << "Computation time for generating the fracture network = " << (t2-t1)/CLOCKS_PER_SEC << " seconds" << endl;
	bool fourier=true;
	// 1.3. Variables for the source term positions
	double source_value=1.0,source_position_x,source_position_y;
	int nb_elec;
	// 2. Simulations
	// 2.1. Basic linear system (defined the coefficients depending on the porous domain properties and which do not depend on the Fourier variable, boundary conditions, and source terms)
	/*DDPModel ddp_model_init(param,domain,net_mesh,fourier,singularity),ddp_model;
	ddp_model_init.DetermineBasicLinearSystemDDP();*/
	// 2.2. Loop for the simulations
	t1=clock();

	cout << "Experiment computation - begin" << endl;
	//map<int,pair<double,double> >::iterator it1=source_term_positions.begin();	
	
	vector<pair<double,double> > source_term_pos_vect;
	for (map<int,pair<double,double> >::iterator it1=source_term_positions.begin();it1!=source_term_positions.end();it1++){
		source_term_pos_vect.push_back(it1->second);
	}

	//#pragma omp parallel
	#pragma omp parallel for num_threads(4)
	for (int i=0;i<source_term_pos_vect.size();i++){
		printf("Experiment %d computed with thread %d \n",i,omp_get_thread_num());
		nb_elec=i+1;source_position_x=source_term_pos_vect[i].first;source_position_y=source_term_pos_vect[i].second;
		//2.1. Variables
		DDPModel ddp_model_init(param,domain,net_mesh,fourier,singularity);
        	/*ddp_model_init.DetermineBasicLinearSystemDDP();
		SourceTermsDef current_source_term_def;current_source_term_def.push_back(SourceTerm(source_position_x,source_position_y,source_value));
		map<double,ublas_matrix> Potential_Matrices;map<double,ublas_vector> Potential_Fractures;
		// 2.2. Simulation for each wave number
		WaveNumberComputation(param,domain,net_mesh,ddp_model_init,current_source_term_def,Potential_Matrices,Potential_Fractures);
		// 2.3. Write the results
		stringstream ss;string str_cpt;ss<<nb_elec;str_cpt=ss.str();
		string output_file_matrix="Matrix_Potential"+str_cpt;
		output_file_matrix=output_file_matrix+".txt";
		string output_file_fract="Fracture_Potential"+str_cpt;
		output_file_fract=output_file_fract+".txt";
		Results results(param.files_param.output_path+output_file_matrix,Potential_Matrices);
		Results results1(param.files_param.output_path+output_file_fract,Potential_Fractures);
	*/
	}
        t2=clock();
        cout << "Computation time for simulating all the experiments = " << (t2-t1)/CLOCKS_PER_SEC << " seconds" << endl;
}

void WaveNumberComputation(Parameters param,Domain domain,NetworkMeshes net_mesh,DDPModel ddp_model,SourceTermsDef current_source_term_def,map<double,ublas_matrix>& Potential_Matrices,map<double,ublas_vector>& Potential_Fractures){
	vector<double> w_vect=param.simu_param.DDP_param.w_vect;
	ublas_vector DFN_potential;ublas_matrix matrix_potential;
	for (vector<double>::iterator it1=w_vect.begin();it1!=w_vect.end();it1++){
                // 2.2. Define the boundary conditions for ddp model (depend on w)
                ddp_model.DDPModelFourierBC(param,domain,net_mesh.domain,current_source_term_def,*it1,net_mesh.border_map,net_mesh.nodes_map);
                // 2.3. Update the linear system with the coefficients which depend on the Fourier variables, boundary conditions, and source terms
                ddp_model.UpdateLinearSystemFourier(*it1);
                // 2.4. Compute electrical potential
                ddp_model.EvaluatePotential(DFN_potential,matrix_potential);
                // 2.4. Store the computed values
                Potential_Matrices[*it1]=matrix_potential;Potential_Fractures[*it1]=DFN_potential;
	}
}

void FourierPotentialDDPSourceTerm(Parameters param,Domain domain,bool singularity){
	// 1. Global variables and fracture network generation
	// 1.1. Variables
	vector<double> w_vect=param.simu_param.DDP_param.w_vect;
	if (param.simu_param.source_terms.size()!=1){
		cout << "WARNING in main (PERFORM.cpp): configuration not implemented" << endl;
	}
	// 1.2 Fracture network generation
	NetworkMeshes net_mesh(domain,param,param.simu_param.simu_option);
	/*cout << "NETWORK RESTRICTRED TO ITS BACKBONE AND UPPER VALUES" << endl;
	// Compute the backbone
	net_mesh=net_mesh.return_backbone(param,EPSILON);
	// Cut the domain
	cout << "NETWORK RESTRICTRED TO THE UPPER FRACTURES" << endl;
	DFNCutting(net_mesh,25);*/
	bool fourier=true;
	// 2. Simulation for each source position
	// 2.1. Parameters for the simulations
	// standard parameters and variables
	double source_value=param.simu_param.source_terms[0].source_value;
	double source_position_y=param.simu_param.source_terms[0].y_position;
	double Lsimu=domain.domain_size_x();
	int cpt=0;
	// definition of the x-coordinate of the source term
	double source_spacing=param.simu_param.source_terms[0].source_spacing;
	vector<double> vect_source_pos_x,vect_values_I;ublas_matrix elec_positions;
	// several positions for the source terms defined in a text file (electrical resistivity survey)
	if (source_spacing==-1){vect_source_pos_x=ReadSourceTermsCoordX(param);}
	// one position for the source terms
	else if (source_spacing==-2){
		double source_position_x=param.simu_param.source_terms[0].x_position;
		vect_source_pos_x.push_back(source_position_x);
	}
	else{
		double source_position_x=param.simu_param.source_terms[0].x_position;
		vect_source_pos_x.push_back(source_position_x);
		//while (source_position_x<=Lsimu-4*source_spacing){
		while (source_position_x<Lsimu-source_spacing){
			source_position_x+=source_spacing;
			vect_source_pos_x.push_back(source_position_x);
		}
	}
	// Basic linear system (defined the coefficients depending on the porous domain properties and which do not depend on the Fourier variable, boundary conditions, and source terms)
	DDPModel ddp_model_init(param,domain,net_mesh,fourier,singularity),ddp_model;
	ddp_model_init.DetermineBasicLinearSystemDDP();
	ublas_vector potential;
	// 2.2. Loop for the simulations
	for (vector<double>::iterator it1=vect_source_pos_x.begin();it1!=vect_source_pos_x.end();it1++){
		cpt++;
		cout << "Experiment " << cpt << " and x_coord=" << *it1 << endl;
		//2.1. Variables
		param.simu_param.source_terms[0]=SourceTerm(*it1,source_position_y,source_value);
		map<double,ublas_matrix> Potential_Matrices;map<double,ublas_vector> Potential_Fractures;
		ublas_vector DFN_potential;ublas_matrix matrix_potential;
		// 2.2. Simulation for each wave number
		for (vector<double>::iterator it1=w_vect.begin();it1!=w_vect.end();it1++){
			ddp_model=ddp_model_init;
			// 2.1. Update the current wave number
			//param.simu_param.DDP_param.w_current=*it1;
			//ddp_model.current_w=*it1;
			//cout << "w=" << ddp_model.current_w << endl;
			// 2.2. Define the boundary conditions for ddp model (depend on w)
			ddp_model.DDPModelFourierBC(param,domain,net_mesh.domain,param.simu_param.source_terms,*it1,net_mesh.border_map,net_mesh.nodes_map);
			// 2.3. Update the linear system with the coefficients which depend on the Fourier variables, boundary conditions, and source terms
			ddp_model.UpdateLinearSystemFourier(*it1);
			// 2.4. Compute electrical potential
			ddp_model.EvaluatePotential(DFN_potential,matrix_potential);
			// 2.4. Store the computed values
			Potential_Matrices[*it1]=matrix_potential;
			Potential_Fractures[*it1]=DFN_potential;
		}
		stringstream ss;string str_cpt;
		ss<<cpt;
		str_cpt=ss.str();
		string output_file_matrix="Matrix_Potential"+str_cpt;
		output_file_matrix=output_file_matrix+".txt";
		string output_file_fract="Fracture_Potential"+str_cpt;
		output_file_fract=output_file_fract+".txt";
		Results results(param.files_param.output_path+output_file_matrix,Potential_Matrices);
		Results results1(param.files_param.output_path+output_file_fract,Potential_Fractures);
	}
}
// Function to read the x-coordinates of the source term
vector<double> ReadSourceTermsCoordX(Parameters param){
        vector<double> vect_source_pos_x;
        //string file_name_coord=param.files_param.input_path+param.files_param.file_name_source_xcoord;
        string file_name_coord=param.files_param.input_path+param.files_param.file_name_source_positions;
        ifstream file(file_name_coord.c_str());
        if (file.is_open()){
                int Nb_coord;double x_coord;
                file >> Nb_coord;
                for (int i=0;i<Nb_coord;i++){
                        file >> x_coord;
                        vect_source_pos_x.push_back(x_coord);
                }
        }
        else{
                cout << "WARNING in ReadSourceTermsCoordX (FourierModelFunctions.cpp): file not found" << endl;
        }
        return vect_source_pos_x;
}

// Function to read the source term positions
map<int,pair<double,double> > ReadSourceTermsCoordinates(Parameters param){
        map<int,pair<double,double> > source_positions;
        string file_name_coord=param.files_param.input_path+param.files_param.file_name_source_positions;
        ifstream file(file_name_coord.c_str());
        int elec_nb;double coord_x,coord_y;
        if (file.is_open()){
                while (!file.eof()){
                        file >> elec_nb >> coord_x >> coord_y;
                        source_positions[elec_nb]=make_pair(coord_x,coord_y);
                }
        }
        else{
                cout << "WARNING in ReadSourceTermsCoordX (FourierModelFunctions.cpp): file not found" << endl;
        }
        return source_positions;
}


