/*
 * DDPStudies.cpp
 *
 *  Created on: Apr 29, 2013
 *      Author: delphineroubinet
 */

#include "DDPStudies.h"
#include "DDPModel.h"
#include "../Utilitaries/Geometry/scale.h"
#include "../InputOutput/Results.h"
#include "../Utilitaries/Storage/Structures.h"
#include "../Utilitaries/Storage/UtilitaryStorage.h"
#include "../Utilitaries/Math/Operations.h"
#include "DDPDFNFunctions.h"
#include "../InputOutput/Results.h"
#include "ExchangeFunctions.h"
#include "../DFNModel/DFNGeneration.h"


using namespace std;


void PotentialComputationDDP(Parameters param,Domain domain,double min_visu,double max_visu){
	// DFN definition
	NetworkMeshes net_mesh(domain,param,param.simu_param.simu_option);
	bool singularity=param.simu_param.singularity;

	// Potential computation
	ublas_vector DFN_potential;ublas_matrix matrix_potential;
	DDPModel ddp_model(param,domain,net_mesh,singularity);
	// Steady-state simulation
	ddp_model.ComputePotential(DFN_potential,matrix_potential);

	cout << "Fracture potential" << endl;
	print_vector(DFN_potential);

	cout << "Matrix potential:" << endl;
	print_matrix(matrix_potential);

	double output_fract_surf,flow_fract,flow_matrix;
	flow_fract=TotalElectricCurrentDFN(param,net_mesh,DFN_potential,ddp_model.bc_map,output_fract_surf);

	cout << "flow fracture = " << flow_fract << endl;

	flow_matrix=TotalElectricCurrentEPM(param,matrix_potential,ddp_model.porous_cond,RIGHT_BORDER);

	cout << "flow matrix = " << flow_matrix << endl;

}


void ConductivityPlotDDP(Parameters param,Domain domain){
	// logarithm discretization of the matrix conductivity
	int N_sigma=11;
	double sigma_min=1e-12,sigma_max=1e-2,current_sigma;
	scale sigma_m_scale(sigma_min,sigma_max,N_sigma,1);
	vector<double> sigma_m_values=sigma_m_scale.scale_vector();
	ublas_matrix matrix_potential;ublas_vector potential_dfn;
	double output_fract_surf,flow,Lx=param.domain_param.Lx,Ly=param.domain_param.Ly,pot_diff=1;
	map<double,double> eq_cond;
	NetworkMeshes net_mesh(domain,param,param.simu_param.simu_option);
	double sigma_f=net_mesh.meshes[0].conductivity;
	for (int i=0;i<N_sigma;i++){
		cout << i << "/" << N_sigma << endl;
		current_sigma=sigma_m_values[i];
		param.domain_param.mat_cond=current_sigma;
		// Discrete dual-porosity model definition
		DDPModel ddp_model(param,domain,net_mesh);
		// Potential computation
		ddp_model.ComputePotential(potential_dfn,matrix_potential);
		// Flow computation
		flow=TotalElectricCurrentDFN(param,net_mesh,potential_dfn,ddp_model.bc_map,output_fract_surf);
		// TO MODIFZY: assume a homogeneous conductivity in the matrix part
		flow+=TotalElectricCurrentEPM(param,matrix_potential,ddp_model.porous_cond,RIGHT_BORDER)/Lx*(Lx-output_fract_surf);
		//flow+=TotalElectricCurrentEPM(param,matrix_potential,RIGHT_BORDER);
		//eq_cond[current_sigma/sigma_f]=flow/Lx*Ly/pot_diff;
		eq_cond[current_sigma]=flow/Lx*Ly/pot_diff;
	}
	string output_file=param.files_param.output_path+param.files_param.output_file;
	Results results(output_file,eq_cond);
}

void ExchangeCoeffBlockScale(EPMSystem epm_model,ublas_matrix potential,double alpha,double& FMExchange_num,double& FMExchange_ddp){

	// Loop on the meshes corresponding to projected fractures
	vector<pair<int,int> > fract_indices=epm_model.fract_indices;

	if (fract_indices.size()==0){
		//cout << "NO  FRACTURES" << endl;
		FMExchange_num=0;
		FMExchange_ddp=0;
		return;
	}

	// Computation of the averaged matrix potential in each block
	double MatrixPotBlock=0;
	int nb_fract_pot=fract_indices.size();
	int nb_matrix_pot=potential.size1()*potential.size2()-nb_fract_pot;
	for (int i=0;i<potential.size1();i++){
		for (int j=0;j<potential.size2();j++){
			if (!exist_in_vector(fract_indices,make_pair(i,j))){
				MatrixPotBlock+=potential(i,j);
			}
		}
	}
	MatrixPotBlock=MatrixPotBlock/nb_matrix_pot;

	double delta_x=epm_model.delta_x,delta_y=epm_model.delta_y;
	int Nx=epm_model.Nx,Ny=epm_model.Ny;double fm_flow=0,fm_flow_ddp=0;double cumul_length=0;
	FMExchange_num=0;FMExchange_ddp=0;

	for (vector<pair<int,int> >::iterator it=fract_indices.begin();it!=fract_indices.end();it++){
		// 1.3.1. Computation of the exchange with the surrounding matrix meshes
		// Check if the surrounding meshes are matrix meshes
		if (it->first>0&&!exist_in_vector(fract_indices,make_pair(it->first-1,it->second))){
			fm_flow+=ComputeFMFlowMesh(epm_model,potential,it->first,it->second,it->first-1,it->second,delta_x)*delta_y;
			cumul_length+=delta_y;
			fm_flow_ddp+=alpha*(potential(it->first,it->second)-MatrixPotBlock)*delta_y;

		}
		if (it->first<Nx&&!exist_in_vector(fract_indices,make_pair(it->first+1,it->second))){
			fm_flow+=ComputeFMFlowMesh(epm_model,potential,it->first,it->second,it->first+1,it->second,delta_x)*delta_y;
			cumul_length+=delta_y;
			fm_flow_ddp+=alpha*(potential(it->first,it->second)-MatrixPotBlock)*delta_y;

		}
		if (it->second>0&&!exist_in_vector(fract_indices,make_pair(it->first,it->second-1))){
			fm_flow+=ComputeFMFlowMesh(epm_model,potential,it->first,it->second,it->first,it->second-1,delta_y)*delta_x;
			cumul_length+=delta_x;
			fm_flow_ddp+=alpha*(potential(it->first,it->second)-MatrixPotBlock)*delta_x;

		}
		if (it->second<Ny&&!exist_in_vector(fract_indices,make_pair(it->first,it->second+1))){
			fm_flow+=ComputeFMFlowMesh(epm_model,potential,it->first,it->second,it->first,it->second+1,delta_y)*delta_x;
			cumul_length+=delta_x;
			fm_flow_ddp+=alpha*(potential(it->first,it->second)-MatrixPotBlock)*delta_x;

		}
	}
	FMExchange_num=fm_flow/cumul_length;
	FMExchange_ddp=fm_flow_ddp/cumul_length;

}

double ComputeFMFlowMesh(EPMSystem epd_model,ublas_matrix potential,int indexx1,int indexy1,int indexx2,int indexy2,double delta){
	double cond1=epd_model.porous_cond(indexx1,indexy1),cond2=epd_model.porous_cond(indexx2,indexy2);
	double pot1=potential(indexx1,indexy1),pot2=potential(indexx2,indexy2);
	return Harmonic_Mean(cond1,cond2)*(pot1-pot2)/delta;
}

void ExchangeStudy(Parameters param,Domain domain,int nb_div,int nb_seed){

		// 1. Parameters and variables of the study
		int nb_blocks=2;	// number of blocks for each division
		double alpha_ddp;
		double cond_mat=param.domain_param.mat_cond;
		map<int,double> flow_num_div_init,flow_ddp_div_init;
		map<int,double> flow_map_num,flow_map_ddp,nb_block_map;
		int Nx_total=param.simu_param.EPM_param.Nx,Ny_total=param.simu_param.EPM_param.Ny;
		// Initialization
		for (int num_div=1;num_div<=nb_div;num_div++){
			flow_map_num[num_div]=0;
			flow_map_ddp[num_div]=0;
			if (num_div==1){
				nb_block_map[num_div]=num_div;
			}
			else{
				nb_block_map[num_div]=nb_block_map[num_div-1]*nb_blocks*nb_blocks;
			}
		}

		// 2. Loop on the number of generated domain (stochastic study)
		for (int num_simu=0;num_simu<nb_seed;num_simu++){

			cout << "num_simu=" << num_simu << endl;

			// 2.1 Definition of the initial domain and first computation
			// Update the seed for generation
			param.simu_param.DFN_param.seed=param.simu_param.DFN_param.seed+num_simu;
			// Definition of the domain
			map<int,Domain> init_domain; init_domain[0]=domain;
			// Define the fracture network
			NetworkMeshes net_mesh_init(domain,param,param.simu_param.simu_option);
			FracturesIntersectionGrid(Nx_total,Ny_total,net_mesh_init);
			SubNetworkMap init_net_mesh; init_net_mesh[0]=net_mesh_init;
			// Define the EPM system
			map<int,EPMSystem> init_epm_model;init_epm_model[0]=EPMSystem(param,domain,init_net_mesh[0],NO_INTERSECTION);

			// Potential computation on the whole domain
			init_epm_model[0].PotentialComputation();
			// Exchanged flow computation
			ComputationSumFMExchange(init_domain,init_net_mesh,init_epm_model,cond_mat,flow_num_div_init,flow_ddp_div_init);
			flow_map_num[1]=average_map(flow_num_div_init);
			flow_map_ddp[1]=average_map(flow_ddp_div_init);

			// 2.2. Loop on the division of the domain
			for (int num_div=2;num_div<=nb_div;num_div++){

				cout << "num_div=" << num_div << endl;

				// 2.2.1. New variables for each division
				map<int,Domain> new_domain_final,new_domain_current;
				SubNetworkMap new_net_mesh_final,new_net_mesh_current;
				map<int,EPMSystem> new_epm_model_final,new_epm_model_current;
				map<int,double> flow_num_div,flow_ddp_div,flow_ddp_div_bis;

				// 2.2.2. Loop on the initial subdomain
				for (SubNetworkMap::iterator it=init_net_mesh.begin();it!=init_net_mesh.end();it++){
					NetworkMeshes current_net_mesh=it->second;
					// 2.2.2.1. Division (or extraction) of the domain, network and potential
					new_domain_current=DefineSubDomains(init_domain[it->first],nb_blocks,nb_blocks);
					new_net_mesh_current=DefineSubNetwork(new_domain_current,current_net_mesh);
					new_epm_model_current=init_epm_model[it->first].DefineSubEPMFract(new_domain_current,new_net_mesh_current,cond_mat,NO_INTERSECTION);
					// 2.2.2.2. Computation of the averaged exchanged flow at the block scale
					ComputationSumFMExchange(new_domain_current,new_net_mesh_current,new_epm_model_current,cond_mat,flow_num_div,flow_ddp_div);
					// 2.2.2.3. Storage of the new domains for next iteration
					insert_map(new_domain_final,new_domain_current);
					insert_map(new_net_mesh_final,new_net_mesh_current);
					insert_map(new_epm_model_final,new_epm_model_current);
				}

				// 2.2.3. Update initial values
				init_domain=new_domain_final;
				init_net_mesh=new_net_mesh_final;
				init_epm_model=new_epm_model_final;

				// 2.2.4. Compute and store the averaged flow
				flow_map_num[num_div]=average_map(flow_num_div);
				flow_map_ddp[num_div]=average_map(flow_ddp_div);
			}
		}

		// Average over the number of simulation
		for (int num_div=1;num_div<=nb_div;num_div++){
			flow_map_num[num_div]=flow_map_num[num_div]/nb_seed;
			flow_map_ddp[num_div]=flow_map_ddp[num_div]/nb_seed;
		}

		// Write results in output files
		string output_file=param.files_param.output_path+param.files_param.output_file,output_file1=param.files_param.output_path+param.files_param.output_file1;
		Results res(output_file,output_file1,flow_map_num,flow_map_ddp,nb_block_map);
}

/*void ExchangeStudy2(Parameters param,Domain domain){

	// 1. Parameters and variables of the study
	if (param.simu_param.DDP_param.Nx!=param.simu_param.DDP_param.Ny){
		cout << "WARNING in ExchangeStudy2: configuration not implemented" << endl;
	}
	int nb_block_max=param.simu_param.DDP_param.Nx;
	int nb_seed=param.simu_param.DDP_param.nb_seed;
	DDPModel ddp_model;SubNetworkMap subnetwork;map<int,Domain> subdomains;
	ublas_vector dfn_potential;ublas_matrix matrix_potential;
	double fm_exchange;

	// 2. Loop on the number of generated domain (stochastic study)
	for (int num_simu=0;num_simu<nb_seed;num_simu++){

		cout << "num_simu=" << num_simu << endl;

		// 2.1 Definition of the initial domain and first computation
		// Update the seed for generation
		param.simu_param.DFN_param.seed=param.simu_param.DFN_param.seed+num_simu;
		// Define the fracture network
		NetworkMeshes net_mesh(domain,param,param.simu_param.simu_option);

		// 2.2. Loop on the division of the domain
		for (int nb_blocks=2;nb_blocks<=nb_block_max;nb_blocks++){

			cout << "nb_block=" << nb_blocks << endl;

			// 2.2.1. Build the DDP model
			param.simu_param.DDP_param.Nx=nb_blocks;
			param.simu_param.DDP_param.Ny=nb_blocks;
			ddp_model=DDPModel(param,domain,net_mesh);
			ddp_model.ComputePotential(dfn_potential,matrix_potential);

			// 2.2.2. Define the subdomains at the block scale
			subdomains=DefineSubDomains(domain,nb_blocks,nb_blocks);
			subnetwork=DefineSubNetwork(subdomains,net_mesh);

			// 2.2.2. Computation of the fracture-matrix exchange at the block scale
			fm_exchange=0;
			for (vector<FractureMesh>::iterator it=net_mesh.meshes.begin();it!=net_mesh.meshes.end();it++){

			}

			for (SubNetworkMap::iterator it=subnetwork.begin();it!=subnetwork.end();it++){
				fm_exchange+=
			}
		}
	}

	// Average over the number of simulation
	for (int num_div=1;num_div<=nb_div;num_div++){
		flow_map_num[num_div]=flow_map_num[num_div]/nb_seed;
		flow_map_ddp[num_div]=flow_map_ddp[num_div]/nb_seed;
	}

	// Write results in output files
	string output_file=param.files_param.output_path+param.files_param.output_file,output_file1=param.files_param.output_path+param.files_param.output_file1;
	Results res(output_file,output_file1,flow_map_num,flow_map_ddp,nb_block_map);
}*/

void ComputationSumFMExchange(map<int,Domain> domain_map,SubNetworkMap network_map,map<int,EPMSystem> epm_map,double cond_mat,map<int,double>& flow_num_div,map<int,double>& flow_ddp_div){

	// 1. Variables
	double alpha_ddp;
	double flow_num,flow_ddp;
	// 2. Loop on the different domains
	for (map<int,Domain>::iterator it=domain_map.begin();it!=domain_map.end();it++){
		alpha_ddp=ReturnAlphaBlock(cond_mat,it->second.domain_size_x(),it->second.domain_size_y(),network_map[it->first]);
		ExchangeCoeffBlockScale(epm_map[it->first],epm_map[it->first].potential,alpha_ddp,flow_num,flow_ddp);
		flow_num_div[flow_num_div.size()]=flow_num/cond_mat;
		flow_ddp_div[flow_ddp_div.size()]=flow_ddp/cond_mat;
	}
}

void ConductivityPlotEPMDFN(Parameters param,Domain domain){
	// logarithm discretization of the matrix conductivity
	int N_sigma=11;
	double sigma_min=1e-12,sigma_max=1e-2,current_sigma;
	scale sigma_m_scale(sigma_min,sigma_max,N_sigma,1);
	vector<double> sigma_m_values=sigma_m_scale.scale_vector();
	ublas_matrix matrix_potential;
	double output_fract_surf,flow,Lx=param.domain_param.Lx,Ly=param.domain_param.Ly,pot_diff=1;
	map<double,double> eq_cond;
	NetworkMeshes net_mesh(domain,param,param.simu_param.simu_option);
	cout << "FracturesIntersectionGrid" << endl;
	FracturesIntersectionGrid(param.simu_param.EPM_param.Nx,param.simu_param.EPM_param.Ny,net_mesh);
	double sigma_f=net_mesh.meshes[0].conductivity;
	for (int i=0;i<N_sigma;i++){
		cout << i << "/" << N_sigma << endl;
		current_sigma=sigma_m_values[i];
		param.domain_param.mat_cond=current_sigma;
		// EPM model with fractures
		cout << "epd_model" << endl;
		EPMSystem epd_model(param,domain,net_mesh,"NO_INTERSECTION");
		// Potential computation
		cout << "PotentialComputation" << endl;
		matrix_potential=epd_model.PotentialComputation();
		// Flow computation
		flow=TotalElectricCurrentEPM(param,matrix_potential,epd_model.porous_cond,RIGHT_BORDER);
		//eq_cond[current_sigma/sigma_f]=flow/Lx*Ly/pot_diff;
		eq_cond[current_sigma]=flow/Lx*Ly/pot_diff;
	}
	string output_file=param.files_param.output_path+param.files_param.output_file;
	Results results(output_file,eq_cond);
}

void PotentialStudy(Parameters param,Domain domain,int nb_div,int nb_seed){

		// 1. Parameters and variables of the study
		int nb_blocks;	// number of blocks for each division
		double cond_mat=param.domain_param.mat_cond;
		int Nx_total=param.simu_param.EPM_param.Nx,Ny_total=param.simu_param.EPM_param.Ny;
		ublas_vector DFN_potential;ublas_matrix matrix_potential_ddp;
		DDPModel ddp_model;
		map<int,Domain> extracted_domains;
		map<int,ublas_matrix> extracted_potential;
		map<int,EPMSystem> extracted_epm_model;
		map<int,double> averaged_potential,map_potential_ddp;
		map<int,double> error_div;
		for (int num_div=2;num_div<=nb_div;num_div++){
			nb_blocks=pow(2.,num_div-1);
			error_div[nb_blocks]=0;
		}

		// 2. Loop on the number of generated domain (stochastic study)
		for (int num_simu=0;num_simu<nb_seed;num_simu++){

			cout << "num_simu=" << num_simu << endl;

			// 2.1 Network generation
			// Update the seed for generation
			param.simu_param.DFN_param.seed=param.simu_param.DFN_param.seed+num_simu;
			// Define the fracture network
			NetworkMeshes net_mesh(domain,param,param.simu_param.simu_option);
			FracturesIntersectionGrid(Nx_total,Ny_total,net_mesh);

			// 2.2. Computation of the reference potential
			EPMSystem epmdfn_model=EPMSystem(param,domain,net_mesh);
			epmdfn_model.PotentialComputation();

			// 2.3. Loop on the number of division for the ddp model
			for (int num_div=2;num_div<=nb_div;num_div++){

				cout << "num_div=" << num_div << endl;

				// 2.3.1. Computation of the potential with ddp simulation
				nb_blocks=pow(2.,num_div-1);
				param.simu_param.DDP_param.Nx=nb_blocks;
				param.simu_param.DDP_param.Ny=nb_blocks;
				ddp_model=DDPModel(param,domain,net_mesh);
				ddp_model.ComputePotential(DFN_potential,matrix_potential_ddp);

				// 2.3.2. Extraction of the reference solution at the block scale
				extracted_domains=DefineSubDomains(domain,nb_blocks,nb_blocks);
				extracted_epm_model=epmdfn_model.DefineSubEPM(extracted_domains,cond_mat);
				averaged_potential=AveragedMatrices(extracted_epm_model);

				// 2.2.3.2. Error computation
				map_potential_ddp=MatrixToMap(matrix_potential_ddp,nb_blocks);
				//error_div[nb_blocks]+=RelativeErrorMap(map_potential_ddp,averaged_potential)/nb_seed;
				error_div[nb_blocks]+=ErrorMapL2(map_potential_ddp,averaged_potential)/nb_seed;
			}
		}

		// Write results in output files
		string output_file=param.files_param.output_path+param.files_param.output_file,output_file1=param.files_param.output_path+param.files_param.output_file1;
		Results res(output_file,error_div);
}



