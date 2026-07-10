/*
 * EllipseComputation.cpp
 *
 *  Created on: May 1, 2013
 *      Author: delphineroubinet
 */

#include "EllipseComputation.h"
#include "../DFNModel/DFNGeneration.h"
#include "../DFNModel/DFNComputation.h"
#include "../EPMModel/EPMDefinition.h"
#include "../DDPModel/DDPModel.h"
#include "../InputOutput/Results.h"

using namespace std;


/**********************************************/
// Evaluation of Ellipse Permeability //
/**********************************************/

// Compute the equivalent permeability of a DFN structure for several angles of the boundary conditions
void EllipseComputation(Parameters param,Domain domain,string model_option){
	// 0. Fracture network definition
	NetworkMeshes init_net_mesh(domain,param,param.simu_param.simu_option);
	// 1. Variables
	map<double,double> perm_ellipse;
	int Nangle=param.simu_param.nb_angle_step;
	double current_angle,extract_length=0.5*init_net_mesh.ReturnDomainLength(),flow,pot_diff=1,
			Lx=extract_length,Ly=extract_length,output_fract_surf;
	NetworkMeshes current_net_mesh,new_net_mesh;ublas_vector potential_dfn;BoundaryConditionsDFN bc_map_dfn;
	//BoundaryConditionsDef bc_map_classic=DefineBoundaryConditions(CONFIG1);
	//BoundaryConditionsDef bc_map_classic=DefineBoundaryConditions(CONFIG2);
	BoundaryConditionsDef bc_map_classic=DefineBoundaryConditions(CONFIG3);
	string face=RIGHT_BORDER;
	param.simu_param.bc_map=bc_map_classic;
	ublas_matrix matrix_potential;
	// 2. Angular loop
	for (int i=0;i<Nangle;i++){
		cout << "Angle " << i << endl;
		// 2.1. Current angle definition
		current_angle=i*2*PI/(Nangle-1);
		// 2.2. Fracture network rotation
		current_net_mesh=DFNRotation(init_net_mesh,current_angle);
		// 2.3. Centered fracture network extraction
		if (!DFNExtraction(current_net_mesh,extract_length)){cout << "WARNING in EllipseComputationRegular (EllipseComputation.cpp): the extracted network is empty" << endl;}
		param.domain_param.Lx=extract_length;param.domain_param.Ly=extract_length;
		// 2.4. Equivalent permeability computation
		if (model_option==DFN){
			// Boundary condition definition
			bc_map_dfn=ReturnBoundCondDFN(current_net_mesh.domain,bc_map_classic,current_net_mesh.border_map,current_net_mesh.nodes_map);
			// Electric potential computation
			potential_dfn=DFNComputationClassic(current_net_mesh,bc_map_dfn);
			// Flow evaluation
			flow=TotalElectricCurrentDFN(param,current_net_mesh,potential_dfn,bc_map_dfn,output_fract_surf);
		}
		else if (model_option==EPM){
			EPMSystem epd_model(param,domain);
			matrix_potential=epd_model.PotentialComputation();
			flow=TotalElectricCurrentEPM(param,matrix_potential,epd_model.porous_cond,face);
		}
		else if (model_option==DDP){
			bool singularity=param.simu_param.singularity;
			// Discrete dual-porosity model definition
			DDPModel ddp_model(param,domain,current_net_mesh,singularity);
			// Potential computation
			ddp_model.ComputePotential(potential_dfn,matrix_potential);
			// Flow computation
			flow=TotalElectricCurrentDFN(param,current_net_mesh,potential_dfn,ddp_model.bc_map,output_fract_surf);
			flow+=TotalElectricCurrentEPM(param,matrix_potential,ddp_model.porous_cond,face)/Lx*(Lx-output_fract_surf);
		}
		else{cout << "WARNING in EllipseComputationRegular (EllipseComputation.cpp): model_option not defined" << endl;}
		// 2.5. Equivalent permeability affectation
		perm_ellipse[current_angle]=flow/Lx*Ly/pot_diff;
	}
	// Polar permeability
	Results results(param.files_param.output_path+param.files_param.output_file,perm_ellipse);
}


