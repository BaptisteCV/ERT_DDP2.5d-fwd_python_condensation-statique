//============================================================================
// Name        : MAIN.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : C++ code for solute transport in heterogeneous fractured porous media
//============================================================================

#include "InputOutput/Parameters.h"
#include "DomainDefinition/Domain.h"
#include "DFNModel/NetworkMeshes.h"
#include "InputOutput/Results.h"
#include "Utilitaries/Storage/Structures.h"
#include "DFNModel/DFNComputation.h"
#include "DFNModel/DFNGeneration.h"
#include "Utilitaries/Visualization/VisuStructures.h"
#include "Utilitaries/Visualization/DisplayResults.h"
#include "EPMModel/EPMDefinition.h"
#include "EllipseComputation/EllipseComputation.h"
#include "DDPModel/DDPModel.h"
#include "DDPModel/DDPStudies.h"
#include "Utilitaries/Math/Operations.h"
#include "DDPModel/ExchangeFunctions.h"
#include "FourierModel/FourierModelFunctions.h"
#include "FourierModel/FourierModelFunctionsCSEigen.h"

#include <iostream>
#include <string>
#include <omp.h>
#include <map>

using namespace std;

int main(int argc, char **argv) {

	clock_t t_begin=clock(),t_end;
	cout << "Electrical Resistivity is starting..." << endl;    

	bool singularity=false,bc_on_border=false;
	// 1. Parameter and domain definition reading
	// 1.1. Parameter reading
	Parameters param(singularity);
	// 1.2. Domain definition
	Domain domain(param.domain_param.Lx,param.domain_param.Ly);
	// 1.3. Display parameters
	double min_visu=-1,max_visu=-1;
    // 1.4. Output results (using ublas sparse matrices)
    std::map<int, ublas_matrix > res_map;

	// 2. Simulation
	// 2.1. DFN generation
	if (param.simu_param.simu_option==DFN_GENERATION){
		//DFNVisu DFN_visu=DFNGeneration(param,domain); // uncomment if DFN_GENERATION option
        //NetworkMeshes net_mesh_inter=DFN_visu.net_mesh; // uncomment if DFN_GENERATION option
	// Compute the backbone
		//NetworkMeshes net_mesh1=DFN_visu.net_mesh,backbone;
		//backbone=net_mesh1.return_backbone(param,EPSILON);
		//DFN_visu.net_mesh=backbone;

		// Cut the domain
		//NetworkMeshes net_mesh2=DFN_visu.net_mesh;
		//DFNCutting(net_mesh2,26.6);
		//DFNCutting(net_mesh2,25);
		//DFN_visu.net_mesh=net_mesh2;

		// DFN writing
		//Results results(param.files_param.output_path+param.files_param.output_file,DFN_visu.net_mesh); // uncomment if DFN_GENERATION option
		// TODO remove display function
		//DisplayResults(argc,argv,DFN_visu);
	}

	else if (param.simu_param.simu_option==FOURIER_POTENTIAL_COMPUTATION){
		// Equivalent Porous Media
		if (param.simu_param.model_option==EPM){
			// Read the source term positions
			map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);
			// Run the simulations
			if (source_term_positions.size()>0){
			    FourierPotentialEPMSourceTerm(param,domain,singularity,bc_on_border,source_term_positions);
			 }
			else if (source_term_positions.size()==0){FourierPotentialEPM(param,domain,singularity);}
			else{cout << "WARNING in main (MAIN.cpp): option not implemented" << endl;}
		}
		else if (param.simu_param.model_option==EPM_CS){
			map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);
			if (source_term_positions.size()>0){
			    eigen_cs::FourierPotentialEPMSourceTermCSEigen(param,domain,singularity,bc_on_border,source_term_positions);
			}
			else if (source_term_positions.size()==0){
				cout << "WARNING in main (MAIN.cpp): EPM_CS without source positions falls back to standard EPM workflow" << endl;
				FourierPotentialEPM(param,domain,singularity);
			}
			else{cout << "WARNING in main (MAIN.cpp): option not implemented" << endl;}
		}
		// Discrete-Dual-Porosity Eigen natif (pas de condensation statique)
		else if (param.simu_param.model_option==DDP_EIGEN){
			map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);
			if (source_term_positions.size()>0){
				double t_wall = omp_get_wtime();
				FourierPotentialDDPSourceTermEigen(param,domain,singularity,source_term_positions);
				cout << "DDP_EIGEN Wall Time = " << omp_get_wtime() - t_wall << " sec" << endl;
			}
			else{ cout << "WARNING in main (MAIN.cpp): DDP_EIGEN without source positions not implemented" << endl; }
		}
		// Discrete-Dual-Porosity with static condensation (Eigen)
		else if (param.simu_param.model_option==DDP_CS){
			map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);
			if (source_term_positions.size()>0){
				double t_wall_ddpcs = omp_get_wtime();
				eigen_cs::FourierPotentialDDPSourceTermCSEigen(param,domain,singularity,bc_on_border,source_term_positions);
				cout << "DDP_CS Wall Time = " << omp_get_wtime() - t_wall_ddpcs << " sec" << endl;
			}
			else if (source_term_positions.size()==0){
				cout << "WARNING in main (MAIN.cpp): DDP_CS without source positions not yet implemented" << endl;
			}
			else{cout << "WARNING in main (MAIN.cpp): option not implemented" << endl;}
		}
		// Discrete-Dual-Porosity model for fractured porous domains
		else if (param.simu_param.model_option==DDP){

		    // Read the source term positions
			map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);

            //----------------------- BEGIN TEST -------------------------------------
			#ifdef TEST_SMALL_MATRIX
            // 1000 matrices of size 8x2
            cout << "TEST SMALL MATRIX 1000 matrices of size 8x2" << endl;
            for (int i=129;i<=1000;i++) {
				cout << "i = " << i << endl;
                source_term_positions[i] = make_pair((double)rand()/RAND_MAX*13,(double)rand()/RAND_MAX*16);
            }
            param.simu_param.DDP_param.Nx=2;param.simu_param.DDP_param.Ny=2;
            #endif

            #ifdef TEST_BIG_MATRIX
            cout << "TEST BIG MATRIX 10 matrices of size 8000x2000" << endl;
            // Changes for test2: 10 matrices of size 8000x2000
            map<int, pair<double, double> > source_term_positions_bis;
            for (int i=1;i<=10;i++) {
                source_term_positions_bis[i] = make_pair((double)rand()/RAND_MAX*13,(double)rand()/RAND_MAX*16);
            }
            source_term_positions=source_term_positions_bis;
            param.simu_param.DDP_param.Nx=2000;param.simu_param.DDP_param.Ny=2000;
            #endif
            //---------------------- END TEST -----------------------------------------

			// Run the simulations
			if (source_term_positions.size()>0){
				double t_wall_ddp = omp_get_wtime();
				FourierPotentialDDPSourceTerm(param,domain,singularity,source_term_positions);
				cout << "DDP Wall Time = " << omp_get_wtime() - t_wall_ddp << " sec" << endl;
			}
			else if (source_term_positions.size()==0){FourierPotentialDDP(param,domain,singularity,bc_on_border);}
			else{cout << "WARNING in main (MAIN.cpp): option not implemented" << endl;}

		}
		else{cout << "WARNING in main (MAIN.cpp): solving method not implemented" << endl;}
	}


	// 2.2. Potential computation and visualization
	else if (param.simu_param.simu_option==POTENTIAL_COMPUTATION){
		// Discrete Fracture Network representation
		if (param.simu_param.model_option==DFN){
			//DFNVisu fracture_network=PotentialComputationDFN(param,domain);
			//DisplayResults(argc,argv,fracture_network);
		}
		// Equivalent Porous Media
		else if (param.simu_param.model_option==EPM){
			EPMSystem epd_model(param,domain);
			ublas_matrix EPM_Potential=epd_model.PotentialComputation();
			Results results(param.files_param.output_path+param.files_param.output_file,EPM_Potential);
		}
		else if (param.simu_param.model_option==EPM_CS){
			cout << "WARNING in main (MAIN.cpp): EPM_CS is only used in Fourier source workflow, fallback to standard EPM for POTENTIAL_COMPUTATION" << endl;
			EPMSystem epd_model(param,domain);
			ublas_matrix EPM_Potential=epd_model.PotentialComputation();
			Results results(param.files_param.output_path+param.files_param.output_file,EPM_Potential);
		}
		// Equivalent Porous Media with a projected DFN
		else if (param.simu_param.model_option==EPMDFN){
			NetworkMeshes net_mesh(domain,param,param.simu_param.simu_option);
			//EPMVisu EPM_Potential=PotentialComputationEPMDFN(param,domain,net_mesh,min_visu,max_visu);
			//DisplayResults(argc,argv,EPM_Potential);
		}
		// Discrete Dual Porosity representation
		else if (param.simu_param.model_option==DDP){
			PotentialComputationDDP(param,domain,min_visu,max_visu);
			//DisplayResults(argc,argv,EPMDFN_Potential,"Matrix");	// "Matrix/Current" to display electric potential or current
		}
		else{cout << "WARNING in main (MAIN.cpp): solving method not implemented" << endl;}
	}

	// 2.3. Permeability ellipse computation for DFN domain
	else if (param.simu_param.simu_option==PERMEABILITY_ELLIPSE){
		EllipseComputation(param,domain,param.simu_param.model_option);
	}

	else if (param.simu_param.simu_option==CONDUCTIVITY_PLOT){
		if (param.simu_param.model_option==DDP){
			ConductivityPlotDDP(param,domain);
		}
		else if (param.simu_param.model_option==EPMDFN){
			ConductivityPlotEPMDFN(param,domain);
		}
	}

	// 2.4. Study of the flow exchanged between fractures and matrix at the matrix-block scale
	else if (param.simu_param.simu_option==EXCHANGE_STUDY){
		//int nb_div=6,nb_seed=1; // nb_div:number of times that we divide the domain and nb_seed:number of generated domain
		//ExchangeStudy(param,domain,nb_div,nb_seed);
		//ExchangeStudy2(param,domain);
	}
	// 2.5. Comparison of the potential with DDP and EPMDFN
	else if (param.simu_param.simu_option==POTENTIAL_STUDY){
		int nb_div=param.simu_param.DDP_param.nb_div,nb_seed=param.simu_param.DDP_param.nb_seed;
		PotentialStudy(param,domain,nb_div,nb_seed);
	}
	else if (param.simu_param.simu_option=="TEST"){
		cout << "TEST" << endl;
		int size1=1e5,size2=size1;
		ublas_matrix matrix_syst(size1,size2);
		//print_matrix(matrix_syst);
	}
	else{cout << "WARNING in main (MAIN.cpp): simulation option not implemented" << endl;}

	t_end=clock();
	cout << "CPU Time = " << (t_end-t_begin)/CLOCKS_PER_SEC << endl;
	//cout << "CPU Time = " << (t_end-t_begin) << endl;

    return 0;
}
