//============================================================================
// Name        : PERFORM.cpp
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

// Tests
//#include </home/leal/Desktop/Doctorat/Scripts/MontpInv_Python/transdimension/refonte_matlab_python/montp_inv/Library/Eigen/IterativeLinearSolvers>
//#include </usr/include/eigen3/Eigen/IterativeLinearSolvers>
//#include </usr/local/lib/python3.10/dist-packages/pybind11/include/pybind11/pybind11.h>
//#include </usr/local/lib/python3.10/dist-packages/pybind11/include/pybind11/stl.h>
//#include </usr/local/lib/python3.10/dist-packages/pybind11/include/pybind11/eigen.h>

// original versions
#include <Eigen/IterativeLinearSolvers>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <map>

using Eigen::SparseMatrix;
using namespace std;
namespace py = pybind11;

typedef Eigen::Triplet<double> T;

static Eigen::SparseMatrix<double> UblasToEigenSparse(const ublas_matrix& in) {
	Eigen::SparseMatrix<double> out(in.size1(), in.size2());
	std::vector<T> triplets;
	for (int i = 0; i < (int)in.size1(); ++i) {
		for (int j = 0; j < (int)in.size2(); ++j) {
			double v = in(i, j);
			if (v != 0.0) {
				triplets.push_back(T(i, j, v));
			}
		}
	}
	out.setFromTriplets(triplets.begin(), triplets.end());
	return out;
}

static std::map<int, Eigen::SparseMatrix<double> > UblasMapToEigenMap(const std::map<int, ublas_matrix>& in) {
	std::map<int, Eigen::SparseMatrix<double> > out;
	for (std::map<int, ublas_matrix>::const_iterator it = in.begin(); it != in.end(); ++it) {
		out[it->first] = UblasToEigenSparse(it->second);
	}
	return out;
}

std::map<int, Eigen::SparseMatrix<double>> perform() {
//std::map<int, Eigen::SparseMatrix<double>, std::less<int>,
//         Eigen::aligned_allocator<std::pair<const int, Eigen::SparseMatrix<double>> > > perform() {

	clock_t t_begin=clock(),t_end;

	//cout << "PERFORM begin" << endl;

	bool singularity=false,bc_on_border=false;
	// 1. Parameter and domain definition reading
	// 1.1. Parameter reading
	Parameters param(singularity);
	// 1.2. Domain definition
	Domain domain(param.domain_param.Lx,param.domain_param.Ly);
	// 1.3. Display parameters
	double min_visu=-1,max_visu=-1;
    // 1.4. Output results
    std::map<int, Eigen::SparseMatrix<double> > res_map;
    //std::map<int, Eigen::SparseMatrix<double>, std::less<int>,
    //     Eigen::aligned_allocator<std::pair<const int, Eigen::SparseMatrix<double>> > > res_map;
    if (param.simu_param.debug_mode){cout << "~~~~~~~~~~ DEBUG mode is activated in C++ forward model ~~~~~~~~~~" << endl;}

	// 2. Simulation
	// 2.1. DFN generation
	if (param.simu_param.simu_option==DFN_GENERATION){
	//	DFNVisu DFN_visu=DFNGeneration(param,domain);

	// Compute the backbone
		/*NetworkMeshes net_mesh1=DFN_visu.net_mesh,backbone;
		backbone=net_mesh1.return_backbone(param,EPSILON);
		DFN_visu.net_mesh=backbone;*/

		// Cut the domain
		/*NetworkMeshes net_mesh2=DFN_visu.net_mesh;
		//DFNCutting(net_mesh2,26.6);
		//DFNCutting(net_mesh2,25);
		DFN_visu.net_mesh=net_mesh2;*/

		// DFN writing
//		Results results(param.files_param.output_path+param.files_param.output_file,DFN_visu.net_mesh);
		// TODO remove display function
		//DisplayResults(argc,argv,DFN_visu);
	}

	else if (param.simu_param.simu_option==FOURIER_POTENTIAL_COMPUTATION){
		// Equivalent Porous Media
		if (param.simu_param.model_option==EPM){
		    if (param.simu_param.debug_mode){cout << "DEBUG: EPM" << endl;}
			// Read the source term positions
			map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);
			// Run the simulations
			if (source_term_positions.size()>0){
			    //FourierPotentialEPMSourceTerm(param,domain,singularity,bc_on_border,source_term_positions);
			    std::map<int, ublas_matrix> res_map_ublas = FourierPotentialEPMSourceTerm(param,domain,singularity,bc_on_border,source_term_positions,true);
			    res_map = UblasMapToEigenMap(res_map_ublas);
			 }
			else if (source_term_positions.size()==0){FourierPotentialEPM(param,domain,singularity);}
			else{cout << "WARNING in main (PERFORM.cpp): option not implemented" << endl;}
		}
		else if (param.simu_param.model_option==EPM_CS){
		    if (param.simu_param.debug_mode){cout << "DEBUG: EPM_CS" << endl;}
		    map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);
		    if (source_term_positions.size()>0){
		        res_map = eigen_cs::FourierPotentialEPMSourceTermCSEigen(param,domain,singularity,bc_on_border,source_term_positions,true);
		    }
		    else if (source_term_positions.size()==0){
		        cout << "WARNING in main (PERFORM.cpp): EPM_CS without source positions falls back to standard EPM workflow" << endl;
		        FourierPotentialEPM(param,domain,singularity);
		    }
		    else{cout << "WARNING in main (PERFORM.cpp): option not implemented" << endl;}
		}
		// Discrete-Dual-Porosity model for fractured porous domains
		else if (param.simu_param.model_option==DDP){
		    if (param.simu_param.debug_mode){cout << "DEBUG: DDP" << endl;}
			// Read the source term positions
			map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);

			//----------------------- BEGIN TEST -------------------------------------
			#ifdef LIB_TEST_SMALL_MATRIX
			// 1000 matrices of size 8x2
			cout << "TEST SMALL MATRIX 1000 matrices of size 8x2" << endl;
			for (int i=129;i<=1000;i++) {
				source_term_positions[i] = make_pair((double)rand()/RAND_MAX*13,(double)rand()/RAND_MAX*16);
			}
			param.simu_param.DDP_param.Nx=2;param.simu_param.DDP_param.Ny=2;
			#endif

			#ifdef LIB_TEST_BIG_MATRIX
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
				//FourierPotentialDDPSourceTerm(param,domain,singularity,source_term_positions);
				std::map<int, ublas_matrix> res_map_ublas = FourierPotentialDDPSourceTerm(param,domain,singularity,source_term_positions,true);
				res_map = UblasMapToEigenMap(res_map_ublas);
				//cout << "after FourierPotentialDDPSourceTerm" << endl;
			}
			else if (source_term_positions.size()==0){FourierPotentialDDP(param,domain,singularity,bc_on_border);}
			else{cout << "WARNING in main (PERFORM.cpp): option not implemented" << endl;}

		}
		else if (param.simu_param.model_option==DDP_EIGEN){
		    if (param.simu_param.debug_mode){cout << "DEBUG: DDP_EIGEN" << endl;}
		    map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);
		    if (source_term_positions.size()>0){
		        res_map = FourierPotentialDDPSourceTermEigen(param,domain,singularity,source_term_positions,true);
		    }
		    else if (source_term_positions.size()==0){FourierPotentialDDP(param,domain,singularity,bc_on_border);}
		    else{cout << "WARNING in main (PERFORM.cpp): option not implemented" << endl;}
		}
		else if (param.simu_param.model_option==DDP_CS){
		    if (param.simu_param.debug_mode){cout << "DEBUG: DDP_CS" << endl;}
		    map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);
		    if (source_term_positions.size()>0){
		        res_map = eigen_cs::FourierPotentialDDPSourceTermCSEigen(param,domain,singularity,bc_on_border,source_term_positions,true);
		    }
		    else if (source_term_positions.size()==0){cout << "WARNING in main (PERFORM.cpp): DDP_CS without source positions not implemented" << endl;}
		    else{cout << "WARNING in main (PERFORM.cpp): option not implemented" << endl;}
		}
		else{cout << "WARNING in main (PERFORM.cpp): solving method not implemented" << endl;}
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
			cout << "WARNING in main (PERFORM.cpp): EPM_CS is only used in Fourier source workflow, fallback to standard EPM for POTENTIAL_COMPUTATION" << endl;
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
		else{cout << "WARNING in main (PERFORM.cpp): solving method not implemented" << endl;}
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
		/*int nb_div=6,nb_seed=1; // nb_div:number of times that we divide the domain and nb_seed:number of generated domain
		ExchangeStudy(param,domain,nb_div,nb_seed);*/
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
	else{cout << "WARNING in main (PERFORM.cpp): simulation option not implemented" << endl;}

	t_end=clock();
	//cout << "CPU Time = " << (t_end-t_begin)/CLOCKS_PER_SEC << endl;

    /*// Test eigen matrix
    int n = 5;
    SparseMatrix<double> A(n,n);
    std::vector<T> coefficients;

	for (int i=0;i<n;i++){
	    for (int j=0;j<n;j++){
            	coefficients.push_back(T(i,j,10.));
                        }
                }
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    std::map<int,Eigen::SparseMatrix<double>> map_matrices;
    map_matrices[0] = A;
	return map_matrices;*/
	if (param.simu_param.debug_mode){
		cout << "DEBUG: check results at the end of perform.cpp" << endl;
		cout << "	res_map size = " << res_map.size() << endl;
		Eigen::SparseMatrix<double> A=res_map[0]; // if vector
		//Eigen::SparseMatrix<double> A=res_map[1]; // if map
		cout << "	res_map[0] = " << A << endl; // if vector
		//cout << "	res_map[1] = " << A << endl; // if map
		cout << "~~~~~~~~~~ End of DEBUG test ~~~~~~~~~~" << endl;
	}
	return res_map;
}

PYBIND11_MODULE(PERFORM, m) {
    m.doc() = "pybind11 electrical resistivity"; // optional module docstring
    m.def("perform", &perform, py::return_value_policy::reference);
}
