/*
 * EPMDefinition.cpp
 *
 *  Created on: Apr 25, 2013
 *      Author: delphineroubinet
 */

#include "EPMDefinition.h"
#include "../DDPModel/DDPDFNFunctions.h"
#include "../Utilitaries/Storage/Structures.h"
#include <boost/math/special_functions/bessel.hpp>
#include "../Utilitaries/Math/Operations.h"
#include <utility>

using namespace std;

// return the coefficients required for the linear system where coeff1 is the coeff for the cell (i,j), coeff2 for its neighbor and coeff3 for the second member
void EPMSystem::DefinitionLinearSystemStandard(int i,int j,string border){
	// 0. Variables
	double coeff1,coeff2,coeff3;
	int index_center=return_index(i,j,Ny),index_neigh=-1;
	CoefficientsLinearSystem(i,j,Nx,Ny,border,bc_def_epm,porous_cond,delta_x,delta_y,coeff1,coeff2,coeff3,index_neigh);
	// 2. Matrix and vector affectation
	matrix_syst(index_center,index_center)+=coeff1;
	if (index_neigh!=-1){matrix_syst(index_center,index_neigh)+=coeff2;}
	vector_syst(index_center)+=coeff3;
}


// return the coefficients required for the linear system where coeff1 is the coeff for the cell (i,j), coeff2 for its neighbor and coeff3 for the second member
void EPMSystem::DefinitionLinearSystemBCOnBorder(int i,int j,string border){
	// 0. Variables
	int index_center=return_index(i,j,Ny),index_neigh=-1;
	double coeff1,coeff2;
	// 1. Evaluation of the standard coefficient
	CoefficientsLinearSystemBCOnBorder(i,j,Nx,Ny,border,bc_def_epm,porous_cond,delta_x,delta_y,coeff1,coeff2,index_neigh);
	// 2. Matrix and vector affectation
	matrix_syst(index_center,index_center)+=coeff1;
	if (index_neigh!=-1){matrix_syst(index_center,index_neigh)+=coeff2;}
	// 3. Implementation for boundary conditions (nodes on the domain border)
	int index_border;bool bc_border=false;
	double coeff1_border,coeff2_border,coeff3_border;
	if (i==0&&border==LEFT_BORDER){
		map<std::string,BoundaryConditions> bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffsBCOnBorder(border,bc_cell[border],delta_x,delta_y,coeff1_border,coeff2_border,coeff3_border);
		index_border=Nx*Ny+j;
		bc_border=true;
	}
	else if (i==Nx-1&&border==RIGHT_BORDER){
		map<std::string,BoundaryConditions> bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffsBCOnBorder(border,bc_cell[border],delta_x,delta_y,coeff1_border,coeff2_border,coeff3_border);
		index_border=Nx*Ny+Ny+j;
		bc_border=true;
	}
	else if (j==0&&border==BOTTOM_BORDER){
		map<std::string,BoundaryConditions> bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffsBCOnBorder(border,bc_cell[border],delta_x,delta_y,coeff1_border,coeff2_border,coeff3_border);
		index_border=Nx*Ny+2*Ny+i;
		bc_border=true;
	}
	else if (j==Ny-1&&border==TOP_BORDER){
		map<std::string,BoundaryConditions> bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffsBCOnBorder(border,bc_cell[border],delta_x,delta_y,coeff1_border,coeff2_border,coeff3_border);
		index_border=Nx*Ny+2*Ny+Nx+i;
		bc_border=true;
	}
	if (bc_border){
		matrix_syst(index_center,index_border)+=coeff2;
		matrix_syst(index_border,index_border)=coeff1_border;
		matrix_syst(index_border,index_center)=coeff2_border;
		vector_syst(index_border)=coeff3_border;
	}
}

void EPMSystem::DefinitionLinearSystem(int i,int j,string border,bool bc_on_border){
	if (bc_on_border){
		DefinitionLinearSystemBCOnBorder(i,j,border);
	}
	else{DefinitionLinearSystemStandard(i,j,border);}
}

// Solve the electric potential
ublas_matrix EPMSystem::PotentialComputation(){
	// Linear system solution
	ublas_vector solution=LinearSystemSolving(matrix_syst,vector_syst);

	// 3. Conversion from vector to matrix
	ublas_matrix potential(Nx,Ny);
	for (int i=0;i<Nx;i++){
		for (int j=0;j<Ny;j++){
			potential(i,j)=solution(return_index(i,j,Ny));
		}
	}
	this->potential=potential;
	return potential;
}

void EPMSystem::print_bc_def(){
	for (BoundaryConditionsEPM::iterator it1=bc_def_epm.begin();it1!=bc_def_epm.end();it1++){
		print(it1->first);
		for (map<std::string,BoundaryConditions>::iterator it2=it1->second.begin();it2!=it1->second.end();it2++){
			cout << it2->first;
			it2->second.print();
		}
	}
}

void EPMSystem::SourceTermsImplementation(bool fourier){
	for (SourceTermsEPM::iterator it=source_terms_epm.begin();it!=source_terms_epm.end();it++){
		int index_cell=return_index(it->first.first,it->first.second,Ny);
		if (!fourier){
			vector_syst(index_cell)+=it->second;
		}
		else{
			vector_syst(index_cell)+=it->second*0.5;
		}
	}
}

double EPMSystem::phi_function(double sigma0,int i,int j){
	if (source_terms_epm.size()!=1){
		cout << "WARNING in phi_function (EPMDefinition.cpp): case not implemented" << endl;
		return-1;
	}
	double I=source_terms_epm.begin()->second;
	CgalPoint2D pt_current=CenterPosition(make_pair(i,j),make_pair(delta_x,delta_y));
	CgalPoint2D pt_term=CenterPosition(make_pair(source_terms_epm.begin()->first.first,source_terms_epm.begin()->first.second),make_pair(delta_x,delta_y));
	double xi=pt_current.x(),zj=pt_current.y();
	double x0=pt_term.x(),z0=Ny*delta_y;
	double a1=w*std::sqrt((xi-x0)*(xi-x0)+(zj-z0)*(zj-z0));
	return I/(2*PI*sigma0)*boost::math::cyl_bessel_k(0,a1);
}

EPMSystem::EPMSystem(Parameters param,NetworkMeshes net_mesh){
	// 0. Parameters and variables
	// Simulation parameters
	Nx=param.simu_param.EPM_param.Nx,Ny=param.simu_param.EPM_param.Ny;
	delta_x=param.domain_param.Lx/Nx,delta_y=param.domain_param.Ly/Ny;
	porous_cond=InitializationMatrix(Nx,Ny,param.domain_param.mat_cond);
	EPMConductivityDefinition(param.domain_param.mat_cond,pointcpp<double>(0,0),net_mesh);
}

EPMSystem::EPMSystem(Parameters param){
	// 0. Parameters and variables
	// Simulation parameters
	Nx=param.simu_param.EPM_param.Nx,Ny=param.simu_param.EPM_param.Ny;
	delta_x=param.domain_param.Lx/Nx,delta_y=param.domain_param.Ly/Ny;
	if (param.domain_param.mat_cond!=-1){
		porous_cond=InitializationMatrix(Nx,Ny,param.domain_param.mat_cond);
	}
	else{
		porous_cond=param.domain_param.mat_cond_dist;
		if (porous_cond.size1()!=Nx||porous_cond.size2()!=Ny){
			cout << "WARNING in EPMSystem (EPMDefinition.cpp): matrix sizes are different" << endl;
		}
	}
	w=0.0;
}

// Standard linear system for porous media
EPMSystem::EPMSystem(Parameters param,Domain domain,bool fourier_syst,bool fourier_bc,bool singularity,bool bc_on_border,bool fourier_update,bool source_term){
	// 0. Parameters and variables
	// Simulation parameters
	Nx=param.simu_param.EPM_param.Nx,Ny=param.simu_param.EPM_param.Ny;
	delta_x=param.domain_param.Lx/Nx,delta_y=param.domain_param.Ly/Ny;
	if (param.domain_param.mat_cond!=-1){
		porous_cond=InitializationMatrix(Nx,Ny,param.domain_param.mat_cond);
	}
	else{
		porous_cond=param.domain_param.mat_cond_dist;
		if (porous_cond.size1()!=Nx||porous_cond.size2()!=Ny){
			cout << "WARNING in EPMSystem (EPMDefinition.cpp): matrix sizes are different" << endl;
		}
	}
	w=param.simu_param.EPM_param.w_current;
	SourceTermsEPM source_terms_epm_;
	cout << "test1a" << endl;
	bc_def_epm=BoundaryConditionsDefinition(domain,Nx,Ny,param.simu_param.bc_map,param.simu_param.source_terms,source_terms_epm_,fourier_bc,w);
	cout << "test1b" << endl;
	source_terms_epm=source_terms_epm_;

	if (singularity){
		SourceTermsEPM source_terms_epm_sing_;
		bc_def_epm_sing=BoundaryConditionsDefinition(domain,Nx,Ny,param.simu_param.bc_map_sing,param.simu_param.source_terms,source_terms_epm_sing_,fourier_bc,w);
		source_terms_epm_sing=source_terms_epm;
	}
	DetermineLinearSystem(fourier_syst,singularity,bc_on_border,fourier_update,source_term);
}

void EPMSystem::SingularityMethodBCOnBorder(){
	// Definition of sigma0
	double sigma_cst;
	ublas_matrix sigma0=ReturnSigma0(porous_cond,sigma_cst);
	int index_center,index_neigh;
	double coeff1,coeff2,final_coeff1;

	for (int i=0;i<Nx;i++){
		for (int j=0;j<Ny;j++){
			final_coeff1=0;
			index_center=return_index(i,j,Ny);
			CoefficientsLinearSystemBCOnBorder(i,j,Nx,Ny,LEFT_BORDER,bc_def_epm_sing,sigma0,delta_x,delta_y,coeff1,coeff2,index_neigh);
			if (i!=0){
				vector_syst(index_center)+=-(-coeff2+matrix_syst(index_center,index_neigh))*phi_function(sigma_cst,i-1,j);
			}
			final_coeff1+=coeff1;
			CoefficientsLinearSystemBCOnBorder(i,j,Nx,Ny,RIGHT_BORDER,bc_def_epm_sing,sigma0,delta_x,delta_y,coeff1,coeff2,index_neigh);
			if (i!=Nx-1){
				vector_syst(index_center)+=-(-coeff2+matrix_syst(index_center,index_neigh))*phi_function(sigma_cst,i+1,j);
			}
			final_coeff1+=coeff1;
			CoefficientsLinearSystemBCOnBorder(i,j,Nx,Ny,BOTTOM_BORDER,bc_def_epm_sing,sigma0,delta_x,delta_y,coeff1,coeff2,index_neigh);
			if (j!=0){
				vector_syst(index_center)+=-(-coeff2+matrix_syst(index_center,index_neigh))*phi_function(sigma_cst,i,j-1);
			}
			final_coeff1+=coeff1;
			CoefficientsLinearSystemBCOnBorder(i,j,Nx,Ny,TOP_BORDER,bc_def_epm_sing,sigma0,delta_x,delta_y,coeff1,coeff2,index_neigh);
			if (j!=Ny-1){
				vector_syst(index_center)+=-(-coeff2+matrix_syst(index_center,index_neigh))*phi_function(sigma_cst,i,j+1);
			}
			final_coeff1+=coeff1;
			vector_syst(index_center)+=(final_coeff1+sigma0(i,j)*w*w*delta_x*delta_y-matrix_syst(index_center,index_center))*phi_function(sigma_cst,i,j);
		}
	}
}


void EPMSystem::SingularityMethod(){
	// Definition of sigma0
	double sigma_cst;
	ublas_matrix sigma0=ReturnSigma0(porous_cond,sigma_cst);
	int index_center,index_neigh;
	double coeff1,coeff2,coeff3,final_coeff1;

	for (int i=0;i<Nx;i++){
		for (int j=0;j<Ny;j++){
			final_coeff1=0;
			index_center=return_index(i,j,Ny);
			CoefficientsLinearSystem(i,j,Nx,Ny,LEFT_BORDER,bc_def_epm_sing,sigma0,delta_x,delta_y,coeff1,coeff2,coeff3,index_neigh);
			if (i!=0){
				vector_syst(index_center)+=-(-coeff2+matrix_syst(index_center,index_neigh))*phi_function(sigma_cst,i-1,j);
			}
			final_coeff1+=coeff1;
			CoefficientsLinearSystem(i,j,Nx,Ny,RIGHT_BORDER,bc_def_epm_sing,sigma0,delta_x,delta_y,coeff1,coeff2,coeff3,index_neigh);
			if (i!=Nx-1){
				vector_syst(index_center)+=-(-coeff2+matrix_syst(index_center,index_neigh))*phi_function(sigma_cst,i+1,j);
			}
			final_coeff1+=coeff1;
			CoefficientsLinearSystem(i,j,Nx,Ny,BOTTOM_BORDER,bc_def_epm_sing,sigma0,delta_x,delta_y,coeff1,coeff2,coeff3,index_neigh);
			if (j!=0){
				vector_syst(index_center)+=-(-coeff2+matrix_syst(index_center,index_neigh))*phi_function(sigma_cst,i,j-1);
			}
			final_coeff1+=coeff1;
			CoefficientsLinearSystem(i,j,Nx,Ny,TOP_BORDER,bc_def_epm_sing,sigma0,delta_x,delta_y,coeff1,coeff2,coeff3,index_neigh);
			if (j!=Ny-1){
				vector_syst(index_center)+=-(-coeff2+matrix_syst(index_center,index_neigh))*phi_function(sigma_cst,i,j+1);
			}
			final_coeff1+=coeff1;
			vector_syst(index_center)+=(final_coeff1+sigma0(i,j)*w*w*delta_x*delta_y-matrix_syst(index_center,index_center))*phi_function(sigma_cst,i,j);
		}
	}
}

// Determine the coefficients of the matrix that are independent on the Fourier variable, boundary conditions, and source terms (no bc_border option)
void EPMSystem::DetermineBasicLinearSystem(){
	// 1. Variables and initialization
	int node_number=Nx*Ny;		// standard nodes + additional nodes for boundary conditions on the domain border
	matrix_syst=InitializationMatrix(node_number,node_number);vector_syst=InitializationVector(node_number);
	int index_center,index_neigh;double coeff;
	// 2. Determination of the coefficients related to the inside meshes
	for (int i=0;i<Nx;i++){
		for (int j=0;j<Ny;j++){
			index_center=return_index(i,j,Ny);
			// Coefficients affectation for the left border
			if (i!=0){
				index_neigh=return_index(i-1,j,Ny);
				coeff=Harmonic_Mean(porous_cond(i,j),porous_cond(i-1,j))*delta_y/delta_x;
				matrix_syst(index_center,index_center)+=coeff;
				matrix_syst(index_center,index_neigh)-=coeff;
			}
			// Coefficients affectation for the right border
			if (i!=Nx-1){
				index_neigh=return_index(i+1,j,Ny);
				coeff=Harmonic_Mean(porous_cond(i,j),porous_cond(i+1,j))*delta_y/delta_x;
				matrix_syst(index_center,index_center)+=coeff;
				matrix_syst(index_center,index_neigh)-=coeff;
			}
			// Coefficients affectation for the bottom border
			if (j!=0){
				index_neigh=return_index(i,j-1,Ny);
				coeff=Harmonic_Mean(porous_cond(i,j),porous_cond(i,j-1))*delta_x/delta_y;
				matrix_syst(index_center,index_center)+=coeff;
				matrix_syst(index_center,index_neigh)-=coeff;
			}
			// Coefficients affectation for the bottom border
			if (j!=Ny-1){
				index_neigh=return_index(i,j+1,Ny);
				coeff=Harmonic_Mean(porous_cond(i,j),porous_cond(i,j+1))*delta_x/delta_y;
				matrix_syst(index_center,index_center)+=coeff;
				matrix_syst(index_center,index_neigh)-=coeff;
			}
		}
	}
}

// Add the fourier term in the linear system (for which the basic definition has been alreday done)
void EPMSystem::UpdateLinearSystemFourier(){
	int index_center;
	for (int i=0;i<Nx;i++){
		for (int j=0;j<Ny;j++){
			index_center=return_index(i,j,Ny);
			matrix_syst(index_center,index_center)+=porous_cond(i,j)*w*w*delta_x*delta_y;
		}
	}
}

void EPMSystem::UpdateLinearSystemBoundaryConditionsAndSourceTerms(bool fourier){
	// 1. Variables
	map<string,BoundaryConditions> bc_cell;
	int i,j,index_center;double coeff1,coeff3;
	// 2. Update the coefficients of the linear system that are related to the boundary conditions
	// Left border
	i=0;
	for (j=0;j<Ny;j++){
		index_center=return_index(i,j,Ny);bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffs(LEFT_BORDER,bc_cell[LEFT_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);
		// 2. Matrix and vector affectation
		matrix_syst(index_center,index_center)+=coeff1;
		vector_syst(index_center)+=coeff3;
	}
	// Right border
	i=Nx-1;
	for (j=0;j<Ny;j++){
		index_center=return_index(i,j,Ny);bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffs(RIGHT_BORDER,bc_cell[RIGHT_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);
		matrix_syst(index_center,index_center)+=coeff1;
		vector_syst(index_center)+=coeff3;
	}
	// Bottom border
	j=0;
	for (i=0;i<Nx;i++){
		index_center=return_index(i,j,Ny);bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffs(BOTTOM_BORDER,bc_cell[BOTTOM_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);
		matrix_syst(index_center,index_center)+=coeff1;
		vector_syst(index_center)+=coeff3;
	}
	// Bottom border
	j=Ny-1;
	for (i=0;i<Nx;i++){
		index_center=return_index(i,j,Ny);bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffs(TOP_BORDER,bc_cell[TOP_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);
		matrix_syst(index_center,index_center)+=coeff1;
		vector_syst(index_center)+=coeff3;
	}
	// 3. Update the coefficients related to the source terms
	SourceTermsImplementation(fourier);
}

void EPMSystem::DetermineLinearSystem(bool fourier,bool singularity,bool bc_on_border,bool fourier_update,bool source_term){
	// Variables
	int node_number=0;
	if (bc_on_border){
		node_number=Nx*Ny+2*Nx+2*Ny;		// standard nodes + additional nodes for boundary conditions on the domain border
	}
	else{node_number=Nx*Ny;}

	// The whole system needs to be determined
	if (!fourier_update){
		matrix_syst=InitializationMatrix(node_number,node_number);vector_syst=InitializationVector(node_number);
		// 1. Linear system definition
		for (int i=0;i<Nx;i++){
			for (int j=0;j<Ny;j++){
				// Coefficients affectation for each border
				DefinitionLinearSystem(i,j,LEFT_BORDER,bc_on_border);
				DefinitionLinearSystem(i,j,RIGHT_BORDER,bc_on_border);
				DefinitionLinearSystem(i,j,BOTTOM_BORDER,bc_on_border);
				DefinitionLinearSystem(i,j,TOP_BORDER,bc_on_border);
				if (fourier){
					int index_center=return_index(i,j,Ny);
					double coeff=porous_cond(i,j)*w*w*delta_x*delta_y;
					matrix_syst(index_center,index_center)+=coeff;
				}
			}
		}
	}
	// the system is already defined and only needs to be update
	else if (fourier_update&&fourier){
		int index_center;
		// 1. Linear system definition
		for (int i=0;i<Nx;i++){
			for (int j=0;j<Ny;j++){
				index_center=return_index(i,j,Ny);
				matrix_syst(index_center,index_center)+=porous_cond(i,j)*w*w*delta_x*delta_y;
			}
		}
	}

	// Linear system modification to remove singularities
	if (singularity){
		SingularityMethod();
	}
	// 2. Definition of the source terms
	if (!singularity&&source_term){
		SourceTermsImplementation(fourier);
	}
}

ublas_matrix ReturnSigma0(ublas_matrix porous_cond,double& sigma_cst){

	sigma_cst=AveragedMatrix(porous_cond);
	ublas_matrix sigma0=InitializationMatrix(porous_cond.size1(),porous_cond.size2(),sigma_cst);
	return sigma0;

	/*ublas_matrix extract_matrix=ExtractMatrix(porous_cond,(porous_cond.size1()-1)*0.375,(porous_cond.size1()-1)*0.625,(porous_cond.size2()-1)*0.75,porous_cond.size2()-1);
	sigma_cst=AveragedMatrix(extract_matrix);
	ublas_matrix sigma0=InitializationMatrix(extract_matrix.size1(),extract_matrix.size2(),sigma_cst);
	return sigma0;*/

	//ublas_matrix extract_matrix=ExtractMatrix(porous_cond,(porous_cond.size1()-1)*0.25,(porous_cond.size1()-1)*0.75,(porous_cond.size2()-1)*0.75,porous_cond.size2()-1);
	/*double value;
	if (!ConstantValue(extract_matrix,value)){
		cout << "WARNING in ReturnSigma0 (EPMDefinition.cpp): the extracted matrix is not constant" << endl;
	}
	else{
		sigma0=InitializationMatrix(porous_cond.size1(),porous_cond.size2(),value);
		sigma_cst=value;
	}
	return sigma0;*/
}

EPMSystem::EPMSystem(Parameters param,Domain domain,NetworkMeshes net_mesh,string option){
	// 0. Parameters and variables
	// Simulation parameters
	Nx=param.simu_param.EPM_param.Nx,Ny=param.simu_param.EPM_param.Ny;
	delta_x=param.domain_param.Lx/Nx,delta_y=param.domain_param.Ly/Ny;
	EPMConductivityDefinition(param.domain_param.mat_cond,domain.min_pt,net_mesh,option);
	SourceTermsEPM source_terms_epm_;
	bc_def_epm=BoundaryConditionsDefinition(domain,Nx,Ny,param.simu_param.bc_map,param.simu_param.source_terms,source_terms_epm_);
	source_terms_epm=source_terms_epm_;
	// Variables
	cout << "Begin InitializationMatrix" << endl;
	matrix_syst=InitializationMatrix(Nx*Ny,Nx*Ny);
	vector_syst=InitializationVector(Nx*Ny);
	// 1. Linear system definition
	for (int i=0;i<Nx;i++){
		for (int j=0;j<Ny;j++){
			// Coefficients affectation for each border
			DefinitionLinearSystem(i,j,LEFT_BORDER);
			DefinitionLinearSystem(i,j,RIGHT_BORDER);
			DefinitionLinearSystem(i,j,BOTTOM_BORDER);
			DefinitionLinearSystem(i,j,TOP_BORDER);
		}
	}
	cout << "End linear system" << endl;
	// 2. Definition of the source terms
	SourceTermsImplementation();
}

// return the matrix of conductivity taking into the porous domain conductivity (matrix) and the fractures (by projection)
void EPMSystem::EPMConductivityDefinition(double porous_cond_scalar,pointcpp<double> min_domain,NetworkMeshes net_mesh,string option){

	// 1. Conductivity of the matrix part
	porous_cond=InitializationMatrix(Nx,Ny,porous_cond_scalar);

	// 2. Conductivity of the fracture part by projection of the fractures
	if (net_mesh.meshes.size()!=0 ){
		// 2.1. Decomposition of the fractures into small fractures included in meshes
		if (option==INTERSECTION){
			FracturesIntersectionGrid(Nx,Ny,net_mesh);
		}

		// 2.2. Determine the meshes containing a bit of fracture
		CgalPoint2D center;pair<int,int> indices;
		for (vector<FractureMesh>::iterator it=net_mesh.meshes.begin();it!=net_mesh.meshes.end();it++){
			center=it->ReturnCenter();
			indices=return_indices((double)center.x(),(double)center.y(),delta_x,delta_y,Nx,Ny,min_domain.i,min_domain.j);
			porous_cond(indices.first,indices.second)=porous_cond(indices.first,indices.second)+it->aperture*(it->conductivity-porous_cond(indices.first,indices.second))/delta_x;
			if (!exist_in_vector(fract_indices,indices)){
				fract_indices.push_back(indices);
			}
		}
	}
}

map<int,EPMSystem> EPMSystem::DefineSubEPMFract(map<int,Domain> new_domains,SubNetworkMap network_map,double cond_mat,string option){
	map<int,EPMSystem> new_epm_models;

	// Definition of the subdomains
	for (map<int,Domain>::iterator it=new_domains.begin();it!=new_domains.end();it++){
		EPMSystem new_epm;
		// Discretization parameters
		new_epm.delta_x=this->delta_x;
		new_epm.delta_y=this->delta_y;
		new_epm.Nx=it->second.domain_size_x()/new_epm.delta_x;
		new_epm.Ny=it->second.domain_size_y()/new_epm.delta_y;
		// Conductivity and fracture position
		new_epm.EPMConductivityDefinition(cond_mat,it->second.min_pt,network_map[it->first],option);
		// Affectation in the map
		new_epm_models[it->first]=new_epm;
	}

	// Extraction of the conductivity and potential
	//map<int,ublas_matrix> new_cond_map=ExtractMatrixMap(this->porous_cond,new_domains,new_epm_models);
	map<int,ublas_matrix> new_potential=ExtractMatrixMap(this->potential,new_domains,new_epm_models);

	// Affectation
	for (map<int,EPMSystem>::iterator it=new_epm_models.begin();it!=new_epm_models.end();it++){
		it->second.potential=new_potential[it->first];
	}
	return new_epm_models;
}

map<int,EPMSystem> EPMSystem::DefineSubEPM(map<int,Domain> new_domains,double cond_mat){
	map<int,EPMSystem> new_epm_models;

	// Definition of the subdomains
	for (map<int,Domain>::iterator it=new_domains.begin();it!=new_domains.end();it++){
		EPMSystem new_epm;
		// Discretization parameters
		new_epm.delta_x=this->delta_x;
		new_epm.delta_y=this->delta_y;
		new_epm.Nx=(floor)(it->second.domain_size_x()/new_epm.delta_x);
		new_epm.Ny=(floor)(it->second.domain_size_y()/new_epm.delta_y);
		// Conductivity and fracture position
		new_epm.EPMConductivityDefinition(cond_mat,it->second.min_pt);
		// Affectation in the map
		new_epm_models[it->first]=new_epm;
	}

	// Extraction of the conductivity and potential
	map<int,ublas_matrix> new_cond_map=ExtractMatrixMap(this->porous_cond,new_domains,new_epm_models);
	map<int,ublas_matrix> new_potential=ExtractMatrixMap(this->potential,new_domains,new_epm_models);

	// Affectation
	for (map<int,EPMSystem>::iterator it=new_epm_models.begin();it!=new_epm_models.end();it++){
		it->second.potential=new_potential[it->first];
	}
	return new_epm_models;
}

map<int,ublas_matrix> ExtractMatrixMap(ublas_matrix init_matrix,map<int,Domain> domain_map,map<int,EPMSystem> epm_model){

	map<int,ublas_matrix> extracted_mat;

	pair<int,int> ext1,ext2;
	for (map<int,Domain>::iterator it=domain_map.begin();it!=domain_map.end();it++){
		ext1=return_indices(it->second.min_pt.i+0.5*epm_model[it->first].delta_x,it->second.min_pt.j+0.5*epm_model[it->first].delta_y,
				epm_model[it->first].delta_x,epm_model[it->first].delta_y,init_matrix.size1(),init_matrix.size2(),it->second.min_pt.i,it->second.min_pt.j);
		ext2=return_indices(it->second.max_pt.i-0.5*epm_model[it->first].delta_x,it->second.max_pt.j-0.5*epm_model[it->first].delta_y,
				epm_model[it->first].delta_x,epm_model[it->first].delta_y,init_matrix.size1(),init_matrix.size2(),it->second.min_pt.i,it->second.min_pt.j);
		/*ext1=return_indices(it->second.min_pt.i+0.5*epm_model[it->first].delta_x,it->second.min_pt.j+0.5*epm_model[it->first].delta_y,
				epm_model[it->first].delta_x,epm_model[it->first].delta_y);
		ext2=return_indices(it->second.max_pt.i-0.5*epm_model[it->first].delta_x,it->second.max_pt.j-0.5*epm_model[it->first].delta_y,
				epm_model[it->first].delta_x,epm_model[it->first].delta_y);*/
		extracted_mat[it->first]=ExtractMatrix(init_matrix,ext1.first,ext2.first,ext1.second,ext2.second);
	}

	return extracted_mat;
}

map<int,double> AveragedMatrices(map<int,EPMSystem> epm_models){
	map<int,double> averaged_potential;
	for (map<int,EPMSystem>::iterator it=epm_models.begin();it!=epm_models.end();it++){
		averaged_potential[it->first]=AveragedMatrix(it->second.potential);
	}

	return averaged_potential;
}

