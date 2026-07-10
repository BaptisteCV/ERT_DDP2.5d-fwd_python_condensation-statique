/*
 * DDPModel.cpp
 *
 *  Created on: Apr 29, 2013
 *      Author: delphineroubinet
 */

#include "DDPModel.h"
#include "../DFNModel/DFNComputation.h"
#include "../EPMModel/EPMDefinition.h"
#include "../Utilitaries/Storage/UblasStructures.h"
#include "DDPDFNFunctions.h"
#include "ExchangeFunctions.h"
#include <boost/math/special_functions/bessel.hpp>
#include "../Utilitaries/Math/Operations.h"
#include <utility>
#include <omp.h>

using namespace std;

double DDPModel::phi_function(double sigma0,int i,int j){
	if (source_terms_epm.size()!=1){
		cout << "WARNING in phi_function (EPMDefinition.cpp): case not implemented" << endl;
		return-1;
	}
	double I=source_terms_epm.begin()->second;
	CgalPoint2D pt_current=CenterPosition(make_pair(i,j),std::make_pair(delta_x,delta_y));
	CgalPoint2D pt_term=CenterPosition(make_pair(source_terms_epm.begin()->first.first,source_terms_epm.begin()->first.second),make_pair(delta_x,delta_y));
	double xi=pt_current.x(),zj=pt_current.y();
	double x0=pt_term.x(),z0=Ny*delta_y;
	double a1=current_w*std::sqrt((xi-x0)*(xi-x0)+(zj-z0)*(zj-z0));
	return I/(2*PI*sigma0)*boost::math::cyl_bessel_k(0,a1);
}

void DDPModel::DefinitionLinearSystem(int index_center,pair<int,int> block_indices,string border){
	// 0. Variables
	double coeff1,coeff2,coeff3;
	int index_neigh=-1;
	CoefficientsLinearSystem(block_indices.first,block_indices.second,Nx,Ny,border,bc_def_epm,porous_cond,delta_x,delta_y,coeff1,coeff2,coeff3,index_neigh);
	// 2. Matrix and vector affectation
	matrix_syst(index_center,index_center)+=coeff1;
	if (index_neigh!=-1){matrix_syst(index_center,index_neigh+nb_nodes_fract)+=coeff2;}
	vector_syst(index_center)+=coeff3;
}


void DDPModel::DefinitionLinearSystemBCOnBorder(int index_center,pair<int,int> block_indices,string border){
	// 0. Variables
	double coeff1,coeff2;
	int index_neigh=-1;int i_block=block_indices.first,j_block=block_indices.second;
	// 1. Definition of the standard coefficients
	CoefficientsLinearSystemBCOnBorder(i_block,j_block,Nx,Ny,border,bc_def_epm,porous_cond,delta_x,delta_y,coeff1,coeff2,index_neigh);
	// 2. Matrix and vector affectation
	matrix_syst(index_center,index_center)+=coeff1;
	if (index_neigh!=-1){matrix_syst(index_center,index_neigh+nb_nodes_fract)+=coeff2;}
	// 3. Implementation for boundary conditions (nodes on the domain border)
	int index_border;bool bc_border=false;
	double coeff1_border,coeff2_border,coeff3_border;
	if (i_block==0&&border==LEFT_BORDER){
		map<std::string,BoundaryConditions> bc_cell=bc_def_epm[make_pair(i_block,j_block)];
		BoundCondCoeffsBCOnBorder(border,bc_cell[border],delta_x,delta_y,coeff1_border,coeff2_border,coeff3_border);
		index_border=Nx*Ny+j_block+nb_nodes_fract;
		bc_border=true;
	}
	else if (i_block==Nx-1&&border==RIGHT_BORDER){
		map<std::string,BoundaryConditions> bc_cell=bc_def_epm[make_pair(i_block,j_block)];
		BoundCondCoeffsBCOnBorder(border,bc_cell[border],delta_x,delta_y,coeff1_border,coeff2_border,coeff3_border);
		index_border=Nx*Ny+Ny+j_block+nb_nodes_fract;
		bc_border=true;
	}
	else if (j_block==0&&border==BOTTOM_BORDER){
		map<std::string,BoundaryConditions> bc_cell=bc_def_epm[make_pair(i_block,j_block)];
		BoundCondCoeffsBCOnBorder(border,bc_cell[border],delta_x,delta_y,coeff1_border,coeff2_border,coeff3_border);
		index_border=Nx*Ny+2*Ny+i_block+nb_nodes_fract;
		bc_border=true;
	}
	else if (j_block==Ny-1&&border==TOP_BORDER){
		map<std::string,BoundaryConditions> bc_cell=bc_def_epm[make_pair(i_block,j_block)];
		BoundCondCoeffsBCOnBorder(border,bc_cell[border],delta_x,delta_y,coeff1_border,coeff2_border,coeff3_border);
		index_border=Nx*Ny+2*Ny+Nx+i_block+nb_nodes_fract;
		bc_border=true;
	}
	if (bc_border){
		matrix_syst(index_center,index_border)+=coeff2;
		matrix_syst(index_border,index_border)=coeff1_border;
		matrix_syst(index_border,index_center)=coeff2_border;
		vector_syst(index_border)=coeff3_border;
	}
}

// Affectation of the coefficients in the DDP linear system corresponding to the porous domain definition with exchange with the DFN
void DDPModel::PorousDDPLinearSystem(vector<FractureMesh> subnet_meshes,int index_block,double alpha_block){
	// 0. Variables
	// Block indexes
	int index_center=index_block+nb_nodes_fract;
	pair<int,int> block_indices=return_indices(index_block,Ny);
	double coeff1=0,coeff2=0,coeff3=0,tot_coeff3=0;
	// 1. Regular finite volume method for a porous domain
	// Linear system definition for each border of the studied block
	DefinitionLinearSystem(index_center,block_indices,LEFT_BORDER);
	DefinitionLinearSystem(index_center,block_indices,RIGHT_BORDER);
	DefinitionLinearSystem(index_center,block_indices,BOTTOM_BORDER);
	DefinitionLinearSystem(index_center,block_indices,TOP_BORDER);
	// Source term implementation
	SourceTermsEPM::iterator it=source_terms_epm.find(block_indices);
	if (it!=source_terms_epm.end()){
		vector_syst(index_center)+=it->second;
	}
	// 2. Additional part for fracture-matrix interaction
	for (vector<FractureMesh>::iterator it1=subnet_meshes.begin();it1!=subnet_meshes.end();it1++){
		CumulatePotFractCoeff(it1->ReturnLength(),alpha_block/it1->conductivity,coeff1,coeff2,coeff3);
		matrix_syst(index_center,it1->p_ori.index)-=alpha_block*coeff1;
		matrix_syst(index_center,it1->p_tar.index)-=alpha_block*coeff2;
		tot_coeff3+=coeff3;
		matrix_syst(index_center,index_center)+=alpha_block*it1->ReturnLength();

	}
	matrix_syst(index_center,index_center)-=alpha_block*tot_coeff3;
}

// Affectation of the coefficients in the DDP linear system corresponding to the porous domain definition with exchange with the DFN
void DDPModel::PorousDDPLinearSystemUpdateFourier(vector<FractureMesh> subnet_meshes,int index_block,double alpha_block){
	// 0. Variables
	// Block indexes
	int index_center=index_block+nb_nodes_fract;
	double coeff1=0,coeff2=0,coeff3=0;
	// 1. Update related to the Fourier variable
	pair<int,int> block_indices=return_indices(index_block,Ny);
	matrix_syst(index_center,index_center)+=porous_cond(block_indices.first,block_indices.second)*current_w*current_w*delta_x*delta_y;
	// 2. Additional part for fracture-matrix interaction
	for (vector<FractureMesh>::iterator it1=subnet_meshes.begin();it1!=subnet_meshes.end();it1++){
		// New definition for k_fract
		CumulatePotFractCoeffFourier(it1->ReturnLength(),alpha_block/(it1->conductivity*it1->aperture),current_w,coeff1,coeff2,coeff3);
		matrix_syst(index_center,it1->p_ori.index)-=alpha_block*coeff1;
		matrix_syst(index_center,it1->p_tar.index)-=alpha_block*coeff2;
		matrix_syst(index_center,index_center)+=alpha_block*(it1->ReturnLength()-coeff3);
	}
}

// Affectation of the coefficients in the DDP linear system corresponding to the porous domain definition with exchange with the DFN
void DDPModel::PorousDDPLinearSystemFourier(vector<FractureMesh> subnet_meshes,int index_block,double alpha_block){
	// 0. Variables
	// Block indexes
	int index_center=index_block+nb_nodes_fract;
	//pair<int,int> block_indices=return_indices(index_block,Ny);
	double coeff1=0,coeff2=0,coeff3=0;

	// 1. Regular finite volume method for a porous domain
	DetermineLinearSystemFourier(index_center);

	// 2. Additional part for fracture-matrix interaction
	for (vector<FractureMesh>::iterator it1=subnet_meshes.begin();it1!=subnet_meshes.end();it1++){
		// New definition for k_fract
		CumulatePotFractCoeffFourier(it1->ReturnLength(),alpha_block/(it1->conductivity*it1->aperture),current_w,coeff1,coeff2,coeff3);
		matrix_syst(index_center,it1->p_ori.index)-=alpha_block*coeff1;
		matrix_syst(index_center,it1->p_tar.index)-=alpha_block*coeff2;
		matrix_syst(index_center,index_center)+=alpha_block*(it1->ReturnLength()-coeff3);
	}
}

void DDPModel::DetermineLinearSystemFourier(int index_block){
	// Variables
	int index_center=index_block-nb_nodes_fract;
	pair<int,int> block_indices=return_indices(index_center,Ny);
	// 1. Linear system definition
	// Coefficients affectation for each border
	DefinitionLinearSystem(index_block,block_indices,LEFT_BORDER);
	DefinitionLinearSystem(index_block,block_indices,RIGHT_BORDER);
	DefinitionLinearSystem(index_block,block_indices,BOTTOM_BORDER);
	DefinitionLinearSystem(index_block,block_indices,TOP_BORDER);

	// Fourier component
	double coeff=porous_cond(block_indices.first,block_indices.second)*current_w*current_w*delta_x*delta_y;
	matrix_syst(index_block,index_block)+=coeff;
}

// Affectation of coefficients corresponding to the boundary conditions of the DFN system
void DDPModel::BCDDPLinearSystem(){
	// 1. Boundary conditions related to DFN part
	for (BoundaryConditionsDFN::iterator it1=bc_map.begin();it1!=bc_map.end();it1++){
		// electric potential affectation
		if (it1->second.bc_type==DIRICHLET){
			// cancel previous values
			for (int i=0;i<matrix_syst.size2();i++){matrix_syst(it1->first,i)=0;}
			// new value affectation
			matrix_syst(it1->first,it1->first)=1;
			vector_syst(it1->first)=it1->second.bc_value;
		}
		// electric current affectation
		else if (it1->second.bc_type==NEUMANN){
			// new value affectation
			vector_syst(it1->first)=it1->second.bc_value;
		}
		else{cout << "WARNING in DFNComputation.cpp: boundary condition not implemented" << endl;}
	}
}

void DDPModel::BCDDPLinearSystemFourier(){
	NodeFracturesMap node_fract_map=ReturnNodeFracturesMap(subnetwork_map);
	ElecPropMap elec_prop_map;int index_block;
	for (BoundaryConditionsDFN::iterator it1=bc_map.begin();it1!=bc_map.end();it1++){
		// electric potential affectation
		if (it1->second.bc_type==DIRICHLET){
			// cancel previous values by looping only on the values related to this node
			// fracture nodes
			elec_prop_map=node_fract_map[it1->first];
			for (ElecPropMap::iterator it2=elec_prop_map.begin(); it2!=elec_prop_map.end(); it2++){
				matrix_syst(it1->first,it2->first)=0;
			}
			// matrix nodes
			index_block=ReturnBlockIndex(subnetwork_map,it1->first);
			matrix_syst(it1->first,index_block+nb_nodes_fract)=0;
			// new value affectation
			matrix_syst(it1->first,it1->first)=1;
			vector_syst(it1->first)=it1->second.bc_value;
		}
		// electric current affectation
		else if (it1->second.bc_type==NEUMANN){
			// new value affectation
			vector_syst(it1->first)=it1->second.bc_value;
		}
		else if (it1->second.bc_type==MIXED_BC){
			// addition of the flux normal to the borders and expressed as a function of phi_i
			index_block=it1->second.block_index;
			double alpha=ReturnAlphaBCMixed(it1->first,index_block);
			matrix_syst(it1->first,it1->first)+=2*subnetwork_map[index_block].return_aperture_from_node(it1->first)*subnetwork_map[index_block].return_conductivity_from_node(it1->first)*alpha/(2+alpha*delta_x);
		}
		else{cout << "WARNING in BCDDPLinearSystemFourier (DDPModel.cpp): boundary condition not implemented" << endl;}
	}
}

// Return the coefficient alpha for the mixed boundary conditions
double DDPModel::ReturnAlphaBCMixed(int index_node,int index_block){
	if (source_terms_epm.size()!=1){cout << "WARNING in ReturnAlphaBCMixed (DDPModel.cpp): number of source term not implemented" << endl;}
	CgalPoint2D node_position=subnetwork_map[index_block].ReturnNodeCoord(index_node),source_term=CenterPosition(source_terms_epm.begin()->first,make_pair(delta_x,delta_y));
	string border;
	if (subnetwork_map[index_block].border_map.find(index_node)!=subnetwork_map[index_block].border_map.end()){border=subnetwork_map[index_block].border_map[index_node];}
	else{cout << "WARNING in ReturnAlphaBCMixed (DDPModel.cpp): border not defined" << endl;}
	double alpha,r=distance_2D(source_term,node_position);;
	if (border==LEFT_BORDER){
		// Definition of alpha
		alpha=current_w*boost::math::cyl_bessel_k(1,current_w*r)/boost::math::cyl_bessel_k(0,current_w*r);
		return alpha*abs(node_position.x()-source_term.x())/r;
	}
	else if (border==RIGHT_BORDER){
		// Definition of alpha
		alpha=current_w*boost::math::cyl_bessel_k(1,current_w*r)/boost::math::cyl_bessel_k(0,current_w*r);
		return alpha*abs(node_position.x()-source_term.x())/r;
	}
	else if (border==BOTTOM_BORDER){
		// Definition of alpha
		alpha=current_w*boost::math::cyl_bessel_k(1,current_w*r)/boost::math::cyl_bessel_k(0,current_w*r);
		return alpha*abs(node_position.y()-source_term.y())/r;
	}
	else{cout << "WARNING in ReturnAlphaBCMixed (DDPModel.cpp): border not implemented" << endl; return -1;}
}

// Affect the corresponding coefficient in the linear system for a given subnetwork (network at the block-scale)
void DDPModel::DFNDDPLinearSystem(vector<FractureMesh> subnet_meshes,int index_block,double alpha_block){
	// 0. Variables
	double coeff,L,a1,a2,a3,k_fract;
	int index1,index2;
	// 1. Loop on the fractures for coefficient affectation and data updating
	for (vector<FractureMesh>::iterator it1=subnet_meshes.begin();it1!=subnet_meshes.end();it1++){
		// 1.1. Coefficients definition and affectation
		coeff=it1->aperture*it1->conductivity;
		L=it1->ReturnLength();
		index1=it1->p_ori.index,index2=it1->p_tar.index;
		k_fract=alpha_block/it1->conductivity;
		// fracture value affectation for the node1 (x=0 or phi_i)
		ReturnDerivativeCoeffFull(0,L,k_fract,a1,a2,a3);
		matrix_syst(index1,index1)+=coeff*a1;
		matrix_syst(index1,index2)+=coeff*a2;
		matrix_syst(index1,index_block+nb_nodes_fract)+=coeff*a3;
		// fracture value affectation for the node2 (x=L or phi_(i+1))
		ReturnDerivativeCoeffFull(L,L,k_fract,a1,a2,a3);
		matrix_syst(index2,index1)-=coeff*a1;
		matrix_syst(index2,index2)-=coeff*a2;
		matrix_syst(index2,index_block+nb_nodes_fract)-=coeff*a3;
	}
}

// Fourier term added
void DDPModel::DFNDDPLinearSystemFourier(vector<FractureMesh> subnet_meshes,int index_block,double alpha_block){
	// 0. Variables
	double coeff,L,a1,a2,a3,b1,b2,b3,k_fract;
	int index1,index2;
	// 1. Loop on the fractures for coefficient affectation and data updating
	for (vector<FractureMesh>::iterator it1=subnet_meshes.begin();it1!=subnet_meshes.end();it1++){
		// 1.1. Coefficients definition and affectation
		coeff=it1->aperture*it1->conductivity;
		L=it1->ReturnLength();
		index1=it1->p_ori.index,index2=it1->p_tar.index;
		// New definition for k_fract coefficient
		k_fract=alpha_block/(it1->conductivity*it1->aperture);

		// 1.2. Definition of the linear system
		// fracture value affectation for the node1 (x=0 or phi_i)
		ReturnDerivativeCoeffFullFourier(0,L,k_fract,current_w,a1,a2,a3);
		matrix_syst(index1,index1)+=coeff*it1->aperture*current_w*current_w;	// integration of phi over the intersection (should be done only once)
		matrix_syst(index1,index1)-=coeff*a1;
		matrix_syst(index1,index2)-=coeff*a2;
		matrix_syst(index1,index_block+nb_nodes_fract)-=coeff*a3;
		// fracture value affectation for the node2 (x=L or phi_(i+1))
		ReturnDerivativeCoeffFullFourier(L,L,k_fract,current_w,a1,a2,a3);
		matrix_syst(index2,index1)+=coeff*a1;
		matrix_syst(index2,index2)+=coeff*a2;
		matrix_syst(index2,index_block+nb_nodes_fract)+=coeff*a3;
	}
}

// Distribute the potential in fracture and matrix potential from the global solution
void DDPModel::PotentialFractureMatrixDistribution(ublas_vector potential,ublas_vector& DFNPotential,ublas_matrix& PorousPotential){
	// 1. DFN potential
	for (int i=0;(unsigned)i<DFNPotential.size();i++){
		DFNPotential(i)=potential(i);
	}
	// 2. Porous potential
	for (int i=0;i<PorousPotential.size1();i++){
		for (int j=0;j<PorousPotential.size2();j++){
			PorousPotential(i,j)=potential(return_index(i,j,Ny)+nb_nodes_fract);
		}
	}
}

// Distribute the potential in fracture and matrix potential from the global solution
void DDPModel::PotentialFractureMatrixDistributionFourier(ublas_vector potential,ublas_vector& DFNPotential,ublas_matrix& PorousPotential){
	// 1. DFN potential
	DFNPotential=InitializationVector(nb_nodes_fract,nb_nodes_fract);
	for (int i=0;(unsigned)i<DFNPotential.size();i++){
		DFNPotential(i)=potential(i);
	}
	// 2. Porous potential
	// Standard nodes
	for (int i=0;i<Nx;i++){
		for (int j=0;j<Ny;j++){
			PorousPotential(i,j)=potential(return_index(i,j,Ny)+nb_nodes_fract);
		}
	}
}

// Electric potential computation for DDP model
void DDPModel::ComputePotential(ublas_vector& DFNPotential,ublas_matrix& PorousPotential){
	// 0. Variables and parameters definition
	// Variables
	double alpha_block,porous_cond_block;pair<int,int> block_indices;
	AlphaBlock=InitializationMatrix(Nx,Ny);
	// Potential matrices initialization
	PorousPotential=InitializationMatrix(Nx,Ny);
	DFNPotential=InitializationVector(nb_nodes_fract,nb_nodes_fract);
	// Linear system initialization
	int nb_node_matrix=Nx*Ny;
	matrix_syst=InitializationMatrix(nb_nodes_fract+nb_node_matrix,nb_nodes_fract+nb_node_matrix);
	vector_syst=InitializationVector(nb_nodes_fract+nb_node_matrix);
	// 1. Build linear system
	// Loop on the blocks (or subnetwork)
	for (SubNetworkMap::iterator it1=subnetwork_map.begin();it1!=subnetwork_map.end();it1++){
		// 1.1. Block properties definition
		block_indices=return_indices(it1->first,Ny);
		porous_cond_block=porous_cond(block_indices.first,block_indices.second);
		alpha_block=ReturnAlphaBlock(porous_cond_block,delta_x,delta_y,it1->second);
		AlphaBlock(block_indices.first,block_indices.second)=alpha_block;
		// 1.2. DFN part affectation and fracture data computation
		DFNDDPLinearSystem(it1->second.meshes,it1->first,alpha_block);
		// 1.3. Porous part affectation
		PorousDDPLinearSystem(it1->second.meshes,it1->first,alpha_block);
		// 1.4. Boundary conditions affectation
		BCDDPLinearSystem();
	}
	// 2. Solve linear system
	ublas_vector potential=LinearSystemSolving(matrix_syst,vector_syst);
	cout << "Number of unknowns = " << potential.size() << endl;
	// 3. Distribute solution in DFN and porous potential
	PotentialFractureMatrixDistribution(potential,DFNPotential,PorousPotential);
}

// Determine the coefficients of the matrix that are independent on the Fourier variable, boundary conditions, and source terms (no bc_border option)
void DDPModel::DetermineBasicLinearSystemDDP(){
	//double t_basic_total = omp_get_wtime();
	// 1. Variables and initialization
	int node_number=Nx*Ny;		// standard nodes + additional nodes for boundary conditions on the domain border
	//double t_alloc = omp_get_wtime();
	matrix_syst=InitializationMatrix(nb_nodes_fract+node_number,nb_nodes_fract+node_number);vector_syst=InitializationVector(nb_nodes_fract+node_number);
	//std::cout << "[PROFILE][DetermineBasicLinearSystemDDP] InitializationMatrix+Vector (ublas alloc): "
	//          << omp_get_wtime() - t_alloc << " sec" << std::endl;
	int index_center,index_neigh;double coeff;
	// 2. Determination of the coefficients related to the inside meshes
	for (int i=0;i<Nx;i++){
		for (int j=0;j<Ny;j++){
			index_center=nb_nodes_fract+return_index(i,j,Ny);
			// Coefficients affectation for the left border
			if (i!=0){
				index_neigh=nb_nodes_fract+return_index(i-1,j,Ny);
				coeff=Harmonic_Mean(porous_cond(i,j),porous_cond(i-1,j))*delta_y/delta_x;
				matrix_syst(index_center,index_center)+=coeff;
				matrix_syst(index_center,index_neigh)-=coeff;
			}
			// Coefficients affectation for the right border
			if (i!=Nx-1){
				index_neigh=nb_nodes_fract+return_index(i+1,j,Ny);
				coeff=Harmonic_Mean(porous_cond(i,j),porous_cond(i+1,j))*delta_y/delta_x;
				matrix_syst(index_center,index_center)+=coeff;
				matrix_syst(index_center,index_neigh)-=coeff;
			}
			// Coefficients affectation for the bottom border
			if (j!=0){
				index_neigh=nb_nodes_fract+return_index(i,j-1,Ny);
				coeff=Harmonic_Mean(porous_cond(i,j),porous_cond(i,j-1))*delta_x/delta_y;
				matrix_syst(index_center,index_center)+=coeff;
				matrix_syst(index_center,index_neigh)-=coeff;
			}
			// Coefficients affectation for the bottom border
			if (j!=Ny-1){
				index_neigh=nb_nodes_fract+return_index(i,j+1,Ny);
				coeff=Harmonic_Mean(porous_cond(i,j),porous_cond(i,j+1))*delta_x/delta_y;
				matrix_syst(index_center,index_center)+=coeff;
				matrix_syst(index_center,index_neigh)-=coeff;
			}
		}
	}
	//std::cout << "[PROFILE][DetermineBasicLinearSystemDDP] total: "
	//          << omp_get_wtime() - t_basic_total << " sec"
	//          << "  (nb_nodes=" << nb_nodes_fract + Nx*Ny << ")" << std::endl;
}

// Electric potential computation for DDP model
void DDPModel::UpdateLinearSystemFourier(){
	// 0. Variables and initialization
	// Variables
	double alpha_block,porous_cond_block;pair<int,int> block_indices;
	AlphaBlock=InitializationMatrix(Nx,Ny);
	// 1. Update the coefficients which depend on the Fourier variable
	// Loop on the blocks (or subnetwork)
	for (SubNetworkMap::iterator it1=subnetwork_map.begin();it1!=subnetwork_map.end();it1++){
		// 1.1. Block properties definition
		block_indices=return_indices(it1->first,Ny);
		porous_cond_block=porous_cond(block_indices.first,block_indices.second);
		alpha_block=ReturnAlphaBlock(porous_cond_block,delta_x,delta_y,it1->second);
		AlphaBlock(block_indices.first,block_indices.second)=alpha_block;
		// 1.2. DFN part affectation and fracture data computation
		DFNDDPLinearSystemFourier(it1->second.meshes,it1->first,alpha_block);
		// 1.3. Porous part affectation
		PorousDDPLinearSystemUpdateFourier(it1->second.meshes,it1->first,alpha_block);
	}
	// 2. Update the coefficients related to the boundary conditions
	// Boundary conditions related to the DFN part
	BCDDPLinearSystemFourier();
	// Boundary conditions related to the porous part
	UpdateLinearSystemBC_EPM(nb_nodes_fract);
	// 3. Update the coefficients related to the source terms
	SourceTermFourier();
}

// Porous part correspond to the rows and columns after nb_nodes_fract
void DDPModel::UpdateLinearSystemBC_EPM(int nb_nodes_fract){
	// 1. Variables
	map<string,BoundaryConditions> bc_cell;
	int i,j,index_center;double coeff1,coeff3;
	// 2. Update the coefficients of the linear system that are related to the boundary conditions
	// Left border
	i=0;
	for (j=0;j<Ny;j++){
		index_center=nb_nodes_fract+return_index(i,j,Ny);bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffs(LEFT_BORDER,bc_cell[LEFT_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);
		// 2. Matrix and vector affectation
		matrix_syst(index_center,index_center)+=coeff1;
		vector_syst(index_center)+=coeff3;
	}
	// Right border
	i=Nx-1;
	for (j=0;j<Ny;j++){
		index_center=nb_nodes_fract+return_index(i,j,Ny);bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffs(RIGHT_BORDER,bc_cell[RIGHT_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);
		matrix_syst(index_center,index_center)+=coeff1;
		vector_syst(index_center)+=coeff3;
	}
	// Bottom border
	j=0;
	for (i=0;i<Nx;i++){
		index_center=nb_nodes_fract+return_index(i,j,Ny);bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffs(BOTTOM_BORDER,bc_cell[BOTTOM_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);
		matrix_syst(index_center,index_center)+=coeff1;
		vector_syst(index_center)+=coeff3;
	}
	// Bottom border
	j=Ny-1;
	for (i=0;i<Nx;i++){
		index_center=nb_nodes_fract+return_index(i,j,Ny);bc_cell=bc_def_epm[make_pair(i,j)];
		BoundCondCoeffs(TOP_BORDER,bc_cell[TOP_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);
		matrix_syst(index_center,index_center)+=coeff1;
		vector_syst(index_center)+=coeff3;
	}
}

// Solve the liner system and define the fracture and matrix potential
void DDPModel::EvaluatePotential(ublas_vector& DFNPotential,ublas_matrix& PorousPotential){

	// Potential matrices initialization
	PorousPotential=InitializationMatrix(Nx,Ny);
	DFNPotential=InitializationVector(nb_nodes_fract,nb_nodes_fract);

	ublas_vector potential=LinearSystemSolving(matrix_syst,vector_syst);
	PotentialFractureMatrixDistribution(potential,DFNPotential,PorousPotential);
}

// Define the matrix and vector of the linear system
void DDPModel::DefineLinearSystemFourier(ublas_matrix& PorousPotential,bool bc_on_border){
	// 0. Variables and parameters definition
	// Variables
	double alpha_block,porous_cond_block;pair<int,int> block_indices;
	AlphaBlock=InitializationMatrix(Nx,Ny);
	// Linear system and potential matrices initialization
	int nb_node_matrix=0;
	if (bc_on_border){
		nb_node_matrix=Nx*Ny+2*Nx+2*Ny;
		PorousPotential=InitializationMatrix(Nx,Ny+4);
	}
	else{
		nb_node_matrix=Nx*Ny;
		PorousPotential=InitializationMatrix(Nx,Ny);
	}

	matrix_syst=InitializationMatrix(nb_nodes_fract+nb_node_matrix,nb_nodes_fract+nb_node_matrix);
	vector_syst=InitializationVector(nb_nodes_fract+nb_node_matrix);
	// 1. Build linear system
	// Loop on the blocks (or subnetwork)
	for (SubNetworkMap::iterator it1=subnetwork_map.begin();it1!=subnetwork_map.end();it1++){
		// 1.1. Block properties definition
		block_indices=return_indices(it1->first,Ny);
		porous_cond_block=porous_cond(block_indices.first,block_indices.second);
		alpha_block=ReturnAlphaBlock(porous_cond_block,delta_x,delta_y,it1->second);
		AlphaBlock(block_indices.first,block_indices.second)=alpha_block;
		// 1.2. DFN part affectation and fracture data computation
		DFNDDPLinearSystemFourier(it1->second.meshes,it1->first,alpha_block);
		// 1.3. Porous part affectation
		PorousDDPLinearSystemFourier(it1->second.meshes,it1->first,alpha_block);
	}
	if (singularity){
		cout << "WARNING in DDPModel (DDPModel.cpp): This method is not implemented in DDP" << endl;
	}
}

void DDPModel::SourceTermFourier(){
	if (!singularity){
		int index_block;
			for (SourceTermsEPM::iterator it=source_terms_epm.begin();it!=source_terms_epm.end();it++){
				index_block=return_index(it->first.first,it->first.second,Ny)+nb_nodes_fract;
				vector_syst(index_block)+=it->second*0.5;
			}
	}
}

// Electric potential computation for DDP model
void DDPModel::ComputePotentialFourier(ublas_vector& DFNPotential,ublas_matrix& PorousPotential,bool bc_on_border){

	// 1.1. Define the matrix and vector of the linear system
	cout << "DefineLinearSystemFourier" << endl;
	DefineLinearSystemFourier(PorousPotential,bc_on_border);

	// 1.2. Boundary conditions affectation
	cout << "BCDDPLinearSystemFourier" << endl;
	BCDDPLinearSystemFourier();

	// 1.3. Source terms
	cout << "SourceTermFourier" << endl;
	SourceTermFourier();

	// 2. Solve linear system
	cout << "LinearSystemSolving" << endl;
	ublas_vector potential=LinearSystemSolving(matrix_syst,vector_syst);

	// 3. Distribute solution in DFN and porous potential
	cout << "PotentialFractureMatrixDistributionFourier" << endl;
	PotentialFractureMatrixDistributionFourier(potential,DFNPotential,PorousPotential);
}

// DDP model definition
DDPModel::DDPModel(Parameters param,Domain domain,NetworkMeshes net_mesh,bool fourier,bool singularity_){
	current_w=0;
	// 0. Conditions of simulation
	if (net_mesh.domain.domain_size_x()!=param.domain_param.Lx||net_mesh.domain.domain_size_y()!=param.domain_param.Ly){
		cout << "WARNING in DDPModel (DDPModel.cpp): domain lengths are not properly defined" << endl;
		return;
	}
	// 1. Porous domain definition
	Nx=param.simu_param.DDP_param.Nx,Ny=param.simu_param.DDP_param.Ny;
	delta_x=param.simu_param.DDP_param.ReturnDeltax(param.domain_param.Lx);
	delta_y=param.simu_param.DDP_param.ReturnDeltay(param.domain_param.Ly);
	if (param.domain_param.mat_cond!=-1){
		porous_cond=InitializationMatrix(Nx,Ny,param.domain_param.mat_cond);
	}
	else{
		porous_cond=param.domain_param.mat_cond_dist;
		if (porous_cond.size1()!=Nx||porous_cond.size2()!=Ny){
			cout << "WARNING in EPMSystem (EPMDefinition.cpp): matrix sizes are different" << endl;
		}
	}

	singularity=singularity_;
	// 2. Boundary conditions are defined in this function for non-Fourier simulations
	if (!fourier){
		SourceTermsEPM source_terms_epm_;
		bc_def_epm=BoundaryConditionsDefinition(domain,Nx,Ny,param.simu_param.bc_map,param.simu_param.source_terms,source_terms_epm_,fourier,current_w);
		source_terms_epm=source_terms_epm_;
		// 2. DFN domain definition
		// DFN boundary conditions
		bc_map=ReturnBoundCondDFN(net_mesh.domain,param.simu_param.bc_map,net_mesh.border_map,net_mesh.nodes_map);
	}
	// 3. Division of the fractured domain in subnetworks
	//double t_ddpnet = omp_get_wtime();
	subnetwork_map=DDPNetworkMeshes(net_mesh,Nx,Ny);
	//std::cout << "[PROFILE][DDPModel ctor] DDPNetworkMeshes: "
	//          << omp_get_wtime() - t_ddpnet << " sec" << std::endl;
	//double t_nbnodes = omp_get_wtime();
	nb_nodes_fract=ReturnNbNodes(subnetwork_map);
	//std::cout << "[PROFILE][DDPModel ctor] ReturnNbNodes: "
	//          << omp_get_wtime() - t_nbnodes << " sec" << std::endl;
}

void DDPModel::DDPModelFourierBC(Parameters param,Domain domain,Domain network_domain,BordersMap network_border_map,NodesMap network_nodes_map){
	bool fourier=true;
	current_w=param.simu_param.DDP_param.w_current;
	SourceTermsEPM source_terms_epm_;
	bc_def_epm=BoundaryConditionsDefinition(domain,Nx,Ny,param.simu_param.bc_map,param.simu_param.source_terms,source_terms_epm_,fourier,current_w);
	source_terms_epm=source_terms_epm_;
	// 2. DFN domain definition
	// DFN boundary conditions
	if (fourier){
		bc_map=ReturnBoundCondDFNFourier(network_domain,param.simu_param.bc_map,network_border_map,network_nodes_map,bc_def_epm,delta_x,delta_y);
	}
	else{
		bc_map=ReturnBoundCondDFN(network_domain,param.simu_param.bc_map,network_border_map,network_nodes_map);
	}
	if (singularity){
		SourceTermsEPM source_terms_epm_sing_;
		bc_def_epm_sing=BoundaryConditionsDefinition(domain,Nx,Ny,param.simu_param.bc_map_sing,param.simu_param.source_terms,source_terms_epm_sing_,fourier,current_w);
		source_terms_epm_sing=source_terms_epm;
	}
}

void print_subnetwork(SubNetworkMap subnetwork_map){
	for (SubNetworkMap::iterator it=subnetwork_map.begin();it!=subnetwork_map.end();it++){
		it->second.print_DFN();
	}
}
