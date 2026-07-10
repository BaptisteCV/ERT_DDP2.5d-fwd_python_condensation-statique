/*
 * EPMComputation.cpp
 *
 *  Created on: Apr 25, 2013
 *      Author: delphineroubinet
 */

#include "EPMComputation.h"
#include "../Utilitaries/Constantes.h"
#include "../Utilitaries/Math/Operations.h"
#include "../DFNModel/NetworkMeshes.h"
#include <boost/math/special_functions/bessel.hpp>


using namespace std;


double ReturnFVMFactor(string border,double delta_x,double delta_y){
	if (border==LEFT_BORDER||border==RIGHT_BORDER){return delta_y/delta_x;}
	else if (border==BOTTOM_BORDER||border==TOP_BORDER){return delta_x/delta_y;}
	else{cout << "WARNING in ReturnRegularCoeff (EPMDefinition.cpp): border not defined" << endl;}
	return -1;
}

double ReturnFVMRegularCoeff(string border,double current_cond,double neigh_cond,double delta_x,double delta_y){
	double value=ReturnFVMFactor(border,delta_x,delta_y);
	return Harmonic_Mean(current_cond,neigh_cond)*value;
}


void BoundCondCoeffs(string border,BoundaryConditions bound_cond,double current_cond,double delta_x,double delta_y,double& coeff1,double& coeff3){
	coeff1=0,coeff3=0;
	if (bound_cond.bc_type==DIRICHLET){
		double value=2*ReturnFVMFactor(border,delta_x,delta_y);
		coeff1=current_cond*value;
		coeff3=coeff1*bound_cond.bc_value;
	}
	else if (bound_cond.bc_type==NEUMANN){coeff3=ReturnFVMFactor(border,delta_x,delta_y)*current_cond*bound_cond.bc_value;}
	else if (bound_cond.bc_type==MIXED_BC){	// in this case, bc_value is the distance from the source term to the studied boundary cell
		// Mixed conditions for Victor : alpha*cos(theta)*phi+beta*der_phi=gamma with beta=1, gamma=0, alpha*cos(theta) defined in the boundary conditions definition (Besselk...)
		coeff1=2*current_cond*bound_cond.bc_value;
		if (border==LEFT_BORDER||border==RIGHT_BORDER){coeff1=coeff1*delta_y/(2+bound_cond.bc_value*delta_x);}
		else if (border==BOTTOM_BORDER||border==TOP_BORDER){coeff1=coeff1*delta_x/(2+bound_cond.bc_value*delta_y);}
		else{cout << "WARNING in ReturnRegularCoeff (EPMDefinition.cpp): border not defined" << endl;}
	}
	else{cout << "WARNING in ReturnRegularCoeff (EPMDefinition.cpp): boundary condition not implemented" << endl;}
}


void BoundCondCoeffsBCOnBorder(string border,BoundaryConditions bound_cond,double delta_x,double delta_y,double& coeff1,double& coeff2,double& coeff3){
	coeff1=0,coeff2=0,coeff3=0;
	if (bound_cond.bc_type==DIRICHLET){
		coeff1=1;
		coeff3=bound_cond.bc_value;
	}
	else if (bound_cond.bc_type==NEUMANN){
		if (border==LEFT_BORDER||border==RIGHT_BORDER){
			coeff1=2/delta_x;
		}
		else if (border==TOP_BORDER||border==BOTTOM_BORDER){
			coeff1=2/delta_y;
		}
		else{cout << "WARNING in ReturnRegularCoeff (EPMDefinition.cpp): border not defined" << endl;}
		coeff2=-coeff1;
		coeff3=bound_cond.bc_value;
	}
	else if (bound_cond.bc_type==MIXED_BC){	// in this case, bc_value is the distance from the source term to the studied boundary cell
		// Mixed conditions for Victor : alpha*cos(theta)*phi+beta*der_phi=gamma with beta=1, gamma=0, alpha*cos(theta) defined in the boundary conditions definition (Besselk...)
		if (border==LEFT_BORDER||border==RIGHT_BORDER){coeff2=-2/delta_x;}
		else if (border==BOTTOM_BORDER||border==TOP_BORDER){coeff2=-2/delta_y;}
		else{cout << "WARNING in ReturnRegularCoeff (EPMDefinition.cpp): border not defined" << endl;}
		coeff1=bound_cond.bc_value-coeff2;
	}
	else{cout << "WARNING in ReturnRegularCoeff (EPMDefinition.cpp): boundary condition not implemented" << endl;}
}


void CoefficientsLinearSystem(int i,int j,int Nx,int Ny,string border,BoundaryConditionsEPM bc_def_epm,ublas_matrix porous_cond,double delta_x,double delta_y,double& coeff1,double& coeff2,double& coeff3,int& index_neigh){
	coeff1=0,coeff2=0,coeff3=0;
	map<std::string,BoundaryConditions> bc_cell=bc_def_epm[make_pair(i,j)];
	// 1. Neighbor and coefficients definition
	if (border==LEFT_BORDER){
		if (i==0){BoundCondCoeffs(border,bc_cell[LEFT_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);}
		else{
			coeff1=ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i-1,j),delta_x,delta_y);
			coeff2=-coeff1;
			index_neigh=return_index(i-1,j,Ny);
		}
	}
	else if (border==RIGHT_BORDER){
		if (i==Nx-1){BoundCondCoeffs(border,bc_cell[RIGHT_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);}
		else{
			coeff1=ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i+1,j),delta_x,delta_y);
			coeff2=-coeff1;
			index_neigh=return_index(i+1,j,Ny);
		}
	}
	else if (border==BOTTOM_BORDER){
		if (j==0){BoundCondCoeffs(border,bc_cell[BOTTOM_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);}
		else{
			coeff1=ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i,j-1),delta_x,delta_y);
			coeff2=-coeff1;
			index_neigh=return_index(i,j-1,Ny);
		}
	}
	else if (border==TOP_BORDER){
		if (j==Ny-1){BoundCondCoeffs(border,bc_cell[TOP_BORDER],porous_cond(i,j),delta_x,delta_y,coeff1,coeff3);}
		else{
			coeff1=ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i,j+1),delta_x,delta_y);
			coeff2=-coeff1;
			index_neigh=return_index(i,j+1,Ny);
		}
	}
	else{cout << "WARNING in CoefficientsLinearSystem (EPMDefinition.cpp): border not defined" << endl;}
}


void CoefficientsLinearSystemBCOnBorder(int i,int j,int Nx,int Ny,string border,BoundaryConditionsEPM bc_def_epm,ublas_matrix porous_cond,double delta_x,double delta_y,double& coeff1,double& coeff2,int& index_neigh){
	coeff1=0,coeff2=0;
	map<std::string,BoundaryConditions> bc_cell=bc_def_epm[make_pair(i,j)];
	// 1. Neighbor and coefficients definition
	if (border==LEFT_BORDER){
		if (i==0){
			coeff1=2*ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i,j),delta_x,delta_y);
			coeff2=-coeff1;
		}
		else{
			coeff1=ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i-1,j),delta_x,delta_y);
			coeff2=-coeff1;
			index_neigh=return_index(i-1,j,Ny);
		}
	}
	else if (border==RIGHT_BORDER){
		if (i==Nx-1){
			coeff1=2*ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i,j),delta_x,delta_y);
			coeff2=-coeff1;
		}
		else{
			coeff1=ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i+1,j),delta_x,delta_y);
			coeff2=-coeff1;
			index_neigh=return_index(i+1,j,Ny);
		}
	}
	else if (border==BOTTOM_BORDER){
		if (j==0){
			coeff1=2*ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i,j),delta_x,delta_y);
			coeff2=-coeff1;
		}
		else{
			coeff1=ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i,j-1),delta_x,delta_y);
			coeff2=-coeff1;
			index_neigh=return_index(i,j-1,Ny);
		}
	}
	else if (border==TOP_BORDER){
		if (j==Ny-1){
			coeff1=2*ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i,j),delta_x,delta_y);
			coeff2=-coeff1;
		}
		else{
			coeff1=ReturnFVMRegularCoeff(border,porous_cond(i,j),porous_cond(i,j+1),delta_x,delta_y);
			coeff2=-coeff1;
			index_neigh=return_index(i,j+1,Ny);
		}
	}
	else{cout << "WARNING in CoefficientsLinearSystem (EPMDefinition.cpp): border not defined" << endl;}
}

// Return the boundary condition value for mixed conditions = n.r/r^2
// where n is the outward normal of the studied border and r the vector from the point source to the considered border position
double ReturnMixedBoundCondValue1(SourceTermsDef source_terms,int i_cell,int j_cell,int Nx,int Ny,Domain domain,string border){
	// 0. Parameters
	double delta_x=domain.ReturnDeltaX(Nx),delta_y=domain.ReturnDeltaY(Ny);
	// 1. Determine for the source and studied cell location
	// source term position
	if (source_terms.size()!=1){cout << "WARNING in ReturnMixedBoundCondValue (EPMDefinition.cpp): case not implemented" << endl;return -1;}

	cout << "g" << endl;

	double x_source=source_terms.begin()->x_position,y_source=source_terms.begin()->y_position;
	CGAL::Point_2<K_CGAL> pt_source(x_source,y_source);
	// studied cell position
	CgalPoint2D cell_center=CenterPosition(make_pair(i_cell,j_cell),make_pair(delta_x,delta_y));
	double x_cell=cell_center.x(),y_cell=cell_center.y();
	// 2. Vectors definition
	// Unit vector definition and correction of the cell position (from the center of the cell to the border)
	double x_normal,y_normal;
	if (border==LEFT_BORDER){
		x_cell-=0.5*delta_x;
		x_normal=-1;y_normal=0;
	}
	else if (border==RIGHT_BORDER){
		x_cell+=0.5*delta_x;
		x_normal=1;y_normal=0;
	}
	else if (border==BOTTOM_BORDER){
		y_cell-=0.5*delta_y;
		x_normal=0;y_normal=-1;
	}
	else if (border==TOP_BORDER){cout << "WARNING in ReturnMixedBoundCondValue (EPMDefinition.cpp): case not implemented" << endl;return -1;}
	CGAL::Point_2<K_CGAL> pt1_normal(x_cell,y_cell),pt2_normal(x_cell+x_normal,y_cell+y_normal);
	CGAL::Vector_2<K_CGAL> normal_vect(pt1_normal,pt2_normal);
	// Vector from the source to the studied position
	CGAL::Vector_2<K_CGAL> source_vect(pt_source,pt1_normal);
	// 3. Computation of the scalar product
	double scalar_product=source_vect*normal_vect;
	// 4. Return n.r/r^2
	return scalar_product/source_vect.squared_length();
}

// Old boundary conditions
double ReturnMixedBoundCondValue2(SourceTermsDef source_terms,int i_cell,int j_cell,int Nx,int Ny,Domain domain,string border,double w){
	// Definition of the studied boundary cell
	double delta_x=domain.ReturnDeltaX(Nx),delta_y=domain.ReturnDeltaY(Ny);
	CgalPoint2D cell_center=CenterPosition(make_pair(i_cell,j_cell),make_pair(delta_x,delta_y));
	//double x_cell=cell_center.x(),y_cell=cell_center.y();
	// Definition of the position of the source term
	if (source_terms.size()!=1){cout << "WARNING in ReturnMixedBoundCondValue (EPMDefinition.cpp): case not implemented" << endl;return -1;}
	double x_source=source_terms.begin()->x_position,y_source=source_terms.begin()->y_position;
	CgalPoint2D pt_source(x_source,y_source);
	double r,alpha;
	if (border==LEFT_BORDER){
		cell_center=CgalPoint2D(cell_center.x()-0.5*delta_x,cell_center.y());
		// Definition of the distance r
		r=distance_2D(pt_source,cell_center);
		// Definition of alpha
		alpha=w*boost::math::cyl_bessel_k(1,w*r)/boost::math::cyl_bessel_k(0,w*r);
		return alpha*abs(cell_center.x()-x_source)/r;
	}
	else if (border==RIGHT_BORDER){
		cell_center=CgalPoint2D(cell_center.x()+0.5*delta_x,cell_center.y());
		// Definition of the distance r
		r=distance_2D(pt_source,cell_center);
		// Definition of alpha
		alpha=w*boost::math::cyl_bessel_k(1,w*r)/boost::math::cyl_bessel_k(0,w*r);
		return alpha*abs(cell_center.x()-x_source)/r;
	}
	else if (border==BOTTOM_BORDER){
		cell_center=CgalPoint2D(cell_center.x(),cell_center.y()+0.5*delta_y);
		// Definition of the distance r
		r=distance_2D(pt_source,cell_center);
		// Definition of alpha
		alpha=w*boost::math::cyl_bessel_k(1,w*r)/boost::math::cyl_bessel_k(0,w*r);
		return alpha*abs(cell_center.y()-y_source)/r;
	}
	cout << "WARNING in ReturnMixedBoundCondValue (EPMDefinition.cpp): case not implemented" << endl;
	return -1;
}

// Return the boundary conditions
BoundaryConditions ReturnBoundaryConditions(Domain domain,string border,int i,int j,int Nx,int Ny,BoundaryConditionsDef bc_map,SourceTermsDef source_terms,bool fourier,double w){
	// regular Dirichlet and Neumann conditions
	if (bc_map[border].bc_type==DIRICHLET||bc_map[border].bc_type==NEUMANN){
		return bc_map[border];
	}
	// gradient conditions
	else if (bc_map[border].bc_type==DIR_GRADIENT){
		string border1,border2;
		ReturnOppositeBorders(border,border1,border2);
		if (bc_map[border1].bc_type!=DIRICHLET||bc_map[border2].bc_type!=DIRICHLET){
			cout << "WARNING in ReturnBoundaryConditions (EPMDefinition.cpp): option not defined" << endl;
		}
		// for the position k along the border "border" composed of Nk cells
		int k,Nk;
		if (border==LEFT_BORDER||border==RIGHT_BORDER){k=j;Nk=Ny;}
		if (border==TOP_BORDER||border==BOTTOM_BORDER){k=i;Nk=Nx;}
		double bc_value=(bc_map[border2].bc_value-bc_map[border1].bc_value)*(k+0.5)/(double)Nk+bc_map[border1].bc_value;
		return BoundaryConditions(DIRICHLET,bc_value);
	}
	else if (bc_map[border].bc_type==MIXED_BC){
		double bc_value=0;
		if (!fourier){
			bc_value=ReturnMixedBoundCondValue1(source_terms,i,j,Nx,Ny,domain,border);
		}
		else{
			bc_value=ReturnMixedBoundCondValue2(source_terms,i,j,Nx,Ny,domain,border,w);
		}
		// if bc_value = nan -> due to large values of wr, leading to K0(wr)=0 and K1(wr)=0, corresponding to positions that are far enough to observe phi(r)=0
		if (std::isnan(bc_value)){return BoundaryConditions(DIRICHLET,0.0);}
		else {return BoundaryConditions(MIXED_BC,bc_value);}
	}
	else{cout << "WARNING in ReturnBoundaryConditions (EPMDefinition.cpp): boundary conditions not defined" << endl;}
	return BoundaryConditions();
}

// define boundary conditions for the EPM model
BoundaryConditionsEPM BoundaryConditionsDefinition(Domain domain,int Nx,int Ny,BoundaryConditionsDef bc_map,SourceTermsDef source_terms,SourceTermsEPM& source_terms_epm,bool fourier,double w){
	BoundaryConditionsEPM bc_def_epm;
	// Left and right borders
	for (int j=0;j<Ny;j++){
		// Left border
		bc_def_epm[make_pair(0,j)].insert(make_pair(LEFT_BORDER,ReturnBoundaryConditions(domain,LEFT_BORDER,0,j,Nx,Ny,bc_map,source_terms,fourier,w)));
		// Right border
		bc_def_epm[make_pair(Nx-1,j)].insert(make_pair(RIGHT_BORDER,ReturnBoundaryConditions(domain,RIGHT_BORDER,Nx-1,j,Nx,Ny,bc_map,source_terms,fourier,w)));
	}
	// Top and bottom borders
	for (int i=0;i<Nx;i++){
		// Bottom border
		bc_def_epm[make_pair(i,0)].insert(make_pair(BOTTOM_BORDER,ReturnBoundaryConditions(domain,BOTTOM_BORDER,i,0,Nx,Ny,bc_map,source_terms,fourier,w)));
		// Top border
		bc_def_epm[make_pair(i,Ny-1)].insert(make_pair(TOP_BORDER,ReturnBoundaryConditions(domain,TOP_BORDER,i,Ny-1,Nx,Ny,bc_map,source_terms,fourier,w)));
	}
	// Source terms
	for (SourceTermsDef::iterator it1=source_terms.begin();it1!=source_terms.end();it1++){
		pair<int,int> indices=return_indices(it1->x_position,it1->y_position,domain.ReturnDeltaX(Nx),domain.ReturnDeltaY(Ny),Nx,Ny);
		source_terms_epm[indices]=it1->source_value;
	}
	return bc_def_epm;
}

// Electric current computation
void ElectricCurrentComputation(Parameters param,ublas_matrix potential,ublas_matrix mat_cond,ublas_matrix & electric_current_x,ublas_matrix & electric_current_y){
	// 1. Parameters
	int Nx=param.simu_param.EPM_param.Nx,Ny=param.simu_param.EPM_param.Ny;
	double Lx=param.domain_param.Lx,Ly=param.domain_param.Ly;
	double delta_x=param.simu_param.EPM_param.ReturnDeltax(Lx),delta_y=param.simu_param.EPM_param.ReturnDeltay(Ly);
	ublas_matrix elec_current_x(Nx-1,Ny),elec_current_y(Nx,Ny-1);
	for (int i=0;i<Nx-1;i++){
		for (int j=0;j<Ny-1;j++){
			elec_current_x(i,j)=-Harmonic_Mean(mat_cond(i,j),mat_cond(i+1,j))*(potential(i+1,j)-potential(i,j))*delta_y/delta_x;
			elec_current_y(i,j)=-Harmonic_Mean(mat_cond(i,j),mat_cond(i,j+1))*(potential(i,j+1)-potential(i,j))*delta_x/delta_y;
		}
		elec_current_x(i,Ny-1)=-Harmonic_Mean(mat_cond(i,Ny-1),mat_cond(i+1,Ny-1))*(potential(i+1,Ny-1)-potential(i,Ny-1))*delta_y/delta_x;
	}
	for (int j=0;j<Ny-1;j++){
		elec_current_y(Nx-1,j)=-Harmonic_Mean(mat_cond(Nx-1,j),mat_cond(Nx-1,j+1))*(potential(Nx-1,j+1)-potential(Nx-1,j))*delta_x/delta_y;
	}
	electric_current_x=elec_current_x;electric_current_y=elec_current_y;
}

// Electric current average on a face
double TotalElectricCurrentEPM(Parameters param,ublas_matrix potential,ublas_matrix mat_cond,string face){
	double elec_flow=0;int Nx=param.simu_param.EPM_param.Nx,Ny=param.simu_param.EPM_param.Ny;
	ublas_matrix electric_current_x, electric_current_y;
	ElectricCurrentComputation(param,potential,mat_cond,electric_current_x,electric_current_y);
	if (face==BOTTOM_BORDER){
		for (int i=0;i<Nx;i++){
			elec_flow+=-electric_current_y(i,0);
		}
		return elec_flow;
	}
	else if (face==RIGHT_BORDER){
		for (int j=0;j<Ny;j++){
			elec_flow+=electric_current_x(Nx-2,j);
		}
		return elec_flow;
	}
	cout << "WARNING in AverageElectricCurrent (EPMDefinition.cpp): case not implemented" << endl;
	return -1;
}

// return equivalent conductivity
double EquivalentConductivityEPM(Parameters param,ublas_matrix potential,ublas_matrix mat_cond){
	if (param.simu_param.bc_map[TOP_BORDER].bc_type==DIRICHLET&&param.simu_param.bc_map[BOTTOM_BORDER].bc_type==DIRICHLET){
		return TotalElectricCurrentEPM(param,potential,mat_cond,BOTTOM_BORDER)/param.domain_param.Lx*param.domain_param.Ly/(param.simu_param.bc_map[TOP_BORDER].bc_value-param.simu_param.bc_map[BOTTOM_BORDER].bc_value);
	}
	cout << "WARNING in AverageElectricCurrent (EPMDefinition.cpp): case not implemented" << endl;
	return -1;
}

//typedef std::map<std::pair<double,double>,std::map<std::string,BoundaryConditions> > BoundaryConditionsEPM;	// map<cell coordinate,map<border,boundary conditions> >
void print(BoundaryConditionsEPM bc_epm){
	for (BoundaryConditionsEPM::iterator it=bc_epm.begin();it!=bc_epm.end();it++){
		cout << "(" << it->first.first << "," << it->first.second << ")" << endl;
		for (map<std::string,BoundaryConditions>::iterator it1=it->second.begin();it1!=it->second.end();it1++){
			cout << it1->first << " " << it1->second.bc_type << " " << it1->second.bc_value << " " << it1->second.block_index << endl;
		}
	}
}
