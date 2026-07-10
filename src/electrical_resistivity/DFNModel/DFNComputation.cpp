/*
 * DFNComputation.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet
 */

#include "DFNComputation.h"


using namespace std;

BoundaryConditionsDFN ReturnBoundCondDFNFourier(Domain domain,BoundaryConditionsDef bc_map_classic,BordersMap border_map,NodesMap nodes_map,BoundaryConditionsEPM bc_def_epm, double delta_x,double delta_y){
	// 1. Variables
	BoundaryConditionsDFN bc_map;int index;string bc_type,border1,border2;double position,border_dist,bc_value;
	// 2. Loop on each node of the fracture network
	for (map<int,string>::iterator it=border_map.begin();it!=border_map.end();it++){
		index=it->first;bc_type=bc_map_classic[it->second].bc_type;
		// 2.1. Classical affectation (boundary condition value is known)
		if (bc_type==DIRICHLET||bc_type==NEUMANN){
			bc_map[index]=bc_map_classic[it->second];
		}
		// 2.2. Specific affectation (boundary condition value has to be defined)
		else if (bc_type==DIR_GRADIENT){
			// Position of the node on the border
			position=domain.ReturnBorderCoordinate(nodes_map[index],border_dist,border1,border2);
			// Boundary condition value at this position
			if (bc_map_classic[border1].bc_type!=DIRICHLET||bc_map_classic[border2].bc_type!=DIRICHLET){
				cout << "WARNING in ReturnBoundCondDFN (DFNComputation.cpp): case not implemented" << endl;
			}
			else{
				bc_value=(bc_map_classic[border2].bc_value-bc_map_classic[border1].bc_value)/border_dist*position+bc_map_classic[border1].bc_value;
				bc_map[index]=BoundaryConditions(DIRICHLET,bc_value);
			}
		}
		else if (bc_type==MIXED_BC){
			CgalPoint2D node_position=nodes_map[index];
			// little trick for numerical issue...
			double x_position=node_position.x(),y_position=node_position.y();
			if (x_position==domain.max_pt.i){x_position=x_position-EPSILON;}
			if (y_position==domain.max_pt.j){y_position=y_position-EPSILON;}
			int Nx=domain.domain_size_x()/delta_x,Ny=domain.domain_size_y()/delta_y;
			pair<int,int> indices=return_indices(x_position,y_position,delta_x,delta_y,Nx,Ny,domain.min_pt.i,domain.min_pt.j);
			map<string,BoundaryConditions> inter_map=bc_def_epm[indices];
			bc_map[index]=BoundaryConditions(MIXED_BC,inter_map[it->second].bc_value);
			bc_map[index].block_index=return_index(indices.first,indices.second,Ny);
		}
		else{cout << "WARNING in ReturnBoundCondDFN (DFNComputation.cpp): boundary condition type not implemented" << endl;}
	}
	return bc_map;
}

BoundaryConditionsDFN ReturnBoundCondDFN(Domain domain,BoundaryConditionsDef bc_map_classic,BordersMap border_map,NodesMap nodes_map){
	// 1. Variables
	BoundaryConditionsDFN bc_map;int index;string bc_type,border1,border2;double position,border_dist,bc_value;
	// 2. Loop on each node of the fracture network
	for (map<int,string>::iterator it=border_map.begin();it!=border_map.end();it++){
		index=it->first;bc_type=bc_map_classic[it->second].bc_type;
		// 2.1. Classical affectation (boundary condition value is known)
		if (bc_type==DIRICHLET||bc_type==NEUMANN){
			bc_map[index]=bc_map_classic[it->second];
		}
		// 2.2. Specific affectation (boundary condition value has to be defined)
		else if (bc_type==DIR_GRADIENT){
			// Position of the node on the border
			position=domain.ReturnBorderCoordinate(nodes_map[index],border_dist,border1,border2);
			// Boundary condition value at this position
			if (bc_map_classic[border1].bc_type!=DIRICHLET||bc_map_classic[border2].bc_type!=DIRICHLET){
				cout << "WARNING in ReturnBoundCondDFN (DFNComputation.cpp): case not implemented" << endl;
			}
			else{
				bc_value=(bc_map_classic[border2].bc_value-bc_map_classic[border1].bc_value)/border_dist*position+bc_map_classic[border1].bc_value;
				bc_map[index]=BoundaryConditions(DIRICHLET,bc_value);
			}
		}
		else if (bc_type==MIXED_BC){
			cout << "WARNING in ReturnBoundCondDFN (DFNComputation.cpp): mixed bc implemented as no flow in DFN" << endl;
		}
		else{cout << "WARNING in ReturnBoundCondDFN (DFNComputation.cpp): boundary condition type not implemented" << endl;}
	}
	return bc_map;
}

ublas_vector DFNComputationClassic(NetworkMeshes net_mesh,BoundaryConditionsDFN bc_map){
	SourceTermDFN source_term_DFN;
	return DFNComputationClassic(net_mesh,bc_map,source_term_DFN);
}
ublas_vector DFNComputationClassic(NetworkMeshes net_mesh,BoundaryConditionsDFN bc_map,SourceTermDFN source_term_DFN){
	ublas_vector DFN_potential=InitializationVector(net_mesh.cpt_inter);
	return DFNComputationClassic(net_mesh,bc_map,source_term_DFN,DFN_potential);
}
ublas_vector DFNComputationClassic(NetworkMeshes net_mesh,BoundaryConditionsDFN bc_map,SourceTermDFN source_term_DFN,
		ublas_vector DFN_potential){
	// 1. Nodes definition
	NodeFracturesMap node_fract_map=net_mesh.ReturnNodeFracturesMap();
	int nb_nodes=node_fract_map.size();
	// 2. Linear system definition
	ublas_matrix mat_syst(nb_nodes,nb_nodes); ublas_vector vect_syst(nb_nodes);
	// 2.1. Loop on each node
	double value;
	for (NodeFracturesMap::iterator it1=node_fract_map.begin();it1!=node_fract_map.end();it1++){
		//if (transient){mat_syst(it1->first,it1->first)+=1;}
		// Matrix affectation - Loop on corresponding fractures
		for (ElecPropMap::iterator it2=it1->second.begin();it2!=it1->second.end();it2++){
			value=it2->second.aperture*it2->second.conductivity/it2->second.length;
			mat_syst(it1->first,it2->first)=-value;
			mat_syst(it1->first,it1->first)+=value;
		}
		// Vector affectation
		vect_syst(it1->first)=0;
	}
	// 3. Boundary conditions affectation
	for (BoundaryConditionsDFN::iterator it1=bc_map.begin();it1!=bc_map.end();it1++){
		// electric potential affectation
		if (it1->second.bc_type==DIRICHLET){
			// cancel previous values
			for (int i=0;i<mat_syst.size2();i++){mat_syst(it1->first,i)=0;}
			// new value affectation
			mat_syst(it1->first,it1->first)=1;
			vect_syst(it1->first)=it1->second.bc_value;
		}
		// electric current affectation
		else if (it1->second.bc_type==NEUMANN){
			// new value affectation
			vect_syst(it1->first)=it1->second.bc_value;
		}
		else{cout << "WARNING in DFNComputation.cpp: boundary condition not implemented" << endl;}
	}
	// 4. Source term affectation
	for (SourceTermDFN::iterator it=source_term_DFN.begin();it!=source_term_DFN.end();it++){
		vect_syst(it->first)+=it->second+DFN_potential(it->first);
	}
	// 5. Linear system solution
	return LinearSystemSolving(mat_syst,vect_syst);
}

double ReturnAlpha(double cond_mat,double dist,double delta_x){
	dist=0.5*delta_x;
	// exchange coefficient for strate (1D exchange)
	double alpha1=PI*PI*cond_mat/(4*dist*dist);
	// exchange coefficient for cylinder (2D exchange)
	double alpha2=5.78*cond_mat/(dist*dist);
	// exchange coefficient for sphere (3D exchange)
	double alpha3=PI*PI*cond_mat/(dist*dist);
	// empirical exchange coefficient
	double aperture=1e-3;dist=16;
	double alpha_star=2*cond_mat/(aperture*dist*dist);
	double alpha4=2*cond_mat/(delta_x/4);
	return alpha4;
}

void ReturnDerivativeCoefficients(double x,double L,double pot_mat,double aperture,double cond_mat,double dist,double delta_x,double& a1,double& a2,double& a3){
	// 1. Variables
	double alpha=ReturnAlpha(cond_mat,dist,delta_x),
	alpha_=sqrt(alpha),q=-alpha*pot_mat,C=exp(alpha_*x)+exp(-alpha_*x),
	gamma=exp(-alpha_*L)-exp(alpha_*L);
	// 2. Coefficient computation
	a1=alpha_*exp(alpha_*x)+alpha_*C*exp(alpha_*L)/gamma;
	a2=-alpha_*C/gamma;
	a3=q*exp(alpha_*x)/alpha_-q*C*(1-exp(alpha_*L))/(gamma*alpha_);
}

ublas_vector DFNMatrixComputation(SubNetworkMap subnetwork_map,BoundaryConditionsDFN bc_map,ublas_matrix matrix_potential,double mat_cond,double delta_x){
	// 1. Variables
	double L,dist,a1,a2,a3;
	int nb_nodes=ReturnNbNodes(subnetwork_map),index1,index2,Ny=matrix_potential.size2();
	ublas_matrix mat_syst=InitializationMatrix(nb_nodes,nb_nodes); ublas_vector vect_syst=InitializationVector(nb_nodes);
	pair<int,int> indices;
	// 2. System definition
	for (SubNetworkMap::iterator it1=subnetwork_map.begin();it1!=subnetwork_map.end();it1++){
		indices=return_indices(it1->first,Ny);
		for (vector<FractureMesh>::iterator it2=it1->second.meshes.begin();it2!=it1->second.meshes.end();it2++){
			L=it2->ReturnLength();
			// averaged distance between the fracture and the block center
			dist=AveragedDistance(it2->ReturnSegment(),CenterPosition(indices,make_pair(delta_x,delta_x)));
			index1=it2->p_ori.index,index2=it2->p_tar.index;
			// fracture value affectation for the node1 (x=0 or phi_i)
			ReturnDerivativeCoefficients(0,L,matrix_potential(indices.first,indices.second),it2->aperture,mat_cond,dist,delta_x,a1,a2,a3);
			mat_syst(index1,index1)+=a1;
			mat_syst(index1,index2)+=a2;
			vect_syst(index1)+=a3;
			// fracture value affectation for the node2 (x=L or phi_(i+1))
			ReturnDerivativeCoefficients(L,L,matrix_potential(indices.first,indices.second),it2->aperture,mat_cond,dist,delta_x,a1,a2,a3);
			mat_syst(index2,index1)-=a1;
			mat_syst(index2,index2)-=a2;
			vect_syst(index2)-=a3;
		}
	}
	// 3. Boundary conditions affectation
	for (BoundaryConditionsDFN::iterator it1=bc_map.begin();it1!=bc_map.end();it1++){
		// electric potential affectation
		if (it1->second.bc_type==DIRICHLET){
			// cancel previous values
			for (int i=0;i<mat_syst.size2();i++){mat_syst(it1->first,i)=0;}
			// new value affectation
			mat_syst(it1->first,it1->first)=1;
			vect_syst(it1->first)=it1->second.bc_value;
		}
		// electric current affectation
		else if (it1->second.bc_type==NEUMANN){
			// new value affectation
			vect_syst(it1->first)=it1->second.bc_value;
		}
		else{cout << "WARNING in DFNComputation.cpp: boundary condition not implemented" << endl;}
	}
	// 5. Linear system solution
	return LinearSystemSolving(mat_syst,vect_syst);
}

/*void DFNDDPLinearSystem(SubNetworkMap subnetwork_map,BoundaryConditionsDFN bc_map,double mat_cond,double delta_x,ublas_matrix& mat_syst,ublas_vector& vect_syst){
	// 1. Variables
	double L,a1,a2,a3,coeff;
	int nb_nodes_fract=ReturnNbNodes(subnetwork_map),index1,index2;
	// 2. System definition
	// 2.1. DFN linear system
	for (SubNetworkMap::iterator it1=subnetwork_map.begin();it1!=subnetwork_map.end();it1++){
		for (vector<FractureMesh>::iterator it2=it1->second.meshes.begin();it2!=it1->second.meshes.end();it2++){
			coeff=it2->aperture*it2->conductivity;
			L=it2->ReturnLength();
			index1=it2->p_ori.index,index2=it2->p_tar.index;
			// fracture value affectation for the node1 (x=0 or phi_i)
			ReturnDerivativeCoeffFull(0,L,it2->conductivity,mat_cond,delta_x,a1,a2,a3);
			mat_syst(index1,index1)+=coeff*a1;
			mat_syst(index1,index2)+=coeff*a2;
			mat_syst(index1,it1->first+nb_nodes_fract)+=coeff*a3;
			// fracture value affectation for the node2 (x=L or phi_(i+1))
			ReturnDerivativeCoeffFull(L,L,it2->conductivity,mat_cond,delta_x,a1,a2,a3);
			mat_syst(index2,index1)-=coeff*a1;
			mat_syst(index2,index2)-=coeff*a2;
			mat_syst(index2,it1->first+nb_nodes_fract)-=coeff*a3;
		}
	}
	// 3. Boundary conditions affectation
	// 3.1. DFN boundary conditions
	for (BoundaryConditionsDFN::iterator it1=bc_map.begin();it1!=bc_map.end();it1++){
		// electric potential affectation
		if (it1->second.bc_type==DIRICHLET){
			// cancel previous values
			for (int i=0;i<mat_syst.size2();i++){mat_syst(it1->first,i)=0;}
			// new value affectation
			mat_syst(it1->first,it1->first)=1;
			vect_syst(it1->first)=it1->second.bc_value;
		}
		// electric current affectation
		else if (it1->second.bc_type==NEUMANN){
			// new value affectation
			vect_syst(it1->first)=it1->second.bc_value;
		}
		else{cout << "WARNING in DFNComputation.cpp: boundary condition not implemented" << endl;}
	}
}*/

ublas_vector DFNComputationRadialDiffusion(NetworkMeshes net_mesh,BoundaryConditionsDFN bc_map,double cond_fract,double conductivity_matrix){
	// 1. Nodes definition
	NodeFracturesMap node_fract_map=net_mesh.ReturnNodeFracturesMap();
	int nb_nodes=node_fract_map.size();
	// 2. Linear system definition
	ublas_matrix mat_syst(nb_nodes,nb_nodes); ublas_vector vect_syst(nb_nodes);
	// 2.1. Loop on each node
	double value1,value2,coeff,length,k;
	for (NodeFracturesMap::iterator it1=node_fract_map.begin();it1!=node_fract_map.end();it1++){
		// Matrix affectation - Loop on corresponding fractures
		for (ElecPropMap::iterator it2=it1->second.begin();it2!=it1->second.end();it2++){
			coeff=it2->second.aperture*it2->second.conductivity;
			length=it2->second.length;
			k=4*conductivity_matrix/(it2->second.aperture*0.5*it2->second.aperture*0.5*it2->second.conductivity);
			value1=coeff*2*sqrt(k)/(1-exp(-2*sqrt(k)*length));
			value2=coeff*2*sqrt(k)*(1+exp(-2*sqrt(k)*length))/(1-exp(-2*sqrt(k)*length));
			mat_syst(it1->first,it2->first)=-value1;
			mat_syst(it1->first,it1->first)+=value2;
		}
		// Vector affectation
		vect_syst(it1->first)=0;
	}
	// 3. Boundary conditions affectation
	for (BoundaryConditionsDFN::iterator it1=bc_map.begin();it1!=bc_map.end();it1++){
		// electric potential affectation
		if (it1->second.bc_type==DIRICHLET){
			// cancel previous values
			for (int i=0;i<mat_syst.size2();i++){mat_syst(it1->first,i)=0;}
			// new value affectation
			mat_syst(it1->first,it1->first)=1;
			vect_syst(it1->first)=it1->second.bc_value;
		}
		// electric current affectation
		else if (it1->second.bc_type==NEUMANN){
			// new value affectation
			vect_syst(it1->first)=it1->second.bc_value;
		}
		else{cout << "WARNING in DFNComputation.cpp: boundary condition not implemented" << endl;}
	}
	// 4. Linear system solution
	return LinearSystemSolving(mat_syst,vect_syst);
}

double AnalyticalSolutionFractureClassic(double x,double x0,double x1,double y0,double y1){
	return (y1-y0)/(x1-x0)*x+(y0*x1-y1*x0)/(x1-x0);
}

double AnalyticalSolutionFractureMatrix(double x,double L,double y0,double y1,double k,double pot_matrix){
	double beta,gammax,gammaL,k_;
	k_=sqrt(k);
	gammax=exp(-k_*x)-exp(k_*x);
	gammaL=exp(-k_*L)-exp(k_*L);
	beta=exp(k_*x)-gammax/gammaL*exp(k_*L);
	return beta*y0+gammax/gammaL*y1+(1-gammax/gammaL-beta)*pot_matrix;
}

double AnalyticalSolutionFractureRadialDiffusion(double x,double x0,double x1,double y0,double y1,double cond_matrix,double aperture,double cond_fract){
	double length=x1-x0;
	double k=4*cond_matrix/(aperture*0.5*aperture*0.5*cond_fract);
	double A=(y0*exp(-sqrt(k)*length)-y1)/(exp(sqrt(k)*(x0-length))-exp(sqrt(k)*x1));
	double B=(y1-y0*exp(sqrt(k)*length))/(exp(-sqrt(k)*x1)-exp(sqrt(k)*(length-x0)));
	return A*exp(sqrt(k)*x)+B*exp(-sqrt(k)*x);
}

double MetricElecFlow(double fract_cond,double pot2,double pot1,double L){
	return -fract_cond*(pot2-pot1)/L;
}

// return the total electric current going out of extremities with a potential enforced to 0
double TotalElectricCurrentDFN(Parameters param,NetworkMeshes net_mesh,ublas_vector potential,BoundaryConditionsDFN bc_map,double& output_fract_surf){
	double total_elec_flow=0;output_fract_surf=0;
	NodeFracturesMap node_fract_map=net_mesh.ReturnNodeFracturesMap();
	for (BoundaryConditionsDFN::iterator it1=bc_map.begin();it1!=bc_map.end();it1++){
		if (it1->second.bc_type==DIRICHLET&&it1->second.bc_value==0){
			for (ElecPropMap::iterator it2=node_fract_map[it1->first].begin();it2!=node_fract_map[it1->first].end();it2++){
				total_elec_flow+=it2->second.aperture*MetricElecFlow(it2->second.conductivity,potential[it1->first],potential[it2->first],it2->second.length);
				output_fract_surf+=it2->second.aperture;
			}
		}
	}
	return total_elec_flow;
}

ublas_vector ComputeFlowVelocities(NetworkMeshes& net_mesh,Parameters param,Domain domain,string option){
	cout << "MODIF FLOW VELOCITY CONDITIONS" << endl;
	param.simu_param.bc_map[LEFT_BORDER]=BoundaryConditions(DIRICHLET,0);
	param.simu_param.bc_map[RIGHT_BORDER]=BoundaryConditions(DIRICHLET,1);
	param.simu_param.bc_map[TOP_BORDER]=BoundaryConditions(DIRICHLET,2);
	param.simu_param.bc_map[BOTTOM_BORDER]=BoundaryConditions(DIRICHLET,3);
	BoundaryConditionsDFN bc_map_dfn=ReturnBoundCondDFN(domain,param.simu_param.bc_map,net_mesh.border_map,net_mesh.nodes_map);
	// hydraulic head computation
	SourceTermDFN source_term_DFN;
	ublas_vector DFN_potential=DFNComputationClassic(net_mesh,bc_map_dfn,source_term_DFN);
	// flow velocity computation
	net_mesh.EvaluateFlowVelocities(DFN_potential);
	return DFN_potential;
}



