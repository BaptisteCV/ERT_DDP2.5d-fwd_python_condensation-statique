/*
 * DDPDFNFunctions.cpp
 *
 *  Created on: May 27, 2013
 *      Author: delphineroubinet
 */

#include "DDPDFNFunctions.h"
#include "../Utilitaries/Constantes.h"


using namespace std;


// Redefine the fracture network with its intersection with the grid
void FracturesIntersectionGrid(int Nx,int Ny,NetworkMeshes& net_mesh){
	// 1. Variables
	double i_min=net_mesh.domain.min_pt.i,j_min=net_mesh.domain.min_pt.j,
			i_max=net_mesh.domain.max_pt.i,j_max=net_mesh.domain.max_pt.j;
	double delta_x=net_mesh.domain.ReturnDeltaX(Nx),delta_y=net_mesh.domain.ReturnDeltaY(Ny);
	// 2. Grid lines definition
	vector<Segment2D> grid_lines;
	// 2.1. Horizontal lines
	for (int j=1;j<Ny;j++){grid_lines.push_back(Segment2D(CgalPoint2D(i_min,j_min+j*delta_y),CgalPoint2D(i_max,j_min+j*delta_y)));}
	// 2.2. Vertical lines
	for (int i=1;i<Nx;i++){grid_lines.push_back(Segment2D(CgalPoint2D(i_min+i*delta_x,j_min),CgalPoint2D(i_min+i*delta_x,j_max)));}
	// 3. Look for intersection for each fracture
	CgalPoint2D inter;int cpt_node=net_mesh.cpt_inter;bool intersection=true;
	FractureMesh current_fract,fract1,fract2;Segment2D seg_fract;CGAL::Object result;FluxPoint2D inter_node;vector<Segment2D>::iterator it2;
	vector<FractureMesh>::iterator it1;
	int syst_size;
	while(intersection){
		syst_size=net_mesh.meshes.size();
		for (size_t i=0;(unsigned)i<syst_size;i++){
			current_fract=net_mesh.meshes[i];
			seg_fract=Segment2D(current_fract.p_ori.p,current_fract.p_tar.p);intersection=true;
			for (it2=grid_lines.begin();it2!=grid_lines.end();it2++){
				intersection=false;
				result=CGAL::intersection(seg_fract,*it2);
				// if there is intersection -> the studied fracture is divided in 2 fractures
				if (CGAL::assign(inter,result)){
					inter_node=FluxPoint2D(inter,cpt_node+1);
					fract1=FractureMesh(current_fract.p_ori,inter_node,current_fract.aperture,current_fract.conductivity);
					fract2=FractureMesh(inter_node,current_fract.p_tar,current_fract.aperture,current_fract.conductivity);
					if (fract1.ReturnLength()>EPSILON&&fract2.ReturnLength()>EPSILON){
						intersection=true;
						net_mesh.meshes[i]=fract1;net_mesh.meshes.push_back(fract2);
						cpt_node++;
						break;
					}
				}
			}
			if (intersection){break;}
		}
	}
	net_mesh.cpt_inter=cpt_node;
}

SubNetworkMap DefineSubNetwork(DomainMap sub_domains,NetworkMeshes net_mesh){
	SubNetworkMap subnetwork;bool in_domain;
	// 1. Subnetwork initialization
	for (DomainMap::iterator it=sub_domains.begin();it!=sub_domains.end();it++){
		subnetwork[it->first].domain=it->second;
	}
	// 2. For each fracture: check that it is in the domain and check if its extremities is on the big domain borders
	for (vector<FractureMesh>::iterator it1=net_mesh.meshes.begin();it1!=net_mesh.meshes.end();it1++){
		in_domain=false;
		for (DomainMap::iterator it2=sub_domains.begin();it2!=sub_domains.end();it2++){
			// Check if each fracture is in each domain
			if (it2->second.IsInDomain(it1->ReturnSegment())){
				subnetwork[it2->first].meshes.push_back(*it1);
				if (subnetwork[it2->first].inter_list.find(it1->p_ori.index)==subnetwork[it2->first].inter_list.end()){
					subnetwork[it2->first].cpt_inter++;
				}
				if (subnetwork[it2->first].inter_list.find(it1->p_tar.index)==subnetwork[it2->first].inter_list.end()){
					subnetwork[it2->first].cpt_inter++;
				}
				in_domain=true;
				// Create the border map where the borders are the BIG domain borders
				if (net_mesh.border_map.find(it1->p_ori.index)!=net_mesh.border_map.end()){
					subnetwork[it2->first].border_map[it1->p_ori.index]=net_mesh.border_map[it1->p_ori.index];
				}
				if (net_mesh.border_map.find(it1->p_tar.index)!=net_mesh.border_map.end()){
					subnetwork[it2->first].border_map[it1->p_tar.index]=net_mesh.border_map[it1->p_tar.index];
				}
				break;
			}
		}
		if (!in_domain){
			cout << "WARNING in DefineSubNetwork (DDPDFNFunctions.cpp): fracture not inside a subdomain" << endl;
			print(it1->p_ori.p);print(it1->p_tar.p);
		}
	}
	return subnetwork;
}

// Define the DFN for the DDP configuration (addition of subfractures corresponding to fractures included in matrix blocks)
SubNetworkMap DDPNetworkMeshes(NetworkMeshes& net_mesh,int Nx,int Ny){
	// 1. Variables
	// 2. Definition of the new fracture network considering intersections with grid lines (if not already computed)
	if (net_mesh.meshes.size()>0&&!net_mesh.inter_grid){FracturesIntersectionGrid(Nx,Ny,net_mesh);}
	// 3. Blocks definition
	DomainMap sub_domains=DefineSubDomains(net_mesh.domain,Nx,Ny);
	// 4. Determination of fractures inside each matrix block
	return DefineSubNetwork(sub_domains,net_mesh);
}

// Return the linear system coefficients to express the derivative at a fracture node taking into account exchange with the matrix
void ReturnDerivativeCoeffFull(double x,double L,double k,double& a1,double& a2,double& a3){
	// 1. Variables
	double k_=sqrt(k),delta=exp(k_*x)+exp(-k_*x),
	gamma=exp(-k_*L)-exp(k_*L);
	// 2. Coefficient computation
	a1=k_*exp(k_*x)+k_*delta*exp(k_*L)/gamma;
	a2=-k_*delta/gamma;
	a3=-k_*exp(k_*x)+k_*delta*(1-exp(k_*L))/gamma;
}

void ReturnDerivativeCoeffFullFourier(double x,double L,double k,double w,double& a1,double& a2,double& a3){
	// 1. Variables
	double k_=sqrt(k+w*w),lambda=exp(k_*x)+exp(-k_*x),gammaL=exp(-k_*L)-exp(k_*L);
	// 2. Coefficient computation
	a1=k_*(exp(k_*x)+lambda*exp(k_*L)/gammaL);
	a2=-k_*lambda/gammaL;
	a3=-k*(a1+a2)/(k_*k_);
}

// return the cumulated potential value along a fracture taking into account interaction with the porous domain
void CumulatePotFractCoeff(double L,double k,double& coeff1,double& coeff2,double& coeff3){
	double k_=sqrt(k),gammaL=exp(-k_*L)-exp(k_*L);
	coeff1=(exp(k_*L)-1)/k_+exp(k_*L)*(exp(-k_*L)+exp(k_*L)-2)/(gammaL*k_);
	coeff2=-(exp(-k_*L)+exp(k_*L)-2)/(gammaL*k_);
	coeff3=L-coeff1-coeff2;
}

void CumulatePotFractCoeffFourier(double L,double k,double w,double& coeff1,double& coeff2,double& coeff3){
	double k_=sqrt(k+w*w),gammaL=exp(-k_*L)-exp(k_*L);
	coeff1=(exp(k_*L)-1)/k_+exp(k_*L)*(exp(-k_*L)+exp(k_*L)-2)/(gammaL*k_);
	coeff2=-(exp(-k_*L)+exp(k_*L)-2)/(gammaL*k_);
	coeff3=(L-coeff1-coeff2)*k/(k_*k_);
}
