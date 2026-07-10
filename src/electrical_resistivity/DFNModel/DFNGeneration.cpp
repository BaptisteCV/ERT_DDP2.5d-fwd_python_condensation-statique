/*
 * DFNGeneration.cpp
 *
 *  Created on: Apr 23, 2013
 *      Author: delphineroubinet
 */

#include <set>
#include "DFNGeneration.h"
#include "../Utilitaries/Storage/UtilitaryStorage.h"
#include "../Utilitaries/Constantes.h"
#include "../RngStream/rngstream.h"
#include "../InputOutput/Results.h"

using namespace std;
/*
DFNVisu DFNGeneration(Parameters param,Domain domain){
	// DFN definition
	NetworkMeshes net_mesh(domain,param,param.simu_param.simu_option,"no_intersection");
	// DFN visualization
	return DFNVisu(param.domain_param.Lx,net_mesh,DFN);
}
*/
void Border_Definition(FluxPoint2D pt_inter,NetworkMeshes& net_mesh){
	string border=net_mesh.domain.ReturnBorder(pt_inter.p);
	if (border!=NO_BORDER){net_mesh.border_map[pt_inter.index]=border;}
	net_mesh.nodes_map[pt_inter.index]=pt_inter.p;
}

void AddFracture(FractureMesh frac_mesh,NetworkMeshes& net_mesh){
	// indices determination
	int cpt_index=net_mesh.cpt_inter;
	// definition of extremity numbering
	int cpt1=net_mesh.ReturnIndex(frac_mesh.p_ori.p),cpt2=net_mesh.ReturnIndex(frac_mesh.p_tar.p);
	if (cpt1==-1){cpt1=cpt_index+1;cpt_index++;}
	if (cpt2==-1){cpt2=cpt_index+1;cpt_index++;}
	// fracture definition
	FluxPoint2D pt1(frac_mesh.p_ori.p,cpt1),pt2(frac_mesh.p_tar.p,cpt2);
	FractureMesh new_fract(pt1,pt2,frac_mesh.aperture,frac_mesh.conductivity);
	net_mesh.meshes.push_back(new_fract);
	// boundary conditions update
	net_mesh.cpt_inter=cpt_index;
	Border_Definition(pt1,net_mesh);
	Border_Definition(pt2,net_mesh);
}


void DivideDomain(pointcpp<double> pt_min,pointcpp<double> pt_max,int nb_division,double aperture,double conductivity,NetworkMeshes& net_mesh){
	// variables
	FractureMesh frac_mesh;
	double x_current=pt_min.i,y_current=pt_min.j,x_min=pt_min.i,y_min=pt_min.j,x_max=pt_max.i,y_max=pt_max.j;
	double fract_spacing_x=(pt_max.i-pt_min.i)/nb_division,fract_spacing_y=(pt_max.j-pt_min.j)/nb_division;
	// creation of the horizontal fractures
	for (int i=0;i<nb_division-1;i++){
		y_current+=fract_spacing_y;
		frac_mesh=FractureMesh(FluxPoint2D(x_min,y_current),FluxPoint2D(x_max,y_current),aperture,conductivity);
		AddFracture(frac_mesh,net_mesh);
	}
	// addition of the vertical fractures
	for (int i=0;i<nb_division-1;i++){
		x_current+=fract_spacing_x;
		frac_mesh=FractureMesh(FluxPoint2D(x_current,y_min),FluxPoint2D(x_current,y_max),aperture,conductivity);
		AddFracture(frac_mesh,net_mesh);
	}
}

void ComputeIntersections(NetworkMeshes& net_mesh){
	// 1. Variables
	bool intersection=true;
	FractureMesh fract1,fract2,new_fract1,new_fract2,new_fract3,new_fract4;
	Segment2D seg1,seg2;CgalPoint2D inter;FluxPoint2D pt_inter;
	int cpt_inter=net_mesh.cpt_inter,cpt_current;
	// 2. Loop on fracture mesh until there is no undefined intersection
	while (intersection){
		intersection=false;
		for (size_t i=0;i<net_mesh.meshes.size();i++){
			// 2.1. First segment definition
			fract1=net_mesh.meshes[i];
			seg1=fract1.ReturnSegment();
			for (size_t j=0;j<net_mesh.meshes.size();j++){
				bool new_node=false;
				if (i!=j){
					// 2.2. Second segment definition
					fract2=net_mesh.meshes[j];
					seg2=fract2.ReturnSegment();
					// 2.3. Check if there is intersection
					CGAL::Object result=CGAL::intersection(seg1,seg2);
					// if there is intersection
					//if (CGAL::assign(inter,result)&&!seg1.IdenticExtremities(inter)&&!seg2.IdenticExtremities(inter)){
					if (!seg1.IdenticExtremities(seg2.source())&&!seg1.IdenticExtremities(seg2.target())&&CGAL::assign(inter,result)){
					//if (seg1.intersection(seg2,inter)){
						// intersection definition
						cpt_current=net_mesh.ReturnIndex(inter);
						if (cpt_current==-1){cpt_current=cpt_inter+1;new_node=true;}
						pt_inter=FluxPoint2D(inter,cpt_current);
						// define the potential new subsegments from the first and second fracture
						new_fract1=FractureMesh(fract1.p_ori,pt_inter,fract1.aperture,fract1.conductivity);
						new_fract2=FractureMesh(pt_inter,fract1.p_tar,fract1.aperture,fract1.conductivity);
						new_fract3=FractureMesh(fract2.p_ori,pt_inter,fract2.aperture,fract1.conductivity);
						new_fract4=FractureMesh(pt_inter,fract2.p_tar,fract2.aperture,fract1.conductivity);
						// check if it is a real new intersection
						int cpt_true_inter=0;
						if (new_fract1.ReturnLength()>EPSILON){cpt_true_inter++;}
						if (new_fract2.ReturnLength()>EPSILON){cpt_true_inter++;}
						if (new_fract3.ReturnLength()>EPSILON){cpt_true_inter++;}
						if (new_fract4.ReturnLength()>EPSILON){cpt_true_inter++;}
						if (cpt_true_inter>=3){intersection=true;}
						// if it is -> addition of the new subsegments
						if (intersection){
							if (new_fract1.ReturnLength()>EPSILON){
								net_mesh.meshes[i]=new_fract1;
								if (new_fract2.ReturnLength()>EPSILON){net_mesh.meshes.push_back(new_fract2);}
							}
							else{net_mesh.meshes[i]=new_fract2;}
							if (new_fract3.ReturnLength()>EPSILON){
								net_mesh.meshes[j]=new_fract3;
								if (new_fract4.ReturnLength()>EPSILON){net_mesh.meshes.push_back(new_fract4);}
							}
							else{net_mesh.meshes[j]=new_fract4;}
							if (new_node){cpt_inter++;}
							//net_mesh.print_DFN();
							break;
						}
					}
				}
			}
			if (intersection){break;}
		}
	}
	// 3. Update network values
	net_mesh.cpt_inter=cpt_inter;
}

// Divide radomly few blocks of the domain and return the list of the divided blocks
vector<Domain> DivideSierpinskiDomain(pointcpp<double> domain_min,int nb_division,double fract_spacing_x,double fract_spacing_y,double aperture,double conductivity,NetworkMeshes& net_mesh){
	// Parameters
	set<int> index_domain_list;int index;pair<int,int> indices;
	vector<Domain> Domain_list;
	// for each domain to divide
	for (int j=0;j<nb_division;j++){
		// Pick the domain to divide
		index=(int)(rand()/(double)RAND_MAX*nb_division*nb_division);
		while (index_domain_list.find(index)!=index_domain_list.end()){
			index=(int)(rand()/(double)RAND_MAX*nb_division*nb_division);
		}
		index_domain_list.insert(index);
		// determine its extremities
		indices=return_indices(index,nb_division);
		pointcpp<double> pt_min_current(domain_min.i+indices.first*fract_spacing_x,domain_min.j+indices.second*fract_spacing_y),
				pt_max_current(domain_min.i+(indices.first+1)*fract_spacing_x,domain_min.j+(indices.second+1)*fract_spacing_y);
		Domain_list.push_back(Domain(pt_min_current,pt_max_current));
		// divide it
		DivideDomain(pt_min_current,pt_max_current,nb_division,aperture,conductivity,net_mesh);
	}
	return Domain_list;
}

// Define Sierpinski structures (whole horizontal and vertical fractures)
NetworkMeshes DFNGenerationSierpinski(Domain domain,double aperture,double conductivity,int nb_division,int level_division,int seed,string option,string intersection_option){
	// 0. Variables
	NetworkMeshes net_mesh;net_mesh.domain=domain;
	pointcpp<double> pt_min_init=domain.min_pt,pt_max_init=domain.max_pt;
	srand(seed);
	// 1. First domain division
	DivideDomain(pt_min_init,pt_max_init,nb_division,aperture,conductivity,net_mesh);
	// 2. Next domain division
	vector<Domain> divided_domains;divided_domains.push_back(domain);
	for (int i=1;i<level_division;i++){
		// division of the domain
		vector<Domain> new_divided_domains,current_divided_domains;
		for (vector<Domain>::iterator it1=divided_domains.begin();it1!=divided_domains.end();it1++){
			double fract_spacing_x=(it1->max_pt.i-it1->min_pt.i)/nb_division,fract_spacing_y=(it1->max_pt.j-it1->min_pt.j)/nb_division;
			current_divided_domains=DivideSierpinskiDomain(it1->min_pt,nb_division,fract_spacing_x,fract_spacing_y,aperture,conductivity,net_mesh);
			new_divided_domains.insert(new_divided_domains.end(),current_divided_domains.begin(),current_divided_domains.end());
		}
		divided_domains=new_divided_domains;
	}
	// 3. Intersections determination
	if (option==POTENTIAL_COMPUTATION&&intersection_option=="intersection"){
		ComputeIntersections(net_mesh);
	}
	return net_mesh;
}

NetworkMeshes DFNGenerationParallel(Domain domain,FractureFamilyMap fract_family_map,string intersection_option){
	// 1. Variables
	NetworkMeshes net_mesh;net_mesh.domain=domain;
	FractureMesh frac_mesh;
	pointcpp<double> pt_min=domain.min_pt,pt_max=domain.max_pt;
	double x_min=pt_min.i,x_max=pt_max.i,aperture,conductivity,Lx=domain.domain_size_x(),Ly=domain.domain_size_y();
	int fracture_nb;
	// 2. Creation of fracture network for each family
	for (FractureFamilyMap::iterator it1=fract_family_map.begin();it1!=fract_family_map.end();it1++){
		NetworkMeshes sub_net_mesh;sub_net_mesh.domain=domain;
		fracture_nb=it1->second.fracture_nb;aperture=it1->second.aperture;conductivity=it1->second.conductivity;
		double y_current,fract_spacing;
		// 2.1. Fracture equally distributed in the domain
		if (it1->second.fract_spacing==-1){
			y_current=pt_min.j;
			fract_spacing=(pt_max.j-pt_min.j)/(fracture_nb+1);
		}
		else{
			fract_spacing=it1->second.fract_spacing;
			y_current=(Ly-fract_spacing*(fracture_nb-1))*0.5-fract_spacing;
		}
		// 2.1. Addition of the given number of fractures
		for (int i=0;i<fracture_nb;i++){
			y_current+=fract_spacing;
			frac_mesh=FractureMesh(FluxPoint2D(x_min,y_current),FluxPoint2D(x_max,y_current),aperture,conductivity);
			//frac_mesh=FractureMesh(FluxPoint2D(x_min-Lx,y_current),FluxPoint2D(x_max+Lx,y_current),aperture,conductivity);
			AddFracture(frac_mesh,sub_net_mesh);
		}
		if (intersection_option=="intersection"){
			ComputeIntersections(sub_net_mesh);
		}
		// 2.2. Rotation of the domain depending of the fracture angle
		sub_net_mesh=DFNRotation(sub_net_mesh,it1->second.ReturnAngleRadian());
		// 2.3. Addition to the global network
		for (vector<FractureMesh>::iterator it2=sub_net_mesh.meshes.begin();it2!=sub_net_mesh.meshes.end();it2++){
			AddFracture(*it2,net_mesh);
		}
		if (intersection_option=="intersection"){
			ComputeIntersections(net_mesh);
		}
	}
	return net_mesh;
}

NetworkMeshes DFNGenerationParallelBorder(Domain domain,double aperture,double conductivity,string intersection_option){
	// 1. Variables
	NetworkMeshes net_mesh;net_mesh.domain=domain;
	FractureMesh frac_mesh;
	pointcpp<double> pt_min=domain.min_pt,pt_max=domain.max_pt;
	double x_min=pt_min.i+EPSILON,x_max=pt_max.i-EPSILON,y_min=pt_min.j+EPSILON,y_max=pt_max.j-EPSILON;
	// 2. Creation of fracture network for each family
	NetworkMeshes sub_net_mesh;sub_net_mesh.domain=domain;
	// 2.1. Left border
	//x_min=x_min+domain.domain_size_x()*0.5;x_max=x_max+domain.domain_size_x()*0.5;
	frac_mesh=FractureMesh(FluxPoint2D(x_min,y_min),FluxPoint2D(x_min,y_max),aperture,conductivity);
	AddFracture(frac_mesh,sub_net_mesh);
	// 2.2. Right border
	/*frac_mesh=FractureMesh(FluxPoint2D(x_max,y_min),FluxPoint2D(x_max,y_max),aperture,conductivity);
	AddFracture(frac_mesh,sub_net_mesh);*/

	if (intersection_option=="intersection"){
		ComputeIntersections(sub_net_mesh);
	}
	// 2.3. Addition to the global network
	for (vector<FractureMesh>::iterator it2=sub_net_mesh.meshes.begin();it2!=sub_net_mesh.meshes.end();it2++){
		AddFracture(*it2,net_mesh);
	}
	if (intersection_option=="intersection"){
		ComputeIntersections(net_mesh);
	}
	return net_mesh;
}

NetworkMeshes DFNGenerationRandomSingle(Domain domain,double fract_aperture,double fract_conductivity,double fract_length, int seed, int nb_fract, string intersection_option){
	NetworkMeshes net_mesh;net_mesh.domain=domain;FractureMesh frac_mesh;
	Segment2D seg;double angle;CgalPoint2D fract_center;
	double Lx=domain.domain_size_x(),Ly=domain.domain_size_y();
	srand(seed);
	for (size_t i=0;(unsigned)i<nb_fract;i++){
		fract_center=CgalPoint2D(Uniform()*Lx,Uniform()*Ly);
		angle=Uniform()*2*PI;
		seg=Segment2D(fract_center,fract_length,angle);
		domain.SegmentIntersectDomain(seg,seg);
		frac_mesh=FractureMesh(seg,fract_aperture,fract_conductivity);
		AddFracture(frac_mesh,net_mesh);
	}
	if (intersection_option=="intersection"){
		ComputeIntersections(net_mesh);
	}
	return net_mesh;
}

// generates a DFN composed of parallel fractures spacing out by the distance d
NetworkMeshes DFNGenerationSugarBox(Domain domain,double fract_spacing,double aperture,double conductivity){
	NetworkMeshes net_mesh;net_mesh.domain=domain;
	FractureMesh frac_mesh;
	vector<FractureMesh> new_meshes;
	// 0. Parameters
	double i_min=domain.min_pt.i,j_min=domain.min_pt.j,i_current,j_current;
	double L=domain.max_pt.i-i_min;
	// number of fractures
	int N=(int)(L/fract_spacing)+1;
	// 1. Fracture definition
	int index=-1;
	// creation of the horizontal fractures
	for (size_t j=0;(unsigned)j<N;j++){
		j_current=j_min+j*fract_spacing;
		index++;
		for (size_t i=0;(unsigned)i<N-1;i++){
			i_current=i_min+i*fract_spacing;
			frac_mesh=FractureMesh(FluxPoint2D(i_current,j_current,index),FluxPoint2D(i_current+fract_spacing,j_current,index+1),aperture,conductivity);
			new_meshes.push_back(frac_mesh);
			index+=1;
		}
	}
	net_mesh.meshes=new_meshes;
	// addition of the vertical fractures
	CgalPoint2D p1,p2;int index1,index2;
	for (size_t i=0;(unsigned)i<N;i++){
		i_current=i_min+i*fract_spacing;
		for (size_t j=0;j<N-1;j++){
			j_current=j_min+j*fract_spacing;
			// determination of the corresponding point
			p1=CgalPoint2D(i_current,j_current),p2=CgalPoint2D(i_current,j_current+fract_spacing);
			// look for the point in the existing fractures
			index1=net_mesh.ReturnIndex(p1);index2=net_mesh.ReturnIndex(p2);
			// addition of the corresponding fracture
			frac_mesh=FractureMesh(FluxPoint2D(p1,index1),FluxPoint2D(p2,index2),aperture,conductivity);
			new_meshes.push_back(frac_mesh);
		}
	}
	net_mesh.meshes=new_meshes;
	return net_mesh;
}

// generate DFN with where fracture length is defined by a power law and fracture density by a percolation parameter
NetworkMeshes DFNGenerationRandom(Domain domain,double power_law_exp,double percolation_param,double fract_aperture,double fract_conductivity,int seed,string intersection_option){
	NetworkMeshes net_mesh;net_mesh.domain=domain;FractureMesh frac_mesh;
	Segment2D seg;double density=0,fract_length,angle;CgalPoint2D fract_center;
	double Lx=domain.domain_size_x(),Ly=domain.domain_size_y();
	//int seed=0;
	srand(seed);
	while (density<percolation_param){
		fract_length=pow(Uniform(),1./(1-power_law_exp));
		fract_center=CgalPoint2D(Uniform()*Lx,Uniform()*Ly);
		angle=Uniform()*2*PI;
		seg=Segment2D(fract_center,fract_length,angle);
		domain.SegmentIntersectDomain(seg,seg);
		frac_mesh=FractureMesh(seg,fract_aperture,fract_conductivity);
		AddFracture(frac_mesh,net_mesh);
		density+=seg.length*seg.length/(Lx*Ly);
	}
	if (intersection_option=="intersection"){
		ComputeIntersections(net_mesh);
	}
	return net_mesh;
}

/********************************************************/
// Functions to manipulate and modify fracture networks //
/********************************************************/

// Return the fracture network rotated of the angle current_angle relatively to the domain center
NetworkMeshes DFNRotation(NetworkMeshes init_net_mesh,double current_angle){
	// 1. Variables
	NetworkMeshes new_net_mesh(init_net_mesh.domain,init_net_mesh.cpt_inter,init_net_mesh.inter_list);
	vector<FractureMesh> init_meshes=init_net_mesh.meshes,new_meshes;FractureMesh current_fract,new_fract;
	CgalPoint2D pt1,pt2;
	double Lx=init_net_mesh.domain.domain_size_x(),Ly=init_net_mesh.domain.domain_size_y();
	// 2. Loop on the fractures of the initial fracture network
	for (int i=0;i<init_meshes.size();i++){
		// 2.1. Current fracture definition
		current_fract=init_meshes[i];new_fract=current_fract;
		// 2.2. New fracture determination
		// Extremities definition and translation (relatively to the domain center instead of the left bottom corner)
		pt1=CgalPoint2D(current_fract.p_ori.p.x()-Lx*0.5,current_fract.p_ori.p.y()-Ly*0.5);
		pt2=CgalPoint2D(current_fract.p_tar.p.x()-Lx*0.5,current_fract.p_tar.p.y()-Ly*0.5);
		// Extremities rotation
		//Transformation rotate(CGAL::ROTATION,sin(current_angle),cos(current_angle));
		//pt1=rotate(pt1);pt2=rotate(pt2);
		//cout << "WARNING in DFNRotation (DFNGeneration.cpp): lines commented -> no rotation" << endl;
		pt1=rotation(pt1,current_angle);pt2=rotation(pt2,current_angle);
		// Translation to go back to the initial reference
		pt1=CgalPoint2D(pt1.x()+Lx*0.5,pt1.y()+Ly*0.5);
		pt2=CgalPoint2D(pt2.x()+Lx*0.5,pt2.y()+Ly*0.5);
		// Fracture addition to the new network
		new_fract.p_ori.p=pt1;new_fract.p_tar.p=pt2;
		new_meshes.push_back(new_fract);
	}
	new_net_mesh.meshes=new_meshes;
	return new_net_mesh;
}

// Return the fracture network defined from a centered extraction of the initial fracture network
bool DFNExtraction(NetworkMeshes& init_net_mesh,double extract_length){
	// 1. Variables
	NetworkMeshes new_net_mesh;vector<FractureMesh> init_meshes=init_net_mesh.meshes;
	FractureMesh current_fract,new_fract;Segment2D current_seg,new_seg;
	double init_length_x=init_net_mesh.domain.domain_size_x(),init_length_y=init_net_mesh.domain.domain_size_y();
	// 2. Definition of the new domain for a centered extraction
	Domain init_domain=init_net_mesh.domain,selected_domain;
	new_net_mesh.domain=Domain(pointcpp<double>(0,0),pointcpp<double>(extract_length,extract_length));
	double L_trans_x=0.5*init_length_x-0.5*extract_length,L_trans_y=0.5*init_length_y-0.5*extract_length;
	selected_domain=Domain(pointcpp<double>(init_length_x*0.5-extract_length*0.5,init_length_y*0.5-extract_length*0.5),pointcpp<double>(init_length_x*0.5+extract_length*0.5,init_length_y*0.5+extract_length*0.5));
	// 3. Definition of the fracture inside the new domain
	bool fracture=false;
	for (int i=0;(unsigned)i<init_meshes.size();i++){
		// 3.1. Current fracture definition
		current_fract=init_meshes[i];current_seg=current_fract.ReturnSegment();
		if (selected_domain.SegmentIntersectDomain(current_seg,new_seg)){
			new_seg=Segment2D(Translate(new_seg.source(),-L_trans_x,-L_trans_y),Translate(new_seg.target(),-L_trans_x,-L_trans_y));
			new_fract=FractureMesh(new_seg,current_fract.aperture,current_fract.conductivity);
			if (new_fract.ReturnLength()>0){AddFracture(new_fract,new_net_mesh);}
			fracture=true;
		}
	}
	// 4. Intersection computation in the new domain
	ComputeIntersections(new_net_mesh);
	init_net_mesh=new_net_mesh;
	return fracture;
}

// Return the fracture network defined from a centered extraction of the initial fracture network
bool DFNExtractionNotCentered(NetworkMeshes& init_net_mesh,double extract_length){
	// 1. Variables
	NetworkMeshes new_net_mesh;vector<FractureMesh> init_meshes=init_net_mesh.meshes;
	FractureMesh current_fract,new_fract;Segment2D current_seg,new_seg;
	double init_length_x=init_net_mesh.domain.domain_size_x(),init_length_y=init_net_mesh.domain.domain_size_y();
	// 2. Definition of the new domain for a centered extraction
	Domain init_domain=init_net_mesh.domain;
	new_net_mesh.domain=Domain(pointcpp<double>(0,0),pointcpp<double>(extract_length,extract_length));
	// 3. Definition of the fracture inside the new domain
	bool fracture=false;
	for (int i=0;(unsigned)i<init_meshes.size();i++){
		// 3.1. Current fracture definition
		current_fract=init_meshes[i];current_seg=current_fract.ReturnSegment();
		if (new_net_mesh.domain.SegmentIntersectDomain(current_seg,new_seg)){
			new_fract=FractureMesh(new_seg,current_fract.aperture,current_fract.conductivity);
			if (new_fract.ReturnLength()>0){AddFracture(new_fract,new_net_mesh);}
			fracture=true;
		}
	}
	// 4. Intersection computation in the new domain
	ComputeIntersections(new_net_mesh);
	init_net_mesh=new_net_mesh;
	return fracture;
}

// Return the fracture network defined from a centered extraction of the initial fracture network
bool DFNExtractionNotCenteredBis(NetworkMeshes& init_net_mesh,double positionx1, double positionx2,double positiony1, double positiony2){
	// 1. Variables
	NetworkMeshes new_net_mesh;vector<FractureMesh> init_meshes=init_net_mesh.meshes;
	FractureMesh current_fract,new_fract;Segment2D current_seg,new_seg;
	double init_length_x=init_net_mesh.domain.domain_size_x(),init_length_y=init_net_mesh.domain.domain_size_y();
	// 2. Definition of the new domain for a centered extraction
	Domain init_domain=init_net_mesh.domain;
	new_net_mesh.domain=Domain(pointcpp<double>(positionx1,positiony1),pointcpp<double>(positionx2,positiony2));
	double L_trans_x=positionx1,L_trans_y=positiony1;
	// 3. Definition of the fracture inside the new domain
	bool fracture=false;
	for (int i=0;i<init_meshes.size();i++){
		// 3.1. Current fracture definition
		current_fract=init_meshes[i];current_seg=current_fract.ReturnSegment();
		if (new_net_mesh.domain.SegmentIntersectDomain(current_seg,new_seg)){
			new_seg=Segment2D(Translate(new_seg.source(),-L_trans_x,-L_trans_y),Translate(new_seg.target(),-L_trans_x,-L_trans_y));
			new_fract=FractureMesh(new_seg,current_fract.aperture,current_fract.conductivity);
			if (new_fract.ReturnLength()>0){AddFracture(new_fract,new_net_mesh);}
			fracture=true;
		}
	}
	// 4. Intersection computation in the new domain
	ComputeIntersections(new_net_mesh);
	init_net_mesh=new_net_mesh;
	return fracture;
}

// return a fracture network where we only keep the fractures that are located above y_cut
void DFNCutting(NetworkMeshes& net_mesh,double y_cut){
	NetworkMeshes new_net_mesh;new_net_mesh.domain=net_mesh.domain;
	// Add a parallel fracture at the cutting position
	FractureMesh new_fract(FluxPoint2D(net_mesh.domain.min_pt.i,y_cut),FluxPoint2D(net_mesh.domain.max_pt.i,y_cut),-1,-1);
	AddFracture(new_fract,net_mesh);
	ComputeIntersections(net_mesh);
	// Loop on the fractures to remove the fractures above y_cut and the added parallel fractur
	for (vector<FractureMesh>::iterator it=net_mesh.meshes.begin();it!=net_mesh.meshes.end();it++){
		if ((it->p_ori.p.y()>=y_cut-EPSILON&&it->p_tar.p.y()>=y_cut-EPSILON)&&(abs(it->p_ori.p.y()-y_cut)>EPSILON||abs(it->p_tar.p.y()-y_cut)>EPSILON)){
			AddFracture(*it,new_net_mesh);
		}
	}
	net_mesh=new_net_mesh;
}
