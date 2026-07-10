/*
 * NetworkMeshes.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#include "NetworkMeshes.h"
#include "DFNGeneration.h"
#include "DFNComputation.h"
#include "../Utilitaries/Geometry/Critical_Operations.h"
#include "../DDPModel/DDPDFNFunctions.h"
#include <fstream>

using namespace std;


/************************************/
// Functions for FractureFamilyDesc //
/************************************/
FractureFamilyDesc::FractureFamilyDesc(){}
FractureFamilyDesc::~FractureFamilyDesc(){}
FractureFamilyDesc::FractureFamilyDesc(int fracture_nb_,double fracture_angle_,double aperture_,double conductivity_,double fract_spacing_){
	fracture_nb=fracture_nb_;
	fracture_angle=fracture_angle_;
	aperture=aperture_;
	conductivity=conductivity_;
	fract_spacing=fract_spacing_;
}
double FractureFamilyDesc::ReturnAngleRadian(){
	return PI*fracture_angle/180;
}

/*******************************/
// Functions for NetworkMeshes //
/*******************************/
NetworkMeshes::NetworkMeshes(){cpt_inter=-1;inter_grid=false;}
NetworkMeshes::~NetworkMeshes(){}
NetworkMeshes::NetworkMeshes(Domain domain_,int cpt_inter_,BordersMap border_map_){
	cpt_inter=cpt_inter_;
	border_map=border_map_;
	domain=domain_;
	inter_grid=false;
}
NetworkMeshes::NetworkMeshes(Domain domain_,int cpt_inter_,set<int> inter_list_){
	domain=domain_;
	cpt_inter=cpt_inter_;
	inter_list=inter_list_;
	inter_grid=false;
}
NetworkMeshes::NetworkMeshes(Domain domain_,Parameters& param,string option,string intersection_option){
	domain=domain_;inter_grid=false;
	string file_name=param.files_param.input_path+param.files_param.file_name_DFN;
	ifstream fichier(file_name.c_str());
	string generation_option;
	fichier >> generation_option;
	param.simu_param.DFN_param.generation_option=generation_option;
	if (generation_option==PARALLEL){
		int family_nb;
		fichier >> family_nb;
		FractureFamilyMap fract_family_map;
		for (int i=0;i<family_nb;i++){
			int fracture_nb;double angle,aperture,conductivity,fract_spacing;
			fichier >> fracture_nb >> angle >> aperture >> conductivity >> fract_spacing;
			fract_family_map[i]=FractureFamilyDesc(fracture_nb,angle,aperture,conductivity,fract_spacing);
		}
		NetworkMeshes net_meshes=DFNGenerationParallel(domain,fract_family_map,intersection_option);
		*this=net_meshes;
	}
	else if (generation_option==PARALLEL_BORDER){
		double aperture,conductivity;
		fichier >> aperture >> conductivity;
		NetworkMeshes net_meshes=DFNGenerationParallelBorder(domain,aperture,conductivity,intersection_option);
		param.simu_param.DFN_param.fract_aperture=aperture;
		param.simu_param.DFN_param.fract_conductivity=conductivity;
		*this=net_meshes;
	}
	else if (generation_option==SINGLE){
		double aperture,conductivity,length;int seed,nb_fract;
		fichier >> aperture >> conductivity >> length >> seed >> nb_fract;
		param.simu_param.DFN_param.fract_aperture=aperture;
		param.simu_param.DFN_param.fract_conductivity=conductivity;
		param.simu_param.DFN_param.fract_length=length;
		param.simu_param.DFN_param.seed=seed;
		param.simu_param.DFN_param.nb_fract=nb_fract;
		NetworkMeshes net_meshes=DFNGenerationRandomSingle(domain,aperture,conductivity,length,seed,nb_fract,intersection_option);
		*this=net_meshes;
	}
	else if (generation_option==SIERPINSKI){
		int nb_division,level_division,seed;double fract_aperture,fract_conductivity;
		fichier >> nb_division >> level_division >> seed >> fract_aperture >> fract_conductivity;
		NetworkMeshes net_meshes=DFNGenerationSierpinski(domain,fract_aperture,fract_conductivity,nb_division,level_division,seed,option,intersection_option);
		*this=net_meshes;
	}
	else if (generation_option==SUGAR_BOX){
		double fract_spacing,fract_aperture,fract_conductivity;
		fichier >> fract_spacing >> fract_aperture >> fract_conductivity;
		NetworkMeshes net_meshes=DFNGenerationSugarBox(domain,fract_spacing,fract_aperture,fract_conductivity);
		*this=net_meshes;
	}
	else if (generation_option==RANDOM){
		double power_law_exp,percolation_param,fract_aperture,fract_conductivity;int seed;
		fichier >> power_law_exp >> percolation_param >> fract_aperture >> fract_conductivity >> seed;
		param.simu_param.DFN_param.power_law_exp=power_law_exp;
		param.simu_param.DFN_param.percolation_param=percolation_param;
		param.simu_param.DFN_param.fract_aperture=fract_aperture;
		param.simu_param.DFN_param.fract_conductivity=fract_conductivity;
		param.simu_param.DFN_param.seed=seed;
		NetworkMeshes net_meshes=DFNGenerationRandom(domain,power_law_exp,percolation_param,fract_aperture,fract_conductivity,seed,intersection_option);
		*this=net_meshes;
	}
	else if (generation_option=="POWER_LAW_ANGLE"){
		double power_law_exp,percolation_param,fract_aperture,fract_conductivity,coeff_theta1,coeff_theta2,l_min;int seed;
		fichier >> power_law_exp >> percolation_param >> coeff_theta1 >> coeff_theta2 >> fract_aperture >> fract_conductivity >> seed >> l_min;
		param.simu_param.DFN_param.power_law_exp=power_law_exp;
		param.simu_param.DFN_param.percolation_param=percolation_param;
		param.simu_param.DFN_param.fract_aperture=fract_aperture;
		param.simu_param.DFN_param.fract_conductivity=fract_conductivity;
		param.simu_param.DFN_param.seed=seed;
		param.simu_param.DFN_param.seed=l_min;

		NetworkMeshes net_mesh;net_mesh.domain=domain;FractureMesh frac_mesh;
		Segment2D seg;double density=0,fract_length,angle;CgalPoint2D fract_center;
		double Lx=domain.domain_size_x(),Ly=domain.domain_size_y();
		srand(seed);
		int cpt_fract=0;
		while (density<percolation_param){
			cpt_fract++;
			fract_length=l_min*pow(Uniform(),1./(1-power_law_exp));
			fract_center=CgalPoint2D(Uniform()*Lx,Uniform()*Ly);
 			//angle of the fracture
            if ((double)cpt_fract/2.==(int)cpt_fract/2){angle=coeff_theta1*PI/180;}
            else {angle=coeff_theta2*PI/180;}   			
			seg=Segment2D(fract_center,fract_length,angle);
			domain.SegmentIntersectDomain(seg,seg);
			frac_mesh=FractureMesh(seg,fract_aperture,fract_conductivity);
			AddFracture(frac_mesh,net_mesh);
			density+=seg.length*seg.length/(Lx*Ly);
		}
		if (intersection_option=="intersection"){
			ComputeIntersections(net_mesh);
		}

		*this=net_mesh;
	}
    else if (generation_option=="FRACTAL"){
        FractureMesh frac_mesh; Segment2D seg1,seg2;double fract_length,angle;CgalPoint2D fract_center;
        double Lx=domain.domain_size_x(),Ly=domain.domain_size_y();
        //if (Lx!=Ly){cout << "PROBLEM in NetworkMeshes (NetworkMeshes.cpp): Lx and Ly must be equal" << endl; exit(0);}

        double r_min,C_density,D_dim,fract_conductivity,fract_aperture,coeff_theta1,coeff_theta2;
		int seed;
        fichier >> r_min >> C_density >> D_dim >> coeff_theta1 >> coeff_theta2 >> fract_aperture >> fract_conductivity >> seed;
        
        int nb_fract=std::ceil(C_density*std::pow(r_min/Lx,-(double)D_dim));
        int cpt_fract=0;

        NetworkMeshes net_mesh=*this;
        srand(seed);
        while (cpt_fract<nb_fract){
                cpt_fract++;
                // center of the fracture
                //fract_center=CgalPoint2D(-0.5*Lx+Uniform()*Lx,-0.5*Ly+Uniform()*Ly);
				fract_center=CgalPoint2D(Uniform()*Lx,Uniform()*Ly);
                // length of the fracture
                fract_length=std::pow(C_density/(double)cpt_fract,1./(double)D_dim)*Lx;
                //angle of the fracture
                if ((double)cpt_fract/2.==(int)cpt_fract/2){angle=coeff_theta1*PI/180;}
                else {angle=coeff_theta2*PI/180;}                // definition of a segment
                seg1=Segment2D(fract_center,fract_length,angle);
				domain.SegmentIntersectDomain(seg1,seg2);
                frac_mesh=FractureMesh(seg2,fract_aperture,fract_conductivity);
                AddFracture(frac_mesh,net_mesh);
        }
        // 4. Intersections determination
        ComputeIntersections(net_mesh);
        *this=net_mesh;
    }
	else if (generation_option==DETERMINISTIC1){
		NetworkMeshes net_mesh;net_mesh.domain=domain;
		int nb_fracture,pos1,pos2;
		string str_line, x_str, y_str;
		double fract_aperture,fract_conductivity,x1,y1,x2,y2;
		Segment2D seg1,seg2;
		FractureMesh frac_mesh;
		fichier >> nb_fracture;
		for (int i=0;i<nb_fracture;i++){
			fichier >> x1 >> y1 >> x2 >> y2 >> fract_aperture >> fract_conductivity;
			FluxPoint2D p_origin(x1,y1),p_target(x2,y2);
			Segment2D seg1(p_origin.p,p_target.p);
			//cout << seg1.source().x() << " " << seg1.source().y() << " " << seg1.target().x() << " " << seg1.target().y() << endl;
			if (!domain.SegmentIntersectDomain(seg1,seg2)){
				cout << "WARNING in NetworkMeshes.cpp: fracture out of the domain" << endl;
			}
			else{
				//cout << seg2.source().x() << " " << seg2.source().y() << " " << seg2.target().x() << " " << seg2.target().y() << endl;
				frac_mesh=FractureMesh(p_origin,p_target,fract_aperture,fract_conductivity);
				AddFracture(frac_mesh,net_mesh);
			}
		}
		if (intersection_option=="intersection"){
			ComputeIntersections(net_mesh);
		}
		*this=net_mesh;
	}
	else if (generation_option==DETERMINISTIC2){
		NetworkMeshes net_mesh;net_mesh.domain=domain;
		int nb_fracture,pos1,pos2;
		string str_line, x_str, y_str;
		double fract_aperture,fract_angle,fract_length,fract_cond;
		FractureMesh frac_mesh;
		fichier >> nb_fracture;
		for (int i=0;i<nb_fracture;i++){
			// definition of the segement's center
			fichier >> str_line;
			pos1 = 0;
			pos1 = str_line.find("(");
			pos2 = str_line.find(",",pos1);
			x_str = str_line.substr(pos1+1,pos2-1);
			pos1 = pos2;
			pos2 = str_line.find(")",pos1);
			y_str = str_line.substr(pos1+1,pos2-pos1-1);
			CgalPoint2D center(atof(x_str.c_str()),atof(y_str.c_str()));

			//cout.precision(30);
			//cout << center.x() << endl;

			// angle
			fichier >> fract_angle;
			// aperture
			fichier >> fract_aperture;
			// length
			fichier >> fract_length;
			// conductivity
			fichier >> fract_cond;
			// fracture definition
			Segment2D seg(center,fract_length,fract_angle);
			frac_mesh=FractureMesh(seg,fract_aperture,fract_cond);
			AddFracture(frac_mesh,net_mesh);
		}
		if (intersection_option=="intersection"){
			ComputeIntersections(net_mesh);
		}
		*this=net_mesh;
	}
	else if (generation_option==DETERMINISTIC3){
                NetworkMeshes net_mesh;net_mesh.domain=domain;
                int nb_fracture,pos1,pos2,nb_fract_cond;
                string str_line, x_str, y_str;
                double fract_aperture,x1,y1,x2,y2,value;
                FractureMesh frac_mesh;
                fichier >> nb_fracture;
                for (int i=0;i<nb_fracture;i++){
			NetworkMeshes net_mesh_temp;net_mesh_temp.domain=domain;
			vector<double> fract_cond_vec;
                        fichier >> x1 >> y1 >> x2 >> y2 >> fract_aperture >> nb_fract_cond;
			for (int j=0;j<nb_fract_cond;j++){
				fichier >> value;
				fract_cond_vec.push_back(value);
			}
                        FluxPoint2D p_origin(x1,y1),p_target(x2,y2);
                        frac_mesh=FractureMesh(p_origin,p_target,fract_aperture,-1);
			AddFracture(frac_mesh,net_mesh_temp);
			// compute the intersections between the fracture and the grid to define the subfractures in net_mesh_temp
			FracturesIntersectionGrid(param.simu_param.DDP_param.Nx,param.simu_param.DDP_param.Ny,net_mesh_temp);
			// Add the subfractures to the main network with heterogeneous conductivity properties
			if (fract_cond_vec.size()!=net_mesh_temp.meshes.size()){cout << "WARNING in NetworkMeshes (NetworkMeshes.cpp): problem in the heterogeneous conductivity definition" << endl; cout << "fract_cond_vec.size() = " << fract_cond_vec.size() << ", net_mesh_temp.meshes.size() = " << net_mesh_temp.meshes.size() << endl;exit(0);}
			for (int j=0;j<nb_fract_cond;j++){
				frac_mesh=FractureMesh(net_mesh_temp.meshes[j].p_ori,net_mesh_temp.meshes[j].p_tar,fract_aperture,fract_cond_vec[j]);
				AddFracture(frac_mesh,net_mesh);
			}
                }
                if (intersection_option=="intersection"){
                        ComputeIntersections(net_mesh);
                }
		inter_grid=true;
                *this=net_mesh;
        }
	else{cout << "WARNING in NetworkMeshes (NetworkMeshes.cpp): generation option not implemented" << endl;}
	fichier.close();
}

void NetworkMeshes::NetworkMeshesFromFile(string code_path, string file_name_domain){
	// local variables
	string str_line, x_str, y_str;
	int pos1, pos2, mesh_index, nb_segment, index;
	double velocity, aperture;
	// Domain definition from files
	string file_name = code_path + "/Input/DFNFiles/";
	file_name += file_name_domain;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		// number of segments to describe the fracture network
		fichier >> nb_segment;
		meshes = vector<FractureMesh>(nb_segment);
		// each line describes a segment (or mesh) composing the fracture network
		fichier >> str_line;
		for (int i=0;i<nb_segment;i++){
			// definition of the segement's first extremity (coordinate and index)
			pos1 = 0;
			pos1 = str_line.find("(");
			pos2 = str_line.find(",",pos1);
			x_str = str_line.substr(pos1+1,pos2-1);
			pos1 = pos2;
			pos2 = str_line.find(")",pos1);
			y_str = str_line.substr(pos1+1,pos2-pos1-1);
			fichier >> index;
			FluxPoint2D p_origin(atof(x_str.c_str()),atof(y_str.c_str()),index);
			// definition of the segement's second extremity
			fichier >> str_line;
			pos1 = 0;
			pos1 = str_line.find("(");
			pos2 = str_line.find(",",pos1);
			x_str = str_line.substr(pos1+1,pos2-1);
			pos1 = pos2;
			pos2 = str_line.find(")",pos1);
			y_str = str_line.substr(pos1+1,pos2-pos1-1);
			fichier >> index;
			FluxPoint2D p_target(atof(x_str.c_str()),atof(y_str.c_str()),index);
			// flow velocity of the studied segment
			fichier >> velocity;
			// segment aperture
			fichier >> aperture;
			// studied segment numbering
			fichier >> mesh_index;
			// numbering of the neighboring meshes
			vector<int> neigh_meshes;
			fichier >> str_line;
			while ((str_line.find("(")==str_line.npos)&&(!fichier.eof())){
				neigh_meshes.push_back(atoi(str_line.c_str()));
				fichier >> str_line;
			}
			// mesh affectation
			FractureMesh mesh(p_origin,p_target,velocity,aperture,mesh_index,neigh_meshes);
			this->meshes[mesh_index] = mesh;
		}}
	else{cout << "WARNING in Parameters (Parameters.cpp) : unknown file" << endl;}
	fichier.close();
}

FractureMesh NetworkMeshes::return_mesh(int mesh_ind){
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		if (it->mesh_index == mesh_ind){
			return *it;
		}
	}
	cout << "WARNING in return_mesh (NetworkMeshes.cpp): mesh number " << mesh_ind << " not found" << endl;
	return FractureMesh();
}

double NetworkMeshes::return_conductivity_from_node(int node_index){
	double conductivity;bool cond_bool=false;
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		if (it->p_ori.index==node_index||it->p_tar.index==node_index){
			if (!cond_bool){
				conductivity=it->conductivity;cond_bool=true;
			}
			else{
				if (conductivity==it->conductivity){return conductivity;}
				else{cout << "WARNING in return_conductivity_from_node (NetworkMeshes.cpp): conductivity not constant" << endl;}
			}
		}
	}
	if (cond_bool){return conductivity;}
	cout << "WARNING in return_conductivity_from_node (NetworkMeshes.cpp): node_index not found" << endl;
	return -1;
}

double NetworkMeshes::return_aperture_from_node(int node_index){
	double aperture;bool ap_bool=false;
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		if (it->p_ori.index==node_index||it->p_tar.index==node_index){
			if (!ap_bool){
				aperture=it->aperture;ap_bool=true;
			}
			else{
				if (aperture==it->aperture){return aperture;}
				else{cout << "WARNING in return_aperture_from_node (NetworkMeshes.cpp): aperture not constant" << endl;}
			}
		}
	}
	if (ap_bool){return aperture;}
	cout << "WARNING in return_conductivity_from_node (NetworkMeshes.cpp): node_index not found" << endl;
	return -1;
}

// return the map of node index and a list of their corresponding fractures
NodeFracturesMap NetworkMeshes::ReturnNodeFracturesMap(){
	// 0. Variables
	NodeFracturesMap node_fract_map;
	// 1. Loop on the fractures (or meshes)
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		// 1.1. Study of the first extremity
		NodeInsertionMap(it->p_ori.index,it->p_tar.index,*it,node_fract_map,it->conductivity);
		// 1.2. Study of the second extremity  (if not defined already)
		NodeInsertionMap(it->p_tar.index,it->p_ori.index,*it,node_fract_map,it->conductivity);
	}
	return node_fract_map;
}

void NetworkMeshes::Translate(double trans_value,double y_min,double y_max){
	FractureMesh fract_mesh;
	FluxPoint2D p_origin,p_target;
	double value;
	for (size_t i=0;i<meshes.size();i++){
		value=meshes[i].p_ori.p.y()+trans_value;
		if (value>y_max){value-=y_max-y_min;}
		p_origin=FluxPoint2D(meshes[i].p_ori.p.x(),value,meshes[i].p_ori.index);
		value=meshes[i].p_tar.p.y()+trans_value;
		if (value>y_max){value-=y_max-y_min;}
		p_target=FluxPoint2D(meshes[i].p_tar.p.x(),value,meshes[i].p_tar.index);
		meshes[i]=FractureMesh(p_origin,p_target,meshes[i].velocity,meshes[i].aperture,meshes[i].mesh_index,meshes[i].neigh_meshes);
	}
}

int NetworkMeshes::ReturnIndex(CgalPoint2D pt){
	for (size_t i=0;i<meshes.size();i++){
		if (Points_Equal(pt,meshes[i].p_ori.p)){return meshes[i].p_ori.index;}
		if (Points_Equal(pt,meshes[i].p_tar.p)){return meshes[i].p_tar.index;}
	}
	return -1;
}

CgalPoint2D NetworkMeshes::ReturnNodeCoord(int index){
	for (size_t i=0;i<meshes.size();i++){
		if (meshes[i].p_ori.index==index){return meshes[i].p_ori.p;}
		if (meshes[i].p_tar.index==index){return meshes[i].p_tar.p;}
	}
	return CgalPoint2D(-1,-1);
}

void NetworkMeshes::print_DFN(){
	int cpt_mesh=0;
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		cout << "Mesh " << cpt_mesh << ":" << endl;
		cout << "Index origin: " << it->p_ori.index << endl;
		print(it->p_ori.p);
		cout << "Index target: " << it->p_tar.index << endl;
		print(it->p_tar.p);
		cpt_mesh++;
	}
}

void NetworkMeshes::print_DFN_conductivity(){
	int cpt_mesh=0;
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		cout << "Mesh " << cpt_mesh << ":" << endl;
		cout << "Conductivity: " << it->conductivity << endl;
		cpt_mesh++;
	}
}

void NetworkMeshes::print_border_map(){
	for (BordersMap::iterator it=border_map.begin();it!=border_map.end();it++){
		cout << it->first << " " << it->second << endl;
	}
}

double NetworkMeshes::ReturnDomainLength(){
	if (domain.domain_size_x()==domain.domain_size_y()){return domain.domain_size_x();}
	cout << "WARNING in ReturnDomainLength (NetworkMeshes.cpp): the domain is not a square" << endl;
	return -1;
}

/*****************************/
// Functions for subnetworks //
/*****************************/
int ReturnNbNodes(SubNetworkMap subnetwork_map){
	set<int> node_index;
	for (SubNetworkMap::iterator it1=subnetwork_map.begin();it1!=subnetwork_map.end();it1++){
		for (vector<FractureMesh>::iterator it2=it1->second.meshes.begin();it2!=it1->second.meshes.end();it2++){
			node_index.insert(it2->p_ori.index);node_index.insert(it2->p_tar.index);
		}
	}
	return node_index.size();
}

// return the map of node index and a list of their corresponding fractures
NodeFracturesMap ReturnNodeFracturesMap(SubNetworkMap subnetwork_map){
	// 0. Variables
	NodeFracturesMap node_fract_map;
	// 1. Loop on the fractures (or meshes)
	for (SubNetworkMap::iterator it1=subnetwork_map.begin();it1!=subnetwork_map.end();it1++){
		for (vector<FractureMesh>::iterator it2=it1->second.meshes.begin();it2!=it1->second.meshes.end();it2++){
			// 1.1. Study of the first extremity
			NodeInsertionMap(it2->p_ori.index,it2->p_tar.index,*it2,node_fract_map,it2->conductivity);
			// 1.2. Study of the second extremity  (if not defined already)
			NodeInsertionMap(it2->p_tar.index,it2->p_ori.index,*it2,node_fract_map,it2->conductivity);
		}
	}
	return node_fract_map;
}

int ReturnBlockIndex(SubNetworkMap subnetwork_map,int fract_node_index){
	for (SubNetworkMap::iterator it1=subnetwork_map.begin();it1!=subnetwork_map.end();it1++){
		for (vector<FractureMesh>::iterator it2=it1->second.meshes.begin();it2!=it1->second.meshes.end();it2++){
			if (it2->p_ori.index==fract_node_index||it2->p_tar.index==fract_node_index){
				return it1->first;
			}
		}
	}
	cout << "WARNING in ReturnBlockIndex (NetworkMeshes.cpp): matrix block not found" << endl;
	return -1;
}

/************************/
// Additional functions //
/************************/
// addition of the node in the numbering map
void NodeInsertionMap(int node1,int node2,FractureMesh fract_mesh,NodeFracturesMap& node_fract_map,double cond_fract){
	ElecPropMap prop_map;
	// if already defined
	NodeFracturesMap::iterator it1=node_fract_map.find(node1);
	if (it1!=node_fract_map.end()){prop_map=it1->second;}
	// addition of the studied segment information in the submap
	prop_map[node2]=ElectricProperties(fract_mesh,cond_fract);
	// addition of the information in the main map
	node_fract_map[node1]=prop_map;
}

// return the list of nodes connected to the domain borders
set<int> NetworkMeshes::return_connected_nodes(){
	set<int> connected_nodes;NodeFracturesMap node_fract_map=ReturnNodeFracturesMap();
	for (BordersMap::iterator it=border_map.begin();it!=border_map.end();it++){
		connected_nodes.insert(it->first);
		add_neigh_nodes(connected_nodes,node_fract_map,it->first);
	}
	return connected_nodes;
}

// returns connected_nodes which contains the nodes connected to node_index (neighbors and neighbors of neighbors...)
void add_neigh_nodes(set<int> & connected_nodes,NodeFracturesMap node_fract_map,int node_index){
	set<int> nodes_to_study;nodes_to_study.insert(node_index);
	while (!nodes_to_study.empty()){
		set<int> new_nodes_to_study;
		for (set<int>::iterator it=nodes_to_study.begin();it!=nodes_to_study.end();it++){
			for (ElecPropMap::iterator it1=node_fract_map[*it].begin();it1!=node_fract_map[*it].end();it1++){
				if (connected_nodes.find(it1->first)==connected_nodes.end()){
					connected_nodes.insert(it1->first);
					new_nodes_to_study.insert(it1->first);
				}
			}
		}
		nodes_to_study=new_nodes_to_study;
	}
}


// return a new networkmeshes from a list of nodes
NetworkMeshes NetworkMeshes::return_network_from_nodes(set<int> connected_nodes){
	NetworkMeshes new_net_mesh;
	// Domain
	new_net_mesh.domain=domain;
	// Connected node numbering
	map<int,int> node_numbering;//<new_index,old_index>
	int cpt_node=-1;

	for (set<int>::iterator it1=connected_nodes.begin();it1!=connected_nodes.end();it1++){
		cpt_node++;
		node_numbering[cpt_node]=*it1;
	}
	// nodes and mesh selection
	vector<FractureMesh> mesh_list;FractureMesh mesh;int old_index,new_index;
	int cpt_mesh=0;map<int,int>::iterator it_map;
	for (map<int,int>::iterator it1=node_numbering.begin();it1!=node_numbering.end();it1++){
		new_index=it1->first;old_index=it1->second;
		// Border map
		if (border_map.find(old_index)!=border_map.end()){new_net_mesh.border_map[new_index]=border_map[old_index];}
		// Node map
		new_net_mesh.nodes_map[new_index]=nodes_map[old_index];
		// Inter list
		if (inter_list.find(old_index)!=inter_list.end()){new_net_mesh.inter_list.insert(new_index);}
		// Meshes
		mesh_list=return_mesh_from_node_ori(old_index);
		for (vector<FractureMesh>::iterator it2=mesh_list.begin();it2!=mesh_list.end();it2++){
			if (connected_nodes.find(it2->p_ori.index)!=connected_nodes.end()&&connected_nodes.find(it2->p_tar.index)!=connected_nodes.end()){
				mesh=*it2;
				for (map<int,int>::iterator it3=node_numbering.begin();it3!=node_numbering.end();it3++){
					if (it3->second==mesh.p_ori.index){mesh.p_ori.index=it3->first;}
				}
				for (map<int,int>::iterator it3=node_numbering.begin();it3!=node_numbering.end();it3++){
					if (it3->second==mesh.p_tar.index){mesh.p_tar.index=it3->first;}
				}
				mesh.mesh_index=cpt_mesh;
				new_net_mesh.meshes.push_back(mesh);
				cpt_mesh++;
			}
		}
	}

	new_net_mesh.cpt_inter=cpt_node;
	return new_net_mesh;
}

vector<FractureMesh> NetworkMeshes::return_mesh_from_node_ori(int node_index){
	vector<FractureMesh> mesh_list;
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		if (it->p_ori.index==node_index){
			mesh_list.push_back(*it);
		}
	}
	return mesh_list;
}

// return the nodes corresponding to meshes with a flow velocity larger than a given threshold
set<int> NetworkMeshes::return_nodes_from_velocity(Parameters param,double velocity_threshold,string option){
	set<int> list_nodes;
	// compute flow velocities
	ComputeFlowVelocities(*this,param,domain,option);
	// selection of the meshes with a flow velocity larger than the velocity threshold
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		if (fabs(it->velocity)>velocity_threshold){
			list_nodes.insert(it->p_ori.index);
			list_nodes.insert(it->p_tar.index);
		}
	}
	return list_nodes;
}

NetworkMeshes NetworkMeshes::return_backbone(Parameters param,double velocity_threshold,string option){
	NetworkMeshes new_network_meshes;
	// selection of the network with nodes connected to the borders
	set<int> connected_nodes=return_connected_nodes();
	new_network_meshes=return_network_from_nodes(connected_nodes);
	// selection of the meshes with flow velocities larger than a given threshold
	set<int> nodes_from_velocity=new_network_meshes.return_nodes_from_velocity(param,velocity_threshold,option);
	new_network_meshes=new_network_meshes.return_network_from_nodes(nodes_from_velocity);
	return new_network_meshes;
}

// return the flow velocity normal to the node
double NetworkMeshes::ReturnFlowVelocityNormal(int node_index){
	if (border_map.find(node_index)==border_map.end()){
		cout << "WARNING in ReturnFlowVelocity (NetworkMEshes.cpp): the node in argument is not on the domain border" << endl;
		return -1;
	}
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		if (it->p_ori.index==node_index){return it->velocity;}
		if (it->p_tar.index==node_index){return -it->velocity;}
	}
	cout << "WARNING in ReturnFlowVelocity (NetworkMEshes.cpp): node not found" << endl;
	return -1;
}

void NetworkMeshes::EvaluateFlowVelocities(ublas_vector DFN_potential){
	FractureMesh fract_mesh;double velocity;
	for (int i=0;(unsigned)i<meshes.size();i++){
		fract_mesh=meshes[i];
		velocity=meshes[i].ReturnConductance()*(DFN_potential(fract_mesh.p_ori.index)-DFN_potential(fract_mesh.p_tar.index))/fract_mesh.ReturnLength();
		meshes[i].velocity=velocity;
	}
}
