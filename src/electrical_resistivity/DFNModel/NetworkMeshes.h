/*
 * NetworkMeshes.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef NETWORKMESHES_H_
#define NETWORKMESHES_H_

#include <map>
#include <set>
#include "FractureMesh.h"
#include "../DomainDefinition/Domain.h"
#include "../DomainDefinition/ElectricProperties.h"
#include "../InputOutput/Parameters.h"


//map<node index,map<neighboring node,electric properties of the corresponding segment> >
typedef std::map<int,ElectricProperties> ElecPropMap;
typedef std::map<int,ElecPropMap> NodeFracturesMap;
typedef std::map<int,std::string> BordersMap;
typedef std::map<int,CgalPoint2D> NodesMap;

/************************************/
// Functions for FractureFamilyDesc //
/************************************/
class FractureFamilyDesc{
public:
	int fracture_nb;
	double fracture_angle;
	double aperture;
	double conductivity;
	double fract_spacing;
public:
	FractureFamilyDesc();
	virtual ~FractureFamilyDesc();
	FractureFamilyDesc(int,double,double,double,double);
	double ReturnAngleRadian();
};

typedef std::map<int,FractureFamilyDesc> FractureFamilyMap;


/*******************************/
// Functions for NetworkMeshes //
/*******************************/
class NetworkMeshes{
public:
	Domain domain;
	std::vector<FractureMesh> meshes;	// list of meshes composing the fracture network
	int cpt_inter;
	BordersMap border_map;	// map<node_index,domain_border>
	NodesMap nodes_map;	// map<node_index,node_coordinate>
	std::set<int> inter_list;
	bool inter_grid;
public:
	NetworkMeshes();
	virtual ~NetworkMeshes();
	NetworkMeshes(Domain,Parameters&,std::string,std::string option="intersection");
	NetworkMeshes(Domain,int,BordersMap);
	NetworkMeshes(Domain,int,std::set<int>);
	void NetworkMeshesFromFile(std::string,std::string);
	FractureMesh return_mesh(int);
	double return_conductivity_from_node(int);
	double return_aperture_from_node(int);
	FractureMesh operator[](int index){return meshes[index];}
	NodeFracturesMap ReturnNodeFracturesMap();
	void Translate(double,double,double);
	int ReturnIndex(CgalPoint2D);
	CgalPoint2D ReturnNodeCoord(int);
	void print_DFN();
	void print_DFN_conductivity();
	void print_border_map();
	double ReturnDomainLength();
	double ComputeFMExchange();
	// Functions for "backbone" computation
	NetworkMeshes return_backbone(Parameters,double,std::string option="gradient");
	std::set<int> return_connected_nodes();
	std::set<int> return_nodes_from_velocity(Parameters,double,std::string option="gradient");
	NetworkMeshes return_network_from_nodes(std::set<int>);
	double ReturnFlowVelocityNormal(int);
	void EvaluateFlowVelocities(ublas_vector);
	std::vector<FractureMesh> return_mesh_from_node_ori(int);
};

void NodeInsertionMap(int,int,FractureMesh,NodeFracturesMap&,double);
void add_neigh_nodes(std::set<int>&,NodeFracturesMap,int);

/*****************************/
// Functions for subnetworks //
/*****************************/
typedef std::map<int,NetworkMeshes> SubNetworkMap;
NodeFracturesMap ReturnNodeFracturesMap(SubNetworkMap);
int ReturnNbNodes(SubNetworkMap);
int ReturnBlockIndex(SubNetworkMap,int);

#endif /* NETWORKMESHES_H_ */
