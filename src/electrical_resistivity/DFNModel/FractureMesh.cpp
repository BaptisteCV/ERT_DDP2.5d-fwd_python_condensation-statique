/*
 * FractureMesh.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */


#include "FractureMesh.h"

using namespace std;


FractureMesh::FractureMesh(FluxPoint2D p_ori_,FluxPoint2D p_tar_,double velocity_,double aperture_,int mesh_index_,std::vector<int> neigh_meshes_){
	p_ori = p_ori_;
	p_tar = p_tar_;
	velocity = velocity_;
	aperture = aperture_;
	mesh_index = mesh_index_;
	neigh_meshes = neigh_meshes_;
}

FractureMesh::FractureMesh(FluxPoint2D p_ori_,FluxPoint2D p_tar_,double aperture_,double cond_fract_){
	p_ori = p_ori_;
	p_tar = p_tar_;
	aperture = aperture_;
	conductivity=cond_fract_;
}

FractureMesh::FractureMesh(Segment2D seg,double aperture_,double cond_fract_){
	p_ori=FluxPoint2D(seg.source());
	p_tar=FluxPoint2D(seg.target());
	aperture = aperture_;
	conductivity=cond_fract_;
}

// fracture length
double FractureMesh::ReturnLength() const{
	return CGAL::sqrt(CGAL::squared_distance(this->p_ori.p,this->p_tar.p));
}

CgalPoint2D FractureMesh::ReturnCenter(){
	Segment2D seg_mesh(p_ori.p,p_tar.p);
	return seg_mesh.center;
}

double FractureMesh::ReturnAngle(){
	Segment2D seg_mesh(p_ori.p,p_tar.p);
	return seg_mesh.get_orientation();
}

Segment2D FractureMesh::ReturnSegment(){
	Segment2D seg_mesh(p_ori.p,p_tar.p);
	return seg_mesh;
}

double FractureMesh::ReturnConductance(){
	return RHO*G*std::pow(aperture,2)/(12*MU);
}


