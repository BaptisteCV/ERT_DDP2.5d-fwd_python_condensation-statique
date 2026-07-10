/*
 * ExchangeFunctions.cpp
 *
 *  Created on: May 27, 2013
 *      Author: delphineroubinet
 */

#include "ExchangeFunctions.h"
#include "../Utilitaries/Constantes.h"
#include "../Utilitaries/Math/Operations.h"

using namespace std;

// Compute the exchange coefficient considering the average distance from the fractures to any position in the block [Li,2008]
// Required the subnetwork included in the block
double ReturnAlphaBlock(double cond_mat,double delta_x,double delta_y,NetworkMeshes subnetwork){

	if (subnetwork.meshes.size()==0){return 0;}

	double av_cond_fract,cum_fract_length;
	double av_dist=ReturnAvDist(cond_mat,delta_x,delta_y,subnetwork,av_cond_fract,cum_fract_length);
	if (av_dist==0){return 0;}
	av_cond_fract=av_cond_fract/cum_fract_length;
	double alpha=min(cond_mat,av_cond_fract)/av_dist;
	//double alpha=0;
	return alpha;
}

double ReturnAvDist(double cond_mat,double delta_x,double delta_y,NetworkMeshes subnetwork,double& av_cond_fract,double& cum_fract_length){
	// 0. Variables
	double av_dist=0;
	av_cond_fract=0,cum_fract_length=0;
	int N_seg=100;CgalVector2D vect;CgalPoint2D pt_ori_fract,pt_current,inter1,inter2;
	//double fract_angle,L=subnetwork.domain.domain_size_x(),dist1,dist2;
	//if (fabs(subnetwork.domain.domain_size_x()-subnetwork.domain.domain_size_y())>EPSILON){cout << "WARNING in ReturnAlphaBlock (ExchangeFunctions.cpp): option not implemented" << endl;}
	double fract_angle,dist1,dist2;
	double L_max=max(subnetwork.domain.domain_size_x(),subnetwork.domain.domain_size_y());	// use only to get a distance large enough to define an orthogonal line which intersects the subnetwork borders
	// 1. Computation of the average distance for each fracture
	for (vector<FractureMesh>::iterator it1=subnetwork.meshes.begin();it1!=subnetwork.meshes.end();it1++){
		// 1.1. Each fracture is discretized in small segments
		pt_ori_fract=it1->p_ori.p;
		vect=CgalVector2D(pt_ori_fract,it1->p_tar.p);
		fract_angle=it1->ReturnSegment().get_orientation();	//angle of the current fracture
		fract_angle+=PI/2;	// modification to determine the orthogonal segments
		// 1.2. Computation of the averaged distance for each segment
		for (int i=0;i<N_seg;i++){
			// definition of the middle of the studied segment
			pt_current=CgalPoint2D(pt_ori_fract.x()+vect.x()/N_seg*(i+0.5),pt_ori_fract.y()+vect.y()/N_seg*(i+0.5));
			// definition of the orthogonal segment
			CgalPoint2D M_end1 = CgalPoint2D(pt_current.x()+2*L_max*std::cos(fract_angle),pt_current.y()+2*L_max*std::sin(fract_angle));
			CgalPoint2D M_end2 = CgalPoint2D(pt_current.x()-2*L_max*std::cos(fract_angle),pt_current.y()-2*L_max*std::sin(fract_angle));
			Segment2D seg_ortho(M_end1,M_end2);
			// definition of the intersections with the domain
			if (!subnetwork.domain.IntersectionBorders(seg_ortho,inter1,inter2)){cout << "WARNING in ReturnAlphaBlock (ExchangeFunctions.cpp): intersections with domain not found" << endl;}
			// definition of the maximum distances from the studied segment to the domain borders
			dist1=distance_2D(pt_current,inter1);dist2=distance_2D(pt_current,inter2);
			av_dist+=0.5*(dist1*dist1+dist2*dist2)/(dist1+dist2)*it1->ReturnLength()/(double)(N_seg);

			//print(it1->p_ori.p);
			//print(it1->p_tar.p);

			//cout << it1->ReturnLength() << endl;
		}
		//av_dist=av_dist/it1->ReturnLength();
		// 1.3. Computation of the averaged fracture conductivity
		av_cond_fract+=it1->conductivity*it1->ReturnLength();
		cum_fract_length+=it1->ReturnLength();
	}
	av_dist=av_dist/(cum_fract_length);
	return av_dist;
}
