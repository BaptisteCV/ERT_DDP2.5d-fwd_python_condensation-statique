/*
 * Domain.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#include "Domain.h"
#include "../Utilitaries/Constantes.h"
#include "../Utilitaries/Storage/UtilitaryStorage.h"

using namespace std;


//////////////////////////
/* Functions for Domain */
//////////////////////////
Domain::Domain(){}
Domain::~Domain(){}
Domain::Domain(double Lx, double Ly){
	min_pt = pointcpp<double>(0,0);
	max_pt = pointcpp<double>(Lx,Ly);
}
Domain::Domain(pointcpp<double> min_,pointcpp<double> max_){
	min_pt=min_;
	max_pt=max_;
}
bool Domain::on_input_limit(pointcpp<double> pt){
	if ((pt.j>max_pt.j)||(pt.j<min_pt.j)){
		cout << "WARNING in on_input_limit (Domain.cpp): point out of the domain" << endl;
		return false;
	}
	if (abs(pt.i-min_pt.i)<EPSILON){return true;}
	return false;
}

bool Domain::on_output_limit(pointcpp<double> pt){
	if ((pt.j>max_pt.j)||(pt.j<min_pt.j)){
		cout << "WARNING in on_input_limit (Domain.cpp): point out of the domain" << endl;
		return false;
	}
	if (abs(pt.i-max_pt.i)<EPSILON){return true;}
	return false;
}

string Domain::ReturnBorder(CgalPoint2D pt){
	if (fabs(pt.x()-min_pt.i)<EPSILON){return LEFT_BORDER;}
	if (fabs(pt.x()-max_pt.i)<EPSILON){return RIGHT_BORDER;}
	if (fabs(pt.y()-min_pt.j)<EPSILON){return BOTTOM_BORDER;}
	if (fabs(pt.y()-max_pt.j)<EPSILON){return TOP_BORDER;}
	return NO_BORDER;
}

double Domain::ReturnDeltaX(int Nx){return domain_size_x()/Nx;}
double Domain::ReturnDeltaY(int Ny){return domain_size_y()/Ny;}

bool Domain::IsInDomain(CgalPoint2D pt){
	if (pt.x()+EPSILON>min_pt.i&&pt.x()-EPSILON<max_pt.i&&pt.y()+EPSILON>min_pt.j&&pt.y()-EPSILON<max_pt.j){
		return true;
	}
	return false;
}

bool Domain::IsInDomain(Segment2D seg){
	CgalPoint2D pt1=seg.source(),pt2=seg.target();
	if (IsInDomain(pt1)&&IsInDomain(pt2)){return true;}
	return false;
}

Segment2D Domain::ReturnBorder(string border){
	if (border==LEFT_BORDER){return Segment2D(CgalPoint2D(min_pt.i,min_pt.j),CgalPoint2D(min_pt.i,max_pt.j));}
	if (border==RIGHT_BORDER){return Segment2D(CgalPoint2D(max_pt.i,min_pt.j),CgalPoint2D(max_pt.i,max_pt.j));}
	if (border==BOTTOM_BORDER){return Segment2D(CgalPoint2D(min_pt.i,min_pt.j),CgalPoint2D(max_pt.i,min_pt.j));}
	if (border==TOP_BORDER){return Segment2D(CgalPoint2D(min_pt.i,max_pt.j),CgalPoint2D(max_pt.i,max_pt.j));}
	cout << "WARNING in ReturnBorder (Domain.cpp): border not defined" << endl;
	return Segment2D();
}

// Return true if seg intersects one of domain borders
// and return the intersection inter
bool Domain::IntersectionBorders(Segment2D seg,CgalPoint2D& inter1,CgalPoint2D& inter2){
	bool intersection=false;
	inter1=CgalPoint2D(NOT_DEFINED,NOT_DEFINED),inter2=CgalPoint2D(NOT_DEFINED,NOT_DEFINED);CgalPoint2D inter;
	// study intersection for each border
	if (Intersect(seg,ReturnBorder(LEFT_BORDER),inter1)){intersection=true;}
	if (Intersect(seg,ReturnBorder(RIGHT_BORDER),inter)){
		if (intersection){inter2=inter;}
		else{inter1=inter;intersection=true;}
	}
	if (Intersect(seg,ReturnBorder(BOTTOM_BORDER),inter)){
		if (intersection){inter2=inter;}
		else{inter1=inter;intersection=true;}
	}
	if (Intersect(seg,ReturnBorder(TOP_BORDER),inter)){
		if (intersection){inter2=inter;}
		else{inter1=inter;intersection=true;}
	}
	return intersection;
}

bool Domain::SegmentIntersectDomain(Segment2D seg1,Segment2D& seg2){
	// 1. Segment is totally inside the domain
	if (IsInDomain(seg1)){seg2=seg1;return true;}
	// 2. Segment intersects the domain
	CgalPoint2D inter1,inter2;
	if (IntersectionBorders(seg1,inter1,inter2)){
		// 2.1. Segment goes through the domain
		if (inter2.x()!=NOT_DEFINED&&inter2.y()!=NOT_DEFINED){
			seg2=Segment2D(inter1,inter2);return true;
		}
		// 2.2. Segment is partially inside the domain
		if (IsInDomain(seg1.source())){seg2=Segment2D(seg1.source(),inter1);return true;}
		if (IsInDomain(seg1.target())){seg2=Segment2D(seg1.target(),inter1);return true;}
	}
	// 3. Segment is not in the domain
	return false;
}

// Return node position from its global position pt to its border position (distance from the bottom or the left border)
// Where the returned position is from border1 to border2
double Domain::ReturnBorderCoordinate(CgalPoint2D pt,double& border_dist,string& border1,string& border2){
	string border=ReturnBorder(pt);
	if (border==BOTTOM_BORDER||border==TOP_BORDER){
		border1=LEFT_BORDER;border2=RIGHT_BORDER;
		border_dist=domain_size_x();
		return pt.x();}
	else if (border==LEFT_BORDER||border==RIGHT_BORDER){
		border1=BOTTOM_BORDER;border2=TOP_BORDER;
		border_dist=domain_size_y();
		return pt.y();
	}
	else {cout << "WARNING in ReturnCurvCoordinate (EllipseComputation.cpp): undefined border" << endl;}
	return -1;
}

void Domain::PrintDomain(){
	print(min_pt);
	print(max_pt);
}

pointcpp<double> Translate(pointcpp<double> pt,double Lx,double Ly){
	return pointcpp<double>(pt.i+Lx,pt.j+Ly);
}

void ReturnOppositeBorders(string border,string& border1,string& border2){
	if (border==LEFT_BORDER||border==RIGHT_BORDER){border1=BOTTOM_BORDER;border2=TOP_BORDER;}
	else if (border==TOP_BORDER||border==BOTTOM_BORDER){border1=LEFT_BORDER;border2=RIGHT_BORDER;}
	else{cout << "WARNING in ReturnOppositeBorders (Domain.cpp): border not defined" << endl;}
}

//////////////////////////////
/* Functions for SubDomains */
//////////////////////////////
DomainMap DefineSubDomains(Domain domain,int Nx,int Ny){
	DomainMap SubDomains;
	double Lx=domain.max_pt.i-domain.min_pt.i,Ly=domain.max_pt.j-domain.min_pt.j;
	double delta_x=Lx/Nx,delta_y=Ly/Ny;
	for (int i=0;i<Nx;i++){
		for (int j=0;j<Ny;j++){
			SubDomains[return_index(i,j,Ny)]=Domain(pointcpp<double>(domain.min_pt.i+i*delta_x,domain.min_pt.j+j*delta_y),pointcpp<double>(domain.min_pt.i+(i+1)*delta_x,domain.min_pt.j+(j+1)*delta_y));
		}
	}
	return SubDomains;
}

