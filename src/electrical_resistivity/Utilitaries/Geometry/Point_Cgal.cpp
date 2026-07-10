/*
 * Point_Cgal.cpp
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */


#include "Point_Cgal.h"
#include "../Constantes.h"

using namespace std;


/**determine if two points are identic*/
int identic(const CgalPoint2D p1, const CgalPoint2D p2, const double eps){
	double x1 = p1.x(), x2 = p2.x(), y1 = p1.y(), y2 = p2.y();
	if((fabs(x1-x2)<eps)&&(fabs(y1-y2)<eps))  return 1; else return 0;
}

/**distance between two points*/
double distance_2D(const CgalPoint2D p1, const CgalPoint2D p2){
	double value = (p1.x()-p2.x())*(p1.x()-p2.x())+(p1.y()-p2.y())*(p1.y()-p2.y());
	return std::sqrt(value);
}

void print(CgalPoint2D pt){
	cout << "(" << pt.x() << "," << pt.y() << ")" << endl;
}


CgalPoint2D rotation(CgalPoint2D pt,double angle){
	double r=sqrt(pt.x()*pt.x()+pt.y()*pt.y());
	double theta=0;
	if (fabs(pt.x())<EPSILON){
		if (pt.y()>0){theta=0.5*PI;}
		else if (pt.y()<0){theta=-0.5*PI;}
	}
	else if (pt.x()<0){
		if (pt.y()>=0){theta=atan(pt.y()/pt.x())+PI;}
		else if (pt.y()<0){theta=atan(pt.y()/pt.x())-PI;}
	}
	else if (pt.x()>0){theta=atan(pt.y()/pt.x());}
	else{cout << "WARNING in rotation (Point_Cgal.cpp): case not implemented" << endl;}
	//return CgalPoint2D(r*cos(theta),r*sin(theta));
	return CgalPoint2D(r*cos(theta+angle),r*sin(theta+angle));
}

CgalPoint2D Translate(CgalPoint2D pt,double L_trans_x,double L_trans_y){
	return CgalPoint2D(pt.x()+L_trans_x,pt.y()+L_trans_y);
}
