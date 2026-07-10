/*
 * UtilitaryStorage.cpp
 *
 *  Created on: Apr 29, 2013
 *      Author: delphineroubinet
 */


#include "UtilitaryStorage.h"
#include <math.h>

using namespace std;


int return_index(int i,int j,int Ny){return i*Ny+j;}

std::pair<int,int> return_indices(int index,int Ny){return std::make_pair(floor(index/Ny),index-floor(index/Ny)*Ny);}
std::pair<int,int> return_indices(double x_pos,double y_pos,double delta_x,double delta_y,int Nx,int Ny,double x_min,double y_min){
	return std::make_pair(min(int((x_pos-x_min)/delta_x),Nx-1),min(int((y_pos-y_min)/delta_y),Ny-1));
}

CgalPoint2D CenterPosition(std::pair<int,int> indices,std::pair<double,double> delta){
	return CgalPoint2D(indices.first*delta.first+0.5*delta.first,indices.second*delta.second+0.5*delta.second);
}

void print(std::pair<int,int> pt){cout << "(" << pt.first << "," << pt.second << ")" << endl;}

double min_vector(vector<double> vect){
	double min_value=vect[0];
	for (int i=0;(unsigned)i<vect.size();i++){
		if (min_value>vect[i]){min_value=vect[i];}
	}
	return min_value;
}

double max_vector(vector<double> vect){
	double max_value=vect[0];
	for (int i=0;(unsigned)i<vect.size();i++){
		if (max_value<vect[i]){max_value=vect[i];}
	}
	return max_value;
}
