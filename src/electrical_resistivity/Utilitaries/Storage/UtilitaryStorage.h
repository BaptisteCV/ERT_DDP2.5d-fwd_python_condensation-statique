/*
 * UtilitaryStorage.h
 *
 *  Created on: Apr 29, 2013
 *      Author: delphineroubinet
 */

#ifndef UTILITARYSTORAGE_H_
#define UTILITARYSTORAGE_H_

#include <utility>
#include "../Geometry/Point_Cgal.h"


int return_index(int,int,int);
std::pair<int,int> return_indices(int,int);
std::pair<int,int> return_indices(double,double,double,double,int,int,double x_min=0,double y_min=0);
CgalPoint2D CenterPosition(std::pair<int,int>,std::pair<double,double>);
void print(std::pair<int,int>);
double min_vector(std::vector<double>);
double max_vector(std::vector<double>);


#endif /* UTILITARYSTORAGE_H_ */
