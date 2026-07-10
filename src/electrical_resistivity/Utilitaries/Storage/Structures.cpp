/*
 * Structures.cpp
 *
 *  Created on: Jun 29, 2012
 *      Author: delphine
 */


#include "Structures.h"
#include <iostream>

using namespace std;


double get_min_value(std::map<int,double> values){
	double min_value = values.begin()->second;
	for (std::map<int,double>::iterator it=values.begin();it!=values.end();it++){
		if (it->second<min_value) min_value = it->second;
	}
	return min_value;
}

double get_max_value(std::map<int,double> values){
	double max_value = values.begin()->second;
	for (std::map<int,double>::iterator it=values.begin();it!=values.end();it++){
		if (it->second>max_value) max_value = it->second;
	}
	return max_value;
}

void print(std::map<int,double> values){
	for (std::map<int,double>::iterator it=values.begin(); it!=values.end(); it++){
		cout << it->first << " " << it->second << endl;
	}
}

void print(std::vector<double> values){
	int N = values.size();
	for (int i=0;i<N;i++){
		cout << values[i] << endl;
	}
}

void print(std::vector<std::pair<int,int> > values){
	for (int i=0;(unsigned)i<values.size();i++){
		cout << "(" << values[i].first << "," << values[i].second << ")" << endl;
	}
}

void print(Distribution values){
	for (Distribution::iterator it=values.begin();it!=values.end();it++){
		cout << "(" << it->first.first << "," << it->first.second << "):" << endl;
		cout << "vector1:" << endl;
		int size1 = it->second.first.size();
		for (int i=0;i<size1;i++){cout << it->second.first[i] << endl;}
		cout << "vector2:" << endl;
		int size2 = it->second.second.size();
		for (int i=0;i<size2;i++){cout << it->second.second[i] << endl;}
	}
}

double define_delta(double min, double max, int N, int lin_log){
	if (N==1)
		return 0.;
	else if (lin_log==SCALE_LOG){
		if (min==0) min = 1;
		return pow(max/min,1./(N-1));
	}
	else if (lin_log==SCALE_LIN){
		return (max-min)/(N-1);
	}
	else{cout << "WARNING in define_delta (Results.cpp): scale option not defined" << endl;}
	return 0;
}

vector<double> define_scale(double min, double max, int N, int lin_log){

	vector<double> scale(N);
	double delta = define_delta(min,max,N,lin_log);
	for (int i=0;i<N;i++){
		if (lin_log==SCALE_LIN){
			scale[i] = min+delta*i;
		}
		else if (lin_log==SCALE_LOG){
			scale[i] = min*pow(delta,i);
		}
		else{cout << "WARNING in define_scale (Results.cpp): scale option not defined" << endl;}
	}
	return scale;
}

bool find_bound_in_vector(const std::vector<double> & vect, const double value, int & index){
	index = -1;
	for (size_t i = 0; i<vect.size()-1; i++){
		if ((vect[i]<=value)&&(vect[i+1]>=value)){
			index = i;
			return true;
		}
	}
	std::cout<<"\n PB IN find_bound_in_vector (Vector_Utilitary.cpp) : value not found";
	return  false;
}

double linear_interpolation_2points(const double xmin, const double xmax, const double ymin, const double ymax, const double x){
	return ymin + (x - xmin) / (xmax - xmin) * (ymax - ymin);
}

/**
* return true if value is found as the first value of one of the pairs in projection_map
* return the corresponding point
*/
bool find_projected_value(map<pointcpp<double>,pair<int,double> > proj_map, int value, pointcpp<double> & point){
	for (map<pointcpp<double>,pair<int,double> >::iterator it = proj_map.begin(); it!=proj_map.end(); it++){
		if (it->second.first == value){
			point = it->first;
			return true;
		}
	}
	cout << "WARNING in find_projected_value (Structures.cpp): projected value not found" << endl;
	return false;
}

bool exist_in_vector(vector<pair<int,int> > vect,pair<int,int> value){
	for (vector<pair<int,int> >::iterator it=vect.begin();it!=vect.end();it++){
		if (*it==value){
			return true;
		}
	}
	return false;
}

double average_map(map<int,double> values){
	double av_value=0;
	for (map<int,double>::iterator it=values.begin();it!=values.end();it++){
		av_value+=fabs(it->second);
	}
	av_value=av_value/values.size();
	return av_value;
}

double ErrorMapL2(map<int,double> map1,map<int,double> map2){
	double res=0;
	if (map1.size()!=map2.size()){
		cout << "WARNING in ErrorMapL2 (Structures.cpp): maps' sizes are different" << endl;
		return res;
	}
	for (map<int,double>::iterator it=map1.begin();it!=map1.end();it++){
		res+=(it->second-map2[it->first])*(it->second-map2[it->first]);
	}
	return res/map1.size();
}

double RelativeErrorMap(map<int,double> map1,map<int,double> map_ref){
	double res=0;
	if (map1.size()!=map_ref.size()){
		cout << "WARNING in ErrorMapL2 (Structures.cpp): maps' sizes are different" << endl;
		return res;
	}
	for (map<int,double>::iterator it=map1.begin();it!=map1.end();it++){
		res+=fabs(it->second-map_ref[it->first])/map_ref[it->first];
	}
	return res/map1.size();
}
