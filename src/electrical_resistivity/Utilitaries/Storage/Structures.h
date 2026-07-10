/*
 * Structures.h
 *
 *  Created on: Jun 29, 2012
 *      Author: delphine
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include "../Constantes.h"
#include "../Geometry/Pointcpp.h"
#include <map>
#include <vector>

typedef std::map<std::pair<double,double>,std::pair<std::vector<double>, std::vector<double> > > Distribution;


double get_min_value(std::map<int,double>);
double get_max_value(std::map<int,double>);
void print(std::map<int,double>);
void print(std::vector<double>);
void print(std::vector<std::pair<int,int> >);
void print(Distribution);
double define_delta(double,double,int,int lin_log = SCALE_LIN);
std::vector<double> define_scale(double,double,int,int lin_log = SCALE_LIN);
bool find_bound_in_vector(const std::vector<double> &,const double,int&);
double linear_interpolation_2points(const double,const double,const double,const double,const double);
bool find_projected_value(std::map<pointcpp<double>,std::pair<int,double> >,int,pointcpp<double>&);
bool exist_in_vector(std::vector<std::pair<int,int> >,std::pair<int,int>);
double average_map(std::map<int,double>);
double ErrorMapL2(std::map<int,double>,std::map<int,double>);
double RelativeErrorMap(std::map<int,double>,std::map<int,double>);

/**returns true if the key 'key' is stored in the map with the precision epsilon and returns the corresponding value*/
template<typename T1,typename T2>
bool exists_key(const std::map<T1,T2> & M, const T1 & key, T2 & value, const double epsilon = EPSILON, bool option_reverse = false){
	for(typename std::map<T1,T2>::const_iterator it=M.begin(); it!=M.end(); it++)
		if (fabs((double)key-it->first)<EPSILON){value = it->second; return true;}
	return false;
}

template<> inline
bool exists_key(const std::map<std::pair<double,double>,std::pair<std::vector<double>,std::vector<double> > > & M, const std::pair<double,double> & key, std::pair<std::vector<double>,std::vector<double> > & value, const double epsilon, bool option_reverse){

	std::map<std::pair<double,double>,std::pair<std::vector<double>,std::vector<double> > >::const_iterator it;
	for(it=M.begin(); it!=M.end(); it++){
		if ((std::fabs(it->first.first-key.first)<epsilon)&&(std::fabs(it->first.second-key.second)<epsilon)){
			value = it->second;
			return true;
		}
	}
	if (option_reverse) exists_key(M, std::make_pair(key.second,key.first), value, epsilon);

	return false;
}

/**returns true if the key 'key' is stored in the map with the precision epsilon and returns the corresponding iterator*/
template<typename T1,typename T2>
bool remove_from_key(std::map<T1,T2> & M, const T1 & key, const double epsilon = EPSILON){
	for(typename std::map<T1,T2>::iterator it=M.begin(); it!=M.end(); it++)
		if (fabs((double)key-it->first)<EPSILON){M.erase(it); ; return true;}
	return false;
}

template<> inline
bool remove_from_key(std::map<std::pair<double,double>,std::pair<std::vector<double>,std::vector<double> > > & M, const std::pair<double,double> & key,
					 const double epsilon){

	std::map<std::pair<double,double>,std::pair<std::vector<double>,std::vector<double> > >::iterator it;
	for(it=M.begin(); it!=M.end(); it++){
		if ((std::fabs(it->first.first-key.first)<epsilon)&&(std::fabs(it->first.second-key.second)<epsilon)){
			M.erase(it);
			return true;
		}
	}
	return false;
}

template<typename T>
void insert_map(std::map<int,T>& init_map,std::map<int,T> new_map){
	int cpt=init_map.size();
	for (typename std::map<int,T>::iterator it=new_map.begin();it!=new_map.end();it++){
		init_map[cpt]=it->second;
		cpt++;
	}
}


#endif /* STRUCTURES_H_ */
