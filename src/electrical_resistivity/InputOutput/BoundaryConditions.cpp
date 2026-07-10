/*
 * BoundaryConditions.cpp
 *
 *  Created on: May 1, 2013
 *      Author: delphineroubinet
 */

#include "BoundaryConditions.h"
#include "../Utilitaries/Constantes.h"
#include <iostream>
#include <fstream>

using namespace std;


void BoundaryConditions::print(){cout << bc_type << "," << bc_value << endl;}

BoundaryConditionsDef ReadBoundaryConditions(std::string file_name,SourceTermsDef& source_terms_def){
	BoundaryConditionsDef bc_map;
	string bc_type;double bc_value,x_coord,y_coord,source_value,source_spacing;
	int nb_source_terms;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		// left border
		fichier >> bc_type >> bc_value;
		bc_map[LEFT_BORDER]=BoundaryConditions(bc_type,bc_value);
		// right border
		fichier >> bc_type >> bc_value;
		bc_map[RIGHT_BORDER]=BoundaryConditions(bc_type,bc_value);
		// top border
		fichier >> bc_type >> bc_value;
		bc_map[TOP_BORDER]=BoundaryConditions(bc_type,bc_value);
		// bottom border
		fichier >> bc_type >> bc_value;
		bc_map[BOTTOM_BORDER]=BoundaryConditions(bc_type,bc_value);
		// source terms
		/*fichier >> nb_source_terms;
		for (int i=0;i<nb_source_terms;i++){
			fichier >> x_coord >> y_coord >> source_value >> source_spacing;
			source_terms_def.push_back(SourceTerm(x_coord,y_coord,source_value,source_spacing));
		}*/
	}
	else{cout << "WARNING in ReadBoundaryConditions (BoundaryConditions.cpp): Boundary conditions parameter file not found" << endl;}
	return bc_map;
}

BoundaryConditionsDef DefineBoundaryConditions(std::string option){
	BoundaryConditionsDef bc_map;
	if (option==CONFIG1){
		bc_map[TOP_BORDER]=BoundaryConditions(DIRICHLET,1);
		bc_map[BOTTOM_BORDER]=BoundaryConditions(DIRICHLET,0);
		bc_map[LEFT_BORDER]=BoundaryConditions(NEUMANN,0);
		bc_map[RIGHT_BORDER]=BoundaryConditions(NEUMANN,0);
	}
	else if (option==CONFIG2){
		bc_map[TOP_BORDER]=BoundaryConditions(DIRICHLET,1);
		bc_map[BOTTOM_BORDER]=BoundaryConditions(DIRICHLET,0);
		bc_map[LEFT_BORDER]=BoundaryConditions(DIR_GRADIENT);
		bc_map[RIGHT_BORDER]=BoundaryConditions(DIR_GRADIENT);
	}
	else if (option==CONFIG3){
		bc_map[LEFT_BORDER]=BoundaryConditions(DIRICHLET,1);
		bc_map[RIGHT_BORDER]=BoundaryConditions(DIRICHLET,0);
		bc_map[TOP_BORDER]=BoundaryConditions(DIR_GRADIENT);
		bc_map[BOTTOM_BORDER]=BoundaryConditions(DIR_GRADIENT);
	}
	else{cout << "WARNING in DefineBoundaryConditions (BoundaryConditions.cpp): option not defined" << endl;}
	return bc_map;
}


void print_BC(BoundaryConditionsDFN bc_map){
	for (BoundaryConditionsDFN::iterator it=bc_map.begin();it!=bc_map.end();it++){
		cout << "Boundary conditions:" << endl;
		cout << it->first << " " << it->second.bc_type << " " << it->second.bc_value << endl;
	}
}

void print_BC(BoundaryConditionsMapCurv bc_map){
	cout << "Boundary conditions:" << endl;
	for (BoundaryConditionsMapCurv::iterator it=bc_map.begin();it!=bc_map.end();it++){
		cout << it->first << " " << it->second.bc_type << " " << it->second.bc_value << endl;
	}
}
