/*
 * BoundaryConditions.h
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet
 */

#ifndef BOUNDARYCONDITIONS_H_
#define BOUNDARYCONDITIONS_H_

#include <string>
#include <map>
#include <vector>

/************************************/
/* Functions for BoundaryConditions */
/************************************/
class BoundaryConditions{
public:
	std::string bc_type;
	double bc_value;
	int block_index;
public:
	BoundaryConditions(){};
	BoundaryConditions(std::string bc_type_,double bc_value_=-1){bc_type=bc_type_;bc_value=bc_value_;block_index=-9999;}
	~BoundaryConditions(){};
	void print();
};

/************************************/
/* Functions for SourceTerms */
/************************************/
class SourceTerm{
public:
	double x_position;
	double y_position;
	double source_value;
	double source_spacing; // to generate several experiments
public:
	SourceTerm(){};
	SourceTerm(double x_position_,double y_position_,double source_value_,double source_spacing_=-1){x_position=x_position_;y_position=y_position_;source_value=source_value_;source_spacing=source_spacing_;}
	~SourceTerm(){};
};

typedef std::map<int,BoundaryConditions> BoundaryConditionsDFN;	// list of the fracture nodes located on the domain borders: <node index,associated boundary conditions>
typedef std::map<std::string,BoundaryConditions> BoundaryConditionsDef;
typedef std::vector<SourceTerm> SourceTermsDef;
typedef std::map<double,BoundaryConditions> BoundaryConditionsMapCurv;


BoundaryConditionsDef ReadBoundaryConditions(std::string,SourceTermsDef&);
BoundaryConditionsDef DefineBoundaryConditions(std::string);

void print_BC(BoundaryConditionsDFN bc_map);
void print_BC(BoundaryConditionsMapCurv bc_map);



#endif /* BOUNDARYCONDITIONS_H_ */
