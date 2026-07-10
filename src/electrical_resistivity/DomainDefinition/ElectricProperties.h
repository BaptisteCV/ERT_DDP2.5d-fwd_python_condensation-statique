/*
 * ElectricProperties.h
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet
 */

#ifndef ELECTRICPROPERTIES_H_
#define ELECTRICPROPERTIES_H_

#include "../DFNModel/FractureMesh.h"

/************************************/
// Functions for ElectricProperties //
/************************************/
class ElectricProperties{
public:
	double length;
	double aperture;
	double conductivity;
public:
	ElectricProperties();
	~ ElectricProperties();
	ElectricProperties(FractureMesh,double);
};



#endif /* ELECTRICPROPERTIES_H_ */
