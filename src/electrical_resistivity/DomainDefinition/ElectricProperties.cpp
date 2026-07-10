/*
 * ElectricProperties.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet
 */

#include "ElectricProperties.h"

/************************************/
// Functions for ElectricProperties //
/************************************/
ElectricProperties::ElectricProperties(){}
ElectricProperties::~ElectricProperties(){}
ElectricProperties::ElectricProperties(FractureMesh fract_mesh,double cond_fract){
	length=fract_mesh.ReturnLength();
	aperture=fract_mesh.aperture;
	conductivity=cond_fract;
}



