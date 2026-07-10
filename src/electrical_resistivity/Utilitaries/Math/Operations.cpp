/*
 * Operations.cpp
 *
 *  Created on: May 27, 2013
 *      Author: delphineroubinet
 */

#include "Operations.h"
#include <iostream>

using namespace std;

double Harmonic_Mean(double value1,double value2){
	return 2./(1./value1+1./value2);
}

double Arithmic_Mean(double value1,double value2){
	return (value1+value2)*0.5;
}
