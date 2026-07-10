/*
 * DDPDFNFunctions.h
 *
 *  Created on: May 27, 2013
 *      Author: delphineroubinet
 */

#ifndef DDPDFNFUNCTIONS_H_
#define DDPDFNFUNCTIONS_H_

#include "../DFNModel/NetworkMeshes.h"


void FracturesIntersectionGrid(int,int,NetworkMeshes&);
SubNetworkMap DefineSubNetwork(DomainMap,NetworkMeshes);
SubNetworkMap DDPNetworkMeshes(NetworkMeshes&,int,int);
void ReturnDerivativeCoeffFull(double,double,double,double&,double&,double&);
void CumulatePotFractCoeff(double,double,double&,double&,double&);

void ReturnDerivativeCoeffFullFourier(double,double,double,double,double&,double&,double&);
void CumulatePotFractCoeffFourier(double,double,double,double&,double&,double&);

#endif /* DDPDFNFUNCTIONS_H_ */
