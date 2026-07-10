/*
 * DFNGeneration.h
 *
 *  Created on: Apr 23, 2013
 *      Author: delphineroubinet
 */

#ifndef DFNGENERATION_H_
#define DFNGENERATION_H_

#include "../DFNModel/NetworkMeshes.h"
#include "../Utilitaries/Visualization/VisuStructures.h"


//DFNVisu DFNGeneration(Parameters,Domain);
NetworkMeshes DFNGenerationParallel(Domain,FractureFamilyMap,std::string intersection_option="intersection");
void DivideDomain(pointcpp<double>,pointcpp<double>,int,double,double,NetworkMeshes&);
NetworkMeshes DFNGenerationSierpinski(Domain,double,double,int,int,int,std::string,std::string intersection_option="intersection");
std::vector<Domain> DivideSierpinskiDomain(pointcpp<double>,int,double,double,double,double,NetworkMeshes&);
void AddFracture(FractureMesh,NetworkMeshes&);
void Border_Definition(FluxPoint2D,NetworkMeshes&);
NetworkMeshes DFNGenerationSugarBox(Domain,double,double,double);
NetworkMeshes DFNGenerationRandom(Domain,double,double,double,double,int,std::string intersection_option="intersection");
NetworkMeshes DFNGenerationParallelBorder(Domain,double aperture,double conductivity,std::string intersection_option="intersection");
NetworkMeshes DFNGenerationRandomSingle(Domain,double,double,double,int,int,std::string intersection_option="intersection");

NetworkMeshes DFNRotation(NetworkMeshes,double);
bool DFNExtraction(NetworkMeshes&,double);
bool DFNExtractionNotCentered(NetworkMeshes&,double);
bool DFNExtractionNotCenteredBis(NetworkMeshes&,double,double,double,double);
void ComputeIntersections(NetworkMeshes&);
void DFNCutting(NetworkMeshes&,double);

#endif /* DFNGENERATION_H_ */
