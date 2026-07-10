/*
 * FractureMesh.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef FRACTUREMESH_H_
#define FRACTUREMESH_H_

#include "../Utilitaries/Geometry/FluxPoint2D.h"
#include "../Utilitaries/Geometry/Pointcpp.h"
#include "../Utilitaries/Math/RandomNumber.h"
#include "../Utilitaries/Geometry/Segment.h"
#include "../Utilitaries/Constantes.h"


// definition of a mesh
class FractureMesh{
public:
	FluxPoint2D p_ori; /**< origin of the mesh. */
	FluxPoint2D p_tar;	/**< target of the mesh. */
	double velocity;	/**< the velocty in the mesh */
	double aperture;	// mesh aperture
	int mesh_index; /**< the index of the mesh .*/
	std::vector<int> neigh_meshes;	/**< the neighboring meshes.*/
	double conductivity;
public:
	FractureMesh(){};
	virtual ~FractureMesh(){};
	FractureMesh(FluxPoint2D,FluxPoint2D,double,double,int,std::vector<int>);
	FractureMesh(FluxPoint2D,FluxPoint2D,double,double);
	FractureMesh(Segment2D,double,double);
	double ReturnLength() const;
	CgalPoint2D ReturnCenter();
	double ReturnAngle();
	Segment2D ReturnSegment();
	double ReturnConductance();
};


#endif /* FRACTUREMESH_H_ */
