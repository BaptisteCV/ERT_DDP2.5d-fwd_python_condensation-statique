/*
 * Constantes.h
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */

#ifndef CONSTANTES_H_
#define CONSTANTES_H_

#include "Geometry/Pointcpp.h"
#include <math.h>

#ifndef EPSILON
//#define EPSILON 1E-6		//epsilon value
#define EPSILON 1E-10          //epsilon value
#endif

static const double NOT_DEFINED=-999;

static const double PI=4.*atan(1.0);

#ifndef RHO
#define RHO 1E3		//fluid density
#endif

#ifndef G
#define G 9.81		//gravity
#endif

#ifndef MU
#define MU 1E-3		//dynamic fluid viscosity
#endif

typedef std::map<pointcpp<double>,std::pair<int,double> > Projection_Map;



// Scale of representation
const int SCALE_LIN = 0;
const int SCALE_LOG = 1;

// Boundary conditions definition
const std::string DIRICHLET="DIRICHLET";
const std::string NEUMANN="NEUMANN";
const std::string DIR_GRADIENT="DIR_GRADIENT";	// Dirichlet conditions where the value along the border follows a gradient
const std::string MIXED_BC="MIXED_BC";
const std::string SOURCE_TERM="SOURCE_TERM";

// Domain borders
const std::string LEFT_BORDER="Left_Border";
const std::string RIGHT_BORDER="Right_Border";
const std::string BOTTOM_BORDER="Bottom_Border";
const std::string TOP_BORDER="Top_Border";
const std::string NO_BORDER="No_Border";

// Simulation option
const std::string DFN_GENERATION="DFN_GENERATION";
const std::string POTENTIAL_COMPUTATION="POTENTIAL_COMPUTATION";
const std::string FOURIER_POTENTIAL_COMPUTATION="FOURIER_POTENTIAL_COMPUTATION";
const std::string PERMEABILITY_ELLIPSE="PERMEABILITY_ELLIPSE";
const std::string CONDUCTIVITY_PLOT="CONDUCTIVITY_PLOT";
const std::string EXCHANGE_STUDY="EXCHANGE_STUDY";
const std::string POTENTIAL_STUDY="POTENTIAL_STUDY";

// DFN Generation option
const std::string PARALLEL="PARALLEL";
const std::string PARALLEL_BORDER="PARALLEL_BORDER";
const std::string SINGLE="SINGLE";
const std::string SIERPINSKI="SIERPINSKI";
const std::string SUGAR_BOX="SUGAR_BOX";
const std::string RANDOM="RANDOM";
const std::string DETERMINISTIC1="DETERMINISTIC1";
const std::string DETERMINISTIC2="DETERMINISTIC2";
const std::string DETERMINISTIC3="DETERMINISTIC3";

// Model option
const std::string DDP="DDP";	// Discrete Dual-Porosity
const std::string DDP_CS="DDP_CS";	// Discrete Dual-Porosity with static condensation (Eigen)
const std::string DDP_EIGEN="DDP_EIGEN";	// Discrete Dual-Porosity with native Eigen containers (no static condensation)
const std::string EPM="EPM";	// Equivalent Porous Media
const std::string EPM_Eigen="EPM_Eigen";	// Equivalent Porous Media (Eigen, no static condensation)
const std::string EPM_CS="EPM_CS";	// Equivalent Porous Media with static condensation
const std::string DFN="DFN";	// Discrete Fracture Network
const std::string EPMDFN="EPMDFN";

// Boundary conditions configuration
const std::string CONFIG1="CONFIG1";	// Configuration for ellipse permeability evaluation with impervious conditions for the lateral borders
const std::string CONFIG2="CONFIG2";	// Configuration for ellipse permeability evaluation with gradient potential on the lateral borders
const std::string CONFIG3="CONFIG3";	// Configuration for ellipse permeability evaluation with gradient potential on the transverse borders

// Option to compute intersections between the fracture network and the grid
const std::string INTERSECTION="intersection";
const std::string NO_INTERSECTION="no_intersection";


const std::string DEFAULT="DEFAULT";




#endif /* CONTANTES_H_ */
