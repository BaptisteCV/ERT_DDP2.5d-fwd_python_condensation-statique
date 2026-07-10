/*
 * Domain.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef DOMAIN_H_
#define DOMAIN_H_

#include "../Utilitaries/Geometry/Pointcpp.h"
#include "../Utilitaries/Geometry/Segment.h"


//////////////////////////
/* Functions for Domain */
//////////////////////////
class Domain {
public:
	pointcpp<double> min_pt;
	pointcpp<double> max_pt;
public:
	Domain();
	virtual ~Domain();
	Domain(double,double);
	Domain(pointcpp<double>,pointcpp<double>);
	bool on_input_limit(pointcpp<double>);
	bool on_output_limit(pointcpp<double>);
	double domain_size_x(){return max_pt.i-min_pt.i;}
	double domain_size_y(){return max_pt.j-min_pt.j;}
	void get_extremities(pointcpp<double> & pmin, pointcpp<double> & pmax){pmin=min_pt; pmax=max_pt;}
	std::string ReturnBorder(CgalPoint2D);
	double ReturnDeltaX(int);
	double ReturnDeltaY(int);
	bool IsInDomain(CgalPoint2D);
	bool IsInDomain(Segment2D);
	Segment2D ReturnBorder(std::string);
	bool SegmentIntersectDomain(Segment2D,Segment2D&);
	bool IntersectionBorders(Segment2D,CgalPoint2D&,CgalPoint2D&);
	double ReturnBorderCoordinate(CgalPoint2D,double&,std::string&,std::string&);
	void PrintDomain();
};

pointcpp<double> Translate(pointcpp<double>,double,double);
void ReturnOppositeBorders(std::string,std::string&,std::string&);

//////////////////////////////
/* Functions for SubDomains */
//////////////////////////////
typedef std::map<int,Domain> DomainMap; // <domain index,domain definition>

DomainMap DefineSubDomains(Domain,int,int);



#endif /* DOMAIN_H_ */
