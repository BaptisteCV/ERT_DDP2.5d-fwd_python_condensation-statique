/*

 * VisuTools.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet*/


/*#include "VisuTools.h"
#include <gl.h>
#include <iostream>
#include <glut.h>

using namespace std;

// set a color corresponding to x value for 0<x<1
void SetColor(double x){

	double r=0,g=0,b=0;

	if (x==0){r=0; g=0; b=0;}
	else if((x>=0)&&(x<0.125)) {r=0; g=0; b=0.5+4*x;}
	else if((x>0.125)&&(x<=0.375)) {r=0; g=4*x-0.5; b=1;}
	else if((x>0.375)&&(x<=0.625)) {r=4*x-1.5; g=1; b=2.5-4*x;}
	else if((x>0.625)&&(x<=0.875)) {r=1; g=3.5-4*x; b=0;}
	else if((x>0.875)&&(x<=1)) {r=4.5-4*x; g=0; b=0;}
	else if (x<0){r=1; g=1; b=1;}
	else{
		cout << "WARNING in SetColor (VisualizationUtilitary.cpp): color not defined" << endl;
		cout << "x = " << x << endl;
		if (x>1){x=1; r=4.5-4*x; g=0; b=0;}
	}

	glColor3f(r,g,b);

}

// set a color corresponding to x value for x_min<x<x_max
void SetColor(double x, double x_min, double x_max){
	if (x_min==x_max) SetColor(0.);
	else if (x==0) SetColor(0.);
	else if (x<0) SetColor(x);
	else SetColor((x-x_min)/(x_max-x_min));
	SetColor((x-x_min)/(x_max-x_min));
}

bool draw_legend_scale( double Lx, double Ly, double min_, double max_, int nb_ticks, int lin_log,double legend_step){
	if (nb_ticks==0){std::cout<<"\n PB IN Visualisation_Glut_D2_Flow.h : nb_ticks is not defined"; return false;}
	//parameters
	char c_temp[50];
	std::string s_temp;
	//void * font = (int *) ((void*)7);
	double bar_thick = 0.05*Lx,step=2;

	for (int i = 0; i<nb_ticks; i++){
		double val = i/(double)(nb_ticks-1);
		double h = val*(max_-min_)+min_;
		Gl_color color = val;
		color.set_color();
		SetColor(val,min_,max_);
		SetColor(h,min_,max_);

		glPolygonMode(GL_FRONT,GL_RGB);
		glBegin(GL_POLYGON);
		glVertex2d(Lx+step,i*Ly/nb_ticks+0.5*legend_step);
		glVertex2d(Lx+step+bar_thick,i*Ly/nb_ticks+0.5*legend_step);
		glVertex2d(Lx+step+bar_thick,(i+1)*Ly/nb_ticks+0.5*legend_step);
		glVertex2d(Lx+step,(i+1)*Ly/nb_ticks+0.5*legend_step);
		glEnd();

		//SetColor(NOIR);
		SetColor(0.);
		//sprintf(c_temp,"%.0f", h);
		sprintf(c_temp,"%.1f", h);
		s_temp = std::string(c_temp);

		double x = Lx+step+bar_thick+step; double y = (i+0.25)*Ly/nb_ticks+0.5*legend_step;
		glRasterPos2f(x, y);                             // set position to start drawing fonts
		for (int c=0; c < s_temp.length(); c++){
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,s_temp.at(c));
		}
	}
	return true;*/
/*	return false;
}
*/

