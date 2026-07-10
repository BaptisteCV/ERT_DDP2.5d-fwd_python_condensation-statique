/*

 * VisuStructures.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet
*/
/*
#include "VisuStructures.h"
#include "VisuTools.h"
#include "../../DFNModel/DFNComputation.h"
#include "../../DDPModel/ExchangeFunctions.h"
#include "../../EPMModel/EPMDefinition.h"

#include <gl.h>
#include <glut.h>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

DFNVisu::DFNVisu(){}
DFNVisu::DFNVisu(Domain domain_,NetworkMeshes net_mesh_,ublas_vector DFN_potential_,double cond_fract_,double cond_matrix_,double fract_spacing_){
	domain=domain_;net_mesh=net_mesh_;DFN_potential=DFN_potential_;cond_fract=cond_fract_;cond_matrix=cond_matrix_;fract_spacing=fract_spacing_;
}
DFNVisu::DFNVisu(double L_,NetworkMeshes net_mesh_,string option_){
	L=L_;net_mesh=net_mesh_;option=option_;
}
DFNVisu::DFNVisu(double L_,NetworkMeshes net_mesh_,ublas_vector DFN_potential_,string option_){
	L=L_;net_mesh=net_mesh_;DFN_potential=DFN_potential_;option=option_;
}
*/
/*
void DFNVisu::display_DFN(){
	// 1. window setting
	glClear(GL_COLOR_BUFFER_BIT);	// clear the window with current clearing color
	glTranslatef(-1,-1,-1);	// change coordinate reference (translate it from the center to the bottom left of the domain)
	//double Lx=domain.domain_size_x(),Ly=domain.domain_size_y();
	double Lx=L,Ly=L;
	glScalef(2./Lx,2./Ly,1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	// 2. value display
	// draw a polygon for each cell with a color representing its value
	// white polygon for the back
	SetColor(-1.);
	glBegin(GL_POLYGON);
	glVertex2f(0,0);
	glVertex2f(0,Ly);
	glVertex2f(Lx,Ly);
	glVertex2f(Lx,0);
	glEnd();
	// plot each fracture
	vector<FractureMesh> meshes=net_mesh.meshes;
	double factor_ap=1,factor_x=0,factor_y=0;
	for (int i=0;(unsigned)i<meshes.size();i++){
		glLineWidth(meshes[i].aperture*factor_ap);
		glBegin(GL_LINES);
		//SetColor(current_mat.get(i,j),min_visu,max_visu);	// line color
		//SetColor(0.);
		glColor3f(0.,0.,0.);
		glVertex2f(meshes[i].p_ori.p.x()+factor_x,meshes[i].p_ori.p.y()+factor_y);
		glVertex2f(meshes[i].p_tar.p.x()+factor_x,meshes[i].p_tar.p.y()+factor_y);
		glEnd();
	}

	// polygon contour
	double contour_width=2;
	SetColor(1.);
	glLineWidth(contour_width);
	glBegin(GL_LINES);
	glVertex2f(0,0);
	glVertex2f(0,Ly);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(Lx,0);
	glVertex2f(Lx,Ly);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(0,0);
	glVertex2f(Lx,0);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(0,Ly);
	glVertex2f(Lx,Ly);
	glEnd();


	glFlush();
}

void DFNVisu::display_Potential_Classic(){
	// 1. window setting
	glClear(GL_COLOR_BUFFER_BIT);	// clear the window with current clearing color
	glTranslatef(-1,-1,-1);	// change coordinate reference (translate it from the center to the bottom left of the domain)
	//double Lx=domain.domain_size_x(),Ly=domain.domain_size_y();
	double Lx=L,Ly=L;
	glScalef(2./Lx,2./Ly,1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	// 2. value display
	// draw a polygon for each cell with a color representing its value
	glPolygonMode(GL_FRONT,GL_FILL);	// polygon with a color
	// white polygon for the back
	glBegin(GL_POLYGON);
	SetColor(-1.);
	glVertex2i(0,0);
	glVertex2i(0,Ly);
	glVertex2i(Lx,Ly);
	glVertex2i(Lx,0);
	glEnd();
	// plot each fracture
	vector<FractureMesh> meshes=net_mesh.meshes;
	double delta_seg=1e-1;
	//double factor_ap=1e5,factor_x=Lx/2,factor_y=Ly/2,value,length;
	double factor_ap=1,factor_x=0,factor_y=0,value,length;
	pointcpp<double> p_ori,p1,p2,p_current;int Nx;
	double min_visu=min_vector(DFN_potential),max_visu=max_vector(DFN_potential);
	for (int i=0;(unsigned)i<meshes.size();i++){
		glLineWidth(meshes[i].aperture*factor_ap);
		CgalVector2D vect(meshes[i].p_ori.p, meshes[i].p_tar.p);
		p_ori=CgalPoint2D(meshes[i].p_ori.p);
		length=CGAL::sqrt(vect.squared_length());
		Nx=length/delta_seg;
		for (int j=0;j<Nx;j++){
			value=AnalyticalSolutionFractureClassic(delta_seg*(j+0.5),0,length,DFN_potential[meshes[i].p_ori.index],DFN_potential[meshes[i].p_tar.index]);
			p1.i=p_ori.i+vect.x()*delta_seg/length*j;
			p1.j=p_ori.j+vect.y()*delta_seg/length*j;
			p2.i=p_ori.i+vect.x()*delta_seg/length*(j+1);
			p2.j=p_ori.j+vect.y()*delta_seg/length*(j+1);
			glBegin(GL_LINES);
			SetColor(value,min_visu,max_visu);	// line color
			glVertex2f(p1.i+factor_x,p1.j+factor_y);
			glVertex2f(p2.i+factor_x,p2.j+factor_y);
			glEnd();
		}
	}
	glFlush();
}

void DFNVisu::display_Potential_Radial(){
	// 1. window setting
	glClear(GL_COLOR_BUFFER_BIT);	// clear the window with current clearing color
	glTranslatef(-1,-1,-1);	// change coordinate reference (translate it from the center to the bottom left of the domain)
	//double Lx=domain.domain_size_x(),Ly=domain.domain_size_y();
	double Lx=L,Ly=L;
	glScalef(2./Lx,2./Ly,1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	// 2. value display
	// draw a polygon for each cell with a color representing its value
	glPolygonMode(GL_FRONT,GL_FILL);	// polygon with a color
	// white polygon for the back
	glBegin(GL_POLYGON);
	SetColor(-1.);
	glVertex2i(0,0);
	glVertex2i(0,Ly);
	glVertex2i(Lx,Ly);
	glVertex2i(Lx,0);
	glEnd();
	// plot each fracture
	vector<FractureMesh> meshes=net_mesh.meshes;
	double delta_seg=1e-1;
	double factor_ap=1e5,factor_x=Lx/2,factor_y=Ly/2,value,length;
	pointcpp<double> p_ori,p1,p2,p_current;int Nx;
	double min_visu=0;//min_vector(DFN_potential);
	double max_visu=1;//max_vector(DFN_potential);

	print_vector(DFN_potential);

	for (int i=0;i<meshes.size();i++){
		glLineWidth(meshes[i].aperture*factor_ap);
		CgalVector2D vect(meshes[i].p_ori.p, meshes[i].p_tar.p);
		p_ori=CgalPoint2D(meshes[i].p_ori.p);
		length=CGAL::sqrt(vect.squared_length());
		Nx=length/delta_seg;
		// fracture value
		for (int j=0;j<Nx;j++){
			value=AnalyticalSolutionFractureRadialDiffusion(delta_seg*(j+0.5),0,length,DFN_potential[meshes[i].p_ori.index],DFN_potential[meshes[i].p_tar.index],
			cond_matrix,meshes[i].aperture,cond_fract);
			p1.i=p_ori.i+vect.x()*delta_seg/length*j;
			p1.j=p_ori.j+vect.y()*delta_seg/length*j;
			p2.i=p_ori.i+vect.x()*delta_seg/length*(j+1);
			p2.j=p_ori.j+vect.y()*delta_seg/length*(j+1);
			glBegin(GL_LINES);
			SetColor(value,min_visu,max_visu);	// line color
			glVertex2f(p1.i+factor_x,p1.j+factor_y);
			glVertex2f(p2.i+factor_x,p2.j+factor_y);
			glEnd();
		}
		// matrix value
		Nx=0.5*fract_spacing/delta_seg;
		for (int j=0;j<Nx;j++){
			value=AnalyticalSolutionMatrixRadialDiffusion(delta_seg*(j+0.5),0,length,DFN_potential[meshes[i].p_ori.index],DFN_potential[meshes[i].p_tar.index],
			cond_matrix,meshes[i].aperture,cond_fract);
			p1.i=p_ori.i+vect.x()*delta_seg/length*j;
			p1.j=p_ori.j+vect.y()*delta_seg/length*j;
			p2.i=p_ori.i+vect.x()*delta_seg/length*(j+1);
			p2.j=p_ori.j+vect.y()*delta_seg/length*(j+1);
			glBegin(GL_LINES);
			SetColor(value,min_visu,max_visu);	// line color
			glVertex2f(p1.i+factor_x,p1.j+factor_y);
			glVertex2f(p2.i+factor_x,p2.j+factor_y);
			glEnd();
		}
	}
	glFlush();
}


EPMVisu::EPMVisu(){}
EPMVisu::EPMVisu(double Lx_,double Ly_,ublas_matrix EPM_potential_,double min_visu_,double max_visu_){
	Lx=Lx_;Ly=Ly_;EPM_potential=EPM_potential_;
	if (min_visu_==-1&&max_visu_==-1){
		min_visu=min_matrix(EPM_potential),max_visu=max_matrix(EPM_potential);
	}
	else{min_visu=min_visu_;max_visu=max_visu_;}
}
EPMVisu::EPMVisu(double Lx_,double Ly_,ublas_matrix EPM_potential_,ublas_matrix current_x_,ublas_matrix current_y_){
	Lx=Lx_;Ly=Ly_;EPM_potential=EPM_potential_;Current_x=current_x_;Current_y=current_y_;
}

void EPMVisu::display_EPM_Potential(){
	double legend_step=25;
	// 1. window setting
	glClear(GL_COLOR_BUFFER_BIT);	// clear the window with current clearing color
	glTranslatef(-1,-1,-1);	// change coordinate reference (translate it from the center to the bottom left of the domain)
	int Nx=EPM_potential.size1(),Ny=EPM_potential.size2();
	glScalef(2./(Lx+legend_step),2./(Ly+legend_step),1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	double delta_x=Lx/Nx,delta_y=Ly/Ny;

	//delta_x=0.1;delta_y=0.1;

	// 2. value display
	// draw a polygon for each cell with a color representing its value
	glPolygonMode(GL_FRONT,GL_FILL);	// polygon with a color
	// white polygon for the back
	glBegin(GL_POLYGON);
	SetColor(-1.);
	glVertex2i(0,0);
	glVertex2i(0,Ny*delta_y+legend_step);
	glVertex2i(Nx*delta_x+legend_step,Ny*delta_y+legend_step);
	glVertex2i(Nx*delta_x+legend_step,0);
	glEnd();

	cout << "min_visu=" << min_visu << " and max_visu=" << max_visu << endl;

	// plot each cell
	for (int i=0;i<EPM_potential.size1();i++){
		for (int j=0;j<EPM_potential.size2();j++){
			glBegin(GL_POLYGON);
			SetColor(EPM_potential(i,j),min_visu,max_visu);	// cell color
			glVertex2f(i*delta_x,j*delta_y+0.5*legend_step);
			glVertex2f(i*delta_x,(j+1)*delta_y+0.5*legend_step);
			glVertex2f((i+1)*delta_x,(j+1)*delta_y+0.5*legend_step);
			glVertex2f((i+1)*delta_x,j*delta_y+0.5*legend_step);
			glEnd();
		}
	}

	int nb_ticks=10;int lin_log=0;
	draw_legend_scale(Lx,Ly,min_visu,max_visu,nb_ticks,lin_log,legend_step);

	glFlush();

}

void EPMVisu::display_EPM_Current_x(){
	// 1. window setting
	glClear(GL_COLOR_BUFFER_BIT);	// clear the window with current clearing color
	glTranslatef(-1,-1,-1);	// change coordinate reference (translate it from the center to the bottom left of the domain)
	int Nx=Current_x.size1(),Ny=Current_x.size2();
	glScalef(2./Nx,2./Ny,1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	// 2. value display
	// draw a polygon for each cell with a color representing its value
	glPolygonMode(GL_FRONT,GL_FILL);	// polygon with a color
	// white polygon for the back
	glBegin(GL_POLYGON);
	SetColor(-1.);
	glVertex2i(0,0);
	glVertex2i(0,Ny);
	glVertex2i(Nx,Ny);
	glVertex2i(Nx,0);
	glEnd();
	// plot each cell
	double min_visu=min_matrix(Current_x),max_visu=max_matrix(Current_x);
	cout << "min_visu =" << min_visu << "and max_visu =" << max_visu << endl;
	for (int i=0;i<Current_x.size1();i++){
		for (int j=0;j<Current_x.size2();j++){
			glBegin(GL_POLYGON);
			SetColor(Current_x(i,j),min_visu,max_visu);	// cell color
			glVertex2i(i,j);
			glVertex2i(i,j+1);
			glVertex2i(i+1,j+1);
			glVertex2i(i+1,j);
			glEnd();
		}
	}
	glFlush();
}

void EPMVisu::display_EPM_Current_y(){
	// 1. window setting
	glClear(GL_COLOR_BUFFER_BIT);	// clear the window with current clearing color
	glTranslatef(-1,-1,-1);	// change coordinate reference (translate it from the center to the bottom left of the domain)
	int Nx=Current_y.size1(),Ny=Current_y.size2();
	glScalef(2./Nx,2./Ny,1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	// 2. value display
	// draw a polygon for each cell with a color representing its value
	glPolygonMode(GL_FRONT,GL_FILL);	// polygon with a color
	// white polygon for the back
	glBegin(GL_POLYGON);
	SetColor(-1.);
	glVertex2i(0,0);
	glVertex2i(0,Ny);
	glVertex2i(Nx,Ny);
	glVertex2i(Nx,0);
	glEnd();
	// plot each cell
	double min_visu=min_matrix(Current_y),max_visu=max_matrix(Current_y);
	cout << "min_visu =" << min_visu << "and max_visu =" << max_visu << endl;
	for (int i=0;i<Current_y.size1();i++){
		for (int j=0;j<Current_y.size2();j++){
			glBegin(GL_POLYGON);
			SetColor(Current_y(i,j),min_visu,max_visu);	// cell color
			glVertex2i(i,j);
			glVertex2i(i,j+1);
			glVertex2i(i+1,j+1);
			glVertex2i(i+1,j);
			glEnd();
		}
	}
	glFlush();
}


EPMDFNVisu::EPMDFNVisu(){}
EPMDFNVisu::EPMDFNVisu(double L_,NetworkMeshes net_mesh_,ublas_vector DFN_potential_,ublas_matrix EPM_potential_){
	L=L_;net_mesh=net_mesh_;DFN_potential=DFN_potential_;EPM_potential=EPM_potential_;
}
EPMDFNVisu::EPMDFNVisu(double L_,SubNetworkMap subnetwork_map_,ublas_vector DFN_potential_,ublas_matrix EPM_potential_,double cond_matrix_,ublas_matrix AlphaBlock_,double min_visu_,double max_visu_){
	L=L_;subnetwork_map=subnetwork_map_;DFN_potential=DFN_potential_;EPM_potential=EPM_potential_;cond_matrix=cond_matrix_;AlphaBlock=AlphaBlock_;
	if (min_visu_==-1&&max_visu_==-1){
		min_visu=min_matrix(EPM_potential),max_visu=max_matrix(EPM_potential);
		min_visu=min(min_matrix(EPM_potential),min_vector(DFN_potential));
		max_visu=max(max_matrix(EPM_potential),max_vector(DFN_potential));
	}
	else{min_visu=min_visu_;max_visu=max_visu_;}
}
EPMDFNVisu::EPMDFNVisu(Parameters param_,SubNetworkMap subnetwork_map_,ublas_vector DFN_potential_,ublas_matrix EPM_potential_){
	param=param_;subnetwork_map=subnetwork_map_;DFN_potential=DFN_potential_;EPM_potential=EPM_potential_;
}

void EPMDFNVisu::display_EPMDFN_Potential(){
	// 1. window setting
	double legend_step=25;
	glClear(GL_COLOR_BUFFER_BIT);	// clear the window with current clearing color
	glTranslatef(-1,-1,-1);	// change coordinate reference (translate it from the center to the bottom left of the domain)

	int Nx=EPM_potential.size1();
        //int Ny=EPM_potential.size2();
	//glScalef(2./Nx,2./Ny,1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	double Lx=L,Ly=L;
	glScalef(2./(Lx+legend_step),2./(Ly+legend_step),1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	double block_size=L/Nx;

	//block_size=10;

	// 2. value display
	double min_visu=min(min_matrix(EPM_potential),min_vector(DFN_potential)),
			max_visu=max(max_matrix(EPM_potential),max_vector(DFN_potential));
	//double min_visu=0,max_visu=1;
	cout << "min_visu=" << min_visu << " and max_visu=" << max_visu << endl;

	// draw a polygon for each cell with a color representing its value
	glPolygonMode(GL_FRONT,GL_FILL);	// polygon with a color
	// white polygon for the back
	glBegin(GL_POLYGON);
	SetColor(-1.);
	glVertex2i(0,0);
	glVertex2i(0,Ly+legend_step);
	glVertex2i(Lx+legend_step,Ly+legend_step);
	glVertex2i(Lx+legend_step,0);
	glEnd();

	// plot each cell
	for (int i=0;i<EPM_potential.size1();i++){
		for (int j=0;j<EPM_potential.size2();j++){
			glBegin(GL_POLYGON);
			SetColor(EPM_potential(i,j),min_visu,max_visu);	// cell color
			glVertex2i(i*block_size,j*block_size);
			glVertex2i(i*block_size,(j+1)*block_size);
			glVertex2i((i+1)*block_size,(j+1)*block_size);
			glVertex2i((i+1)*block_size,j*block_size);
			glEnd();
		}
	}

	vector<FractureMesh> meshes=net_mesh.meshes;
	double delta_seg=1e-1;
	//double factor_ap=1e5,factor_x=Lx/2,factor_y=Ly/2,value,length;
	double factor_ap=1,factor_x=0,factor_y=0,value,length;
	pointcpp<double> p_ori,p1,p2,p_current;
	for (int i=0;(unsigned)i<meshes.size();i++){
		glLineWidth(meshes[i].aperture*factor_ap);
		CgalVector2D vect(meshes[i].p_ori.p, meshes[i].p_tar.p);
		p_ori=CgalPoint2D(meshes[i].p_ori.p);
		length=CGAL::sqrt(vect.squared_length());
		Nx=length/delta_seg;
		for (int j=0;j<Nx;j++){
			value=AnalyticalSolutionFractureClassic(delta_seg*(j+0.5),0,length,DFN_potential[meshes[i].p_ori.index],DFN_potential[meshes[i].p_tar.index]);
			p1.i=p_ori.i+vect.x()*delta_seg/length*j;
			p1.j=p_ori.j+vect.y()*delta_seg/length*j;
			p2.i=p_ori.i+vect.x()*delta_seg/length*(j+1);
			p2.j=p_ori.j+vect.y()*delta_seg/length*(j+1);
			glBegin(GL_LINES);
			SetColor(value,min_visu,max_visu);	// line color
			glVertex2f(p1.i+factor_x,p1.j+factor_y);
			glVertex2f(p2.i+factor_x,p2.j+factor_y);
			glEnd();
		}
	}

	int nb_ticks=10;int lin_log=0;
	draw_legend_scale(Lx,Ly,min_visu,max_visu,nb_ticks,lin_log,legend_step);

	glFlush();

}

void EPMDFNVisu::display_EPMDFN_Potential_Matrix(){
	// 1. window setting
	double legend_step=20;
	glClear(GL_COLOR_BUFFER_BIT);	// clear the window with current clearing color
	glTranslatef(-1,-1,-1);	// change coordinate reference (translate it from the center to the bottom left of the domain)

	int Nx=EPM_potential.size1(),Ny=EPM_potential.size2();
	//glScalef(2./Nx,2./Ny,1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	double Lx=L,Ly=L;
	glScalef(2./(Lx+legend_step),2./(Ly+legend_step),1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	double block_size=L/Nx;

	//block_size=0.1;

	// 2. value display

	// draw a polygon for each cell with a color representing its value
	glPolygonMode(GL_FRONT,GL_FILL);	// polygon with a color
	// white polygon for the back
	glBegin(GL_POLYGON);
	SetColor(-1.);
	glVertex2i(0,0);
	glVertex2i(0,Ly+legend_step);
	glVertex2i(Lx+legend_step,Ly+legend_step);
	glVertex2i(Lx+legend_step,0);
	glEnd();

	double min_visu=min(min_matrix(EPM_potential),min_vector(DFN_potential)),
			max_visu=max(max_matrix(EPM_potential),max_vector(DFN_potential));
	//double min_visu=0,max_visu=1;
	cout << "min_visu=" << min_visu << " and max_visu=" << max_visu << endl;
	// plot each cell
	for (int i=0;i<EPM_potential.size1();i++){
		for (int j=0;j<EPM_potential.size2();j++){
			glBegin(GL_POLYGON);
			SetColor(EPM_potential(i,j),min_visu,max_visu);	// cell color
			//SetColor(EPM_potential(i,j)/max_visu,min_visu/max_visu,1);	// cell color
			glVertex2f(i*block_size,j*block_size+0.5*legend_step);
			glVertex2f(i*block_size,(j+1)*block_size+0.5*legend_step);
			glVertex2f((i+1)*block_size,(j+1)*block_size+0.5*legend_step);
			glVertex2f((i+1)*block_size,j*block_size+0.5*legend_step);
			glEnd();
		}
	}

	// fracture plot
	double delta_seg;//=1e-1;
	double alpha_block;
	pair<double,double> indices; int N_seg=10;
	//double factor_ap=1e5,factor_x=Lx/2,factor_y=Ly/2,value,length;
	double factor_ap=1,factor_x=0,factor_y=0,value,length;
	pointcpp<double> p_ori,p1,p2,p_current;
	for (SubNetworkMap::iterator it1=subnetwork_map.begin();it1!=subnetwork_map.end();it1++){
		indices=return_indices(it1->first,Ny);
		//alpha_block=ReturnAlphaBlock(cond_matrix,block_size,block_size);
		alpha_block=AlphaBlock(indices.first,indices.second);
		for (vector<FractureMesh>::iterator it2=it1->second.meshes.begin();it2!=it1->second.meshes.end();it2++){
			glLineWidth(it2->aperture*factor_ap);
			CgalVector2D vect(it2->p_ori.p,it2->p_tar.p);
			p_ori=CgalPoint2D(it2->p_ori.p);
			length=CGAL::sqrt(vect.squared_length());
			delta_seg=length/N_seg;
			//N_seg=length/delta_seg;
			for (int j=0;j<N_seg;j++){
				value=AnalyticalSolutionFractureMatrix(delta_seg*(j+0.5),length,DFN_potential[it2->p_ori.index],DFN_potential[it2->p_tar.index],alpha_block/it2->conductivity,EPM_potential(indices.first,indices.second));
				p1.i=p_ori.i+vect.x()*delta_seg/length*j;
				p1.j=p_ori.j+vect.y()*delta_seg/length*j;
				p2.i=p_ori.i+vect.x()*delta_seg/length*(j+1);
				p2.j=p_ori.j+vect.y()*delta_seg/length*(j+1);
				glBegin(GL_LINES);
				SetColor(value,min_visu,max_visu);	// line color
				//SetColor(value/max_visu,min_visu/max_visu,1);	// line color
				glVertex2f(p1.i+factor_x,p1.j+factor_y+0.5*legend_step);
				glVertex2f(p2.i+factor_x,p2.j+factor_y+0.5*legend_step);
				glEnd();
			}
		}
	}

	int nb_ticks=10;int lin_log=0;
	draw_legend_scale(Lx,Ly,min_visu,max_visu,nb_ticks,lin_log,legend_step);
	//draw_legend_scale(Lx,Ly,0,1,nb_ticks,lin_log,legend_step);

	glFlush();

}

void EPMDFNVisu::display_EPMDFN_Current(){
	// 1. window setting
	glClear(GL_COLOR_BUFFER_BIT);	// clear the window with current clearing color
	glTranslatef(-1,-1,-1);	// change coordinate reference (translate it from the center to the bottom left of the domain)

	int Nx=EPM_potential.size1(),Ny=EPM_potential.size2();
	//glScalef(2./Nx,2./Ny,1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	double Lx=L,Ly=L;
	glScalef(2./Lx,2./Ly,1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	double block_size=L/Nx;
	// 2. value display
	double min_visu=min(min_matrix(EPM_potential),min_vector(DFN_potential)),
			max_visu=max(max_matrix(EPM_potential),max_vector(DFN_potential));
	// plot each cell
	ublas_matrix matrix_current_x,matrix_current_y;
	cout << "WARNING: electric current not computed in display_EPMDFN_Current" << endl;
	//ElectricCurrentComputation(param,EPM_potential,matrix_current_x,matrix_current_y);
	double value,delta_y=param.simu_param.DDP_param.ReturnDeltay(Ly);

	double delta_seg;//=1e-1;
	double alpha_block;
	pair<double,double> indices; int N_seg=10;
	//double factor_ap=1e5,factor_x=Lx/2,factor_y=Ly/2,value,length;
	double factor_ap=1,factor_x=0,factor_y=0,length,pot1,pot2;
	pointcpp<double> p_ori,p1,p2,p_current;

	vector<double> values;
	//double min_visu=0,max_visu=0.2;
	for (int i=0;i<matrix_current_x.size1();i++){
		for (int j=0;j<matrix_current_x.size2();j++){
			glBegin(GL_POLYGON);
			value=fabs(matrix_current_x(i,j))/delta_y;
			SetColor(value,min_visu,max_visu);	// cell color
			values.push_back(value);
			glVertex2f((i+0.5)*block_size,j*block_size);
			glVertex2f((i+0.5)*block_size,(j+1)*block_size);
			glVertex2f((i+1.5)*block_size,(j+1)*block_size);
			glVertex2f((i+1.5)*block_size,j*block_size);
			glEnd();
		}
	}

	// electric current in the y-direction
	for (int i=0;i<matrix_current_y.size1();i++){
		for (int j=0;j<matrix_current_y.size2();j++){
			glBegin(GL_POLYGON);
			value=matrix_current_y(i,j);
			SetColor(value,min_visu,max_visu);	// cell color
			values.push_back(value);
			glVertex2f((i+0.5)*block_size,j*block_size);
			glVertex2f((i+0.5)*block_size,(j+1)*block_size);
			glVertex2f((i+1.5)*block_size,(j+1)*block_size);
			glVertex2f((i+1.5)*block_size,j*block_size);
			glEnd();
		}
	}

	for (SubNetworkMap::iterator it1=subnetwork_map.begin();it1!=subnetwork_map.end();it1++){
		indices=return_indices(it1->first,Ny);
		//alpha_block=ReturnAlphaBlock(cond_matrix,block_size,block_size);
		alpha_block=AlphaBlock(indices.first,indices.second);
		for (vector<FractureMesh>::iterator it2=it1->second.meshes.begin();it2!=it1->second.meshes.end();it2++){
			glLineWidth(it2->aperture*factor_ap);
			CgalVector2D vect(it2->p_ori.p,it2->p_tar.p);
			p_ori=CgalPoint2D(it2->p_ori.p);
			length=CGAL::sqrt(vect.squared_length());
			delta_seg=length/N_seg;
			//N_seg=length/delta_seg;
			for (int j=0;j<N_seg;j++){
				pot1=AnalyticalSolutionFractureMatrix(delta_seg*j,length,DFN_potential[it2->p_ori.index],DFN_potential[it2->p_tar.index],alpha_block/it2->conductivity,EPM_potential(indices.first,indices.second));
				pot2=AnalyticalSolutionFractureMatrix(delta_seg*(j+1),length,DFN_potential[it2->p_ori.index],DFN_potential[it2->p_tar.index],alpha_block/it2->conductivity,EPM_potential(indices.first,indices.second));
				//value=fabs(it2->aperture*MetricElecFlow(it2->conductivity,pot2,pot1,delta_seg));
				value=fabs(MetricElecFlow(it2->conductivity,pot2,pot1,delta_seg));
				values.push_back(value);
				p1.i=p_ori.i+vect.x()*delta_seg/length*j;
				p1.j=p_ori.j+vect.y()*delta_seg/length*j;
				p2.i=p_ori.i+vect.x()*delta_seg/length*(j+1);
				p2.j=p_ori.j+vect.y()*delta_seg/length*(j+1);
				glBegin(GL_LINES);
				SetColor(value,min_visu,max_visu);	// line color
				glVertex2f(p1.i+factor_x,p1.j+factor_y);
				glVertex2f(p2.i+factor_x,p2.j+factor_y);
				glEnd();
			}
		}
	}

	cout << "min_visu=" << min_vector(values) << endl;
	cout << "max_visu=" << max_vector(values) << endl;

	glFlush();
}
*/
