/*

 * DisplayResults.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet*/

/*
#include "DisplayResults.h"
#include <glut.h>
#include <gl.h>
#include <glu.h>

using namespace std;

DisplayResults::DisplayResults() {}
DisplayResults::~DisplayResults() {}

DFNVisu* current_results;
EPMVisu* current_results2;
EPMDFNVisu* current_res;

void glutReshapeFunc(int width, int height){

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

}

void display_DFN(){current_results->display_DFN();}
void display_Potential_Classic(){current_results->display_Potential_Classic();}
//void display_Potential_Radial(){current_results->display_Potential_Radial();}

void SetupRC(){
	glClearColor(1,1,1,1);	// initialize window color to white
}

void keyboardFunc(unsigned char key, int x, int y){
	//current_results->keyboardFunc(key);
	glutPostRedisplay();
}

DisplayResults::DisplayResults(int argc,char **argv,DFNVisu results_){
	// class members affectation
	results = results_;
	current_results = &results_;

	//int factor=5;//20;
	int width=200;//(int)results.L*factor;
	int height=200;//(int)results.L*factor;

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(width,height);
	glutInitWindowPosition(200,200);
	glutInit(&argc, argv);
	if (results.option=="DFN"){glutCreateWindow("Discrete Fracture Network");}
	else if (results.option=="POTENTIAL_CLASSIC"){glutCreateWindow("Electric Potential in DFN");}
	else if (results.option=="POTENTIAL_RADIAL"){glutCreateWindow("Electric Potential in Radial_DFN");}
	else {cout << "WARNING in DisplayResults (DisplayResulys.cpp): option not implemented" << endl;}
	glutKeyboardFunc(keyboardFunc);
	glutReshapeFunc(width,height);
	if (results.option=="DFN"){glutDisplayFunc(display_DFN);}
	else if (results.option=="POTENTIAL_CLASSIC"){glutDisplayFunc(display_Potential_Classic);}
	//else if (results.option=="POTENTIAL_RADIAL"){glutDisplayFunc(display_Potential_Radial);}
	else {cout << "WARNING in DisplayResults (DisplayResulys.cpp): option not implemented" << endl;}
	glutMainLoop();

}

void display_EPM_Potential(){current_results2->display_EPM_Potential();}
void display_EPM_Current_x(){current_results2->display_EPM_Current_x();}
void display_EPM_Current_y(){current_results2->display_EPM_Current_y();}
DisplayResults::DisplayResults(int argc, char **argv, EPMVisu results_,string option) {

	// class members affectation
	results2 = results_;
	current_results2 = &results_;

	int factor=3;
	int width=0,height=0;
	if (option=="Potential"){
		width=300;//results2.Lx*factor;
		height=300;//results2.Ly*factor;
	}
	else if (option=="Current_x"){
		width=results2.Current_x.size1()*factor;
		height=results2.Current_x.size2()*factor;
	}
	else if (option=="Current_y"){
		width=results2.Current_y.size1()*factor;
		height=results2.Current_y.size2()*factor;
	}
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(width,height);
	glutInitWindowPosition(200,200);
	glutInit(&argc, argv);
	if (option=="Potential"){glutCreateWindow("Electric Potential");}
	else if (option=="Current_x"){glutCreateWindow("Electric Current_x");}
	else if (option=="Current_y"){glutCreateWindow("Electric Current_y");}
	glutKeyboardFunc(keyboardFunc);
	glutReshapeFunc(width,height);
	if (option=="Potential"){glutDisplayFunc(display_EPM_Potential);}
	else if (option=="Current_x"){glutDisplayFunc(display_EPM_Current_x);}
	else if (option=="Current_y"){glutDisplayFunc(display_EPM_Current_y);}
	glutMainLoop();
}


void display_EPMDFN_Potential(){current_res->display_EPMDFN_Potential();}
void display_EPMDFN_Potential_Matrix(){current_res->display_EPMDFN_Potential_Matrix();}
void display_EPMDFN_Current(){current_res->display_EPMDFN_Current();}

DisplayResults::DisplayResults(int argc,char **argv,EPMDFNVisu res_EPMDFN_,string option){
	// class members affectation
	res_EPMDFN=res_EPMDFN_;
	current_res=&res_EPMDFN_;

	int factor=1;
	int width=(int)results.L*factor;
	int height=(int)results.L*factor;

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(width,height);
	glutInitWindowPosition(200,200);
	glutInit(&argc, argv);
	glutCreateWindow("Fracture-Matrix Potential");
	glutKeyboardFunc(keyboardFunc);
	glutReshapeFunc(width,height);
	if (option=="Matrix"){glutDisplayFunc(display_EPMDFN_Potential_Matrix);}
	else if (option=="Current"){glutDisplayFunc(display_EPMDFN_Current);}
	else {glutDisplayFunc(display_EPMDFN_Potential);}
	glutMainLoop();

}
*/


