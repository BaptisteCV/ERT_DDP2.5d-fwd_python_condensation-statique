/*
 * Results.cpp
 *
 *  Created on: Jun 29, 2012
 *      Author: delphine
 */


#include "Results.h"
#include "../Utilitaries/Storage/Structures.h"
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;


////////////////////////////////
// Results for DFN Generation //
////////////////////////////////
Results::Results(string output_file_,NetworkMeshes net_meshes_){
	output_file=output_file_;net_meshes=net_meshes_;
	WriteDFNDefinition();
}

void Results::WriteDFNDefinition(){
	std::ofstream output(output_file.c_str(), std::ofstream::out);
	vector<FractureMesh> meshes=net_meshes.meshes;
	//output.precision(30);
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		output << it->p_ori.p.x() << " " << it->p_ori.p.y() << " " << it->p_tar.p.x() << " " << it->p_tar.p.y() << " " << it->aperture << " " << it->conductivity << endl;
	}
	output.close();
}

Results::Results(string output_file_,map<double,double> map_perm_){
	output_file=output_file_;map_perm=map_perm_;
	std::ofstream output(output_file.c_str(), std::ofstream::out);
	output.precision(10);
	for (map<double,double>::iterator it=map_perm.begin();it!=map_perm.end();it++){
		output << it->first << " " << it->second << endl;
	}
	output.close();
}

Results::Results(string output_file_,map<int,double> map_res1_){
	output_file=output_file_;map_res1=map_res1_;
	std::ofstream output(output_file.c_str(), std::ofstream::out);
	//output.precision(10);
	for (map<int,double>::iterator it=map_res1.begin();it!=map_res1.end();it++){
		output << it->first << " " << it->second << endl;
	}
	output.close();
}

Results::Results(string output_file_,map<double,ublas_matrix> map_mat_res_){
	output_file=output_file_;map_mat_res=map_mat_res_;
	std::ofstream output(output_file.c_str(), std::ofstream::out);
	ublas_matrix current_mat;
	for (map<double,ublas_matrix>::iterator it=map_mat_res.begin();it!=map_mat_res.end();it++){
		current_mat=it->second;
		/*for (int i=0;i<current_mat.size1();i++){
			for (int j=0;j<current_mat.size2();j++){
				output << current_mat(i,j) << " ";
			}
			output << endl;
		}*/
		for (int j=current_mat.size2()-1;j>=0;j--){
			for (int i=0;i<current_mat.size1();i++){
				output << current_mat(i,j) << " ";
			}
			output << endl;
		}
	}
	output.close();
}

Results::Results(string output_file_, map<double, eigen_cs::DenseMatrix> map_mat_res_eigen_) {
	output_file = output_file_;
	ofstream output(output_file.c_str(), ofstream::out);
	for (auto& kv : map_mat_res_eigen_) {
		const eigen_cs::DenseMatrix& mat = kv.second;
		for (int j = (int)mat.cols() - 1; j >= 0; j--) {
			for (int i = 0; i < (int)mat.rows(); i++)
				output << mat(i, j) << " ";
			output << "\n";
		}
	}
	output.close();
}

Results::Results(string output_file_,map<double,ublas_vector> map_fract_res_){
	output_file=output_file_;map_fract_res=map_fract_res_;
	std::ofstream output(output_file.c_str(), std::ofstream::out);
	ublas_vector current_vector;
	for (map<double,ublas_vector>::iterator it=map_fract_res.begin();it!=map_fract_res.end();it++){
		current_vector=it->second;
		for (int i=0;i<current_vector.size();i++){
				//output << "(" << i << "," << j << ")" << current_mat(i,j) << endl;
				output << current_vector(i) << endl;
		}
	}
	output.close();
}

Results::Results(string output_file_,map<double,ublas_vector> map_fract_res_,int Nx,int Ny){
        output_file=output_file_;map_fract_res=map_fract_res_;
        std::ofstream output(output_file.c_str(), std::ofstream::out);
        ublas_vector current_vector;
        for (map<double,ublas_vector>::iterator it=map_fract_res.begin();it!=map_fract_res.end();it++){
                current_vector=it->second;
                for (int j=Ny-1;j>=0;j--){
                        for (int i=0;i<Nx;i++){
                                output << current_vector(return_index(i,j,Ny)) << " ";
                        }
                        output << endl;
                }
        }
        output.close();
}

Results::Results(string output_file_,ublas_matrix res1_){
	output_file=output_file_;res1=res1_;
	std::ofstream output(output_file.c_str(), std::ofstream::out);
	for (int i=0;i<res1.size1();i++){
		for (int j=0;j<res1.size2();j++){
			//output << "(" << i << "," << j << ")" << res1(i,j) << endl;
			output << res1(i,j) << " ";
		}
		output << endl;
	}
	output.close();
}

Results::Results(string output_file_,string output_file1_,map<int,double> map_res1_,map<int,double> map_res2_,map<int,double> map_res3_){
	output_file=output_file_;map_res1=map_res1_;
	output_file1=output_file1_;map_res2=map_res2_;
	map_res3=map_res3_;
	std::ofstream output(output_file.c_str(), std::ofstream::out);
	//output.precision(10);
	for (map<int,double>::iterator it=map_res1.begin();it!=map_res1.end();it++){
		output << map_res3[it->first] << " " << it->second << endl;
	}
	output.close();
	std::ofstream output1(output_file1.c_str(), std::ofstream::out);
	//output.precision(10);
	for (map<int,double>::iterator it=map_res2.begin();it!=map_res2.end();it++){
		output1 << map_res3[it->first] << " " << it->second << endl;
	}
	output1.close();
}
