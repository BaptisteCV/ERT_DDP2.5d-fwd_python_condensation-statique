/*
 * UblasStructures.cpp
 *
 *  Created on: Apr 23, 2013
 *      Author: delphineroubinet
 */

#include "UblasStructures.h"
#include "UtilitaryStorage.h"

using namespace std;

void print_matrix(ublas_matrix ublas_mat,int option){
	int Nx=ublas_mat.size1(),Ny=ublas_mat.size2();
	cout << "Matrix values = " << endl;
	if (option==0){
		for (int i=0;i<Nx;i++){
			for (int j=0;j<Ny;j++){
				cout << ublas_mat(i,j) << " ";
			}
			cout << endl;
		}
	}
	else if (option==1){
		for (int j=Ny-1;j>=0;j--){
			for (int i=0;i<Nx;i++){
				cout << ublas_mat(i,j) << " ";
			}
			cout << endl;
		}
	}
	else{cout << "WARNING in print_matrix (UblasStructures.cpp): option not defined" << endl;}
}

void print_vector(ublas_vector ublas_vect){
	cout << "Vector values = " << endl;
	for (int i=0;i<ublas_vect.size();i++){
		cout << ublas_vect(i) << endl;
	}
}

ublas_vector InitializationVector(int size){
	return boost::numeric::ublas::zero_vector<double>(size);
}

ublas_vector InitializationVector(int size,double value){
	//ublas_vector vect(size,value);
	ublas_vector vect(size);
	for (int i=0;i<vect.size();i++){
		vect(i)=value;
	}
	return vect;
}

ublas_matrix InitializationMatrix(int size1,int size2){
	return boost::numeric::ublas::zero_matrix<double>(size1,size2);
}

ublas_matrix InitializationMatrix(int size1,int size2,double value){
	//ublas_matrix mat(size1,size2,value);
	ublas_matrix mat(size1,size2);
	for (int i=0;i<mat.size1();i++){
		for (int j=0;j<mat.size2();j++){
			mat(i,j)=value;
		}
	}
	return mat;
}

ublas_matrix RelativeMatrix(ublas_matrix matrix_potential,ublas_matrix potential_host){
	ublas_matrix relative_potential(matrix_potential.size1(),matrix_potential.size2());
	if (matrix_potential.size1()!=potential_host.size1()||matrix_potential.size2()!=potential_host.size2()){
		cout << "WARNING in RelativeMatrix (UblasStructures.cpp): matrix have different sizes" << endl;
		return relative_potential;
	}
	for (int i=0;i<matrix_potential.size1();i++){
		for (int j=0;j<matrix_potential.size2();j++){
			relative_potential(i,j)=(potential_host(i,j)-matrix_potential(i,j))/potential_host(i,j);
			//relative_potential(i,j)=matrix_potential(i,j)/potential_host(i,j);
		}
	}
	return relative_potential;
}

void LinesReplacementMatrix1(int line_num,ublas_matrix& init_matrix,ublas_matrix block_matrix){
	// 1. Variables
	int Nx_init = init_matrix.size1(), Ny_init = init_matrix.size2();
	int Nx_block = block_matrix.size1(), Ny_block = block_matrix.size2();
	int Nx_total = Nx_init-1+Nx_block, Ny_total = Ny_init-1+Ny_block;
	ublas_matrix new_matrix(Nx_total,Ny_total);	// assumes that the matrix is initialized at 0
	// 2. Affectation of the initial values before the specific node
	for (int i=0;i<line_num;i++){
		for (int j=0;j<Ny_init;j++){
			new_matrix(i,j) = init_matrix(i,j);
		}
	}
	// 3. Insertion of the new block at the specific line
	for (int i=0;i<Nx_block;i++){
		for (int j=0;j<Ny_block;j++){
			new_matrix(line_num+i,line_num+j) = block_matrix(i,j);
		}
	}
	// 4. Affectation of the initial values after the specific node
	for (int i=line_num;i<Nx_init;i++){
		for (int j=0;j<Ny_init;j++){
			new_matrix(line_num+Nx_block+i,line_num+Ny_block+j) = init_matrix(i,j);
		}
	}
	init_matrix = new_matrix;
}

// write block_matrix into init_matrix from line_num to the init_matrix size
void LinesReplacementMatrix2(int line_num,ublas_matrix& init_matrix,ublas_matrix block_matrix){
	for (int i=line_num;i<init_matrix.size1();i++){
		for (int j=line_num;j<init_matrix.size1();j++){
			init_matrix(i,j) = init_matrix(i-line_num,j-line_num);
		}
	}
}

// replace the line_num line by nb_lines with the value value
void LinesReplacementVector(int line_num,ublas_vector& init_vector,ublas_vector block_vector){

	// 1. Variables
	int Nx_init = init_vector.size();
	int Nx_block = block_vector.size();
	ublas_vector new_vector(Nx_init-1+Nx_block);	// assumes that the vector is initialized at 0
	// 2. Affectation of the initial values before the specific line
	for (int i=0;i<line_num;i++){
		new_vector(i) = init_vector(i);
	}
	// 3. Insertion of the new block
	for (int i=0;i<Nx_block;i++){
		new_vector(line_num+i) = block_vector(i);
	}
	// 4. Affectation of the initial values after the specific line
	for (int i=line_num;i<Nx_init;i++){
		new_vector(line_num+Nx_block+i) = init_vector(i);
	}

	init_vector = new_vector;
}

double min_vector(ublas_vector vect){
	ublas_vector::iterator it=vect.begin();
	double min_value=*it;
	for (int i=0;i<vect.size();i++){
		if (min_value>vect[i]){min_value=vect[i];}
	}
	return min_value;
}

double max_vector(ublas_vector vect){
	ublas_vector::iterator it=vect.begin();
	double max_value=*it;
	for (int i=0;i<vect.size();i++){
		if (max_value<vect[i]){max_value=vect[i];}
	}
	return max_value;
}

double min_matrix(ublas_matrix matrix){
	double min_value=matrix(0,0);
	for (int i=0;i<matrix.size1();i++){
		for (int j=0;j<matrix.size2();j++){
			if (min_value>matrix(i,j)){min_value=matrix(i,j);}
		}
	}
	return min_value;
}

double max_matrix(ublas_matrix matrix){
	double max_value=matrix(0,0);
	for (int i=0;i<matrix.size1();i++){
		for (int j=0;j<matrix.size2();j++){
			if (max_value<matrix(i,j)){max_value=matrix(i,j);}
		}
	}
	return max_value;
}

double ErrorComputation(ublas_matrix matrix1,ublas_matrix matrix2){
	double error=0;
	if (matrix1.size1()!=matrix2.size1()||matrix1.size2()!=matrix2.size2()){cout << "WARNING in ErrorComputation (UblasStructures.cpp): matrix has to have the same sizes" << endl; return error;}
	for (int i=0;i<matrix1.size1();i++){
		for (int j=0;j<matrix1.size2();j++){
			error+=fabs(matrix1(i,j)-matrix2(i,j));
		}
	}
	return error/(matrix1.size1()*matrix1.size2());
}

double ErrorComputationL2(ublas_matrix matrix1,ublas_matrix matrix2){
	double error=0;
	if (matrix1.size1()!=matrix2.size1()||matrix1.size2()!=matrix2.size2()){cout << "WARNING in ErrorComputation (UblasStructures.cpp): matrix has to have the same sizes" << endl; return error;}
	for (int i=0;i<matrix1.size1();i++){
		for (int j=0;j<matrix1.size2();j++){
			error+=(matrix1(i,j)-matrix2(i,j))*(matrix1(i,j)-matrix2(i,j));
		}
	}
	return error/(matrix1.size1()*matrix1.size2());
}

bool PositiveVector(ublas_vector vect){
	for (int i=0;i<vect.size();i++){
		if (vect[i]<0){return false;}
	}
	return true;
}

bool PositiveMatrix(ublas_matrix matrix){
	for (int i=0;i<matrix.size1();i++){
		for (int j=0;j<matrix.size2();j++){
			if (matrix(i,j)<0){return false;}
		}
	}
	return true;
}

ublas_matrix ExtractMatrix(ublas_matrix init_matrix,int i1, int i2, int j1, int j2){
	ublas_matrix extract_matrix(i2-i1+1,j2-j1+1);
	if (i1<0||j1<0||i2>init_matrix.size1()-1||j2>init_matrix.size2()-1){
		cout << "WARNING in ExtractMatrix (UblasStructures.cpp): matrix extraction not possible" << endl;
	}
	else{
		int cpt_i,cpt_j;
		cpt_i=0;
		for (int i=i1;i<=i2;i++){
			cpt_j=0;
			for (int j=j1;j<=j2;j++){
				extract_matrix(cpt_i,cpt_j)=init_matrix(i,j);
				cpt_j++;
			}
			cpt_i++;
		}
	}
	return extract_matrix;
}

double AveragedMatrix(ublas_matrix matrix){
	double res=0;
	for (int i=0;i<matrix.size1();i++){
		for (int j=0;j<matrix.size2();j++){
			res+=matrix(i,j);
		}
	}
	return res/(matrix.size1()*matrix.size2());
}

map<int,double> MatrixToMap(ublas_matrix matrix,int Ny){
	map<int,double> res;
	int index;
	for (int i=0;i<matrix.size1();i++){
		for (int j=0;j<matrix.size2();j++){
			index=return_index(i,j,Ny);
			res[index]=matrix(i,j);
		}
	}
	return res;
}

bool ConstantValue(ublas_matrix matrix,double& cst){
	double res=matrix(0,0);
	for (int i=0;i<matrix.size1();i++){
		for (int j=0;j<matrix.size2();j++){
			if (res!=matrix(i,j)){
				return false;
			}
		}
	}
	cst=res;
	return true;
}
