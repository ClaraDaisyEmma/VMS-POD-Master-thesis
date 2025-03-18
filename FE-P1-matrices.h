//this file provides the FE matrices for P1
#include <iostream>
#include "personal-copy-eigen-3.4.0/Eigen/Dense" 

using namespace Eigen;

MatrixXd B_0(int dim){ 
    if (dim < 2) {
        std::cerr << "dimension error" << std::endl;
        exit(1);
    }
    MatrixXd result_B0=MatrixXd::Zero(dim,dim);
    double h=1.0/(dim+1);
   
    result_B0(0,0)=2/h;
    result_B0(0,1)=-1/h;
    for(int l=1; l<=dim-2; l++){
        result_B0(l,l-1)=-1/h;
        result_B0(l,l)=2/h;
        result_B0(l,l+1)=-1/h;
    }
    result_B0(dim-1,dim-2)=-1/h;
    result_B0(dim-1,dim-1)=2/h;
    
    return result_B0; 
}

MatrixXd B_2(int dim){
    if (dim < 2) {
        std::cerr << "dimension error" << std::endl;
        exit(1);
    }
    MatrixXd result_B2=MatrixXd::Zero(dim,dim);
    double h=1.0/(dim+1);

    result_B2(0,0)=(2*h)/3;
    result_B2(0,1)=h/6;
    for(int l=1; l<=dim-2; l++){
        result_B2(l,l-1)=h/6;
        result_B2(l,l)=(2*h)/3;
        result_B2(l,l+1)=h/6;             
    }
    result_B2(dim-1,dim-2)=h/6;
    result_B2(dim-1,dim-1)=(2*h)/3;
    
    return result_B2; 
}