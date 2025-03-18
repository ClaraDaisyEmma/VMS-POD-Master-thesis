//this file provides the FE matrices for P2
#include <iostream>
#include "personal-copy-eigen-3.4.0/Eigen/Dense" 

using namespace Eigen;

MatrixXd B_0(int dim){ 
    if (dim < 2) {
        std::cerr << "dimension error" << std::endl;
        exit(1);
    }
    MatrixXd result_B0=MatrixXd::Zero(dim,dim);
    double h=2.0/(dim+1);

    result_B0(0,0)=16/(3*h);
    result_B0(0,1)=-8/(3*h);

    result_B0(1,0)=-8/(3*h);
    result_B0(1,1)=14/(3*h);
    result_B0(1,2)=-8/(3*h);
    result_B0(1,3)=1/(3*h);

    for(int l=2; l<=dim-3; l++){
        if(l%2==0){      
            result_B0(l,l-1)=-8/(3*h);
            result_B0(l,l)=16/(3*h);
            result_B0(l,l+1)=-8/(3*h);
        }
        else{    
            result_B0(l,l-2)=1/(3*h);
            result_B0(l,l-1)=-8/(3*h);
            result_B0(l,l)=14/(3*h);
            result_B0(l,l+1)=-8/(3*h);
            result_B0(l,l+2)=1/(3*h);
        }       
    }     
    result_B0(dim-2,dim-4)=1/(3*h);
    result_B0(dim-2,dim-3)=-8/(3*h);
    result_B0(dim-2,dim-2)=14/(3*h);
    result_B0(dim-2,dim-1)=-8/(3*h);

    result_B0(dim-1,dim-2)=-8/(3*h);
    result_B0(dim-1,dim-1)=16/(3*h);
   
    return result_B0; 
}

MatrixXd B_2(int dim){
    if (dim < 2) {
        std::cerr << "dimension error" << std::endl;
        exit(1);
    }
    MatrixXd result_B2=MatrixXd::Zero(dim,dim);;
    double h=2.0/(dim+1);

    result_B2(0,0)=(16*h)/30;
    result_B2(0,1)=(2*h)/30;

    result_B2(1,0)=(2*h)/30;
    result_B2(1,1)=(8*h)/30;
    result_B2(1,2)=(2*h)/30;
    result_B2(1,3)=(-1*h)/30;

    for(int l=2; l<=dim-3; l++){
        if(l%2==0){     
            result_B2(l,l-1)=(2*h)/30;
            result_B2(l,l)=(16*h)/30;
            result_B2(l,l+1)=(2*h)/30;
        }
        else{
            result_B2(l,l-2)=(-1*h)/30;
            result_B2(l,l-1)=(2*h)/30;
            result_B2(l,l)=(8*h)/30;
            result_B2(l,l+1)=(2*h)/30;
            result_B2(l,l+2)=(-1*h)/30;
        }        
    }
    result_B2(dim-2,dim-4)=(-1*h)/30;
    result_B2(dim-2,dim-3)=(2*h)/30;
    result_B2(dim-2,dim-2)=(8*h)/30;
    result_B2(dim-2,dim-1)=(2*h)/30;

    result_B2(dim-1,dim-2)=(2*h)/30;
    result_B2(dim-1,dim-1)=(16*h)/30;

    return result_B2; 
}