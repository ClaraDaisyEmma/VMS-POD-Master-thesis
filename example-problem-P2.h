//The file provides the snapshots in the form of the exact solution at the nodes and timesteps, 
//the matrix for the convective field and the matrix of the forcing term values for t and x.
#include <iostream>
#include "personal-copy-eigen-3.4.0/Eigen/Dense" 

using namespace Eigen;

const double z=0.01;            //parameter of the example problem

MatrixXd Snapshots(int dim_space, int dim_time){
    if (dim_space < 1 || dim_time < 1) {
        std::cerr << "dimension error" << std::endl;
        exit(1);
    }
    MatrixXd result_Snap=MatrixXd::Zero(dim_space,dim_time+1);   
    double space_step_length=1.0/(dim_space+1), time_step_length= 1.0/(dim_time),t,x;                                 
   
    for(int l1=0;l1<=dim_space-1; l1++){
        x=(l1+1)*space_step_length;
        for(int l2=0; l2<= dim_time; l2++){
            t=l2*time_step_length;
            result_Snap(l1,l2)= exp(-t+1-100*((x-t)*(x-t)))*(x-(exp(x/z)-1)/(exp(1/z)-1));
        }
    }
    return result_Snap;
}

MatrixXd N_k(int dim){ 
    if (dim < 3) {
        std::cerr << "dimension error" << std::endl;
        exit(1);
    }
    MatrixXd result_Nk=MatrixXd::Zero(dim,dim);
    // Define matrix values:
    result_Nk(0,1)=-2.0/3.0;
    result_Nk(1,0)=2.0/3.0;
    result_Nk(1,2)=-2.0/3.0;
    result_Nk(1,3)=1.0/6.0;
    for(int l=2; l<=dim-3; l++){
        if(l%2==0){      
            result_Nk(l,l-1)=2.0/3.0;
            result_Nk(l,l+1)=-2.0/3.0;
        }
        else{    
            result_Nk(l,l-2)=-1.0/6.0;
            result_Nk(l,l-1)=2.0/3.0;
            result_Nk(l,l+1)=-2.0/3.0;
            result_Nk(l,l+2)=1.0/6.0;
        }   
    }
    result_Nk(dim-2,dim-4)=-1.0/6.0;
    result_Nk(dim-2,dim-3)=2.0/3.0;
    result_Nk(dim-2,dim-1)=-2.0/3.0;
    result_Nk(dim-1,dim-2)=2.0/3.0;
  
    return result_Nk;
}

MatrixXd forcing_term(int dim_space, int dim_time, double epsylon){ 
    if (dim_space < 1 || dim_time < 1) {
        std::cerr << "dimension error" << std::endl;
        exit(1);
    }
    MatrixXd result_f=MatrixXd::Zero(dim_space,dim_time+1);
    double space_step_length=1.0/(dim_space+1), time_step_length= 1.0/(dim_time),t,x;
    
    for(int l1=0;l1<=(dim_space-1); l1++){
        x=(l1+1)*space_step_length;
        for(int l2=0; l2<= dim_time; l2++){
            t=(l2)*time_step_length;
            result_f(l1,l2)=exp(-t+1-100*(x-t)*(x-t))* (-epsylon*( (40000*(x-t)*(x-t)-200)*(x-(exp(x/z)-1)/(exp(1/z)-1))
            -400*(x-t)*(1-exp(x/z)/(z*(exp(1/z)-1)))-exp(x/z)/(z*z*(exp(1/z)-1)))+(1-exp(x/z)/(z*(exp(1/z)-1))));     
        }
    }
    return result_f;
}