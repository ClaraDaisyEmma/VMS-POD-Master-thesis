//this file saves matrices as a csv file
#include <iostream>
#include "personal-copy-eigen-3.4.0/Eigen/Dense" 
#include <fstream> 

using namespace Eigen;

void save(MatrixXd matrix,const std::string& name){
    int dim_time= matrix.cols(), dim_space=matrix.rows();
    double x, t, space_step_length=1.0/(dim_time-1 ), h=1.0/(dim_space+1);
  
    std::ofstream os (name);
    os << "x_position, t_position, value ,"<< std::endl;
    for(int l1=0; l1<=dim_time-1; l1++){
        t = l1*space_step_length;
        x=0;
        os << x<< ", " << t << ", " << 0 << ","<< std::endl; //Randwert
        for (int l2=1; l2<=dim_space; l2++){        
            x = l2*h;
            os << x<< ", " << t << ", " << matrix(l2-1,l1)<< ","<< std::endl;
        }
        x=1;
        os << x << ", " << t << ", " << 0<< "," << std::endl; //Randwert
    }
    os.close();
}