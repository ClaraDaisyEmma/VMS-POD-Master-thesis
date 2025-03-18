#include <iostream>
#include "personal-copy-eigen-3.4.0/Eigen/Dense"                   
#include "save-csv.h" 
#include "example-problem-P2.h"
#include "FE-P2-matrices.h"

using namespace Eigen;

//This function does the proper orthogonal decomposition.
MatrixXd POD(MatrixXd Snapshots, int r, MatrixXd X){
    MatrixXd snapshotproduct=(Snapshots.transpose()*X*Snapshots), result_POD=MatrixXd::Zero(Snapshots.rows(),r);
    int dim_time= Snapshots.cols(); 
    double lamda;                           
 
    SelfAdjointEigenSolver<MatrixXd> es(snapshotproduct);   //Eigenvalue solver for self adjoint matrices, 
                                                            //provides the eigenvalues in ascending order of size;
    //creation of the modes:
    for(int k=0; k<r ; k++){
        lamda=(es.eigenvalues()(dim_time-1-k));
        if (lamda < 1e-12) {  // Avoid division by zero
            std::cerr << "Warning: Small eigenvalue encountered for the POD modes!" << std::endl;
            continue;
        }
        result_POD(all,k)= 1/sqrt(lamda)*Snapshots*es.eigenvectors().col(dim_time-1-k);         
    }
    return result_POD;
}

//This function provides a matrix with an ON basis of L^R. 
MatrixXd resolved_ON(MatrixXd PODbasis,int R, MatrixXd B){
    MatrixXd result_ON=MatrixXd::Zero(PODbasis.rows(),R);
    double divisor;
    for(int i=0; i<R; i++){
        result_ON(all,i)=PODbasis(all,i);
        int j=0;
        while(j<i){
            result_ON(all,i)=result_ON(all,i)-(result_ON(all,j).transpose()*B*PODbasis(all,i))*result_ON(all,j);
            j=j+1;
        }
        divisor=sqrt(((result_ON(all,i)).transpose()*B*result_ON(all,i)));
        if (divisor < 1e-12){   // Avoid division by zero
            std::cerr << "Warning: risk to divide by zero" << std::endl;
            continue;
        }
        result_ON(all,i)=result_ON(all,i)/divisor;
    }
    return result_ON;
}

int main(){
    //M= number of FE cells; N= number of time steps; r= number of POD-modes; R= numer of resolved-scale POD-modes;
    //nodes= numer of FE nodes; Eps= diffusion coefficient; DeltaT= lenght of the time intervall; 
    //Eps_plus= artivicial viscosity coefficient; error= averarge error of the POD solution and the snapshots;
    constexpr int M=250, N=500, r=40, R=20;                                                        
    int nodes=2*M-1;                                                                               
    constexpr  double Eps=1e-10, DeltaT=1.0/N ,Eps_plus=4.0/M;      
    double  error=0 ;                             
                                                       
    //VMS_POD_solution= matrix with VMS-POD solution vectors as columns; snapshots= matrix with snapshots as columns; 
    //B0= matrix with H^1_0 inner product of the FE basis funktions; B2= matrix with L^2 inner product of the FE basis funktions;
    //LHS= left hand side of the linear system; F= matrix with the forcing term solutions f(t,x), row i: for x_i, column j: for t_j; 
    //Nk= matrix with (b_i', b_j)_L^2, whith the FE basis funktions b_i, b_j; POD_modes= matrix with POD modes as columns;
    //Rev_transformation= revers transforamtion of the VMS-POD solution onto the FE ansatz space;
    //unresolved_scales= matrix for the artificial viscosity; resolved_scales= matrix of the ON basis of L^R;
    MatrixXd    VMS_POD_solution=MatrixXd::Zero(r, N+1), snapshots=MatrixXd::Zero(nodes,N+1), 
                B0, B2, LHS, F, Nk, POD_modes, Rev_transformation, unresolved_scales, resolved_scales;               

    //RHS= right hand side of the linear system; AV= averrage mean value of the sfull order solution; 
    //AVRHS= term on the RHS that is influenced by the AV
    VectorXd RHS, AV=VectorXd::Zero(nodes), AVRHS;
           
    if(R==0 || R>r)
    {
        std::cout << "Chose $r>R>0$"<< std::endl;
        return EXIT_FAILURE;
    }
    //The hadder file FE-P2-matrices.h includes the functions B_0(int dim) and B_2(int dim), which provides (nodex x nodes)-matrices.
    B0=B_0(nodes);
    B2=B_2(nodes);
    //The hadder file example-problem-P2.h includes the functions Snapshots(int numbNotes, int timesteps) and 
    //forcing_term(int dim_space, int dim_time, double epsylon), which provides (nodex x N+1)-matrices, 
    //aswell as N_k(int dim), which provides a (nodex x nodes)-matrix.
    snapshots=Snapshots(nodes, N);
    Nk=N_k(nodes);
    F=forcing_term(nodes,N, Eps);
    
    //Computing the mean value and subtract from the snapshots:
    for(int k=0; k<N+1; k++){
        AV=AV+snapshots(all,k);
    }
    AV=AV/((double)N +1);
    for(int k=0; k<N+1; k++){
        snapshots(all, k)=snapshots(all,k)-AV;
    }
    //The function POD(MatrixXd Snapshots, int r, MatrixXd X) provides an (nodes x r)-matrix of the POD basis functions
    POD_modes=POD(snapshots,r,(B0+B2));
    //The function resolved_ON(MatrixXd PODbasis,int R, MatrixXd B) provides an (nodes x R)-matrix
    resolved_scales=resolved_ON(POD_modes,R,B0);
    //Computing the artificial viscosity matrix:
    unresolved_scales=POD_modes- (resolved_scales*(resolved_scales.transpose()*(B0)*POD_modes));
    //Computing the startingpoint u^r_0
    VMS_POD_solution(all,0)=POD_modes.transpose()*(B0+B2)*snapshots(all,0);
   
    AVRHS=DeltaT*POD_modes.transpose()*(B2+Nk.transpose()+Eps*B0)*AV; 
    LHS = ((DeltaT/2)+1)*POD_modes.transpose()*B2*POD_modes + (DeltaT/2)*POD_modes.transpose()*Nk.transpose()*POD_modes+ 
        (DeltaT/2)*Eps*POD_modes.transpose()*B0*POD_modes +(DeltaT/2)*Eps_plus*unresolved_scales.transpose()*B0*unresolved_scales;
    RHS=LHS*VMS_POD_solution(all,0);
    //Loop with linear systems to finde the solution:
    for(int k=1; k<=N; k++){
        RHS=POD_modes.transpose()*B2*((DeltaT/2)*(F(all,k-1)+F(all,k))+2*POD_modes*VMS_POD_solution(all,k-1))-AVRHS-RHS;
        VMS_POD_solution(all,k)=  LHS.fullPivHouseholderQr().solve(RHS); 
     
    }
    //Transform back to finite element space:
    Rev_transformation=POD_modes*VMS_POD_solution;                   
    //Averrage L^2 error:
    for(int k=0; k<=N; k++){
        
        error=error+sqrt((snapshots(all,k)-Rev_transformation(all,k)).transpose()*B2*(snapshots(all,k)-Rev_transformation(all,k)));
    }
    error=error/((double)N+1);
    std::cout<< "Average error: " << error << std::endl;
    //Adding the meanvalue:
    for(int k=0; k<N+1; k++){
        Rev_transformation(all,k)=Rev_transformation(all,k)+AV;
    }

    //the hadder save-csv.h saves the results in a csv file with the function: void save(MatrixXd matrix,const std::string& name)
    const std::string dataName_results = "VMS-POD-P2-results.csv";
    save(Rev_transformation, dataName_results);

    return EXIT_SUCCESS;
}