/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 ROTATION COORDINATE DESCENT
%
%
%  This package contains the source code which implements the
%  Rotation Coordinate Descent (RCD and RCDL) in
%
%                 Rotation Coordinate Descent for 
%             Fast Globally Optimal Rotation Averaging
%            
%
%  The source code and demo are suplied for academic use only.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "rcd.hpp"
#include "utils.hpp"

#include <iostream>
#include <ctime> // time_t
#include <ctype.h> //tolower
#include <stdio.h>
#include <algorithm>
#include <numeric>  // iota
#include <vector>
#include <random>
#include <fstream>
#include <chrono>

// eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

// SuiteSparse
#include "cholmod.h"
#include "SuiteSparseQR.hpp"
#include "umfpack.h"

/*  RCD for rotation averaging : estimate n rotations given m relative rotations.
    Return : 
             Q           : Absolute camera orientations in (Quaternion)
             final_objval: Final objective value
             num_epoch   : Number of epoch 

    Input :
            I            : Camera connections
            R            : Relative camera orientations
            Q            : Initialised absolute camera orientations
            f            : First index
            max_epoch    : Maximum number of epoch
*/

namespace RCD_ROTAVG
{
    void rcd(const I_t &I, const SpMat &R, Mat &Q, const int f, const int max_epoch, 
            int &num_epoch, double &total_runtime, double &final_objval)
    {
  
        const Long m = Q.rows();     // number of absolute rotations
        double obj;
        SpMat W(3*m, 3);
        Mat X = Mat::Zero(3*m,3);   // Store rotation matrices as X with a size of (3m*3)

        for (int i=0; i<m; i++)
        {
            const auto &q_vec = Q.row(i);
            Quat q(q_vec(3), q_vec(0), q_vec(1), q_vec(2));
            X.block(3*i,0,3,3) = q.normalized().toRotationMatrix(); 
        }     
        double oldObj = -0.5*(X.transpose() * R * X).trace();
        std::cout << "Status: Running RCD" << std::endl;
        for (int epoch = 0; epoch < max_epoch; epoch++)        
        {
            auto tic_global = std::chrono::high_resolution_clock::now();          
            for (int k=0; k<m; k++)
            {
                W = R.middleCols(3*k, 3);
                W.makeCompressed();
                Mat BW = X*(X.transpose()*W);
                Mat A = W.transpose()*BW;
                Eigen::JacobiSVD<Eigen::Matrix3d> svd(A, Eigen::ComputeThinU);
                auto S = svd.singularValues();
                auto U = svd.matrixU();
                Vec3 s =  S.array().sqrt().inverse();
                Mat aux = BW * (U*s.asDiagonal()) * U.transpose();
                X.topRows(3*k) = aux.topRows(3*k);
                X.block(3*k, 0, 3, 3) = Mat3::Identity();
                X.bottomRows(3*(m-k-1)) = aux.bottomRows(3*(m-k-1));    
            }       
    
            auto elapsed_global = std::chrono::high_resolution_clock::now() - tic_global;   
            auto time_diag_pass_nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed_global).count();   
            obj = -0.5*(X.transpose() * R * X).trace();    
            
            total_runtime = total_runtime + time_diag_pass_nanoseconds;
            printf("Status: Epoch %d | Current objval %.12e Previous Objval %.12e Epoch time[s] %.2e Total time[s] %.2e\n", epoch+1, obj, oldObj, time_diag_pass_nanoseconds*1e-9, total_runtime*1e-9);

            if ( (oldObj-obj) / fmax(abs(oldObj),1) <= 1e-9)
            {
                final_objval = obj;
                num_epoch = epoch+1;
                for (int i=0; i<m; i++)
                {
                    Mat3 Ri = X.middleRows(3*i, 3);
                    if (Ri.determinant()<0)
                    { 
                        Ri *= -1;
                    }
                    Quat q(Ri);
                    Q.row(i) = Vec4(q.x(), q.y(), q.z(), q.w());
                }
                break;
            }
            oldObj = obj;
        }
    }    
   
} // end namespace
