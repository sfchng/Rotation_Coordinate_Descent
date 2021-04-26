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

#ifndef rcd_hpp
#define rcd_hpp

#include "utils.hpp"
#include <stdio.h>
#include "SuiteSparseQR.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <limits>


namespace RCD_ROTAVG
{
    #define EPS 2.2204e-16


    void rcd(const I_t &I, const SpMat &R, Mat &Q, const int f, const int rcd_max_epoch, int &rcd_num_epoch, double &rcd_runtime, double &final_objval);    
}


#endif /* rcd_hpp */
