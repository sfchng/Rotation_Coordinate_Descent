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

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>
#include "SuiteSparseQR.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <limits>

namespace RCD_ROTAVG
{
    #define EPS 2.2204e-16
    const double DBL_MAX_ = std::numeric_limits<double>::max();
    
    typedef SuiteSparse_long Long;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, Long> SpMat;
    typedef Eigen::MatrixXd Mat;
    typedef Eigen::Matrix3d Mat3;
    typedef Eigen::VectorXd Vec;
    typedef Eigen::Vector3d Vec3;
    typedef Eigen::Vector4d Vec4;
    typedef Eigen::Quaterniond Quat;
    typedef Eigen::Triplet<double> T;
    typedef std::vector< std::pair <int,int> > I_t;
    
    // --------------------------------------------------------------------
    //   Function interfaces
    // --------------------------------------------------------------------
    // n    --  number of unknown rotations.
    // f    --  number of fixed rotations.
    // Q    --  (n+f)x4 matrix with absolute rotations. First f rotations are fixed.
    // QQ   --  mX4 matrix with relative rotations.
    // I    --  Camera connections.
    void init_mst(Mat &Q, const Mat &QQ, const I_t &I, const int f);

    SpMat make_A(const int n, const int f, const I_t &I);
    
    // Relative rotations block matrix
    SpMat make_R(const int n, const Mat &QQ, const I_t &I);
    
    // ---- extra functions -------------------------------------------------
    
    // normalise non-fixed quaternions in Q
    void quat_normalised(Mat &Q, const int f);


    void enforce_directions(Mat &QQ, I_t &I);

    

}


#endif /* utils_hpp */
