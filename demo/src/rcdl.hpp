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


#ifndef rcdl_hpp
#define rcdl_hpp

#include <stdio.h>
#include "SuiteSparseQR.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <limits>


namespace RCDL_ROTAVG
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

    // ----------------------------------------------------------------------
    //     Cost
    // ----------------------------------------------------------------------
    enum Cost {L2, L1, L15, L05, Geman_McClure, Huber, Pseudo_Huber, Andrews,
        Bisquare, Cauchy, Fair, Logistic, Talwar, Welsch};
    
    inline
    std::ostream &operator<<( std::ostream &os, const Cost cost )
    {
        switch (cost) {
            case L2: os << "L2"; break;
            case L1: os << "L1"; break;
            case L15: os << "L1.5"; break;
            case L05: os << "L0.5"; break;
            case Geman_McClure: os << "Geman-McClure"; break;
            case Huber: os << "Huber"; break;
            case Pseudo_Huber: os << "Pseudo-Huber"; break;
            case Andrews: os << "Andrews"; break;
            case Bisquare: os << "Bisquare"; break;
            case Cauchy: os << "Cauchy"; break;
            case Fair:  os << "Fair"; break;
            case Logistic: os << "Logistic"; break;
            case Talwar: os << "Talwar"; break;
            case Welsch: os << "Welsch"; break;
        }
        return os;
    }
    
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

    void l1ra(const Mat &QQ, const I_t &I, const SpMat &A,
              Mat &Q, const int f, const int max_iters, double change_th,
              int &iter, double &runtime);
    
    void irls(const Mat &QQ, const I_t &I, const SpMat &A,
              Cost cost, double sigma, Mat &Q, const int f,
              const int max_iters, double change_th, Vec &weights,
              int &iteration, double &runtime, double &runtime_solve);

     void rcdl(const Mat &QQ, const I_t &I, const SpMat &R, Mat &Q, const int f, const int rcd_max_epoch, double change_th,
           int &rcd_epoch_out, double &rcd_runtime, 
           const SpMat &A, Cost cost, double sigma, const int local_max_iter, Vec &weights, int &local_iters_out, double &local_runtime, int &local_update_iters, double &total_local_runtime,
           double &final_objval, double &local_solve_runtime, int &local_call_iters);     
    
    
    // ---- extra functions -------------------------------------------------
    
    // normalise non-fixed quaternions in Q
    void quat_normalised(Mat &Q, const int f);


    void enforce_directions(Mat &QQ, I_t &I); 

}


#endif /* rcdl_hpp */
