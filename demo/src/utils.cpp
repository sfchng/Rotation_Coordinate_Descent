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

// eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

// SuiteSparse
#include "cholmod.h"
#include "SuiteSparseQR.hpp"
#include "umfpack.h"



namespace RCD_ROTAVG
{
    void enforce_directions(Mat &QQ, I_t &I)
    {
        const Long m = I.size();
        
        // impose i<j
        for (int k=0; k<m; k++)
        {
            const auto &e = I[k];
            const int &i = e.first;
            const int &j = e.second;
            
            if (i>j) // reverse direction
            {
                // swap indices
                auto aux = I[k].first;
                I[k].first = e.second;
                I[k].second = aux;
                QQ.row(k).head(3) *= -1; //QQ.row(k).head(3);
            }
        }
    }


   SpMat make_R(const int n, const Mat &QQ, const I_t &I)
    {
        // I encoded relative rotations as a sparse symmetric matrix R
        assert(n >= 0);
        
        const Long m = I.size();
        SpMat R(3*n, 3*n);
        
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        tripletList.reserve(2*9*m);
        
        int block_i, block_j;
        for (int k=0; k<m; k++)
        {
            const auto &e = I[k];
            
            block_i = 3 * e.first;
            block_j = 3 * e.second;
            
            const auto &q_vec = QQ.row(k);
            Quat q(q_vec(3), q_vec(0), q_vec(1), q_vec(2));
            Mat3 R_ij = q.normalized().toRotationMatrix().transpose();
            
            for (int ii=0; ii<3; ii++)
            {
                for (int jj=0; jj<3; jj++)
                {
                    tripletList.push_back(T(block_i+ii, block_j+jj, R_ij(ii,jj)));
                    tripletList.push_back(T(block_j+jj, block_i+ii, R_ij(ii,jj)));
                }
            }
        }
        
        R.setFromTriplets(tripletList.begin(), tripletList.end());
        
        
        return R;
    }

    Vec4 quat_mult(const Vec4 &q1, const Vec4 &q2)
    {
        Quat a1(q1(3), q1(0), q1(1), q1(2));
        Quat a2(q2(3), q2(0), q2(1), q2(2));
        a1 *= a2;
        return Vec4(a1.x(), a1.y(), a1.z(),a1.w());
    }    
    
    void init_mst(Mat &Q, const Mat &QQ, const I_t &I, const int f)
    {
        assert(f>0); // at least one rotation must be fixed.
        
        const Long &m = QQ.rows();
        const int n = (int) Q.rows();
        
        std::vector<bool> flags(n,false);
        
        flags[0] = true;
        int count = 1;
        
        while (count < n)
        {
            bool span_flag = false;
            
            for (Long j=0; j<m; j++)
            {
                const auto &e = I[j];
                const int &e1 = e.first;
                const int &e2 = e.second;
                
                if ( flags[e1] && !flags[e2] )
                {
                    if (e2 >= f) // do not change known rotations
                    {
                        Q.row(e2) = quat_mult( QQ.row(j), Q.row(e1) );
                    }
                    
                    if (!flags[e2])
                    {
                        count++;
                        flags[e2] = true;
                    }
                    span_flag = true;
                }
                
                if (!flags[e1] && flags[e2])
                {
                    if (e1 >= f)
                    {
                        Vec4 QQj_inv = QQ.row(j);
                        QQj_inv(3) *= -1;
                        Q.row(e1) = quat_mult( QQj_inv, Q.row(e2) );
                    }
                    
                    if (!flags[e1])
                    {
                        count++;
                        flags[e1] = true;
                    }
                    span_flag = true;
                }
            }
            
            if( !span_flag && count < n )
            {
                std::cerr << "Relative rotations DO NOT SPAN all the nodes in the VIEW GRAPH\n";
                std::cerr << "Number of nodes in Spanning Tree = " << count << "\n";
                std::cerr << "Connected Nodes are given as output\n";
                std::cerr << "Remove extra nodes and retry." << std::endl;
                std::exit(-1);
            }
        }
    }
    
    
    void quat_normalised(Mat &Q, const int f)
    {
        const Long &n = Q.rows();
        for (Long i=f; i<n; i++)
        {
            Quat q(Q(i,3), Q(i,0), Q(i,1), Q(i,2));
            q = q.normalized();
            Q.row(i) << q.x(), q.y(), q.z(), q.w();
        }
    }

   
} // end namespace
