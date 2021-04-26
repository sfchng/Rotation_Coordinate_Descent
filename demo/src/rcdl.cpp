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


#include "rcdl.hpp"

#include <iostream>
#include <ctime> // time_t
#include <ctype.h> //tolower
#include <stdio.h>
#include <algorithm>
#include <numeric>  // iota
#include <vector>
#include <random>
#include <chrono>

// eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

// SuiteSparse
#include "cholmod.h"
#include "SuiteSparseQR.hpp"
#include "umfpack.h"


// debugging
#include <fstream>

namespace RCDL_ROTAVG
{
#define VERBOSE false
    
    cholmod_sparse viewAsCholmod(Eigen::Ref<SpMat> mat)
    {
        cholmod_sparse res;
        res.nzmax   = mat.nonZeros();
        res.nrow    = mat.rows();
        res.ncol    = mat.cols();
        res.p       = mat.outerIndexPtr();
        res.i       = mat.innerIndexPtr();
        res.x       = mat.valuePtr();
        res.z       = 0;
        res.sorted  = 1;
        if(mat.isCompressed())
        {
            res.packed  = 1;
            res.nz = 0;
        }
        else
        {
            res.packed  = 0;
            res.nz = mat.innerNonZeroPtr();
        }
        
        res.itype = CHOLMOD_LONG; //CHOLMOD_INTLONG;
        res.xtype = CHOLMOD_REAL;
        res.dtype = CHOLMOD_DOUBLE;
        res.stype = 0;
        
        return res;
    }
    
    
    cholmod_dense viewAsCholmod(Mat& mat)
    {
        cholmod_dense res;
        
        res.nrow   = mat.rows();
        res.ncol   = mat.cols();
        res.nzmax  = res.nrow * res.ncol;
        res.d      = mat.outerStride();
        res.x      = (void*)(mat.derived().data());
        res.z      = 0;
        
        res.xtype = CHOLMOD_REAL;
        res.dtype = CHOLMOD_DOUBLE;
        
        return res;
    }
    
    
    Vec4 quat_mult(const Vec4 &q1, const Vec4 &q2)
    {
        Quat a1(q1(3), q1(0), q1(1), q1(2));
        Quat a2(q2(3), q2(0), q2(1), q2(2));
        a1 *= a2;
        return Vec4(a1.x(), a1.y(), a1.z(),a1.w());
    }
    
    
    // Compute  inv(Qj) Qij Qi for all connections in I
    Mat delta_rel(const I_t &I, const Mat &QQ, const Mat &Q) //i
    {
        const Long &m = I.size();
        
        Mat resp(m, 4);
        Mat Q_inv = Q;
        Q_inv.col(3).array() *= -1;
        
        for (Long k=0; k<m; k++)
        {
            const auto &e = I[k];
            const Long &i = e.first;
            const Long &j = e.second;
            
            resp.row(k) = quat_mult( Q_inv.row(j), quat_mult( QQ.row(k), Q.row(i) ) );
        }
        
        return resp;
    }
    
    
    //solve sparse linear system  (Ax = b) by using umfpack
    Vec linsolve(Long *Ap, Long *Ai,  double *Ax, Vec &b )
    {
        const Long n = b.size();
        Vec x(n);
        
        Long status;
        double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
        
        void *Symbolic, *Numeric ;
        //    cholmod_common Common, *cc ;
        //    cc = &Common;
        //    cc->print = 0;
        //    cholmod_l_start (cc) ;
        //    cc->print = 0;
        
        //status = umfpack_dl_symbolic (n, n, Ap, Ai, Ax, &Symbolic, Control, Info);
        status = umfpack_dl_symbolic (n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);
        
        if (status < 0)
        {
            umfpack_dl_report_info (Control, Info) ;
            umfpack_dl_report_status (Control, status) ;
            std::cerr<< "umfpack_dl_symbolic failed" <<std::endl;
            std::exit(-1);
        }
        
        // produces output in mode release...
        //status = umfpack_dl_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
        status = umfpack_dl_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
        if (status < 0)
        {
            //umfpack_dl_report_info (Control, Info) ;
            //umfpack_dl_report_status (Control, status) ;
            std::cerr<< "umfpack_dl_numeric failed" <<std::endl;
            std::exit(-1);
        }
        
        umfpack_dl_free_symbolic (&Symbolic) ;
        status = umfpack_dl_solve (UMFPACK_A, Ap, Ai, Ax, x.data(), b.data(), Numeric, Control, Info) ;
        
        if (status < 0)
        {
            umfpack_dl_report_info (Control, Info) ;
            umfpack_dl_report_status (Control, status) ;
            std::cerr<< "umfpack_dl_solve failed" << std::endl ;
            std::exit(-1);
        }
        
        umfpack_dl_free_numeric (&Numeric) ;
        
        //    cholmod_l_finish(cc);
        
        return x;
    }
    
    
    void sp_vec_to_squared_mat(const Long n, const Long nz, const Long *Ai, Long *Ap_new, Long *Ai_new)
    {
        long prev_j = 0;
        long nz_j = 0; // count nz in column j;
        
        Ap_new[0] = 0;
        for (long d=0; d<nz; d++)
        {
            long k = Ai[d];
            long j = k / n;
            long i = k % n;
            
            //std::cout<< "(" << i << ", " << j <<std::endl;
            if (prev_j != j) // we moved to a new column
            {
                // upgrade Ai from old_j until j-1
                for (long t = prev_j; t<=j-1; t++)
                {
                    Ap_new[t+1] = Ap_new[prev_j] + nz_j;
                }
                
                // start counting a new column
                nz_j = 1;
            }
            else
            {
                nz_j++;
            }
            
            Ai_new[d] = i;
            
            prev_j = j;
        }
        
        for (long t = prev_j; t<n; t++)
        {
            Ap_new[t+1] = Ap_new[prev_j] + nz_j;
        }
    }
    
    
    Vec l1decode_pd(const Vec &x0, const SpMat &A, const Vec &y,
                    const int pdmaxiter, const SpMat &AtA)
    {
        const double PDTOL = 1e-3;
        
        const long n = x0.size();
        const long m = y.size();
        
        double alpha = 0.01;
        double beta = 0.5;
        double mu = 10;
        
        //gradf0 = [zeros(N,1); ones(M,1)];
        Vec gradf0(n+m);
        gradf0.head(n).setZero();
        gradf0.tail(m).setOnes();
        
        Vec x = x0;
        Mat Ax = A*x;
        
        //u = (0.95)*abs(y-Ax) + (0.10)*max(abs(y-Ax));
        Vec u, y_Ax_abs;
        y_Ax_abs = y-Ax;
        y_Ax_abs = y_Ax_abs.array().abs();
        u  = y_Ax_abs * 0.95;
        u.array() += y_Ax_abs.maxCoeff() * 0.10;
        
        Vec fu1 = Ax - y - u;
        Vec fu2 = -Ax + y - u;
        
        Vec lamu1 = -fu1.array().inverse();
        Vec lamu2 = -fu2.array().inverse();;
        
        Mat Atv = A.transpose()*(lamu1-lamu2);
        
        double sdg = -(fu1.dot(lamu1) + fu2.dot(lamu2));
        double tau = mu*2*m / sdg;
        
        Vec rcent(2*m); //[-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
        rcent.head(m) = -lamu1.array()*fu1.array();
        rcent.tail(m) = -lamu2.array()*fu2.array();
        rcent.array() -= (1.0/tau);
        
        Vec rdual; //  = gradf0 + [Atv; -lamu1-lamu2];
        rdual = gradf0;
        rdual.head(n).array() += Atv.array();
        rdual.tail(m).array() -= lamu1.array();
        rdual.tail(m).array() -= lamu2.array();
        
        Vec rdual_rcent(n+m + 2*m); //= norm([rdual; rcent]);
        rdual_rcent.head(n+m) = rdual;
        rdual_rcent.tail(2*m) = rcent;
        double resnorm = rdual_rcent.norm();
        
        int pditer = 0;
        bool done = (sdg < PDTOL) || (pditer >= pdmaxiter);
        
        Vec xp, up;
        
        while (!done)
        {
            pditer++;
            
            //w2 = -1 - 1/tau*(1./fu1 + 1./fu2);
            Vec w2 = -1 -1.0/tau*(fu1.array().inverse() + fu2.array().inverse());
            
            //   sig1 = -lamu1./fu1 - lamu2./fu2;
            Vec sig1 = -lamu1.array()/fu1.array() - lamu2.array()/fu2.array();
            Vec sig2 =  lamu1.array()/fu1.array() - lamu2.array()/fu2.array();
            Vec sigx = sig1.array() - sig2.array().square().array()/sig1.array();
            
            //w1 = -1/tau*(A'*(-1./fu1 + 1./fu2));
            Vec w1 = -fu1.array().inverse() + fu2.array().inverse();
            w1 = -1.0/tau * ( A.transpose() * w1 );
            
            //w1p = w1 - A'*((sig2./sig1).*w2);
            Vec w1p = ( sig2.array() / sig1.array() ).array() * w2.array();
            w1p = w1 - A.transpose() * w1p;
            
            SpMat sigx_sp = sigx.sparseView();
            SpMat H11p_vec = AtA * sigx_sp;
            
            long nz = H11p_vec.nonZeros();
            double *Hx = H11p_vec.valuePtr();
            Long *Hi = H11p_vec.innerIndexPtr();
            Long *Hi_new = new Long[nz];
            Long *Hp_new = new Long[n+1];
            
            sp_vec_to_squared_mat(n, nz, Hi, Hp_new, Hi_new);
            
            Vec dx = linsolve(Hp_new, Hi_new, Hx, w1p);
            
            delete[] Hp_new;
            delete[] Hi_new;
            
            Vec Adx = A*dx;
            
            // du = (w2 - sig2.*Adx)./sig1;
            Vec du = (w2.array() - sig2.array() * Adx.array()).array()/sig1.array() ;
            
            //dlamu1 = -(lamu1./fu1).*(Adx-du) - lamu1 - (1/tau)*1./fu1;
            Vec dlamu1 =  -lamu1.array() / fu1.array();  //e-6
            dlamu1.array() *=  (Adx-du).array() ; // here the error goes big
            dlamu1.array() -= lamu1.array();
            dlamu1.array() -= (1.0/tau)*fu1.array().inverse().array();
            
            // dlamu2 = (lamu2./fu2).*(Adx + du) -lamu2 - (1/tau)*1./fu2;
            Vec dlamu2 =  lamu2.array() / fu2.array();
            dlamu2.array() *=  (Adx+du).array() ;
            dlamu2.array() -= lamu2.array();
            dlamu2.array() -= (1.0/tau)*fu2.array().inverse().array();
            
            // Atdv = A'*(dlamu1-dlamu2);
            Vec Atdv = A.transpose() * (dlamu1-dlamu2);
            
            //   % make sure that the step is feasible: keeps lamu1,lamu2 > 0, fu1,fu2 < 0
            //   indl = find(dlamu1 < 0);  indu = find(dlamu2 < 0);
            //   s = min([1; -lamu1(indl)./dlamu1(indl); -lamu2(indu)./dlamu2(indu)]);
            double s = 1;
            for (int i=0; i<m; i++)
            {
                double dlamu1_i = dlamu1(i);
                if (dlamu1_i<0)
                {
                    s = fmin(s, -lamu1(i)/dlamu1_i);
                }
                
                double dlamu2_i = dlamu2(i);
                if (dlamu2_i<0)
                {
                    s = fmin(s, -lamu2(i)/dlamu2_i);
                }
            }
            // indl = find((Adx-du) > 0);
            //  indu = find((-Adx-du) > 0);
            // s = (0.99)*min([s; -fu1(indl)./(Adx(indl)-du(indl)); -fu2(indu)./(-Adx(indu)-du(indu))]);
            for (int i=0; i<m; i++)
            {
                //indl = find((Adx-du) > 0);  indu = find((-Adx-du) > 0);
                
                double Adx_du_i = Adx(i)-du(i);
                if (Adx_du_i>0)
                {
                    s = fmin(s, -fu1(i)/Adx_du_i );
                }
                
                double neg_Adx_du_i = -Adx(i)-du(i);
                if (neg_Adx_du_i>0)
                {
                    s = fmin(s, -fu2(i)/neg_Adx_du_i );
                }
            }
            s *= 0.99;
            
            // backtrack
            bool suffdec = false;
            int backiter = 0;
            
            //VecX xp, up;
            Vec Axp, Atvp;
            Vec lamu1p, lamu2p;
            Vec fu1p, fu2p;
            Vec rdp;
            while(!suffdec)
            {
                xp = x + s*dx;
                up = u + s*du;
                
                Axp = Ax + s*Adx;
                Atvp = Atv + s*Atdv;
                // here Atdv is inconsistent so Atvp is too
                
                lamu1p = lamu1 + s*dlamu1;
                lamu2p = lamu2 + s*dlamu2;
                
                fu1p = Axp - y - up;
                fu2p = -Axp + y - up;
                
                // rdp = gradf0 + [Atvp; -lamu1p-lamu2p];
                rdp = gradf0;
                rdp.head(n).array() += Atvp.array();
                rdp.tail(m).array() += (-lamu1p-lamu2p).array();
                
                // rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
                Vec rcp(2*m);
                rcp.head(m) = -lamu1p.array() * fu1p.array();
                rcp.tail(m) = -lamu2p.array() * fu2p.array();
                rcp.array() -= 1.0/tau;
                
                // suffdec = (norm([rdp; rcp]) <= (1-alpha*s)*resnorm);
                suffdec = sqrt(rdp.squaredNorm() + rcp.squaredNorm()) <= (1-alpha*s)*resnorm ;
                
                s *= beta;
                backiter++;
                if (backiter > 32)
                {
                    std::cout<<"Stuck backtracking, returning last iterate."<<std::endl;
                    xp = x;
                    return xp;
                }
            }
            
            // next iteration
            x = xp;
            u = up;
            
            Ax = Axp;
            Atv = Atvp;
            
            lamu1 = lamu1p;
            lamu2 = lamu2p;
            
            fu1 = fu1p;
            fu2 = fu2p;
            
            // surrogate duality gap
            //   sdg = -(fu1'*lamu1 + fu2'*lamu2);
            sdg = -(fu1.dot(lamu1) + fu2.dot(lamu2));
            
            tau = mu*2*m/sdg;
            
            // rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
            rcent.head(m) = -lamu1.array()*fu1.array();
            rcent.tail(m) = -lamu2.array()*fu2.array();
            rcent.array() -= (1.0/tau);
            
            rdual = rdp;
            
            //resnorm = norm([rdual; rcent]);
            resnorm = sqrt( rdual.squaredNorm() + rcent.squaredNorm() );
            
            done = (sdg < PDTOL) || (pditer >= pdmaxiter);
            
            //printf("Iteration = %d, tau = %8.3e, Primal = %8.3e, PDGap = %8.3e, Dual res = %8.3e\n",
            //pditer, tau, u.sum(), sdg, rdual.norm());
            //printf("                  H11p condition number = %8.3e", hcond);
        }
        
        return xp;
    }
    
    
    void exp_map(Mat &W)
    {
        // theta = sqrt(sum(W(:,2:4).*W(:,2:4),2));
        W.col(3) = W.leftCols(3).rowwise().norm();
        const Vec &theta = W.col(3);
        
        Vec ang_coef = (theta / 2.0).array().sin().array() / theta.array();
        
        // W(:,1) = cos(theta/2);
        //W.col(0) = (theta / 2.0).array().cos();
        W.col(3) /= 2.0;
        W.col(3) = W.col(3).array().cos();
        
        // W(:,2:4) = W(:,2:4).*repmat(sin(theta/2)./theta,[1,3]);
        //VecX ang_coef = (theta / 2.0).array().sin().array() / theta.array();
        W.col(0).array() *= ang_coef.array();
        W.col(1).array() *= ang_coef.array();
        W.col(2).array() *= ang_coef.array();
        
        // W(isnan(W)) = 0;
        W = W.unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; });
    }
    
    
    // transform rows of w to [ r1*theta, r2*theta, r3*theta, theta]
    // where theta is the angle of rotation
    //       r = [r1 r2 r3] is the axis of rotation
    void log_map(Mat &w)
    {
        Mat s2 = w.leftCols(3).rowwise().norm();
        
        // w(:,1)=2*atan2(s2,w(:,1));
        const Long &m = w.rows();
        for(int i=0; i<m; ++i)
        {
            double &theta = w.col(3)(i);
            
            theta = 2*atan2(s2(i), theta);
            
            if (theta < -EIGEN_PI)
            {
                theta += 2*EIGEN_PI;
            }
            else if (theta >= EIGEN_PI)
            {
                theta -= 2*EIGEN_PI;
            }
        }
        
        // B=w(:,2:4).*repmat(w(:,1)./s2,[1,3]);
        Mat aux = w.col(3).array() / s2.array();
        w.col(0).array() *= aux.array();
        w.col(1).array() *= aux.array();
        w.col(2).array() *= aux.array();
        
        // B(isnan(B))=0;% This tackles the devide by zero problem.
        for( int i=0; i<m; ++i)
        {
            if (s2(i)< EPS )
                w.row(i).head(3).setZero();
        }
    }
    
    
    //void ls_solve( Mat &X, SpMat &DA, Mat &DB)
    template<typename Derived>
    void ls_solve( Eigen::MatrixBase<Derived> &X, SpMat &DA, Mat &DB)
    {
        // least squares by using suite sparse
        const Long &n = X.rows();
        cholmod_common Common, *cc ;
        cholmod_sparse DA_chol;
        cholmod_dense *X_chol, DB_chol;
        
        cc = &Common ;
        cholmod_l_start (cc) ;
        
        DA_chol = viewAsCholmod(DA);
        DB_chol = viewAsCholmod(DB);
        X_chol = SuiteSparseQR<double>(&DA_chol, &DB_chol, cc) ;
        
        X = Eigen::Map<Mat> (reinterpret_cast<double*> (X_chol->x), n, 3);
        
        cholmod_l_free_dense (&X_chol, cc) ;
        cholmod_l_finish(cc);
    }
    
    
    void irls(const Mat &QQ, const I_t &I, const SpMat &A,
              Cost cost, double sigma, Mat &Q, const int f,
              const int max_iters, double change_th, Vec &weights,
              int &iters, double &runtime, double &runtime_solve)
    {

        clock_t tic, toc;
        
        Long m = QQ.rows();
        Long n = Q.rows() - f;  // number of variables
        
        Mat w = Mat::Zero(m, 4);  // delta rel rotations
        Mat W = Mat::Zero(n, 4);  // delta abs rotations
        Eigen::Map <Mat> W3 (W.data(), n, 3);
        
        double score = DBL_MAX_; //inf
        iters = 0;
        
        weights.setOnes();
        
        
        if (VERBOSE)
         printf("IRLS iteration: %4d; Time: %7.3f; Change: %0.12f\n",
               0, (double)(toc-tic)/CLOCKS_PER_SEC, 0.0);
                
        
        Mat DB(m,3);
        SpMat D(m,m);
        std::vector<T> tripletList;
        tripletList.reserve(m);

        //std::cout <<"Change th" << change_th << std::endl;
        while( ( score > change_th ) && (iters < max_iters) )
        {
            auto tic_local = std::chrono::high_resolution_clock::now();     
            w = delta_rel(I, QQ, Q);
            log_map(w);
            
            // make D
            tripletList.clear();
            for(int i=0; i<m; i++)
            {
                tripletList.push_back( T( i, i, weights(i) ) );
            }
            D.setFromTriplets(tripletList.begin(), tripletList.end());
            

            SpMat DA = D*A;
            DB.col(0) = DB.col(1) = DB.col(2) = weights.array();
            DB.col(0).array() *= w.col(0).array();
            DB.col(1).array() *= w.col(1).array();
            DB.col(2).array() *= w.col(2).array();
      
            auto tic_ls_solve = std::chrono::high_resolution_clock::now();
            ls_solve(W3, DA, DB);
            auto elapsed_ls_solve = std::chrono::high_resolution_clock::now()-tic_ls_solve;
            auto t_nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed_ls_solve).count();
            runtime_solve += t_nanoseconds;

     
            Mat E = A*W3 - w.leftCols(3);
            
            score = W3.rowwise().norm().mean();

            printf("IRLS Solve: %4d; Time[s]: %.2e; Change: %.20f\n",
                   iters+1, t_nanoseconds*1e-9, score);        

            
            exp_map(W);
            
            // Upgrade abs rotations
            for (Long i=0; i<n; i++)
            {
                Q.row(i+f) = quat_mult(Q.row(i+f), W.row(i));

            }
 
            iters++;
            auto elapsed_toc_local = std::chrono::high_resolution_clock::now()-tic_local;
            auto t_nanoseconds_local = std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed_toc_local).count();   
            runtime += t_nanoseconds_local;

            if (VERBOSE)
            {
            toc = clock();
            printf("IRLS iteration: %4d; Time: %0.12f; Change: %.20f\n",
                   iters, (double)(toc-tic)/CLOCKS_PER_SEC, score);   
            }  
        }
        
    }
    
    
    SpMat make_A(const int n, const int f, const I_t &I)
    {
        assert(n >= 0 && f >= 0);
        assert(n-f > 1);
        
        const Long m = I.size();
        SpMat A(m, n-f);
        
        int i, j;
        for (int k=0; k<m; k++)
        {
            const auto &e = I[k];
            const int &e1 = e.first;
            const int &e2 = e.second;
            
            j = e2-f;
            if (j < 0) continue;
            A.coeffRef(k, j) =  1;
            
            i = e1-f;
            if (i < 0) continue;
            A.coeffRef(k, i) = -1;
        }

        /*std::cout << "PRINT A" << std::endl;
        std::cout << A << std::endl;
        */
        A.makeCompressed();
        return A;
    }
    


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



    
    SpMat make_C(const int f, const I_t &I)
    {
        const Long m = I.size();
        SpMat C(m, f);
        
        for (int k=0; k<m; k++)
        {
            const auto &e = I[k];
            const int &i = e.first;
            const int &j = e.second;
            
            if (i<f)
            {
                C.coeffRef(k, i) = 1;
            }
            
            if (j<f)
            {
                C.coeffRef(k, j) = -1;
            }
        }
        C.makeCompressed();
        return C;
    }
    
    
    
    //f: fixed vars
    SpMat make_AtA (const Long n, const Long f, const I_t &I)
    {
        const Long m = I.size();
        SpMat AtA(n*n,m);
        
        for (int k=0; k<m; k++)
        {
            const auto &e = I[k];
            const Long &e1 = e.first;
            const Long &e2 = e.second;
            
            const Long i = e1-f;
            const Long j = e2-f;
            
            if (i>=0)
            {
                const Long r1 = n*i + i;
                AtA.coeffRef(r1,k) = 1;
            }
            
            if (j>=0)
            {
                const Long r2 = n*j + j;
                AtA.coeffRef(r2,k) = 1;
            }
            
            if ( i>=0 && j>=0 )
            {
                const Long r3 = n*i + j;
                const Long r4 = n*j + i;
                AtA.coeffRef(r3,k) = -1;
                AtA.coeffRef(r4,k) = -1;
            }
        }
        AtA.makeCompressed();
        
        return AtA;
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

    void rcdl(const Mat &QQ, const I_t &I, const SpMat &R, Mat &Q, const int f, const int rcd_max_epoch, double change_th,
            int &rcd_epoch_out, double &total_runtime, const SpMat &A, 
            Cost cost, double sigma, const int irls_max_iter, Vec &weights, int &irls_total_iter, double &irls_total_runtime, int &local_update_iters, double &total_local_runtime, 
            double &final_objval, double &irls_total_runtime_solve, int &local_call_iters)
    {
        int irls_iter=0;

        /* Precision */
        int precision = std::numeric_limits<double>::max_digits10;

        int local_update_index = 0;
        int s = 0;
        
        const Long m = Q.rows();     // number of unknown rotations

        /* Initialise matrix */
        Mat Q_bar = Mat::Zero(m,4);   

        // Store rotation matrices as X with a size of (3m*3)
        Mat X = Mat::Zero(3*m,3);

        for (int i=0; i<m; i++)
        {
            const auto &q_vec = Q.row(i);
            Quat q(q_vec(3), q_vec(0), q_vec(1), q_vec(2));
            X.block(3*i,0,3,3) = q.normalized().toRotationMatrix(); 
        }
 
        // Evaluate solution prior to optimisation
        double oldObj = -0.5*(X.transpose() * R * X).trace();
        double obj;
 
        SpMat W(3*m, 3);

        for (int epoch = 0; epoch < rcd_max_epoch; epoch++)
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

            if ( (oldObj-obj) / fmax(abs(oldObj),1) <= 1e-9)
            {
                printf("Status: Epoch %d | Current objval %.12e Previous Objval %.12e Epoch time[s] %.2e Total time[s] %.2e\n", epoch+1, obj, oldObj, time_diag_pass_nanoseconds*1e-9, total_runtime*1e-9); 
            
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
                
                final_objval = obj;
                rcd_epoch_out = epoch +1;
                break;
            } 

            /* Initialise storage array of runtime */
            double irls_runtime = 0; double irls_runtime_solve = 0;        

            local_update_index+=1;
            if ( s ==0 || epoch % s == 0)
            {
                std::cout<< "Status: Calling local solver" << std::endl;
                auto tic_local = std::chrono::high_resolution_clock::now();      
                for (int i=0; i<m; i++)
                {
                    Mat3 Ri = X.middleRows(3*i, 3);
                    if (Ri.determinant()<0)
                    {
                        Ri = -Ri;
                    }       
                    Quat q(Ri);
                    Q_bar.row(i) = Vec4(q.x(), q.y(), q.z(), q.w());   
                }

                // irls
                irls(QQ, I, A, cost, sigma, Q_bar, f, irls_max_iter, change_th, weights, irls_iter, irls_runtime, irls_runtime_solve);
                
                auto elapsed_local = std::chrono::high_resolution_clock::now() - tic_local;   
                auto time_local = std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed_local).count();          

                // Store rotation matrices as X with a size of (3m*3)  
                Mat X_bar = Mat::Zero(3*m,3);
                for (int i=0; i<m; i++)
                {
                    const auto &q_vec = Q_bar.row(i);
                    Quat q(q_vec(3), q_vec(0), q_vec(1), q_vec(2));
                    X_bar.block(3*i,0,3,3) = q.normalized().toRotationMatrix();   
                } 

                // Evaluate new objective value
                double obj_local = -0.5*(X_bar.transpose() * R * X_bar).trace();

                if (obj_local < obj)
                {
                    obj = obj_local;
                    Q = Q_bar;
                    X = X_bar;
                    local_update_iters = local_update_iters + 1;
                }
                else 
                {   s +=2;
                }
                total_local_runtime = total_local_runtime + time_local;
                irls_total_runtime = irls_total_runtime + irls_runtime;
                irls_total_iter = irls_total_iter + irls_iter;
                irls_total_runtime_solve = irls_total_runtime_solve + irls_runtime_solve;
            }  

            printf("Status: Epoch %d | Current objval %.12e Previous Objval %.12e Epoch time[s] %.2e Total time[s] %.2e\n", epoch+1, obj, oldObj, time_diag_pass_nanoseconds*1e-9, total_runtime*1e-9);
            if ( (oldObj-obj) / fmax(abs(oldObj),1) <= 1e-9)
            {
                final_objval = obj;
                rcd_epoch_out = epoch+1;
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
