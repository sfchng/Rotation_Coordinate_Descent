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
#include <iostream>
#include "rcd.hpp"
#include "utils.hpp"
#include <fstream>
#include <set>

using namespace RCD_ROTAVG;


int main(int argc, const char * argv[])
{
    std::string demo_notice =
    "RCD for rotation averaging\n"
    "\n";
    
    std::string usage =
    "Usage:\n"
    "\n"
    "test input_file [output_file [rcd_epoch]\n"
    "\n"
    "  input_file      --  input file with format:\n"
    "      m n f\n"
    "      [relative rotations block]\n"
    "      [absolute rotations block]\n"
    "\n"
    "      where:\n"
    "          m  --  Number of relative rotations.\n"
    "          n  --  Number of absolute rotations.\n"
    "          f  --  Number of fixed absolute rotations. The first f rotations are fixed.\n"
    "                 If f=0, the first rotations is set to I; otherwise, the first f rotations\n"
    "                 are fixed with the first f rotations in the block of absolute rotations.\n"
    "\n"
    "      relative rotation block  --  Each line contains a relative rotation [w x y z] between\n"
    "                                   views i and j (with i,j integers and i<j) in the format:\n"
    "          i j w x y z\n"
    "\n"
    "      absolute rotation block  --  Each line contains an absolute rotation [w x y z] in the\n"
    "                                   order of the views. At least f rotations must be given. \n"
    "                                   Rotations from line (f+1) are used for initialisation if given.\n"
    "                                   Each line has the format:\n"
    "          w x y z\n"
    "\n"
    "  output_file     --  Output file with an absolute rotations [w x y z] and IRLS weights.\n"
    "                      The first n lines are the absolute rotations in the same order\n"
    "                      of the views). The m next lines are the weights. Default output file \n"
    "                      is l1_irls_out.txt. Rotation format:\n"
    "      w x y z\n"
    "\n"
    "  rcd_epoch       --  RCD maximum number of epoch. Default value is 10000\n"
    "\n"
    "  mst initialisation --  Minimum Spanning Tree initialisation. Default value is 0\n"
    "\n";    
    
    
    
    std::cout << demo_notice << std::endl;
    
    if (argc-1<1)
    {
        std::cout << usage <<std::endl;
        std::exit(-1);
    }
    
    //---------------------------------------------------------------------
    // Check input and output size
    //---------------------------------------------------------------------
    if (argc-1<1 || argc-1>4)
    {
        std::cerr << "Invalid number of arguments. Expected at least 1";
        std::cerr << "and at most 3 arguments." << std::endl;
        std::cerr << usage << std::endl;
        std::exit(-1);
    }
    
    //---------------------------------------------------------------------
    //    Check and read input
    //---------------------------------------------------------------------
    const char *input_file = argv[1];
    std::cout<< "Reading input measurements: "<< argv[1] << std::endl;
    
    std::ifstream myfile (input_file);
    if (!myfile.is_open())
    {
        std::cerr << "Unable to open file " << input_file << std::endl;
        std::exit(-1);
    }
    
    int m, n, f; // # rel rots, # abs rots, # fixed abs rots
    myfile >> m >> n >> f;
    
    std::cout<< "# Absolute Rotations : " << n << std::endl;
    std::cout<< "# Relative Rotations : " << m << std::endl;
    
    I_t I;
    I.reserve(m);
    Mat QQ = Mat::Zero(m,4);
    Mat Q = Mat::Zero(n,4);
    
    int i=0, j=0;
    
    // read rel rots
    std::set<int> vertices;
    int e1, e2;
    double x, y, z, w;
    while (j < m)
    {
        if (myfile >> e1 >> e2 >> w >> x >> y >> z)
        {
            I.push_back(std::make_pair(e1, e2) ); // NOT READY TO BE USED
            vertices.insert(e1);
            vertices.insert(e2);
            QQ.row(j++) << x, y, z, w;
        }
        else
        {
            std::cerr << "Corrupt input file: inconsistent number of connections." << std::endl;
            std::exit(-1);
        }
    }
    
    e1 = 0;
    std::map<int,int> vertex_to_idx;
    for (auto const& x : vertices)
    {
        vertex_to_idx[x] = e1++;
    }

    for (auto &c: I)
    {
        c.first = vertex_to_idx[c.first];
        c.second = vertex_to_idx[c.second];
    }
    
    
    // read abs rots
    while (i < n) // f
    {
        if (myfile >> w >> x >> y >> z)
        {
            Q.row(i++) << x, y, z, w;
        }
        else
        {
            break;
        }
    }
    myfile.close();
    
    if (i<f)
    {
        std::cerr << "Insuficient number of absolute rotations. At least "<<f<< " must be given." << std::endl;
        std::exit(-1);
    }
    
    //check input_file
    int max = -1;
    for (auto e: I)
    {
        if (e.second > max)
            max = e.second;
    }
    if ( n != max+1)
    {
        std::cerr << "Corrupt input file: check abs rotations" << std::endl;
        std::exit(-1);
    }

    //---- Output file ----------------------------------------------------
    const char *output_file = (argc-1>1) ? argv[2]: "rcd_out.txt";
    //std::cout << "Output file: " << output_file <<std::endl;
    
    // ---- Maximum num of epoch -----------------------------------------------
    const int rcd_max_epoch = (argc-1>2) ? std::atoi(argv[3]):10000;
    //std::cout<< "Global maximum number of epoch: " << rcd_max_epoch << std::endl;

    // ---- Minimum Spanning Tree Initialisation
    const int FLAG_mst_init = (argc-1>3) ? std::atoi(argv[4]):0;
    //std::cout<< "MST Initialisation: " << FLAG_mst_init << std::endl;


    enforce_directions(QQ, I);
    
    SpMat R = make_R(n, QQ, I);

    
    if (f==0)
    {
        Q.row(0) << 0, 0, 0, 1;
        std::cout << "set first abs rot = I" << std::endl;
        f=1;
    }

    const int init_f = (i>f) ? i : f;
    if (FLAG_mst_init == 1)
    {
        init_mst(Q, QQ, I, init_f);
    }
   
    int rcd_num_epoch = 0; double rcd_runtime = 0; 
    double final_objval = 0;

    rcd(I, R, Q, f, rcd_max_epoch, rcd_num_epoch, rcd_runtime, final_objval);
    
    quat_normalised(Q, f);

    std::cout << std::endl;
    printf("RCD Report: Epoch: %d, Runtime[s]:%.2e, Final cost: %.12e, Termination: CONVERGENCE\n", 
    rcd_num_epoch, rcd_runtime*1e-9, final_objval);

    std::ofstream out(output_file);
    if (out.is_open())
    {
        Mat Q2 = Q;
        Q2.col(0) = Q.col(3);
        Q2.col(1) = Q.col(0);
        Q2.col(2) = Q.col(1);
        Q2.col(3) = Q.col(2);
        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision);
        out << Q2.format(HeavyFmt) << "\n";
        out.close();
    }
    else
    {
        std::cerr << "Error: Unable to save results." << std::endl;
    }
    
    return 0;
}
