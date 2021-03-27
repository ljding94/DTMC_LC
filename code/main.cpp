// Copyright[2020] [Lijie Ding]
#include "dtmc_lc.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char const *argv[])
{
    std::clock_t c_start = std::clock();
    double beta = 1;
    int N;      // number of beads
    int rb = 0; //radius for dis-shape initialization
    int Ne = 1; // number of edges
    int Lr;     // number of bead per line
    int L;      // number of bead on vertical direction, used for cylinder
                // initialization
    // also help adjusting the fixed distance of two beads
    double d0 = 1.5;
    double l0 = 1.73;
    double kar;
    double karg;
    double lam;
    double Kd;
    double q; // Kt = Kd*q
    double Kt;
    double Cn;
    double kard;
    // double delta_s = 0.2;
    double delta_s = 0.1;
    double delta_theta = 0.5;

    std::string folder;

    N = std::atoi(argv[1]);
    Ne = std::atoi(argv[2]);
    L = std::atoi(argv[3]);
    kar = std::atof(argv[4]);
    karg = std::atof(argv[5]);
    lam = std::atof(argv[6]);
    Kd = std::atof(argv[7]);
    q = std::atof(argv[8]);
    Cn = std::atof(argv[9]);
    kard = std::atof(argv[10]);

    Kt = q * Kd;
    // will skip bin measurement if set to 0;
    int bn_r = 0;       // penetration depth bin_num
    int bn_un2 = 0;     //un2 distribution measurement bin_num
    int bn_l = 0;       // extension experiment bin_num
    int bn_g = 0;      // pair distribution related bin_num
    double del_r = 0.2; // pair distribution r increment
    std::vector<int> bin_nums{bn_r, bn_un2, bn_l, bn_g};
    bool fix_bead_on = 0;
    if (L == -1)
    {
        //use rhombus shape initilization
        rb = int(std::sqrt(N / 4));
    }
    else if (L == 0)
    {
        // use rhombus shape initialization
        L = int(std::sqrt(N));
        if (Ne == 2)
        {
            //use cylinder initialization
            L = 10;
            if (N % L != 0)
            {
                std::cout << "N % L != 0 \n";
            }
            Lr = N / L;
        }
    }
    else if (L > 0)
    { // fix two beads, use rhombus shape initilization
        fix_bead_on = 1;
    }
    std::string finfo;
    finfo = "N" + std::string(argv[1]) + "_Ne" + std::string(argv[2]) + "_L" + std::string(argv[3]) + "_kar" + std::string(argv[4]) + "_karg" + std::string(argv[5]) + "_lam" + std::string(argv[6]) + "_Kd" + std::string(argv[7]) + "_q" + std::string(argv[8]) + "_Cn" + std::string(argv[9]) + "_kard" + std::string(argv[10]);
    dtmc_lc membrane(beta, N, rb, Ne, L, d0, l0, kar, karg, lam, Kd, Kt, Cn, kard);
    N = membrane.mesh.size();
    if (argc == 12)
    {
        // use "triangulation kar lam local" for local running
        // used for local running!
        std::cout << "running on local machine\n";
        folder = "../data/scratch_local";
        if (fix_bead_on)
        {
            int fix_bead0 = int(N / (2 * L)) * L - 1;
            if (std::abs((membrane.mesh[fix_bead0 + L - 1].R[0] -
                          membrane.mesh[fix_bead0].R[0]) -
                         (L * d0 - d0)) > 0.1)
            {
                std::cout << "not exact fixed length\n";
            }
            membrane.fixed_beads = {fix_bead0, fix_bead0 + L - 1};
        }
        membrane.State_write(folder + "/State_" + finfo + "_init.txt");
        membrane.Thermal(0, int(N / (delta_s * delta_s)) + 1, 2, delta_s, delta_theta);
        // membrane.O_MC_measure(5, 1, int(N / (ds * ds)) + 1, ds,
        membrane.O_MC_measure(20, 5, int(N / (delta_s * delta_s)) + 1, delta_s, delta_theta, folder, finfo, bin_nums, del_r);
        // membrane.O_MC_measure(2, 1, 0, delta_s, delta_theta, folder,
        // finfo,bin_num_r);
        membrane.State_write(folder + "/State_" + finfo + ".txt");

        return 0;
    }
    else if (argc == 11)
    {
        // ccv running
        folder = "/users/lding3/scratch";

        if (fix_bead_on)
        {
            int fix_bead0 = int(N / (2 * L)) * L - 1;
            if ((membrane.mesh[fix_bead0 + L - 1].R[0] -
                 membrane.mesh[fix_bead0].R[0]) != L * d0 - d0)
            {
                std::cout << "not exact fixed length\n";
            }
            membrane.fixed_beads = {fix_bead0, fix_bead0 + L - 1};
        }
        // membrane.Thermal(500, N / (delta_s * delta_s), 10,
        // delta_s,delta_theta);
        // used 2000, 4000 for manuscript
        membrane.Thermal(2000, N / (delta_s * delta_s), 1, delta_s,
                         delta_theta);
        membrane.O_MC_measure(4000, 50, N / (delta_s * delta_s) + 1, delta_s,
                              delta_theta, folder, finfo, bin_nums, del_r);
        membrane.State_write(folder + "/State_" + finfo + ".txt");

        return 0;
    }
}