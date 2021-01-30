#ifndef _DTMC_LC_H
#define _DTMC_LC_H

#include "Eigen/Dense"
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

// observable
struct observable
{
    // energy
    double E;
    // geometric
    double I2H2;
    std::vector<double> Les;
    // crystalline
    double Tp2uu;
    double Tuuc;
    // coupling
    double Tun2;
    double IKun2;
    // miscellany
    double IdA; // integral of dA
    double I2H; // integral of dA(2H)
    double IK;  // integral of dA(K)
    //TODO: implement this observable
    int Bond_num; // total number of bonds
};
struct vertex
{
    // configuration related
    std::vector<double> R{0, 0, 0}; // position (x,y,z)
    std::vector<double> u{0, 0, 1}; // rod orientation (ux,uy,uz)
    std::vector<double> n{0, 0, 0}; // membrane normal

    std::vector<int> nei; // index of neighbors
    // nei[k+1],nei[k] are connected!!!
    int edge_num; // which edge
    std::vector<int> edge_nei;
    // neighbors form edge with this one (if exist)

    std::vector<int> nei2flip; // index in nei that could be flipped with

    // measurement related
    std::vector<double> dAn2H; // in bulk: (dA, 2H), on edge (0,0)
    // energy related (directly)
    double ds; // beads in bulk: 0 , beads on edge 0.5*(l_{i+}+l_{i-})
    // double dskg; // in bulk: 0, on edge: kg*ds
    double dAK; // in bulk: K*dA, on edge: 0
    double un2; // local un2
};
class dtmc_lc
{
public:
    // eternal parameters
    double beta; // system temperature
    int N;       // number of beads
    int rb;      // radius for disk-shape initialization
    int Ne;      // number of edges
    int L;       // vertical height of cylinder initialization
    // also used to set fixed distance betwen two beads
    double l0; // tether maximum length
    // double l1;
    //  edge tether maximum length, not needed for lc model
    double kar;  // mean curvature bending stiffness
    double karg; // gaussian curvature bending stiffness just for geodesic curvature
    double lam;  // line tension coefficient
    double Kd;   // liquid crystal interaction moduli
    double Kt;   // liquid crystall twist interaction moduli
    double Cn;   // liquid crystal to membrane normal moduli
    double kard; // depletion-like gaussian curvature bending stiffness

    // system configuration

    std::vector<vertex> mesh; // the mesh
    std::vector<std::vector<int>>
        edge_lists; // list of the edge beads sorted, linked
    std::vector<std::pair<int, int>> bond_list;
    // bond_list[i] is a pair of two bead connected in [bulk!}
    // notice only bond in the bulk are included, and there will be on
    // repetition due to symmetry,
    std::vector<int> fixed_beads; // beads can't be moved

    void mesh_bead_info_update(std::vector<int> ind_relate);

    observable Ob_sys;
    observable Ob_m(std::vector<int> ind_relate,
                    std::vector<std::pair<int, int>> bond_relate);
    void Ob_sys_update(observable Ob_new, observable Ob_old);

    // measure observables related to these ind and bond
    double E_m(observable Ob);

    // randomnumber generators
    std::mt19937 gen;
    std::uniform_int_distribution<> rand_pos;  // random position
    std::uniform_real_distribution<> rand_uni; // uniform distribution

    // initialization
    dtmc_lc(double beta_, int N_, int rb_, int Ne_, int L_, double d0_, double l0_,
            double kar_, double karg_, double lam_, double Kd_, double Kt_, double Cn_,
            double kard_);

    // put the beads and bonds in to position accordingly
    void init_rhombus_shape(double d0_);
    void init_disk_shape(int rb, double d0_);
    void init_cylinder_shape(double d0_);

    void reset_config();
    void push_neis_back(int i, std::vector<int> nei_dist);
    void push_eneis_back(int i, std::vector<int> enei_dist);
    void push_bneis_list(int i, std::vector<int> bnei_dist);

    // bond in bulk to bond_list
    // d0_ initial vertex-vertex distance
    // other parameters are the same as in the class

    // local energy-related measurement
    // _m stand for measurement
    double ds_m(int index);                 // length of the local edge index
    std::vector<double> dAn2H_m(int index); // measure and set dA and |2H|
    double dAK_m(int index);                // K*dA measure the gauss curvature
    double dskg_m(int index);               // kg*ds, for gaussian
    std::vector<double> n_m(int index);     // surface normal
    double p2uu_m(int ind_i, int ind_j);    // spin spin interaction of i and j
    double uuc_m(int ind_i, int ind_j);     // spin twist interaction of i and j
    double un2_m(int index);                // spin normal interaction of i
    double dEgeo_m(int index);              // energy of the local vertex

    // Gyration tensor measurement
    std::vector<double> Gij_m();
    // normalized twist from center of mass measurement
    std::vector<double> un2r_m(int bin_num);
    // modified
    // (u(r)*n(r))^2/(1-(u(r)*r)^2), how much director tilt about local membrane
    // normal, disregard the splay parts

    std::vector<double> unu2r_m(int bin_num);
    // (u(r)*nu)^2, how mush director twist about membrane nematic director
    std::vector<double> un2dis_m(int bin_num);
    // distribution of un2 among beads, good indication for pi wall formation
    std::vector<double> nu0nu2l_m(int bin_num);
    //(nu * nu(l))^2, how each slab nematic director
    std::vector<double> nunu2lcov_m(int bin_num);
    //\sum_x (nu(s) * nu(s+l))^2, how each slab nematic director
    // twist, correlation function
    std::vector<double> Qevec_m(std::vector<int> bead_list);

    // useful tools
    double distance2(int ind_1, int ind_2);
    double distance(int ind_1, int ind_2);
    double distancefp(int ind_1, std::vector<double> p);
    double innerproduct(std::vector<double> a, std::vector<double> b);
    void delete_bond_list(int ind_i, int ind_j);
    // delete bond i-j, including both <i,j> and <j,i> in the bulk bond_list

    observable get_related_local_observables(std::vector<int> ind_list);

    // Monte Carlo updates
    int bead_metropolis(double delta_s);
    // execute bead move update in [-ds,ds]^3
    // return 1: accepted, 0: rejected

    int spin_metropolis(double delta_theta);
    // execute spin move update along random axis in [-dtheta,dtheta]
    // return 1: accepted, 0: rejected 0

    int bond_metropolis();
    // execute bond switch update
    // return 1: accepted, 0: rejected

    int edge_metropolis();
    // execute bond remove/add update
    // return 1: accepted, 0: rejected

    // experiment
    void State_write(std::string filename);
    // write system state to file
    void State_write_seq(std::string filename, int MC_sweeps, int step_p_sweep,
                         double delta_s, double delta_theta);

    void State_load(std::string state_file);
    // load system state from file
    void Thermal(int MC_sweeps, int step_p_sweep, int beta_steps,
                 double delta_s, double delta_theta);
    // thermalisation of the system, starting from beta=0 (hot start)
    void O_MC_measure(int MC_sweeps, int sweep_p_G, int step_p_sweep,
                      double delta_s, double delta_theta, std::string folder,
                      std::string finfo, int bin_num_r, int bin_num_un2, int bin_num_l);
    // measure the obserables
    // energy versus lambda curve testing

    // little tool
    int energy_check(); // check if the energy set are correct, for debug use
    int sort_nei(int index);
    int list_a_nei_b(std::vector<int> a, int b);
    int check_nei_connect();
    int check_duplication(int ind_i);
};
#endif
