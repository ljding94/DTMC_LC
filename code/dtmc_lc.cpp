#include "dtmc_lc.h"
#include <cmath>
#include <iostream>
#include <random>
#define PI 3.14159265358979323846

// initialization
dtmc_lc::dtmc_lc(double beta_, int N_, int rb_, int Ne_, int L_, double d0_, double l0_, double kar_, double karg_, double lam_, double Kd_, double Kt_, double Cn_, double kard_)
{
    // system related
    beta = beta_;
    N = N_;
    rb = rb_;
    Ne = Ne_;
    L = L_;
    l0 = l0_; // sigma0 is always 1;

    // energy related
    // geometric
    kar = kar_;
    karg = karg_;
    lam = lam_;
    // orientational
    Kd = Kd_;
    Kt = Kt_;
    // coupling
    Cn = Cn_;
    kard = kard_;

    edge_lists.resize(Ne);
    for (int n = 0; n < Ne; n++)
    {
        edge_lists[n].clear();
    }
    bond_list.clear();

    // new way, assume (N+2)%L != 0
    // if L=sqrt(N), N can't be 700
    if (Ne == 1)
    {
        if (rb != 0)
        {
            init_disk_shape(rb, d0_);
        }
        else
        {
            init_rhombus_shape(d0_);
        }
    }
    else if (Ne == 2)
    {
        init_cylinder_shape(d0_);
    }
    // above shape setting take care of beads position
    // bond_list (excluding edge bond) and edge_list

    // set inital observable value
    Ob_sys.E = 0;
    Ob_sys.I2H2 = 0;
    Ob_sys.Les.resize(Ne);
    for (int e = 0; e < Ne; e++)
    {
        Ob_sys.Les[e] = 0;
    }
    Ob_sys.Tp2uu = 0;
    Ob_sys.Tuuc = 0;
    Ob_sys.Tun2 = 0;
    Ob_sys.IKun2 = 0;
    Ob_sys.IdA = 0;
    Ob_sys.I2H = 0;
    Ob_sys.IK = 0;
    Ob_sys.Bond_num = 0.5 * bond_list.size();
    for (int n = 0; n < Ne; n++)
    {
        Ob_sys.Bond_num += edge_lists[n].size();
    }

    for (int i = 0; i < mesh.size(); i++)
    {
        // vertex info measurement
        mesh[i].n = n_m(i);
        mesh[i].dAn2H = dAn2H_m(i);
        mesh[i].ds = ds_m(i);
        mesh[i].un2 = un2_m(i);
        mesh[i].dAK = dAK_m(i);
    }
    // geometric energy related
    for (int i = 0; i < mesh.size(); i++)
    {
        Ob_sys.I2H2 += mesh[i].dAn2H[0] * mesh[i].dAn2H[1] * mesh[i].dAn2H[1];
        if (mesh[i].edge_num != -1)
        {
            Ob_sys.Les[mesh[i].edge_num] += mesh[i].ds;
        }
        // crystalline energy related
        for (int j = 0; j < mesh[i].nei.size(); j++)
        {
            // factor of 1/2 is for double counting
            Ob_sys.Tp2uu += 0.5 * p2uu_m(i, mesh[i].nei[j]);
            Ob_sys.Tuuc += 0.5 * uuc_m(i, mesh[i].nei[j]);
        }
        // coupling related
        Ob_sys.Tun2 += mesh[i].un2;
        Ob_sys.IKun2 += mesh[i].dAK * mesh[i].un2;
        // miscellany
        Ob_sys.IdA += mesh[i].dAn2H[0];
        Ob_sys.I2H += mesh[i].dAn2H[1] * mesh[i].dAn2H[0];
        Ob_sys.IK += mesh[i].dAK;
    }
    Ob_sys.E = 0.5 * Cn * N;
    Ob_sys.E += E_m(Ob_sys);

    // set random number generators
    std::random_device rd;
    std::mt19937 gen_set(rd());
    std::uniform_int_distribution<> rand_pos_set(0, mesh.size() - 1); // random pos
    std::uniform_real_distribution<> rand_uni_set(0, 1);
    gen = gen_set;
    rand_pos = rand_pos_set;
    rand_uni = rand_uni_set;

} // end of triangulation

void dtmc_lc::init_rhombus_shape(double d0_)
{
    int x_n, y_n; // position of the vertex in the two vector coordinate
    mesh.resize(N);
    if ((N + 2) % L == 1 || (N + 2) % L >= (L - 1))
    {
        std::cout << "invalid N!\n";
    }
    for (int i = 0; i < N; i++)
    {
        // assign position
        if (i < N - N % L - 2)
        {
            x_n = (i + 1) % L;
            y_n = (i + 1) / L;
        }
        else
        {
            // ignore upper-right corner due to #-of-nei limit
            x_n = (i + 2) % L;
            y_n = (i + 2) / L;
        }

        mesh[i].R[0] = d0_ * (x_n + 0.5 * y_n);
        mesh[i].R[1] = d0_ * 0.5 * std::sqrt(3) * y_n;
        mesh[i].R[2] = 0;

        // put bonds
        if (x_n == L - 1)
        {
            // right line
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (y_n == 0)
            {
                // lower-right corner
                push_neis_back(i, {L, L - 1, -1});
                push_eneis_back(i, {L, -1});
                push_bneis_list(i, {L - 1});
            }
            else if (y_n == (N + 1) / L - 2)
            {
                // upper-right corner
                push_neis_back(i, {L - 1, -1, -L});
                push_eneis_back(i, {L - 1, -L});
                push_bneis_list(i, {-1});
            }
            else
            {
                push_neis_back(i, {L, L - 1, -1, -L});
                push_eneis_back(i, {L, -L});
                push_bneis_list(i, {L - 1, -1});
            }
        }
        else if (y_n == (N + 1) / L - 1 && x_n != 0)
        {
            if (x_n > (N + 1) % L)
            {
                edge_lists[0].push_back(i);
                mesh[i].edge_num = 0;
                if (x_n == L - 2)
                {
                    if ((N + 2) % L < (L - 2))
                    {
                        push_neis_back(i, {-1, -L, -L + 1});
                        push_eneis_back(i, {-1, -L + 1});
                        push_bneis_list(i, {-L});
                    }
                    else
                    {
                        push_neis_back(i, {L - 2, -1, -L, -L + 1});
                        push_eneis_back(i, {L - 2, -L + 1});
                        push_bneis_list(i, {-1, -L});
                    }
                }
                else if (x_n == (N + 1) % L + 1)
                {
                    push_neis_back(i, {L - 2, -1, -L, -L + 1, 1});
                    push_eneis_back(i, {L - 2, 1});
                    push_bneis_list(i, {-1, -L, -L + 1});
                }
                else
                {
                    push_neis_back(i, {-1, -L, -L + 1, 1});
                    push_eneis_back(i, {-1, 1});
                    push_bneis_list(i, {-L, -L + 1});
                }
            }
            else
            {
                // in the bulk
                mesh[i].edge_num = -1;
                push_neis_back(i, {1, L - 1, L - 2, -1, -L, -L + 1});
                push_bneis_list(i, {1, L - 1, L - 2, -1, -L, -L + 1});
            }
        }
        else if (y_n == (N + 1) / L)
        {
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (x_n == (N + 1) % L)
            {
                push_neis_back(i, {-1, -L + 1, -L + 2});
                push_eneis_back(i, {-1, -L + 2});
                push_bneis_list(i, {-L + 1});
            }
            else if (x_n == 0)
            {
                push_neis_back(i, {-L + 1, -L + 2, 1});
                push_eneis_back(i, {-L + 1, 1});
                push_bneis_list(i, {-L + 2});
            }
            else
            {
                push_neis_back(i, {-1, -L + 1, -L + 2, 1});
                push_eneis_back(i, {-1, 1});
                push_bneis_list(i, {-L + 1, -L + 2});
            }
        }
        else if (x_n == 0)
        {
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (y_n == 1)
            {
                push_neis_back(i, {-L + 1, 1, L});
                push_eneis_back(i, {-L + 1, L});
                push_bneis_list(i, {1});
            }
            else if (y_n == (N + 1) / L - 1)
            {
                push_neis_back(i, {-L, -L + 1, 1, L - 1});
                push_eneis_back(i, {-L, L - 1});
                push_bneis_list(i, {-L + 1, 1});
            }
            else
            {
                push_neis_back(i, {-L, -L + 1, 1, L});
                push_eneis_back(i, {-L, L});
                push_bneis_list(i, {-L + 1, 1});
            }
        }
        else if (y_n == 0)
        {
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (x_n == 1)
            {
                push_neis_back(i, {1, L, L - 1});
                push_eneis_back(i, {1, L - 1});
                push_bneis_list(i, {L});
            }
            else
            {
                push_neis_back(i, {1, L, L - 1, -1});
                push_eneis_back(i, {1, -1});
                push_bneis_list(i, {L, L - 1});
            }
        }
        else
        {
            mesh[i].edge_num = -1; // in the bulk
            push_neis_back(i, {1, L, L - 1, -1, -L, -L + 1});
            push_bneis_list(i, {1, L, L - 1, -1, -L, -L + 1});
        }
    }
}

void dtmc_lc::init_disk_shape(int rb, double d0_)
{ //TODO: set up dis shape initialization bases old code
    int x_n, y_n;
    double disi2c, dise2c;
    int nei_pos = 0;

    std::vector<vertex> init_mesh; // initial rhombus-shape mesh
    init_mesh.resize((4 * rb + 1) * (4 * rb + 1));
    std::vector<int> init2after_index(init_mesh.size(), -1);

    // f[# in init_mesh] = # in mesh, -1 as default
    mesh.clear();
    std::vector<int> nei_dist = {1, 4 * rb + 1, 4 * rb, -1, -4 * rb - 1, -4 * rb};
    dise2c = std::pow(d0_ * rb, 2);
    for (int i = 0; i < init_mesh.size(); i++)
    {
        x_n = i % (4 * rb + 1) - 2 * rb;
        y_n = i / (4 * rb + 1) - 2 * rb; //-2r is for centralization
        init_mesh[i].R[0] = d0_ * (x_n + 0.5 * y_n);
        init_mesh[i].R[1] = d0_ * 0.5 * std::sqrt(3) * y_n;
        init_mesh[i].R[2] = 0;
        // put neighbors as if every on is in-bulk
        init_mesh[i].edge_num = -1;
        //ignore bead on init_mesh edge
        for (int j = 0; j < 6; j++)
        {
            nei_pos = i + nei_dist[j];
            if (0 <= nei_pos && nei_pos < init_mesh.size())
            {
                //std::cout << "nei_pos=" << nei_pos << "\n";
                init_mesh[i].nei.push_back(nei_pos);
            }
        }
        // find position in real mesh;

        disi2c = std::pow(init_mesh[i].R[0], 2) + std::pow(init_mesh[i].R[1], 2);
        if (disi2c <= dise2c)
        {
            mesh.push_back(init_mesh[i]);
            init2after_index[i] = mesh.size() - 1;
        }
    }
    //corresting mesh nei index
    N = mesh.size(); // assign the system size parameter
    std::vector<int> new_nei_cache;
    std::vector<int> ind_nei_cache;

    for (int i = 0; i < mesh.size(); i++)
    {
        new_nei_cache.clear();
        ind_nei_cache.clear();
        for (int j = 0; j < mesh[i].nei.size(); j++)
        {
            if (init2after_index[mesh[i].nei[j]] != -1)
            {
                new_nei_cache.push_back(init2after_index[mesh[i].nei[j]]);
                ind_nei_cache.push_back(init2after_index[mesh[i].nei[j]] - i);
            }
        } // end for mesh[i].nei

        mesh[i].nei = new_nei_cache;

        // edge bead?
        if (mesh[i].nei.size() < 6)
        {
            mesh[i].edge_num = 0;
            edge_lists[0].push_back(i);
        }
        else
        { // if in-bulk, add all bond to bond list
            push_bneis_list(i, ind_nei_cache);
        }
    } // end of mesh connection
    // time take care of beads on the edge
    //FIXME: fix nei order
    // FIXME:
    // edge_nei[0] must be nei[0] edge_nei[1] must be nei[-1]
    int nei_p, ind;
    std::vector<int> enei_cache; // neighbors on the edge
    std::vector<int> bnei_cache; // neighbors in the bulk
    for (int i = 0; i < edge_lists[0].size(); i++)
    {
        ind = edge_lists[0][i];
        enei_cache.clear();
        bnei_cache.clear();
        for (int j = 0; j < mesh[ind].nei.size(); j++)
        {
            nei_p = mesh[ind].nei[j];
            if (mesh[nei_p].nei.size() < 6)
            { // nei_p is on the edge
                enei_cache.push_back(nei_p - ind);
            }
            else
            { // nei_p is in the bulk
                bnei_cache.push_back(nei_p - ind);
            }
        } // end for ind nei

        push_eneis_back(ind, enei_cache);
        // to be sorted
        push_bneis_list(ind, bnei_cache);
    }
    int eind, eind_next, eind_cache;
    eind = edge_lists[0][0];

    do
    {
        eind_next = mesh[eind].edge_nei[1];
        eind_cache = mesh[eind_next].edge_nei[0];
        if (eind_cache != eind)
        { //fix edge_nei direction of eind_next
            mesh[eind_next].edge_nei[0] = eind;
            mesh[eind_next].edge_nei[1] = eind_cache;
        }
        sort_nei(eind_next); // sort nei list based on edge-nei
        eind = eind_next;
    } while (eind != edge_lists[0][0]);
}

void dtmc_lc::init_cylinder_shape(double d0_)
{
    int x_n, y_n; // position of the vertex in the two vector coordinate
    // cylinder initial shape
    int Lr = N / L;
    double R = d0_ / (2 * std::sin(PI / Lr));
    mesh.resize(N);
    for (int i = 0; i < N; i++)
    {
        // assign position
        x_n = i % Lr;
        y_n = i / Lr;
        mesh[i].R[0] = R * std::cos(2 * PI * (x_n + 0.5 * y_n) / Lr);
        mesh[i].R[1] = R * std::sin(2 * PI * (x_n + 0.5 * y_n) / Lr);
        mesh[i].R[2] = d0_ * 0.5 * std::sqrt(3) * y_n;

        // put bonds
        if (y_n == 0)
        {
            // lower edge
            mesh[i].edge_num = 0;
            edge_lists[0].push_back(i);
            if (x_n == 0)
            {
                push_neis_back(i, {1, Lr, 2 * Lr - 1, Lr - 1});
                push_eneis_back(i, {1, Lr - 1});
                push_bneis_list(i, {Lr, 2 * Lr - 1});
            }
            else if (x_n == Lr - 1)
            {
                push_neis_back(i, {-Lr + 1, Lr, Lr - 1, -1});
                push_eneis_back(i, {-Lr + 1, -1});
                push_bneis_list(i, {Lr, Lr - 1});
            }
            else
            {
                push_neis_back(i, {1, Lr, Lr - 1, -1});
                push_eneis_back(i, {1, -1});
                push_bneis_list(i, {Lr, Lr - 1});
            }
        }
        else if (y_n == L - 1)
        {
            // upper edge
            mesh[i].edge_num = 1;
            edge_lists[1].push_back(i);
            if (x_n == 0)
            {
                push_neis_back(i, {Lr - 1, -Lr, -Lr + 1, 1});
                push_eneis_back(i, {Lr - 1, 1});
                push_bneis_list(i, {-Lr, -Lr + 1});
            }
            else if (x_n == Lr - 1)
            {
                push_neis_back(i, {-1, -Lr, -2 * Lr + 1, -Lr + 1});
                push_eneis_back(i, {-1, -Lr + 1});
                push_bneis_list(i, {-Lr, -2 * Lr + 1});
            }
            else
            {
                push_neis_back(i, {-1, -Lr, -Lr + 1, 1});
                push_eneis_back(i, {-1, 1});
                push_bneis_list(i, {-Lr, -Lr + 1});
            }
        }
        else
        {
            // internal
            mesh[i].edge_num = -1;
            if (x_n == 0)
            {
                push_neis_back(i, {1, Lr, 2 * Lr - 1, Lr - 1, -Lr, -Lr + 1});
                push_bneis_list(i, {1, Lr, 2 * Lr - 1, Lr - 1, -Lr, -Lr + 1});
            }
            else if (x_n == Lr - 1)
            {
                push_neis_back(i, {-Lr + 1, Lr, Lr - 1, -1, -Lr, -2 * Lr + 1});
                push_bneis_list(i, {-Lr + 1, Lr, Lr - 1, -1, -Lr, -2 * Lr + 1});
            }
            else
            {
                push_neis_back(i, {1, Lr, Lr - 1, -1, -Lr, -Lr + 1});
                push_bneis_list(i, {1, Lr, Lr - 1, -1, -Lr, -Lr + 1});
            }
        }
    }
}

void dtmc_lc::push_neis_back(int i, std::vector<int> nei_dist)
{
    for (int j = 0; j < nei_dist.size(); j++)
    {
        mesh[i].nei.push_back(i + nei_dist[j]);
    }
}
void dtmc_lc::push_eneis_back(int i, std::vector<int> enei_dist)
{
    for (int j = 0; j < enei_dist.size(); j++)
    {
        mesh[i].edge_nei.push_back(i + enei_dist[j]);
    }
}
void dtmc_lc::push_bneis_list(int i, std::vector<int> bnei_dist)
{
    std::pair<int, int> bond;
    for (int j = 0; j < bnei_dist.size(); j++)
    {
        bond.first = i;
        bond.second = i + bnei_dist[j];
        bond_list.push_back(bond);
    }
}
void dtmc_lc::delete_bond_list(int ind_i, int ind_j)
{
    std::pair<int, int> bond0, bond1;
    bond0.first = ind_i;
    bond0.second = ind_j;
    bond1.first = ind_j;
    bond1.second = ind_i;
    for (int i = 0; i < bond_list.size(); i++)
    {
        if (bond_list[i] == bond0 || bond_list[i] == bond1)
        {
            bond_list.erase(bond_list.begin() + i);
            i--;
        }
    }
}