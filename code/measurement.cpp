#include "dtmc_lc.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#define PI 3.14159265358979323846

#pragma region : tools
double dtmc_lc::innerproduct(std::vector<double> a, std::vector<double> b)
{
    double result = 0;
    for (int k = 0; k < a.size(); k++)
    {
        result += a[k] * b[k];
    }
    return result;
}
std::vector<double> crossproduct(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> cp{0, 0, 0};
    cp[0] = a[1] * b[2] - a[2] * b[1];
    cp[1] = a[2] * b[0] - a[0] * b[2];
    cp[2] = a[0] * b[1] - a[1] * b[0];
    return cp;
}
double dtmc_lc::distance2(int ind_1, int ind_2)
{
    double result = 0;
    for (int k = 0; k < 3; k++)
    {
        result += std::pow(mesh[ind_2].R[k] - mesh[ind_1].R[k], 2);
    }
    return result;
}
double dtmc_lc::distance(int ind_1, int ind_2)
{
    double result = distance2(ind_1, ind_2);
    result = std::sqrt(result);
    return result;
}
double dtmc_lc::distancefp(int ind_1, std::vector<double> p)
{
    double result = 0;
    for (int k = 0; k < 3; k++)
    {
        result += std::pow(p[k] - mesh[ind_1].R[k], 2);
    }
    result = std::sqrt(result);
    return result;
}

int dtmc_lc::sort_nei(int index)
{
    // nei list can be sorted as long as two edge_nei are next to each other in the nei list
    std::vector<int> index_nei = mesh[index].nei;
    std::vector<int> index_enei = mesh[index].edge_nei;
    if (mesh[index].edge_nei.size() == 0)
    {
        return 0;
    }
    /*
    don't understand why I added this condition
    if (mesh[index].nei.size() == 3)
    {
        return 0;
    }
    */
    while (mesh[index].nei[0] != mesh[index].edge_nei[0])
    {
        mesh[index].nei.push_back(mesh[index].nei[0]);
        mesh[index].nei.erase(mesh[index].nei.begin());
    }
    if (mesh[index].nei[1] == mesh[index].edge_nei[1])
    {
        mesh[index].nei.push_back(mesh[index].nei[0]);
        mesh[index].nei.erase(mesh[index].nei.begin());
    }
    else if (mesh[index].nei.back() != mesh[index].edge_nei[1])
    {
        std::cout << "wrong, can't sort neighbors!\n";
        return 0;
    }
    return 1;
}

int dtmc_lc::list_a_nei_b(std::vector<int> a, int b)
{
    for (int j = 0; j < a.size(); j++)
    {
        if (a[j] == b)
        {
            return j;
        }
    }
    return -1;
}

int dtmc_lc::check_nei_connect()
{
    int j_next;
    for (int i = 0; i < mesh.size(); i++)
    {
        for (int j = 0; j < mesh[i].nei.size(); j++)
        {
            j_next = 0;
        }
    }
    return 0;
}
int dtmc_lc::check_duplication(int ind_i)
{
    for (int j = 0; j < mesh[ind_i].nei.size() - 1; j++)
    {
        for (int k = j + 1; k < mesh[ind_i].nei.size(); k++)
        {
            if (mesh[ind_i].nei[j] == mesh[ind_i].nei[k])
            {
                std::cout << "duplicated nei!\n";
                return 1;
            }
        }
    }
    return 0;
}

#pragma endregion

#pragma region : energy related

void dtmc_lc::mesh_bead_info_update(std::vector<int> ind_relate)
{
    int ind;
    for (int i = 0; i < ind_relate.size(); i++)
    {
        ind = ind_relate[i];
        mesh[ind].n = n_m(ind);
        mesh[ind].dAn2H = dAn2H_m(ind);
        mesh[ind].ds = ds_m(ind);
        mesh[ind].dAK = dAK_m(ind);
        mesh[ind].un2 = un2_m(ind);
    }
}

observable dtmc_lc::Ob_m(std::vector<int> ind_relate,
                         std::vector<std::pair<int, int>> bond_relate)
{

    int ind, ind_i, ind_j;
    observable Ob;
    // set inital observable value
    Ob.E = 0;
    Ob.I2H2 = 0;
    Ob.Les.resize(Ne);
    for (int e = 0; e < Ne; e++)
    {
        Ob.Les[e] = 0;
    }
    Ob.Tp2uu = 0;
    Ob.Tuuc = 0;
    Ob.Tun2 = 0;
    Ob.IKun2 = 0;
    Ob.IdA = 0;
    Ob.I2H = 0;
    Ob.IK = 0;
    Ob.Bond_num = 0;
    for (int i = 0; i < ind_relate.size(); i++)
    {
        ind = ind_relate[i];
        // geometric terms
        Ob.I2H2 += mesh[ind].dAn2H[0] * mesh[ind].dAn2H[1] * mesh[ind].dAn2H[1];
        Ob.Les[mesh[ind].edge_num] += mesh[ind].ds;
        // crystalline terms
        // coupling terms
        Ob.Tun2 += mesh[ind].un2;
        Ob.IKun2 += mesh[ind].dAK * mesh[ind].un2;
        // miscellany terms
        Ob.IdA += mesh[ind].dAn2H[0];
        Ob.I2H += mesh[ind].dAn2H[0] * mesh[ind].dAn2H[1];
        Ob.IK += mesh[ind].dAK;
    }
    // crystalline terms
    for (int i = 0; i < bond_relate.size(); i++)
    {
        ind_i = bond_relate[i].first;
        ind_j = bond_relate[i].second;
        Ob.Tp2uu += p2uu_m(ind_i, ind_j);
        Ob.Tuuc += uuc_m(ind_i, ind_j);
    }
    Ob.Bond_num += bond_relate.size();
    Ob.E = E_m(Ob);
    return Ob;
}
void dtmc_lc::Ob_sys_update(observable Ob_new, observable Ob_old)
{
    Ob_sys.E += Ob_new.E - Ob_old.E;
    Ob_sys.I2H2 += Ob_new.I2H2 - Ob_old.I2H2;
    for (int e = 0; e < Ne; e++)
    {
        Ob_sys.Les[e] += Ob_new.Les[e] - Ob_old.Les[e];
    }
    Ob_sys.Tp2uu += Ob_new.Tp2uu - Ob_old.Tp2uu;
    Ob_sys.Tuuc += Ob_new.Tuuc - Ob_old.Tuuc;
    Ob_sys.Tun2 += Ob_new.Tun2 - Ob_old.Tun2;
    Ob_sys.IKun2 += Ob_new.IKun2 - Ob_old.IKun2;
    Ob_sys.IdA += Ob_new.IdA - Ob_old.IdA;
    Ob_sys.I2H += Ob_new.I2H - Ob_old.I2H;
    Ob_sys.IK += Ob_new.IK - Ob_old.IK;
    Ob_sys.Bond_num += Ob_new.Bond_num - Ob_old.Bond_num;
}
double dtmc_lc::E_m(observable Ob)
{
    double E = 0;
    E += 0.5 * kar * Ob.I2H2 + karg * Ob.IK;
    for (int e = 0; e < Ne; e++)
    {
        E += lam * Ob.Les[e];
    }
    E += -Kd * Ob.Tp2uu - Kt * Ob.Tuuc;
    E += -0.5 * Cn * Ob.Tun2 + kard * Ob.IKun2;
    return E;
}

double dtmc_lc::ds_m(int index)
{
    if (mesh[index].edge_nei.size() == 0)
    {
        return 0;
    }
    double Le_l = 0;
    int nei_ind;
    double distance;
    for (int i = 0; i < mesh[index].edge_nei.size(); i++)
    {
        nei_ind = mesh[index].edge_nei[i];
        distance = 0;
        for (int j = 0; j < 3; j++)
        {
            // Rij[j] = mesh[index].R[j] - mesh[nei_ind].R[j];
            distance += std::pow(mesh[index].R[j] - mesh[nei_ind].R[j], 2);
        }
        distance = std::sqrt(distance);
        Le_l += 0.5 * distance;
    }
    if (Le_l > 2)
    {
        std::cout << "Le_l=" << Le_l << "\n";
    }
    return Le_l;
}

std::vector<double> dtmc_lc::dAn2H_m(int index)
{
    std::vector<double> dAn2H{0.0, 0.0};
    if (mesh[index].edge_nei.size() != 0)
    {
        return dAn2H;
    }
    double l_ij, l_ij0, l_ij2, l_jj0, l_jj2;
    double dot_ij0j, dot_ij2j, cot0, cot2, theta0, theta2;
    double sigma_ij, sigma_i;
    int ind_j0, ind_j, ind_j2;

    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ij0{0, 0, 0};
    std::vector<double> r_ij2{0, 0, 0};
    std::vector<double> r_jj0{0, 0, 0};
    std::vector<double> r_jj2{0, 0, 0};
    /* vertex configuration:
               - (j2)
             -
        (i)- - - (j)
             -
               - (j0)
    */
    sigma_i = 0;
    std::vector<double> dH_i{0, 0, 0}; // 2H_i, depends on r_ij
    for (int j = 0; j < mesh[index].nei.size(); j++)
    {
        // get nei index
        ind_j = mesh[index].nei[j];
        if (j == 0)
        {
            ind_j0 = mesh[index].nei[mesh[index].nei.size() - 1];
            ind_j2 = mesh[index].nei[j + 1];
        }
        else if (j == mesh[index].nei.size() - 1)
        {
            ind_j0 = mesh[index].nei[j - 1];
            ind_j2 = mesh[index].nei[0];
        }
        else
        {
            ind_j0 = mesh[index].nei[j - 1];
            ind_j2 = mesh[index].nei[j + 1];
        }
        // end set nei_ind
        // triangulation site to site vectors
        for (int k = 0; k < 3; k++)
        {
            r_ij[k] = mesh[ind_j].R[k] - mesh[index].R[k];
            r_ij0[k] = mesh[ind_j0].R[k] - mesh[index].R[k];
            r_ij2[k] = mesh[ind_j2].R[k] - mesh[index].R[k];
            r_jj0[k] = mesh[ind_j0].R[k] - mesh[ind_j].R[k];
            r_jj2[k] = mesh[ind_j2].R[k] - mesh[ind_j].R[k];
        }
        // site-site distance
        l_ij = std::sqrt(innerproduct(r_ij, r_ij));
        l_ij0 = std::sqrt(innerproduct(r_ij0, r_ij0));
        l_ij2 = std::sqrt(innerproduct(r_ij2, r_ij2));
        l_jj0 = std::sqrt(innerproduct(r_jj0, r_jj0));
        l_jj2 = std::sqrt(innerproduct(r_jj2, r_jj2));
        // inner product for cot calculation
        dot_ij0j = innerproduct(r_ij0, r_jj0);
        theta0 = std::acos(dot_ij0j / (l_ij0 * l_jj0));
        dot_ij2j = innerproduct(r_ij2, r_jj2);
        theta2 = std::acos(dot_ij2j / (l_ij2 * l_jj2));
        // get cot
        /*
        cot0 = dot_ij0j / std::sqrt(std::pow(l_ij0 * l_jj0, 2) -
                                    std::pow(dot_ij0j, 2));
        cot2 = dot_ij2j / std::sqrt(std::pow(l_ij2 * l_jj2, 2) -
                                    std::pow(dot_ij2j, 2));
        */
        cot0 = std::cos(theta0) / std::sin(theta0);
        cot2 = std::cos(theta2) / std::sin(theta2);
        sigma_ij = 0.5 * l_ij * (cot0 + cot2);
        sigma_i += 0.25 * sigma_ij * l_ij;
        for (int k = 0; k < 3; k++)
        {
            dH_i[k] += sigma_ij / l_ij * r_ij[k];
        }
    } // end for nei
    for (int k = 0; k < 3; k++)
    {
        dH_i[k] = dH_i[k] / sigma_i;
    }
    dAn2H[0] = sigma_i;
    dAn2H[1] = std::sqrt(innerproduct(dH_i, dH_i));
    return dAn2H;
}

double dtmc_lc::dskg_m(int index)
{
    if (mesh[index].edge_nei.size() == 0)
    {
        return 0;
    }
    double kg;
    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ik{0, 0, 0};
    double l_ij, l_ik;
    int ind_j, ind_k;
    double cos_jk;
    /* vertex configuration:
               - (k)
             -
        (i)- - - (j)
    */
    std::vector<int> nei_original;
    // sort neighbors
    nei_original = mesh[index].nei;
    sort_nei(index);
    // get thetai
    // kg = PI - sum(theta_i)
    kg = PI;
    for (int j = 0; j < mesh[index].nei.size() - 1; j++)
    {
        cos_jk = 0;
        l_ij = 0;
        l_ik = 0;
        ind_j = mesh[index].nei[j];
        ind_k = mesh[index].nei[j + 1];
        for (int k = 0; k < 3; k++)
        {
            r_ij[k] = mesh[ind_j].R[k] - mesh[index].R[k];
            l_ij += r_ij[k] * r_ij[k];
            r_ik[k] = mesh[ind_k].R[k] - mesh[index].R[k];
            l_ik += r_ik[k] * r_ik[k];
            cos_jk += r_ij[k] * r_ik[k];
        }
        cos_jk = cos_jk / (std::sqrt(l_ij * l_ik));
        kg -= std::acos(cos_jk);
    }
    mesh[index].nei = nei_original; // put nei sort back;
    return kg;
}

double dtmc_lc::dAK_m(int index)
{
    // account both in-bulk and on-edge beads

    // if (mesh[index].edge_nei.size() != 0) {

    double K;
    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ik{0, 0, 0};
    double l_ij, l_ik;
    int ind_j, ind_k;
    double cos_jk;
    /* vertex configuration:
               - (k)
             -
        (i)- - - (j)
    */
    // get thetai
    // K = 2*PI - sum(theta_i) if in-bulk
    // K = PI - sum(theta_i) if on-edge
    std::vector<int> nei_original;
    // sort neighbors
    nei_original = mesh[index].nei;
    if (mesh[index].edge_nei.size() == 0)
    {
        K = 2 * PI;
    }
    else
    {
        return 0;
        K = PI;
        sort_nei(index);
    }

    for (int j = 0; j < mesh[index].nei.size(); j++)
    {
        cos_jk = 0;
        l_ij = 0;
        l_ik = 0;
        ind_j = mesh[index].nei[j];
        if (j == (mesh[index].nei.size() - 1))
        {
            if (mesh[index].edge_nei.size())
            {
                break;
                mesh[index].nei = nei_original;
            }
            ind_k = mesh[index].nei[0];
        }
        else
        {
            ind_k = mesh[index].nei[j + 1];
        }
        for (int k = 0; k < 3; k++)
        {
            r_ij[k] = mesh[ind_j].R[k] - mesh[index].R[k];
            l_ij += r_ij[k] * r_ij[k];
            r_ik[k] = mesh[ind_k].R[k] - mesh[index].R[k];
            l_ik += r_ik[k] * r_ik[k];
            cos_jk += r_ij[k] * r_ik[k];
        }
        cos_jk = cos_jk / (std::sqrt(l_ij * l_ik));
        K -= std::acos(cos_jk);
    }
    return K;
}

std::vector<double> dtmc_lc::n_m(int index)
{
    // measure normal, work for both in-bulk and on-edge beads
    int on_edge = 0; // flag indicating if it's an edge bead
    std::vector<double> n_ind{0, 0, 0};
    double n_ind_abs;
    int j_next, ind_j, ind_k;
    double l_ij, l_ik;
    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ik{0, 0, 0};
    std::vector<double> n_jk{0, 0, 0};
    double n_jk_abs;
    double theta_jk;
    std::vector<int> nei_original;
    nei_original =
        mesh[index].nei; // thus it has correct size from the beginning
    if (mesh[index].edge_nei.size() != 0)
    {
        on_edge = 1;
        sort_nei(index);
    }
    for (int j = 0; j < (mesh[index].nei.size() - on_edge); j++)
    {
        j_next = (j + 1) % mesh[index].nei.size();
        ind_j = mesh[index].nei[j];
        ind_k = mesh[index].nei[j_next];
        for (int k = 0; k < 3; k++)
        {
            r_ij[k] = mesh[ind_j].R[k] - mesh[index].R[k];
            r_ik[k] = mesh[ind_k].R[k] - mesh[index].R[k];
        }
        l_ij = std::sqrt(innerproduct(r_ij, r_ij));
        l_ik = std::sqrt(innerproduct(r_ik, r_ik));
        n_jk = crossproduct(r_ij, r_ik);
        n_jk_abs = std::sqrt(innerproduct(n_jk, n_jk));
        theta_jk = std::asin(n_jk_abs / (l_ij * l_ik));
        for (int k = 0; k < 3; k++)
        {
            n_jk[k] = n_jk[k] / n_jk_abs;
            n_ind[k] += theta_jk * n_jk[k];
        }
    }
    n_ind_abs = std::sqrt(innerproduct(n_ind, n_ind));
    for (int k = 0; k < 3; k++)
    {
        n_ind[k] = n_ind[k] / n_ind_abs;
    }
    // put nei sort back
    mesh[index].nei = nei_original;

    return n_ind;
}

double dtmc_lc::p2uu_m(int ind_i, int ind_j)
{
    // spin spin interaction of i and j
    double p2uu = std::pow(innerproduct(mesh[ind_i].u, mesh[ind_j].u), 2);
    return 1.5 * p2uu - 0.5;
}
double dtmc_lc::uuc_m(int ind_i, int ind_j)
{
    // spin twist interaction of i and j
    double uuc_local, lij;
    std::vector<double> rij{0, 0, 0};
    for (int k = 0; k < mesh[ind_i].R.size(); k++)
    {
        rij[k] = mesh[ind_j].R[k] - mesh[ind_i].R[k];
    }
    lij = std::sqrt(innerproduct(rij, rij));
    uuc_local = innerproduct(crossproduct(mesh[ind_i].u, mesh[ind_j].u), rij);
    uuc_local = uuc_local * innerproduct(mesh[ind_i].u, mesh[ind_j].u);
    uuc_local = uuc_local / lij;
    return uuc_local;
}

double dtmc_lc::un2_m(int index)
{
    // spin normal interaction of i
    return std::pow(innerproduct(mesh[index].u, mesh[index].n), 2);
}

#pragma endregion

#pragma region : membrane structure related
std::vector<double> dtmc_lc::Gij_m()
{
    // G_ij = 1/(2N)\sum_n\sum_m{(r_i(n)-r_i(m))(r_j(n)-r_j(m))}
    // G_ij = <r_i r_j> - <r_i><r_j>
    std::vector<double> Gij(9, 0);
    double xc = 0, yc = 0, zc = 0;
    for (int n = 0; n < mesh.size(); n++)
    {
        // G_xx
        Gij[0] += mesh[n].R[0] * mesh[n].R[0];
        // G_xy
        Gij[1] += mesh[n].R[0] * mesh[n].R[1];
        // G_xz
        Gij[2] += mesh[n].R[0] * mesh[n].R[2];
        // G_yx Gij[3] = Gij[1]
        // G_yy
        Gij[4] += mesh[n].R[1] * mesh[n].R[1];
        // G_yz
        Gij[5] += mesh[n].R[1] * mesh[n].R[2];
        // G_zx Gij[6] = Gij[2]
        // G_zy Gij[7] = Gij[5]
        // G_zz
        Gij[8] += mesh[n].R[2] * mesh[n].R[2];
        // center
        xc += mesh[n].R[0];
        yc += mesh[n].R[1];
        zc += mesh[n].R[2];
    }
    xc /= N;
    yc /= N;
    zc /= N;
    Gij[0] = Gij[0] / N - xc * xc;
    Gij[1] = Gij[1] / N - xc * yc;
    Gij[2] = Gij[2] / N - xc * zc;
    Gij[3] = Gij[1];
    Gij[4] = Gij[4] / N - yc * yc;
    Gij[5] = Gij[5] / N - yc * zc;
    Gij[6] = Gij[2];
    Gij[7] = Gij[5];
    Gij[8] = Gij[8] / N - zc * zc;
    return Gij;
}

std::vector<double> dtmc_lc::Qevec_m(std::vector<int> bead_list)
{
    int ind;
    std::vector<double> nu_return{0, 0, 0};
    Eigen::Matrix3d Q;
    Q << 0, 0, 0, 0, 0, 0, 0, 0, 0;
    std::vector<double> Rc{0, 0, 0};
    // get Q tensor and center of mass
    for (int i = 0; i < bead_list.size(); i++)
    {
        ind = bead_list[i];

        Q(0, 0) += mesh[ind].u[0] * mesh[ind].u[0];
        Q(0, 1) += mesh[ind].u[0] * mesh[ind].u[1];
        Q(0, 2) += mesh[ind].u[0] * mesh[ind].u[2];
        Q(1, 1) += mesh[ind].u[1] * mesh[ind].u[1];
        Q(1, 2) += mesh[ind].u[1] * mesh[ind].u[2];
        Q(2, 2) += mesh[ind].u[2] * mesh[ind].u[2];
    }
    Q(1, 0) = Q(0, 1);
    Q(2, 0) = Q(0, 1);
    Q(2, 1) = Q(1, 2);
    Q /= mesh.size();
    Q *= 1.5;
    Q(0, 0) -= 0.5;
    Q(1, 1) -= 0.5;
    Q(2, 2) -= 0.5;

    // find largest eigen value and corresponding eigenvector
    Eigen::EigenSolver<Eigen::Matrix3d> es(Q);
    Eigen::Vector3d eval = es.eigenvalues().real();
    Eigen::Matrix3d evec = es.eigenvectors().real();
    Eigen::Vector3d nu;
    double eval_max = eval(0);
    int eval_max_i = 0;
    for (int i = 1; i < 3; i++)
    {
        if (eval(i) - eval_max > 0)
        {
            eval_max = eval(i);
            eval_max_i = i;
        }
    } // bubble it
    nu = evec.col(eval_max_i);
    nu_return[0] = nu(0);
    nu_return[1] = nu(1);
    nu_return[2] = nu(2);
    return nu_return;
}

std::vector<double> dtmc_lc::un2r_m(int bin_num)
{
    // make the range normalized to 1 using average distance to the edge
    std::vector<double> un2r;
    std::vector<int> Countr;
    double r_cache, r_c2e, r2_c2e, r_c2e_std;
    double del_r;
    int bin;
    double r;
    std::vector<double> Rc{0, 0, 0};
    std::vector<double> ri{0, 0, 0}; // position vector from cmo to i
    double ur_cache, ur2_cache;
    // find center of mass
    for (int i = 0; i < mesh.size(); i++)
    {
        Rc[0] += mesh[i].R[0];
        Rc[1] += mesh[i].R[1];
        Rc[2] += mesh[i].R[2];
    }
    Rc[0] /= N;
    Rc[1] /= N;
    Rc[2] /= N;
    // find average edge to center distance
    r_c2e = 0;
    r2_c2e = 0;
    for (int i = 0; i < edge_lists[0].size(); i++)
    {
        r_cache = distancefp(edge_lists[0][i], Rc);
        r_c2e += r_cache;
        r2_c2e += r_cache * r_cache;
    }
    r_c2e /= edge_lists[0].size();
    r2_c2e /= edge_lists[0].size();
    r_c2e_std = std::sqrt(r2_c2e - r_c2e * r_c2e);
    del_r = r_c2e / bin_num;

    // find un2r
    for (int k = 0; k < bin_num; k++)
    {
        un2r.push_back(0);
        Countr.push_back(0);
    }
    for (int i = 0; i < mesh.size(); i++)
    {
        r = distancefp(i, Rc);
        ri[0] = mesh[i].R[0] - Rc[0];
        ri[1] = mesh[i].R[1] - Rc[1];
        ri[2] = mesh[i].R[2] - Rc[2];
        ur_cache = innerproduct(ri, mesh[i].u);
        ur_cache /= r;
        ur2_cache = ur_cache * ur_cache;
        bin = int(r / del_r);
        if (bin < bin_num)
        {
            un2r[bin] += mesh[i].un2 / (1 - ur2_cache);
            Countr[bin] += 1;
        }
    }
    // normalize
    for (int b = 0; b < bin_num; b++)
    {
        if (Countr[b])
        {
            un2r[b] /= Countr[b];
        }
    }
    un2r.insert(un2r.begin(), r_c2e_std);
    un2r.insert(un2r.begin(), r_c2e);
    return un2r;
}

std::vector<double> dtmc_lc::unu2r_m(int bin_num)
{
    std::vector<double> unu2r;
    std::vector<int> Countr;
    std::vector<int> bead_list;
    std::vector<double> nu;
    int bin;
    double r_cache, r_c2e, r2_c2e, r_c2e_std;
    double del_r;
    double r;
    double unu;
    std::vector<double> Rc{0, 0, 0};
    // center of mass
    for (int i = 0; i < mesh.size(); i++)
    {
        Rc[0] += mesh[i].R[0];
        Rc[1] += mesh[i].R[1];
        Rc[2] += mesh[i].R[2];
    }
    Rc[0] /= N;
    Rc[1] /= N;
    Rc[2] /= N;

    // find average edge to center distance
    r_c2e = 0;
    r2_c2e = 0;
    for (int i = 0; i < edge_lists[0].size(); i++)
    {
        r_cache = distancefp(edge_lists[0][i], Rc);
        r_c2e += r_cache;
        r2_c2e += r_cache * r_cache;
    }
    r_c2e /= edge_lists[0].size();
    r2_c2e /= edge_lists[0].size();
    r_c2e_std = std::sqrt(r2_c2e - r_c2e * r_c2e);
    del_r = r_c2e / bin_num;

    // find largest eigen value and corresponding eigenvector
    for (int i = 0; i < N; i++)
    {
        bead_list.push_back(i);
    }
    nu = Qevec_m(bead_list);

    for (int k = 0; k < bin_num; k++)
    {
        unu2r.push_back(0);
        Countr.push_back(0);
    }
    for (int i = 0; i < mesh.size(); i++)
    {
        r = distancefp(i, Rc);
        bin = int(r / del_r);
        if (bin < bin_num)
        {
            unu = innerproduct(nu, mesh[i].u);
            unu2r[bin] += unu * unu;
            Countr[bin] += 1;
        }
    }
    // normalize
    for (int b = 0; b < bin_num; b++)
    {
        if (Countr[b])
        {
            unu2r[b] /= Countr[b];
        }
    }
    unu2r.insert(unu2r.begin(), r_c2e_std);
    unu2r.insert(unu2r.begin(), r_c2e);
    return unu2r;
}

std::vector<double> dtmc_lc::nu0nu2l_m(int bin_num)
{
    // bin_num depends on L*d0
    std::vector<double> nunu2l;
    std::vector<double> nu0, nul;
    double del_l, l, nu0nul;
    int bin;
    std::vector<double> fix_R0;
    std::vector<std::vector<int>> bead_lists;
    // assort bead into different list according to position(per slab)
    bead_lists.resize(bin_num);
    fix_R0 = {mesh[fixed_beads[0]].R[0], mesh[fixed_beads[1]].R[0]};
    del_l = (fix_R0[1] - fix_R0[0]) / bin_num;
    // sort beads
    for (int i = 0; i < mesh.size(); i++)
    {
        l = mesh[i].R[0] - fix_R0[0];
        bin = int(l / del_l);
        if (bin < bin_num && l > 0)
        {
            bead_lists[bin].push_back(i);
        }
    }
    // calculate nu for each slab
    nu0 = Qevec_m(bead_lists[0]);
    nunu2l.push_back(1);
    for (int b = 1; b < bin_num; b++)
    {
        nul = Qevec_m(bead_lists[b]);
        nu0nul = innerproduct(nu0, nul);
        nunu2l.push_back(nu0nul * nu0nul);
    }

    return nunu2l;
}

std::vector<double> dtmc_lc::nunu2lcov_m(int bin_num)
{
    // bin_num depends on L*d0
    // correlation function measurement
    // TODO: finish is measurement for extension experiment
    std::vector<double> nunu2lcov;
    std::vector<std::vector<double>> nul;
    double del_l, l, nulnulpt;
    int bin;
    std::vector<double> fix_R0;
    std::vector<std::vector<int>> bead_lists;
    // assort bead into different list according to position(per slab)
    bead_lists.resize(bin_num);
    fix_R0 = {mesh[fixed_beads[0]].R[0], mesh[fixed_beads[1]].R[0]};
    del_l = (fix_R0[1] - fix_R0[0]) / bin_num;
    // sort beads
    for (int i = 0; i < mesh.size(); i++)
    {
        l = mesh[i].R[0] - fix_R0[0];
        bin = int(l / del_l);
        if (bin < bin_num && l > 0)
        {
            bead_lists[bin].push_back(i);
        }
    }
    // calculate nu for each slab
    for (int b = 0; b < bin_num; b++)
    {
        nul.push_back(Qevec_m(bead_lists[b]));
    }
    // calculated correlation function
    nunu2lcov.resize(bin_num);
    for (int db = 0; db < bin_num; db++)
    {
        nunu2lcov[db] = 0;
    }
    for (int b = 0; b < bin_num; b++)
    {
        for (int c = b; c < bin_num; c++)
        {
            nulnulpt = innerproduct(nul[b], nul[c]);
            nunu2lcov[c - b] += nulnulpt * nulnulpt;
        }
    }
    for (int b = 0; b < bin_num; b++)
    {
        nunu2lcov[b] /= (bin_num - b);
    }
    return nunu2lcov;
}

#pragma endregion
