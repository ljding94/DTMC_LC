#include "dtmc_lc.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>

void dtmc_lc::State_write(std::string filename)
{
    std::ofstream f(filename);
    if (f.is_open())
    {
        // find maxim nei.size()
        int max_nei_size = 0;
        for (int i = 0; i < mesh.size(); i++)
        {
            if (mesh[i].nei.size() > max_nei_size)
            {
                max_nei_size = mesh[i].nei.size();
            }
        }
        f << "beta=" << beta << "\n";
        f << "N=" << N << "\n";
        f << "rb=" << rb << "\n";
        f << "l0=" << l0 << "\n";
        f << "max_nei_size=" << max_nei_size << "\n";
        f << "x,y,z,ux,uy,uz,dA,2H,ds,dAK,un2,edge_num,edge_neibs,neibs";
        for (int i = 0; i < mesh.size(); i++)
        {
            f << "\n"
              << mesh[i].R[0] << "," << mesh[i].R[1] << "," << mesh[i].R[2];
            f << "," << mesh[i].u[0] << "," << mesh[i].u[1] << ","
              << mesh[i].u[2];
            f << "," << mesh[i].dAn2H[0] << "," << mesh[i].dAn2H[1] << ","
              << mesh[i].ds << "," << mesh[i].dAK << "," << mesh[i].un2;
            f << "," << mesh[i].edge_num;
            for (int j = 0; j < mesh[i].edge_nei.size(); j++)
            {
                f << "," << mesh[i].edge_nei[j];
            }
            for (int j = 0; j < 2 - mesh[i].edge_nei.size(); j++)
            {
                f << ",-1";
            }
            for (int j = 0; j < mesh[i].nei.size(); j++)
            {
                f << "," << mesh[i].nei[j];
            }
            for (int j = 0; j < max_nei_size - mesh[i].nei.size(); j++)
            {
                f << ",-1";
            }
        }
    }
    f.close();
} // end State_write

void dtmc_lc::State_load(std::string state_file)
{
    // FIXME: this function is indeed broken now
    // need to be re con struct based on new state write code above
    std::ifstream f(state_file);
    std::string buff;
    int max_nei_size;
    while (f.good())
    {
        // skip 8 lines
        for (int i = 0; i < 13; i++)
        {
            std::getline(f, buff);
            // std::cout << buff << "\n";
        }
        // load geometric obsevables
        /*
        std::getline(f, buff, '=');
        std::getline(f, buff);
        Les[0] = std::stof(buff);
        // kinda don't wanna use it any more
        std::getline(f, buff, '=');
        std::getline(f, buff);
        I2H2 = std::stof(buff);
        std::getline(f, buff, '=');
        std::getline(f, buff);
        IK = std::stof(buff);
        std::getline(f, buff, '=');
        std::getline(f, buff);
        */
        max_nei_size = std::stof(buff);
        // skipp another line
        std::getline(f, buff);
        // std::cout << buff << "\n";

        int mcount = 0;

        for (mcount = 0; mcount < mesh.size(); mcount++)
        {
            // load observables
            std::getline(f, buff, ',');
            mesh[mcount].ds = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dAn2H[0] = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dAn2H[1] = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dAK = std::stof(buff);
            // load pos
            for (int k = 0; k < 3; k++)
            {
                std::getline(f, buff, ',');
                mesh[mcount].R[k] = std::stof(buff);
            }
            // load edge_nei
            mesh[mcount].edge_nei.clear();
            for (int k = 0; k < 2; k++)
            {
                std::getline(f, buff, ',');
                if (std::stof(buff) != -1)
                {
                    mesh[mcount].edge_nei.push_back(std::stof(buff));
                }
            }
            // load neibs
            // should read to the end of line, but how?
            mesh[mcount].nei.clear();
            for (int k = 0; k < max_nei_size - 1; k++)
            {
                std::getline(f, buff, ',');
                if (std::stof(buff) != -1)
                {
                    mesh[mcount].nei.push_back(std::stof(buff));
                }
            }
            std::getline(f, buff);
            if (std::stof(buff) != -1)
            {
                mesh[mcount].nei.push_back(std::stof(buff));
            }
        }
    }
    //E = 0;
}

void dtmc_lc::State_write_seq(std::string filename, int MC_sweeps,
                              int step_p_sweep, double delta_s,
                              double delta_theta)
{
    std::string savename;
    for (int sweep_n = 0; sweep_n < MC_sweeps; sweep_n++)
    {
        Thermal(1, step_p_sweep, 1, delta_s, delta_theta);
        savename = filename;
        savename.insert(savename.size() - 4, "_" + std::to_string(sweep_n));
        State_write(savename);
    }
}

void dtmc_lc::Thermal(int MC_sweeps, int step_p_sweep, int beta_steps,
                      double delta_s, double delta_theta)
{
    beta = 0.0;
    for (int n_beta = 0; n_beta < beta_steps; n_beta++)
    {
        beta += 1.0 / beta_steps;
        for (int sweep_n = 0; sweep_n < MC_sweeps / beta_steps; sweep_n++)
        {
            // std::cout << sweep_n << "/" << MC_sweeps << "\n";
            for (int i = 0; i < step_p_sweep; i++)
            {
                bead_metropolis(delta_s);
                spin_metropolis(delta_theta);
                bond_metropolis();
                bond_metropolis();
                if (i % int(std::sqrt(N)) == 0)
                {
                    edge_metropolis();
                }
            }
            // std::cout << "thermo, beta=" << beta << "," << sweep_n << "/"<<
            // MC_sweeps << "\n";
        }
    }
}
void dtmc_lc::O_MC_measure(int MC_sweeps, int sweep_p_G, int step_p_sweep,
                           double delta_s, double delta_theta,
                           std::string folder, std::string finfo,
                           int bin_num_r, int bin_num_un2)
{
    std::vector<double> E_all;
    std::vector<double> I2H2_all;
    std::vector<std::vector<double>> Les_all;
    Les_all.resize(Ne);
    std::vector<double> Tp2uu_all;
    std::vector<double> Tuuc_all;
    std::vector<double> Tun2_all;
    std::vector<double> IKun2_all;
    std::vector<double> IdA_all;
    std::vector<double> I2H_all;
    std::vector<double> IK_all;
    std::vector<int> Bond_num_all;

    std::vector<std::vector<double>> Gij_all;
    std::vector<std::vector<double>> un2r_all;
    std::vector<std::vector<double>> un2dis_all;
    std::vector<std::vector<double>> unu2r_all;
    std::vector<std::vector<double>> nu0nu2l_all;
    std::vector<std::vector<double>> nunu2lcov_all;
    double bead_accept = 0;
    double spin_accept = 0;
    double bond_accept = 0;
    double edge_accept = 0;

    bool fix_bead;
    if (fixed_beads.size() == 2)
    {
        fix_bead = 1;
    }

    std::clock_t c_start = std::clock();
    for (int sweep_n = 0; sweep_n < MC_sweeps; sweep_n++)
    {
        std::cout << sweep_n << "/" << MC_sweeps << "\n";
        for (int i = 0; i < step_p_sweep; i++)
        {
            bead_accept += bead_metropolis(delta_s);
            spin_accept += spin_metropolis(delta_theta);
            bond_accept += bond_metropolis();
            bond_accept += bond_metropolis();
            if (i % int(std::sqrt(N)) == 0)
            {
                edge_accept += edge_metropolis();
            }
        }
        E_all.push_back(Ob_sys.E);
        I2H2_all.push_back(Ob_sys.I2H2);
        for (int e = 0; e < Ne; e++)
        {
            Les_all[e].push_back(Ob_sys.Les[e]);
        }
        Tp2uu_all.push_back(Ob_sys.Tp2uu);
        Tuuc_all.push_back(Ob_sys.Tuuc);
        Tun2_all.push_back(Ob_sys.Tun2);
        IKun2_all.push_back(Ob_sys.IKun2);
        IdA_all.push_back(Ob_sys.IdA);
        I2H_all.push_back(Ob_sys.I2H);
        IK_all.push_back(Ob_sys.IK);
        Bond_num_all.push_back(Ob_sys.Bond_num);

        if (sweep_n % sweep_p_G == 0)
        {
            if (fix_bead)
            {
                nu0nu2l_all.push_back(nu0nu2l_m(bin_num_r));
                nunu2lcov_all.push_back(nunu2lcov_m(bin_num_r));
            }
            else
            {
                if (bin_num_un2 != 0)
                {
                    un2r_all.push_back(un2r_m(bin_num_r));
                }
                if (bin_num_un2 != 0)
                {
                    un2dis_all.push_back(un2dis_m(bin_num_un2));
                }
                // don't waste time doing the measurement if it's not to be used
                //unu2r_all.push_back(unu2r_m(bin_num_r));
            }
        }
    }

    bead_accept /= MC_sweeps * step_p_sweep;
    spin_accept /= MC_sweeps * step_p_sweep;
    bond_accept /= MC_sweeps * step_p_sweep * 2;
    edge_accept /= MC_sweeps * (step_p_sweep / int(std::sqrt(N)));

    std::clock_t c_end = std::clock();
    std::ofstream f(folder + "/O_MC_" + finfo + ".txt");
    if (f.is_open())
    {
        f << "beta=" << beta << "\n";
        f << "N=" << N << "\n";
        f << "Ne=" << Ne << "\n";
        f << "L=" << L << "\n";
        f << "l0=" << l0 << "\n";
        f << "step_p_sweep=" << step_p_sweep << "\n";
        f << "MC_sweeps=" << MC_sweeps << "\n";
        f << "CPUtime=" << double(c_end - c_start) * 1.0 / CLOCKS_PER_SEC
          << "\n";
        f << "bead_accept," << bead_accept << "\n";
        f << "spin_accept," << spin_accept << "\n";
        f << "bond_accept," << bond_accept << "\n";
        f << "edge_accept," << edge_accept << "\n";
        f << "E,";
        for (int e = 0; e < Ne; e++)
        {
            f << "Les[" << e << "],";
        }
        f << "IdA,I2H,I2H2,IK,Tp2uu,Tuuc,Bond_num,Tun2,IKun2\n";
        for (int i = 0; i < E_all.size(); i++)
        {
            f << E_all[i] << ",";
            for (int e = 0; e < Ne; e++)
            {
                f << Les_all[e][i] << ",";
            }
            f << IdA_all[i] << "," << I2H_all[i] << "," << I2H2_all[i] << ","
              << IK_all[i] << "," << Tp2uu_all[i] << "," << Tuuc_all[i] << ","
              << Bond_num_all[i] << "," << Tun2_all[i] << "," << IKun2_all[i]
              << "\n";
        }
    }
    f.close();

    if (fix_bead)
    {
        std::ofstream f_nu0nu2l(folder + "/nu0nu2l_MC_" + finfo + ".txt");
        if (f_nu0nu2l.is_open())
        {
            f_nu0nu2l << "bin_num=" << bin_num_r << "\n";
            f_nu0nu2l << "(n_u(0)*n_u(l))^2\n";
            for (int i = 0; i < nu0nu2l_all.size(); i++)
            {
                f_nu0nu2l << nu0nu2l_all[i][0];
                for (int j = 1; j < nu0nu2l_all[i].size(); j++)
                {
                    f_nu0nu2l << "," << nu0nu2l_all[i][j];
                }
                f_nu0nu2l << "\n";
            }
        }
        f_nu0nu2l.close();

        std::ofstream f_nunu2lcov(folder + "/nunu2lcov_MC_" + finfo + ".txt");
        if (f_nunu2lcov.is_open())
        {
            f_nunu2lcov << "bin_num=" << bin_num_r << "\n";
            f_nunu2lcov << "ave_s{(n_u(s)*n_u(s+l))^2}\n";
            for (int i = 0; i < nunu2lcov_all.size(); i++)
            {
                f_nunu2lcov << nunu2lcov_all[i][0];
                for (int j = 1; j < nunu2lcov_all[i].size(); j++)
                {
                    f_nunu2lcov << "," << nunu2lcov_all[i][j];
                }
                f_nunu2lcov << "\n";
            }
        }
        f_nunu2lcov.close();
    }
    else
    {
        if (bin_num_r != 0)
        {
            std::ofstream f_un2r(folder + "/un2r_MC_" + finfo + ".txt");
            // std::ofstream f_unu2r(folder + "/unu2r_MC_" + finfo + ".txt");
            if (f_un2r.is_open())
            {
                f_un2r << "bin_num=" << bin_num_r << "\n";
                f_un2r << "r, r_std,(u*n_u)^2(r)\n";
                for (int i = 0; i < un2r_all.size(); i++)
                {
                    f_un2r << un2r_all[i][0];
                    for (int j = 1; j < un2r_all[i].size(); j++)
                    {
                        f_un2r << "," << un2r_all[i][j];
                    }
                    f_un2r << "\n";
                }
            }
            f_un2r.close();
        }
        if (bin_num_un2 != 0)
        {
            std::ofstream f_un2dis(folder + "/un2dis_MC_" + finfo + ".txt");
            if (f_un2dis.is_open())
            {
                f_un2dis << "bin_num=" << bin_num_un2 << "\n";
                f_un2dis << "un2 density/ un2*bin_num_un2\n";
                for (int i = 0; i < un2dis_all.size(); i++)
                {
                    f_un2dis << un2dis_all[i][0];
                    for (int j = 1; j < un2dis_all[i].size(); j++)
                    {
                        f_un2dis << "," << un2dis_all[i][j];
                    }
                    f_un2dis << "\n";
                }
                f_un2dis.close();
            }
        }
        /*
        if (f_unu2r.is_open()) {
            f_unu2r << "bin_num=" << bin_num_r << "\n";
            f_unu2r << "r,r_std,(u*n_u)^2(r)\n";
            for (int i = 0; i < unu2r_all.size(); i++) {
                f_unu2r << unu2r_all[i][0];
                for (int j = 1; j < unu2r_all[i].size(); j++) {
                    f_unu2r << "," << unu2r_all[i][j];
                }
                f_unu2r << "\n";
            }
        }
        f_unu2r.close();
        */
    }
}
