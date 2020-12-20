#!/usr/local/bin/python3
import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *


def main():
    print("hello! dtmc_lc analysis")
    #config_plot_xyz("../data/scratch_local/State_N50_Ne1_L0_kar5_lam3.0_Kd2.0_q0.0_Cn2.0_kargd2.0.txt",mesh=1)
    config_plot3D("../data/scratch_local/State_N400_Ne1_L-1_kar15_lam6.0_Kd0.0_q0.0_Cn0.0_kargd0.0.txt")
    config_plot3D("../data/scratch_local/State_N400_Ne1_L0_kar15_lam6.0_Kd0.0_q0.0_Cn0.0_kargd0.0.txt")
    config_plot3D("../data/scratch_local/State_N400_Ne1_L30_kar15_lam6.0_Kd0.0_q0.0_Cn0.0_kargd0.0.txt")
    config_plot3D("../data/scratch_local/State_N400_Ne2_L0_kar15_lam6.0_Kd0.0_q0.0_Cn0.0_kargd0.0.txt")

    #config_plot3D("../data/scratch_local/State_N0_rd10_Ne1_L0_kar10_lam5.0_Kd0.0_q0.0_Cn0.0_kargd0.0.txt")
    #config_plot3D("../data/Ne1/Dec8_2020/State_N400_Ne1_L0_kar15_lam5.0_Kd4.0_Kt6.0_Cn5.0_kargd1.0.txt")
    #config_plot3D("../data/Ne1/Dec8_2020/State_N400_Ne1_L0_kar15_lam5.0_Kd4.0_Kt6.0_Cn5.0_kargd5.0.txt")
    #config_plot3D("../data/Ne1/Dec17_2020/State_N200_Ne1_L0_kar10_lam4.0_Kd2.0_q0.3_Cn2.0_kargd3.0 .txt")
    return 0
    foldername = "../data/Ne1/Dec19_2020"
    print("analyzing "+foldername)
    pars = []
    colors = []
    alphas = []
    Ns = [100,200,300,400]
    N = 400
    Ne=1
    Nes = [1,2]
    #Ls = [13,14,15, 16, 17,18, 19, 20, 21, 22, 23, 24, 25, 26,27, 28, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40,41,42,43,44,45,46,47,48,49,50]
    #Ls=[20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60]
    Ls=np.linspace(100,200,21)
    L=0
    kars = [1,5,10,15]
    kar = 15
    lams = np.linspace(4.0,50.0,93)
    lam=6.0
    Kds = np.arange(4.0,20.1,2.0)
    Kd=3.0
    qs=np.arange(0.2,2.1,0.2)
    q=0.2
    Cns=np.arange(3.0,20.1,2.0)
    Cn=3.0
    kargds=np.arange(1.0,10.1,1.0)
    kargd=0.0
    Ne1pars = []
    Ne2pars = []
    for q in qs:
        pars.append([N, Ne, L, kar,lam, Kd, q, Cn, kargds])
    par_nm = ["N", "Ne", "L", "kar","lam","Kd","q", "Cn", "kargd"]
    par_dg = [0,0,0,0,1,1,1,1,1] # number of digit for each
    mod="kargd"
    for i in range(len(pars)):
        N, Ne, L, kar,  lam, Kd, q,Cn,kargd = pars[i]
        print(pars[i])
        O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, tau_c=6)
        if(1):
            pass
            #twistl_stat_plot(foldername,pars[i],par_nm,par_dg,mode=mod,d0=1.5,head="nunu2lcov",tag=r"$K_d=%.1f$"%Kd,leg_num=4,bin_num=20)
            #twistl_stat_plot(foldername,pars[i],par_nm,par_dg,mode=mod,d0=1.5,head="nu0nu2l",tag=r"$K_d=%.1f$"%Kd,leg_num=4,bin_num=20)
            pass
            #twistr_stat_plot(foldername,pars[i],par_nm,par_dg,mode=mod,head="un2r",tag=r"$Kd=%.1f$"%Kd,leg_num=4,bin_num=80)

        for kargd in kargds[::2]:
            if(i%2==0):
                filename = foldername + "/State_N%.0f_Ne%.0f_L%.0f_kar%.0f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_kargd%.1f.txt" % (N, Ne, L, kar,lam,Kd,q,Cn, kargd)
                config_plot_xyz(filename, tag=r"$K_d=%.1f,q=%.1f,C_n=%.1f,\kappa_d=%.1f$" % (Kd,q,Cn,kargd),Format="png")

        #print("sleeping...")
        #time.sleep(10)
    colors = None
    alphas = None


    #lamp_pars_plot(foldername, pars,par_nm,par_dg,mode=mod,head="un2r")
    Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)

    # additional test



if __name__ == '__main__':
    main()
