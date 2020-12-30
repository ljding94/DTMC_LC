#!/usr/local/bin/python3
import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *


def main():
    print("hello! dtmc_lc analysis")
    #config_plot_xyz("../data/scratch_local/State_N1200_Ne1_L-1_kar10_lam6.0_Kd0.0_q0.0_Cn0.0_kargd0.0.txt",mesh=1,rod=1)
    #config_plot3D("../data/Ne1/Dec27_2020/State_N1000_Ne1_L-1_kar15_lam5.0_Kd3.0_q2.0_Cn3.0_kargd0.0.txt",mesh=0,rod=1)
    #config_plot3D("../data/scratch_local/State_N200_Ne1_L-1_kar100_karg20_lam3.0_Kd0.0_q0.0_Cn0.0_kard0.0_init.txt",mesh=0,rod=1)
    #config_plot3D("../data/scratch_local/State_N200_Ne1_L-1_kar100_karg20_lam4.0_Kd0.0_q0.0_Cn0.0_kard0.0.txt",mesh=1,rod=0)
    #config_plot3D("../data/scratch_local/State_N200_Ne1_L-1_kar3_lam5.0_Kd2.0_q2.0_Cn2.0_kargd0.0.txt",mesh=1,rod=1)

    #return 0
    foldername = "../data/Ne1/Dec29_2020_1"
    print("analyzing "+foldername)
    pars = []
    colors = []
    alphas = []
    Ns = [100,200,300,400]
    N = 1000
    Ne=1
    Nes = [1,2]
    #Ls = [13,14,15, 16, 17,18, 19, 20, 21, 22, 23, 24, 25, 26,27, 28, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40,41,42,43,44,45,46,47,48,49,50]
    #Ls=[20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60]
    Ls=np.linspace(100,200,21)
    L=-1
    kars = [1,5,10,15]
    kar = 15
    lams = np.arange(2.0,11.1,0.5)
    lam=6.0
    Kds = np.arange(4.0,40.1,4.0)
    Kd=0.0
    qs=np.arange(2.2,4.1,0.2)
    q=0.8
    Cns=np.arange(4.0,40.1,4.0)
    Cn=16.0
    kargds=np.arange(1.0,5.1,1.0)
    kargd=0.0
    Ne1pars = []
    Ne2pars = []
    for Kd in Kds:
        pars.append([N, Ne, L, kar, lam, Kd, q, Cns, kargd])
    par_nm = ["N", "Ne", "L", "kar","lam","Kd","q", "Cn", "kargd"]
    par_dg = [0,0,0,0,1,1,1,1,1] # number of digit for each
    mod="Cn"
    for i in range(len(pars)):
        print("hallo")
        N, Ne, L, kar,  lam, Kd, q,Cn,kargd = pars[i]
        O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, tau_c=6)
        if(1):
            pass
            #twistl_stat_plot(foldername,pars[i],par_nm,par_dg,mode=mod,d0=1.5,head="nunu2lcov",tag=r"$K_d=%.1f$"%Kd,leg_num=4,bin_num=20)
            #twistl_stat_plot(foldername,pars[i],par_nm,par_dg,mode=mod,d0=1.5,head="nu0nu2l",tag=r"$K_d=%.1f$"%Kd,leg_num=4,bin_num=20)
            #pass
            twistr_stat_plot(foldername,pars[i],par_nm,par_dg,mode=mod,head="un2r",tag=r"$Kd=%.1f$"%Kd,leg_num=1,bin_num=60)

        for Cn in Cns[::1]:
            if(i%1==0):
                filename = foldername + "/State_N%.0f_Ne%.0f_L%.0f_kar%.0f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_kargd%.1f.txt" % (N, Ne, L, kar,lam,Kd,q,Cn, kargd)
                #config_plot_xyz(filename, tag=r"$K_d=%.1f,q=%.1f,C_n=%.1f$" % (Kd,q,Cn),Format="png")

        #print("sleeping...")
        #time.sleep(10)
    colors = None
    alphas = None

    lamp_pars_plot(foldername, pars,par_nm,par_dg,mode=mod,head="un2r")
    #Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)

    # additional test

if __name__ == '__main__':
    main()