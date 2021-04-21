#!/usr/local/bin/python3
import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *


def main():
    print("hello! dtmc_lc analysis")
    #config_plot3D("../data/Ne1/Mar23_2021/State_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd7.0_q2.0_Cn7.0_kard0.0.txt",mesh=1,rod=0)
    #config_plot3D("../data/Ne2/Apr7_2021/State_N200_Ne2_L0_kar100_karg0.0_lam6.0_Kd7.0_q1.0_Cn7.0_kard0.0.txt",mesh=1,rod=0)
    #config_nu2_dis("../data/Ne1/Jan14_2021/State_N1000_Ne1_L-1_kar15_karg0.0_lam5.0_Kd5.0_q1.0_Cn5.0_kard0.0.txt")
    #return 0


    #foldername = "../data/Ne1/Apr4_2021"
    foldername = "../data/Ne2/Apr12_2021" # kinda switch to Ne2 simulation for next paper
    #foldername = "../data/Ne1/Mar23_2021" # patching up simulation for softmatter paper
    print("analyzing "+foldername)
    pars = []
    colors = []
    alphas = []
    Ns = [100,200,300,400]
    N = 400
    Ne=2
    Nes = [1,2]
    #Ls = [13,14,15, 16, 17,18, 19, 20, 21, 22, 23, 24, 25, 26,27, 28, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40,41,42,43,44,45,46,47,48,49,50]
    #Ls=[20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60]
    Ls=np.arange(20,70.1,2)
    L=0 # -1 for circular initialization
    kars = np.arange(20,100.1,20)
    kar = 100
    kargs = np.arange(0.0,10.1,1.0)
    karg=0.0
    #lams = np.arange(4.0,14.1,0.5)
    lams = np.arange(4.0,10.1,0.5)
    lam=6.0
    Kds = np.arange(3.0,9.1,2.0)
    #Kds = [0.0,0.4,0.8,1.2,1.6,2.0,3.0,4.0,5.0,6.0]
    Kd=7.0
    qs=np.arange(0.0,2.01,0.1)
    #qs = [1.5,2.0,2.5,3.0,3.5,4.0]
    #qs = [0.0,0.4,0.8,1.0,1.2,1.4,1.6,2.0,2.4,2.8]
    q=6.0
    #Cns = np.arange(3.0,15.1,1.0)
    Cns = np.arange(3.0,9.1,2.0)
    #Cns = np.arange(2.0,10.1,2.0)
    Cn=5.0
    kards=np.arange(3.0,9.1,2.0)
    kard=0.0
    Ne1pars = []
    Ne2pars = []
    for Cn in Cns:
        #Kd=Cn # special case when two parameters are set to the same
        pars.append([N, Ne, L, kar, karg, lam, Kd, qs, Cn, kard])
    par_nm = ["N", "Ne", "L", "kar","karg","lam","Kd","q", "Cn", "kard"]
    par_dg = [0,0,0,0,1,1,1,1,1,1] # nsumber of digit for each
    mod="q"
    for i in range(len(pars)):
        print("analyzing",pars[i])
        N, Ne, L, kar, karg, lam, Kd, q,Cn,kard = pars[i]
        O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, tau_c=10)
        if(1):
            pass
            #tiltr_stat_plot(foldername,pars[i],par_nm,par_dg,mode=mod,head="un2r",tag=r"$Kd=%.1f$"%Kd,leg_num=1)
            #un2dis_stat_plot(foldername,pars[i],par_nm,par_dg,mode=mod,head="un2dis",tag=r"$K_d=%.0f$"%(Kd),leg_num=5)
            #uucr_stat_ana(foldername, pars[i], par_nm, par_dg, mode=mod, head="uucr", dr=0.2, tag="",leg_num=5)

        for q in qs[::2]:
            if(i%1==0):
                pass
                filename = foldername + "/State_N%.0f_Ne%.0f_L%.0f_kar%.0f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_kard%.1f.txt" % (N, Ne, L, kar,karg,lam,Kd,q,Cn, kard)
                config_plot_xyz(filename,mesh=0,rod=1,tag=r"$K_d=%.1f,q=%.1f,C_n=%.1f$" % (Kd,q,Cn),Format="png")
                #config_plot_xyz(filename, mesh=0,rod=0,cvt_map="Mean",cmap_smooth=3,tag=r"$\bar{\kappa}=%.1f,\lambda=%.1f,C_n=%.1f$" % (karg,lam,Cn),Format="png")

                #config_plot_xyz(filename, mesh=0,rod=0,cvt_map="Gaussian",tag=r"$\bar{\kappa}=%.1f,\lambda=%.1f,C_n=%.1f$" % (karg,lam,Cn),Format="png")

        #print("sleeping...")
        #time.sleep(10)
    colors = None
    alphas = None

    #lamp_pars_plot(foldername, pars,par_nm,par_dg,mode=mod,head="un2r")
    #un2below_pars_plot(foldername,pars,par_nm,par_dg,mode=mod,head="un2dis")
    Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)

    # additional test

if __name__ == '__main__':
    main()