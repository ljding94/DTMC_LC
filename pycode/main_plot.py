#!/usr/local/bin/python3
import numpy as np
from cal import *
from plot import *
from numeric_io import *
import time
import sys

def main():
    #xyznENP_plot("../data/pyscratch_local",1,1,np.pi,0)
    # local running
    #folder = "../data/pydata/Feb27_2021"
    #folder = "../data/pydata/Mar8_2021"
    folder = "../data/pydata/Mar13_2021"
    Cn=100
    ms=[1]
    qs=np.arange(1.50,4.51,0.05)
    lams = np.arange(0.00,5.01,0.05)
    qd=2
    Emin_r1op_qthetaop_m1_pcolormesh(folder,Cn,qs,qd,lams)
    #qs=np.arange(2.0,3.21,0.4)
    #qs = [2.5]
    #Emin_lam_plot(folder,Cn,ms,qs)
    #qs=np.arange(2.000,2.501,0.005)
    #lams = np.arange(2.000,2.501,0.005)
    #phase_of_Emin_pcolormesh(folder,Cn,ms,qs,qd,lams)
    #qs=np.arange(1.00,5.01,0.5)
    #Emin_lam_plot(folder,Cn,ms,qs)

    #LineWidth, FontSize, LabelSize = 1,9,8
    #print("(施工中) trying to plot publication-quality figures")
    #phase_E_fill_between_plot(folder, Cn,ms,qs,lams, LineWidth, FontSize, LabelSize)


if __name__ == "__main__":
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))