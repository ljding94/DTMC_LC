#!/usr/local/bin/python3
from cal import *
from plot import *
from numeric_io import *
import time
import sys

def main():

    # local running
    folder = "../data/pyscratch_local"
    par = [1.0,100,2,2.0]
    lams = np.arange(0.0,10.1,0.2)
    Cn=100
    ms=[1]
    qs=[2.0,2.5,3.0]

    Emin_lam_plot(folder,Cn,ms,qs)


if __name__ == "__main__":
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))