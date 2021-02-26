from cal import *
from plot import *
from numeric_io import *
import time
import sys

def main():
    if(len(sys.argv)>2):
        # on ccv
        folder = "/users/lding3/scratch"
        lams = np.arange(0.0,10.01,0.05)
        Cn= float(sys.argv[1])
        m = int(sys.argv[2])
        q = float(sys.argv[3])
        bn_r,bn_phi=2000,2000
    else:
        # local running
        folder = "../data/pyscratch_local"
        lams = np.arange(0.0,10.1,0.2)
        Cn,m,q0=100,1,3.0
        bn_r,bn_phi=500,1000

    E_lams_cal(folder,lams,Cn,m,q,bn_r,bn_phi,method="SLSQP")
    #Emin_lam_plot(folder,Cn=100,ms=[1],qs=[2.0,2.5,3.0])


if __name__ == "__main__":
    start_time = time.time()
    main()
    #print("--- %s seconds ---" % (time.time() - start_time))