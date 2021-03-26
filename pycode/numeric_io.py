import numpy as np
from cal import *


def DelEmin_par_run(par,bn_r,bn_phi,method):

    lam,Cn,m,q0 = par
    #print(method)
    opt = opt_r1_qtheta_fun(obj_DelE,lam,Cn,m,q0,bn_r,bn_phi,method)

    savename = "DelEmin_lam%.1f_Cn%.0f_m%d_q%.1f.txt"%(lam,Cn,m,q0)
    with open (savename,"w") as f:
        f.write("DelEmin,r1min,qthetamin\n")
        f.write("%f,%f,%f"%(opt.fun,opt.x[0],opt.x[1]))

def E_lams_cal(folder,lams,Cn,m,q0,bn_r,bn_phi,method):
    Emins, r1mins, qthetamins = [],[],[]
    for i in range(len(lams)):
        lam=lams[i]
        print(lam,"/",lams[-1])
        opt = opt_r1_qtheta_fun(obj_E_ENP,lam,Cn,m,q0,bn_r,bn_phi,method)
        Emins.append(opt.fun)
        r1mins.append(opt.x[0])
        qthetamins.append(opt.x[1])

    savename = folder+"/Emin_lams_Cn%.0f_m%d_q%.2f.txt"%(Cn,m,q0)
    with open (savename,"w") as f:
        f.write("lam,Emin,r1min,qthetamin\n")
        for i in range(len(Emins)):
            f.write("%f,%f,%f,%f\n"%(lams[i],Emins[i],r1mins[i],qthetamins[i]))

