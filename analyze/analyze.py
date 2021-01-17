import numpy as np
import matplotlib.pyplot as plt
from plot import *
from autocorrelation import *
from scipy.optimize import curve_fit
from scipy import odr

def find_cpar_ind(par_nm,mode):
    cpar_ind = -1
    for i in range(len(par_nm)):
        if par_nm[i]==mode:
            cpar_ind=i
            break
    return cpar_ind

def O_stat_ana(foldername,par,par_nm,par_dg, mode, tau_c=6):
    E_ave, E_tau, E_err = [], [], []
    Ne = par[find_cpar_ind(par_nm,"Ne")]
    Les_ave, Les_tau, Les_err = [[] for i in range(Ne)], [[] for i in range(Ne)], [[] for i in range(Ne)]
    IdA_ave, IdA_tau, IdA_err = [], [], []
    I2H_ave, I2H_tau, I2H_err = [], [], []
    I2H2_ave, I2H2_tau, I2H2_err = [], [], []
    IK_ave, IK_tau, IK_err = [], [], []
    Tp2uu_ave, Tp2uu_tau, Tp2uu_err = [], [], []
    Tuuc_ave, Tuuc_tau, Tuuc_err = [], [], []
    Tun2_ave, Tun2_tau, Tun2_err = [], [], []
    IKun2_ave, IKun2_tau, IKun2_err = [], [], []
    p2uu_ave, p2uu_tau, p2uu_err = [], [], []
    uuc_ave, uuc_tau, uuc_err = [], [], []
    un2_ave,un2_tau,un2_err = [],[],[]
    if(Ne==2):
        Ledif_ave,Ledif_tau,Ledif_err=[],[],[]
    cpar_ind = find_cpar_ind(par_nm,mode)
    cpar = par[cpar_ind]
    for i in range(len(cpar)):
        par_dealing = par[:]
        par_dealing[cpar_ind] = par[cpar_ind][i]
        f2rtail = "MC"
        for j in range(len(par_dealing)):
            f2rtail+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing[j])
        f2rtail+=".txt"
        file2read = foldername + "/O_"+f2rtail
        data = np.loadtxt(file2read, skiprows=13, delimiter=",", unpack=True)
        E = data[0]
        Les = data[1:1+Ne]
        IdA,I2H,I2H2,IK,Tp2uu,Tuuc,Bond_num,Tun2,IKun2 = data[1+Ne:]
        p2uu = Tp2uu/Bond_num
        uuc = Tuuc/Bond_num
        # N in file name is not real N
        N =par[find_cpar_ind(par_nm,"N")]
        Ndict={200:187,400:367,800:721,1000:823,1600:1459}
        if(par[find_cpar_ind(par_nm,"L")]==-1):
            N=Ndict[N]
        un2=Tun2/N
        # Ne2 case, need Ledif for additional info
        if(Ne==2):
            Ledif = np.abs(Les[0]-Les[1])
            Ledif_ave.append(np.average(Ledif))
            rho, cov0 = autocorrelation_function_fft(Ledif)
            tau, tau_err = tau_int_cal_rho(rho,tau_c)
            Ledif_tau.append(tau)
            Ledif_err.append(np.sqrt(2 * tau / len(Ledif) * cov0))


        print("energy slicing E",E)
        print("Les[0]",Les[0])
        Et = E-0.5*par[find_cpar_ind(par_nm,"kar")]*I2H2
        print("E-0.5kar*I2H2",Et)
        Et = Et-par[find_cpar_ind(par_nm,"lam")]*Les[0]
        print("E-0.5kar*I2H2-lam*L",Et)
        Et=Et+par[find_cpar_ind(par_nm,"Kd")]*Tp2uu
        print("E-0.5kar*I2H2-lam*L+Kd*Tp2uu",Et)
        Et=Et+par[find_cpar_ind(par_nm,"Kd")]*par[find_cpar_ind(par_nm,"q")][i]*Tuuc
        print("E-0.5kar*I2H2-lam*L+Kd*Tp2uu+Kd*q*Tuuc",Et)
        Et=Et+0.5*par[find_cpar_ind(par_nm,"Cn")]*(Tun2-N)
        print("E-0.5kar*I2H2-lam*L+Kd*Tp2uu+Kd*q*Tuuc+0.5Cn*(Tun2-N)",Et)


        # E
        E_ave.append(np.average(E))
        rho, cov0 = autocorrelation_function_fft(E)
        #print("cov0-var(E)",cov0-np.var(E))
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoE.pdf")
        E_tau.append(tau)
        E_err.append(np.sqrt(2 * tau / len(E) * cov0))

        # Le
        for e in range(Ne):
            Les_ave[e].append(np.average(Les[e]))
            rho, cov0 = autocorrelation_function_fft(Les[e])
            tau, tau_err = tau_int_cal_rho(rho,tau_c)
            # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoLe.pdf")
            Les_tau[e].append(tau)
            Les_err[e].append(np.sqrt(2 * tau / len(Les[e]) * cov0))

        # IdA
        IdA_ave.append(np.average(IdA))
        rho, cov0 = autocorrelation_function_fft(IdA)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        IdA_tau.append(tau)
        IdA_err.append(np.sqrt(2 * tau / len(IdA) * cov0))

        # I2H
        I2H_ave.append(np.average(I2H))
        rho, cov0 = autocorrelation_function_fft(I2H)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        I2H_tau.append(tau)
        I2H_err.append(np.sqrt(2 * tau / len(I2H) * cov0))

        # I2H2
        I2H2_ave.append(np.average(I2H2))
        rho, cov0 = autocorrelation_function_fft(I2H2)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoI2H2.pdf")
        I2H2_tau.append(tau)
        I2H2_err.append(np.sqrt(2 * tau / len(I2H2) * cov0))

        # Ikg
        IK_ave.append(np.average(IK))
        rho, cov0 = autocorrelation_function_fft(IK)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoIK.pdf")
        IK_tau.append(tau)
        IK_err.append(np.sqrt(2 * tau / len(IK) * cov0))

        # p2uu
        p2uu_ave.append(np.average(p2uu))
        rho, cov0 = autocorrelation_function_fft(p2uu)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autop2uu.pdf")
        p2uu_tau.append(tau)
        p2uu_err.append(np.sqrt(2 * tau / len(p2uu) * cov0))

        # uuc
        uuc_ave.append(np.average(uuc))
        rho, cov0 = autocorrelation_function_fft(uuc)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autouuc.pdf")
        uuc_tau.append(tau)
        uuc_err.append(np.sqrt(2 * tau / len(uuc) * cov0))

        # un2
        un2_ave.append(np.average(un2))
        rho, cov0 = autocorrelation_function_fft(un2)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoun2.pdf")
        un2_tau.append(tau)
        un2_err.append(np.sqrt(2 * tau / len(un2) * cov0))

        # Tun2
        '''
        Tun2_ave.append(np.average(Tun2))
        rho, cov0 = autocorrelation_function_fft(Tun2)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoTun2.pdf")
        Tun2_tau.append(tau)
        Tun2_err.append(np.sqrt(2 * tau / len(Tun2) * cov0))
        '''

        # IKun2
        IKun2_ave.append(np.average(IKun2))
        rho, cov0 = autocorrelation_function_fft(IKun2)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoIKun2.pdf")
        IKun2_tau.append(tau)
        IKun2_err.append(np.sqrt(2 * tau / len(IKun2) * cov0))

    # only changed "lam" and "B" mode here, others waiting for further decision
    # generalize using par_nm list
    f2stail = "MC"
    for j in range(len(par)):
        if(j==cpar_ind):
            f2stail+="_"+par_nm[j]+"s"
        else:
            f2stail+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing[j])
    f2stail+="_ana.txt"
    # only changed "lam" and "B" mode here, others waiting for further decision
    savefile = foldername + "/O_" + f2stail

    with open(savefile, "w") as f:
        f.write(mode+",E_ave,E_tau,E_err")
        for e in range(Ne):
            f.write(",Les_ave[%d],Les_tau[%d],Les_err[%d]"%(e,e,e))
        f.write(",IdA_ave,IdA_tau,IdA_err,I2H_ave,I2H_tau,I2H_err,I2H2_ave,I2H2_tau,I2H2_err,IK_ave,IK_tau,IK_err,p2uu_ave,p2uu_tau,p2uu_err,uuc_ave,uuc_tau,uuc_err,un2_ave,un2_tau,un2_err,IKun2_ave,IKun2_tau,IKun2_err")
        if(Ne==2):
            f.write(",Ledif_ave,Ledif_tau,Ledif_err")
        f.write("\n")
        for i in range(len(cpar)):
            f.write("%f,%f,%f,%f" % (cpar[i], E_ave[i], E_tau[i], E_err[i]))
            for e in range(Ne):
                f.write(",%f,%f,%f"%(Les_ave[e][i],Les_tau[e][i], Les_err[e][i]))
            f.write(",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f"%(IdA_ave[i], IdA_tau[i], IdA_err[i],I2H_ave[i], I2H_tau[i], I2H_err[i],I2H2_ave[i], I2H2_tau[i], I2H2_err[i], IK_ave[i], IK_tau[i], IK_err[i], p2uu_ave[i], p2uu_tau[i], p2uu_err[i], uuc_ave[i], uuc_tau[i], uuc_err[i], un2_ave[i], un2_tau[i], un2_err[i],IKun2_ave[i],IKun2_tau[i],IKun2_err[i]))
            if(Ne==2):
                f.write(",%f,%f,%f"%(Ledif_ave[i], Ledif_tau[i], Ledif_err[i]))
            f.write("\n")

def data_eigen_Iij(filename):
    Iijs = np.loadtxt(filename,delimiter=",",skiprows=1)
    Ieigs = []
    for Iij in Iijs:
        #print("Iij",Iij)
        #print(np.reshape(Iij,(3,3)))
        w,v = np.linalg.eig(np.reshape(Iij,(3,3)))
        Ieigs.append(np.sort(w))
    Ieigs = np.transpose(Ieigs)
    return Ieigs


def O_stat_ana_Ls(foldername, L, kar, lam):
    En_ave, En_std = [], []
    L_ave, L_std = [], []
    for i in range(len(L)):
        config_plot_xyz(foldername + "/State_L" +
                        str(L[i]) + "_kar" + str(kar) + "_lam" + str(lam) + ".txt")
        print("dealing with ", L[i], kar, lam)
        file2read = foldername + "/O_MC_L" + \
            str(L[i]) + "_kar" + str(kar) + "_lam" + str(lam) + ".txt"
        En, L_e = np.loadtxt(
            file2read, skiprows=7, delimiter=",", unpack=True)
        i2H = (En - lam * L_e) / kar
        O_kar_lam_MCstep_plot(L_e, i2H, En, file2read[:-4] + ".pdf")
        L_ave.append(np.average(L_e))
        L_std.append(np.std(L_e))
        En_ave.append(np.average(En))
        En_std.append(np.std(En))
    savefile = foldername + "/O_MC_kar" + \
        str(kar) + "_lam" + str(lam) + "ana.txt"
    with open(savefile, "w") as f:
        f.write("L,En_ave,En_std,L_ave,L_std\n")
        for i in range(len(L)):
            f.write("%f,%f,%f,%f,%f\n" % (
                L[i], En_ave[i], En_std[i], L_ave[i], L_std[i]))

def tan_fit(x,tan0, lam_p,tanc):
    return tan0*np.exp((x-1)/lam_p)+tanc

def sqrt_fit(x,a):
    return a/np.sqrt(x)
def exp_fit(x,a,b):
    return a*np.power(x,b)

def odr_tan_fit(p,x):
    return p[0]*np.exp((x-1)/p[1])+p[2]

def odr_lamp_popt(xdata,ydata,xerr,yerr):
    tan = odr.Model(odr_tan_fit)
    data = odr.RealData(xdata,ydata,sx=xerr,sy=yerr)
    odr_fit = odr.ODR(data,tan,beta0=[0.5,0.2,0.1])
    output=odr_fit.run()
    #print("output.beta",output.beta)
    #print("output.sd_beta",output.sd_beta)
    #output.pprint()
    return (output.beta,output.sd_beta)

def twistr_stat_plot(foldername, par, par_nm, par_dg, mode,head="un2r",tag="",leg_num=5,bin_num=40):
    Ne = par[find_cpar_ind(par_nm,"Ne")]
    cpar_ind = find_cpar_ind(par_nm,mode)
    cpar = par[cpar_ind]
    unu2r_all,unu2rerr_all = [],[]
    r_ave, r_err = [],[]
    r_std_ave, r_std_err = [],[]
    lamp, lamperr = [],[] #penetration depth
    for i in range(len(cpar)):
        par_dealing = par[:]
        par_dealing[cpar_ind] = par[cpar_ind][i]
        f2rtail = "MC"
        for j in range(len(par_dealing)):
            f2rtail+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing[j])
        f2rtail+=".txt"
        file2read = foldername + "/"+head+"_"+f2rtail
        data = np.loadtxt(file2read, skiprows=2, delimiter=",", unpack=True)
        r,r_std=data[:2]
        # analyze unu2
        unu2r_all.append(np.average(data[2:],axis=1))
        unu2rerr_all.append(np.std(data[2:],axis=1)/np.sqrt(len(data[1])))


        # analyze average edge to center distance
        r_ave.append(np.average(r))
        rho, cov0 = autocorrelation_function_fft(r)
        tau, tau_err = tau_int_cal_rho(rho,6)
        r_err.append(np.sqrt(2 * tau / len(r) * cov0))

        r_std_ave.append(np.average(r_std))
        r_std_err.append(0)
        '''
        rho, cov0 = autocorrelation_function_fft(r_std)
        tau, tau_err = tau_int_cal_rho(rho,6)
        r_std_err.append(np.sqrt(2 * tau / len(r_std) * cov0))
        print("tau,r_std,r_std_err",tau,r_std_ave[-1],r_std_err[-1])
        '''



    if(leg_num>len(cpar)):
        leg_num=cpar
    leg_ind = np.linspace(0,len(cpar)-1,leg_num,dtype=np.int)

    ppi = 72
    # LineWidth, FontSize, LabelSize = 1, 9, 8
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r"\usepackage{physics}")
    fig, axs = plt.subplots(2, 2, figsize=(
        246 / ppi*2, 246 / ppi * 2*0.8))  # , sharex=True
    rplot = np.linspace(1/bin_num,1,bin_num)-0.5/bin_num
    rperr = np.ones(len(rplot))/(bin_num*np.sqrt(3))
    # use 1/sqrt(3) for uniform distribution error

    for i in range(len(cpar)):
        pass
        Linestyle=":"
        Label=None
        if(i in leg_ind):
            Linestyle="-"
            Label=mode+"=%.*f"%(par_dg[cpar_ind],cpar[i])
        axs[0,0].errorbar(rplot,unu2r_all[i],yerr=unu2rerr_all[i],linestyle=Linestyle,label=Label)
        cos = np.sqrt(unu2r_all[i])
        tan_half = np.sqrt((1-cos)/(1+cos))
        tan_halferr = unu2rerr_all[i]/(tan_half*2*cos*(1+cos)*(1+cos))
        #print("len(tan_half)",len(tan_half))
        #print("tan_half",tan_half)
        # find lam_p
        fi=int(bin_num*0.5)
        #print("tan_half[fi:]",tan_half[fi:])
        popt,pcov=curve_fit(tan_fit,rplot[fi:],tan_half[fi:],bounds=(0, [1., 1., 0.5]),sigma=tan_halferr[fi:],absolute_sigma=True)
        popterr = np.diag(pcov)**0.5
        print("curve_fit",popt,popterr)
        print("try odr fit")
        popt,popterr = odr_lamp_popt(rplot[fi:],tan_half[fi:],rperr[fi:],tan_halferr[fi:])
        print("odr fit",popt,popterr)

        rp=rplot[fi:]
        if(i in leg_ind):
            axs[1,0].errorbar(rplot,tan_half-popt[2],yerr=tan_halferr,linestyle="None",marker="o",mfc="None",ms=3)
            axs[1,0].plot(rp,tan_fit(rp,*popt)-popt[2],color=axs[1,0].lines[-1].get_color(),linestyle=Linestyle,label=Label)

        #print("popt,popterr",popt,popterr)
        lamp.append(popt[1])
        lamperr.append(popterr[1])

    #record penetration depth data
    f2stail = "MC"
    for j in range(len(par)):
        if(j==cpar_ind):
            f2stail+="_"+par_nm[j]+"s"
        else:
            f2stail+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing[j])
    f2stail+="_ana.txt"
    # only changed "lam" and "B" mode here, others waiting for further decision
    savefile = foldername +"/"+head+"_" + f2stail
    with open(savefile,"w") as f:
        f.write(mode+",lamp,lamperr,r,r_err,r_std,r_std_err\n")
        for i in range(len(cpar)):
            f.write("%f,%f,%f,%f,%f,%f,%f\n"%(cpar[i],lamp[i],lamperr[i],r_ave[i],r_err[i],r_std_ave[i],r_std_err[i]))

    #axs[1,0].plot([rplot[int(bin_num/3)],rplot[int(bin_num/3)]],[0,1],"--k")
    axs[0,0].set_ylabel(r"$\cos^2{\theta(r/\overline{r})}$")
    axs[1,0].set_ylabel(r"$\tan(\theta(r/\overline{r})/2)$")
    axs[1,0].set_yscale("log")
    axs[1,0].set_xlabel(r"$r/\overline{r}$")
    axs[1,0].legend()

    # plot some analysis data
    axs[0,1].errorbar(cpar,r_ave,yerr=r_err)
    axs[1,1].errorbar(cpar,lamp,yerr=lamperr,linestyle="None",marker="o",mfc="None",ms=5)
    #axs[1,1].plot(cpar,lamp)
    #popt,pcov=curve_fit(sqrt_fit,cpar,lamp, sigma=lamperr, absolute_sigma=1)
    #popterr = np.diag(pcov)**0.5
    #axs[1,1].fill_between(cpar,sqrt_fit(cpar,*popt-popterr),sqrt_fit(cpar,*popt+popterr),color=axs[1,1].lines[-1].get_color(),alpha=0.4,label=r"$\sim 1/\sqrt{C_n}$")
    axs[1,1].legend()
    axs[0,1].set_ylabel(r"$\overline{r}$")
    axs[1,1].set_ylabel(r"$\lambda_p$")
    axs[1,1].set_xlabel(mode)
    if(mode=="Cn"):
        axs[1,1].set_xlabel(r"$C_n$")
    axs[1,1].set_ylim(0.0,0.5)
    fig.suptitle(tag)
    plt.tight_layout()
    #plt.show()
    plt.savefig(savefile[:-4]+".pdf",format="pdf")
    plt.close()


def cos2_fit(l,p):
    cos=np.cos(2*np.pi*l/p)
    return cos*cos

def cos_fit(l,p):
    cos = np.cos(4*np.pi*l/p)
    return cos

def twistl_stat_plot(foldername, par, par_nm, par_dg, mode, d0=1.5, head="nunu2l",tag="",leg_num=5,bin_num=40):
    Ne = par[find_cpar_ind(par_nm,"Ne")]
    cpar_ind = find_cpar_ind(par_nm,mode)
    cpar = par[cpar_ind]
    nunu2l_all = []
    nunu2lerr_all=[]
    p,perr=[],[]
    Ls=[]
    if(mode=="L"):
        Ls=cpar
    for i in range(len(cpar)):
        if(mode!="L"):
            Ls.append(par[find_cpar_ind(par_nm,"L")])
        par_dealing = par[:]
        par_dealing[cpar_ind] = par[cpar_ind][i]
        f2rtail = "MC"
        for j in range(len(par_dealing)):
            f2rtail+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing[j])
        f2rtail+=".txt"
        file2read = foldername + "/"+head+"_"+f2rtail
        data = np.loadtxt(file2read, skiprows=2, delimiter=",", unpack=True)
        nunu2l_all.append(np.average(data,axis=1))
        nunu2lerr_all.append(np.std(data,axis=1)/np.sqrt(len(data[0])))

    f2stail = "MC"
    for j in range(len(par)):
        if(j==cpar_ind):
            f2stail+="_"+par_nm[j]+"s"
        else:
            f2stail+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing[j])
    f2stail+="_ana.txt"
    # only changed "lam" and "B" mode here, others waiting for further decision
    savefile = foldername +"/"+head+"_" + f2stail

    if(leg_num>len(cpar)):
        leg_num=cpar
    leg_ind = np.linspace(0,len(cpar)-1,leg_num,dtype=np.int)

    ppi = 72
    # LineWidth, FontSize, LabelSize = 1, 9, 8
    plt.rc('text', usetex=True)
    fig, axs = plt.subplots(2, 2, figsize=(
        246 / ppi*2, 246 / ppi * 2 *0.8))  # , sharex=True
    l_norm = np.linspace(1/bin_num,1,bin_num)-0.5/bin_num
    lf=(np.array(Ls)-1)*d0
    for i in range(len(cpar)):
        pass
        Linestyle=":"
        Label=None
        if(i in leg_ind):
            Linestyle="-"
            Label=mode+"=%.*f"%(par_dg[cpar_ind],cpar[i])
        lplot=lf[i]*l_norm

        cos2l = 2*nunu2l_all[i]-1
        axs[0,0].errorbar(lplot,cos2l+0.1*i,yerr=2*nunu2lerr_all[i],linestyle=Linestyle,label=Label)
        #axs[0,0].errorbar(lplot,nunu2l_all[i]+0.1*i,yerr=nunu2lerr_all[i],linestyle=Linestyle,label=Label)

        # find lam_p by fitting cos2
        #popt,pcov=curve_fit(cos2_fit,lplot,nunu2l_all[i],p0=3*lf[i],bounds=(0.2*lf[i],10*lf[i]))
        popt,pcov=curve_fit(cos_fit,lplot,cos2l,p0=lf[i],bounds=(0.1*lf[i],10*lf[i]),maxfev=1000)
        popterr = np.diag(pcov)**0.5
        p.append(popt[0])
        perr.append(popterr[0])
        if(i in leg_ind):
            Label=mode+"=%.*f"%(par_dg[cpar_ind],cpar[i])
            #axs[1,0].errorbar(lplot,nunu2l_all[i]+0.1*i,yerr=nunu2lerr_all[i],linestyle="None")
            axs[1,0].errorbar(lplot,cos2l+0.1*i,yerr=2*nunu2lerr_all[i],linestyle="None")
            axs[1,0].plot(lplot,cos_fit(lplot,*popt)+0.1*i,linestyle=Linestyle,color=axs[1,0].lines[-1].get_color(),label=Label)
        #print("popt,popterr",popt,popterr)
            # find and plot fft
            cos2l = 2*nunu2l_all[i]-1
            nuqnuq=scipy.fft.dct(cos2l)
            n = len(cos2l)
            freq = np.array(range(n))*np.pi/n*d0
            q= np.array(range(n))*d0/(4*n) #q=1/p ish
            #make freq unit of 1/\sigma0
            axs[0,1].plot(q[:int(n*2/3)],nuqnuq[:int(n*2/3)],linestyle="None",marker="o",mfc="None",ms=3,color=axs[1,0].lines[-1].get_color(),label=Label)

    axs[0,1].set_ylabel(r"$\mathcal{F}[\cos{2\theta(l)}]$")
    axs[0,1].set_xlabel(r"$q(1/\sigma_0)$")
    axs[0,1].legend()
    axs[1,1].errorbar(cpar,p,yerr=perr,linestyle="None",marker="o",mfc="None",ms=3)
    axs[0,0].set_ylabel(r"$\cos{2\theta(l)}$")
    axs[0,0].legend()
    axs[1,1].set_ylabel(r"$\cos^2{\theta(l)}$")
    axs[1,1].set_xlabel(r"$l$")
    axs[1,1].set_ylabel(r"$P$")
    axs[1,1].set_xlabel(mode)


    fig.suptitle(tag)
    plt.tight_layout()
    plt.savefig(savefile[:-4]+".pdf",format="pdf")
    #plt.show()
    plt.close()
