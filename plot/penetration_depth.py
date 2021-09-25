import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from config import *
from scipy import stats


def get_tiltr_conf_data(foldername,Kd,q,Cns):
    print("getting un2r tilt(r) configuration")
    filenames,rots,xlims,ylims = [],[],[],[]
    for i in range(len(Cns)):
        filenames.append(foldername+"/State_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd%.1f_q%.1f_Cn%.1f_kard0.0.txt"%(Kd,q,Cns[i]))
    rots.append([0.6,0])
    xlims.append([8,15])
    ylims.append([-8,8])
    rots.append([-0.4,0])
    xlims.append([8,15])
    ylims.append([-8,8])
    rots.append([-0.6,0])
    xlims.append([8,15])
    ylims.append([-8,8])
    return (filenames,rots,xlims,ylims)

def get_tilt_tan_data(foldername,Kd,q,Cns):
    print("getting un2r tilt(r) data")
    un2rs,un2rerrs = [], []
    for i in range(len(Cns)):
        file2read = foldername+"/un2r_MC_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd%.1f_q%.1f_Cn%.1f_kard0.0.txt"%(Kd,q,Cns[i])
        data = np.loadtxt(file2read, skiprows=2, delimiter=",", unpack=True)
        un2rs.append(np.average(data[2:],axis=1))
        un2rerrs.append(np.std(data[2:],axis=1)/np.sqrt(len(data[1])))
    bin_num = len(data[2:])
    rplot = np.linspace(1/bin_num,1,bin_num)-0.5/bin_num
    #rperr = np.ones(len(rplot))/(bin_num*np.sqrt(3))
    cos = np.sqrt(un2rs)
    tan_half = np.sqrt((1-cos)/(1+cos))
    tan_halferr = un2rerrs/(tan_half*2*cos*(1+cos)*(1+cos))
    return (rplot,tan_half,tan_halferr)




def get_pdepth_related_data(foldername,qs,Kds):
    print("getting penetration depth related data")
    #qs = [0.3,0.5,0.7]
    #Kds = np.arange(3.0,10.1,1.0)
    #Cns = np.arange(5.0,15.1,1.0) # mode variable
    data = []
    for i in range(len(qs)):
        data.append([])
        for j in range(len(Kds)):
            #  lamp data
            filename = foldername +"/un2r_MC_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd%.1f_q%.1f_Cns_kard0.0_ana.txt"%(Kds[j],qs[i])
            data[i].append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))
            #print("data[i][j]",data[i][j])
            data[i][j][0]/=Kds[j]
            #print("data[i][j][0]/Kds[j]",data[i][0])
    data = np.transpose(np.array(data), axes=(2,0,1,3))
    cpar,lamp,lamperr,r,r_err,r_std,r_std_err=data
    return (cpar,lamp,lamperr)

def pdepth_plot(LineWidth, FontSize, LabelSize):
    print("plotting penetration depth figure")
    #foldername="../data/Ne1/Feb14_2021"
    foldername="../data/Ne1/Mar2_2021"
    qs = np.array([0.3,0.5,0.7])
    Kds = np.arange(3.0,10.1,1.0)
    # get related observables

    # plot the figure
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.2))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    axconf = plt.subplot2grid((3,2), (0, 0)) # configuration, cur to show the edge
    axtiltr = plt.subplot2grid((3,2), (0, 1)) # single
    ax = plt.subplot2grid((3,2), (1, 0),colspan=2,rowspan=2)

    # plot configuration
    Kd,q,Cns=10.0,0.7,[15.0,10.0,5.0]
    rn=len(Cns)
    filename,rot,xlims,ylims = get_tiltr_conf_data(foldername,Kd,q,Cns)
    axconf.text(-0.2,0, r"$C=$", fontsize=FontSize)
    axconfsubs = []
    for i in range(rn):
        axconfsubs.append(axconf.inset_axes([i/rn, 0.1, 1/rn, 0.9]))
        ax_config_xy(axconfsubs[i],filename[i],Color="gray",LineWidth=LineWidth*1.65,rotation=rot[i],bead=1,mesh=1,rod=1,d=0.65,xlim=xlims[i],ylim=ylims[i])
        axconf.text(i/rn+0.1,0, r"$%.0f$"%Cns[i], fontsize=FontSize)
    axconf.set_frame_on(False)
    axconf.set_yticks([])
    axconf.set_xticks([])
    axconf.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)
    x1, y1 = -0.2, 0.9
    axconf.text(axconf.get_xlim()[1]*x1+axconf.get_xlim()[0]* (1-x1),  axconf.get_ylim()[1]*y1+axconf.get_ylim()[0]* (1-y1), r"(a)", fontsize=FontSize)

    # plot tilt(r) functions
    rplot,tan_halfs,tan_halferrs = get_tilt_tan_data(foldername,Kd,q,Cns)
    lr=int(len(rplot)/2)
    markers=["v","s","p"]
    msize=4
    ls=["-","--",":"]
    for i in range(rn):
        #axtiltr.errorbar(rplot[lr:],tan_halfs[i][lr:],yerr = tan_halferrs[i][lr:],color="gray",marker=markers[i],ms=msize,mfc="None",linestyle="None",linewidth=LineWidth,label=r"$C_n=%.1f$"%Cns[i])
        axtiltr.plot(rplot[lr:],tan_halfs[i][lr:],color="black",linestyle=ls[i],linewidth=LineWidth,label=r"$%.0f$"%Cns[i])
        #axtiltr.errorbar(rplot[lr:],tan_halfs[i][lr:],yerr = tan_halferrs[i][lr:],color="gray",linestyle=ls[i],linewidth=LineWidth,label=r"$C_n=%.1f$"%Cns[i])
    axtiltr.tick_params(which="both",direction="in", top="on", right="on",bottom="on", labelsize=LabelSize)
    axtiltr.xaxis.set_label_position("top")
    axtiltr.xaxis.tick_top()
    axtiltr.xaxis.set_major_locator(MultipleLocator(0.2))
    axtiltr.xaxis.set_minor_locator(MultipleLocator(0.1))
    axtiltr.set_xlabel(r"$r/\overline{r}$")
    axtiltr.yaxis.set_major_locator(MultipleLocator(0.1))
    axtiltr.yaxis.set_minor_locator(MultipleLocator(0.05))
    axtiltr.set_ylabel(r"$\tan(\theta/2)$")
    axtiltr.legend(title=r"$C=$",ncol=1 ,title_fontsize=LabelSize,fontsize=LabelSize, handlelength=1.2, columnspacing=0.5,handletextpad=0.4, frameon=False)
    x1, y1 = -0.2, 0.9
    axtiltr.text(axtiltr.get_xlim()[1]*x1+axtiltr.get_xlim()[0]* (1-x1),  axtiltr.get_ylim()[1]*y1+axtiltr.get_ylim()[0]* (1-y1), r"(b)", fontsize=FontSize)

    # plot lamp versus Kd/Cn
    cpars,lamps,lamperrs = get_pdepth_related_data(foldername,qs,Kds)
    #cpars[i,j] i for  q, j for Kd cpar = Cn/Kd here
    #print(np.shape(cpars))
    mcolors = ["red","green","blue"]
    markers = ["+","o","s"]
    msize=4
    malpha=0.7
    for i in range(len(qs)):
        cpar,lamp,lamperr = cpars[i].flatten(),lamps[i].flatten(),lamperrs[i].flatten()
        res = stats.linregress(np.power(cpar,-0.5), lamp)
        print("res.rvalue**2",res.rvalue**2)
        ax.errorbar(np.power(cpar,-0.5),lamp,yerr=lamperr,marker=markers[i],mfc="None",linestyle="None",color=mcolors[i],ms=msize,alpha=malpha,label=r"$%.1f$"%qs[i])
    ax.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelsize=LabelSize)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_xlabel(r"$\sqrt{\epsilon_{LL}/C}$",fontsize=FontSize)
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.set_ylabel(r"$\lambda_p$",fontsize=FontSize)
    ax.legend(title=r"$k_c=$",loc="lower right",ncol=1 ,title_fontsize=LabelSize,fontsize=LabelSize, handlelength=1.2, columnspacing=0.5,handletextpad=0.4, frameon=False)
    x1, y1 = 0.1, 0.9
    ax.text(ax.get_xlim()[1]*x1+ax.get_xlim()[0]* (1-x1),  ax.get_ylim()[1]*y1+ax.get_ylim()[0]* (1-y1), r"(c)", fontsize=FontSize)

    plt.tight_layout(pad=0.5)
    plt.savefig("figures/penetration_depth.pdf")
    plt.close()
