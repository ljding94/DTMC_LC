import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from config import *
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from scipy.ndimage.filters import gaussian_filter

def get_negK_conf_data(foldername,q,lamCns,mode="lamCn"):
    print("getting un2r tilt(r) configuration")
    filenames,rots = [],[]
    if(mode=="q"):
        # for negK lam plot
        for i in range(len(q)):
            lam,Cn=lamCns
            Kd=Cn
            #Kd=10.0
            filenames.append(foldername+"/State_N500_Ne1_L-1_kar100_karg0.0_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_kard0.0.txt"%(lam,Kd,q[i],Cn))
        rots=[(0.0,-0.8),(0.0,-0.8),(0.0,-0.8)]
    else:
        # for negK pcolormesh plot
        for i in range(len(lamCns)):
            lam,Cn=lamCns[i]
            Kd=Cn
            #Kd=10.0
            filenames.append(foldername+"/State_N500_Ne1_L-1_kar100_karg0.0_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_kard0.0.txt"%(lam,Kd,q,Cn))
        rots=[(2.8,-0.8),(1.0,-1.1),(0.7,-1.0),(0.7,-0.8)]
    return (filenames,rots)


def get_negK_related_data(foldername,q,Cn,mode="Cn"):
    print("getting negK, rippling related O data")
    #get ana data from the folder
    data = []
    if(mode=="q"):
        for i in range(len(q)):
            # usual O data
            Kd = Cn
            filename = foldername +"/O_MC_N500_Ne1_L-1_kar100_karg0.0_lams_Kd%.1f_q%.1f_Cn%.1f_kard0.0_ana.txt"%(Kd,q[i],Cn)
            data.append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))
    else:
        for i in range(len(Cn)):
            # usual O data
            Kd = Cn[i]
            #Kd=10.0
            filename = foldername +"/O_MC_N500_Ne1_L-1_kar100_karg0.0_lams_Kd%.1f_q%.1f_Cn%.1f_kard0.0_ana.txt"%(Kd,q,Cn[i])
            data.append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))

    data = np.transpose(np.array(data), axes=(1, 0, 2))
    print("negK np.shape(data)",np.shape(data))
    cpar,IK_ave,IK_err,uuc_ave,uuc_err = data[0],data[16],data[18],data[22],data[24]
    return (cpar,IK_ave,IK_err,uuc_ave,uuc_err)


def negK_lam_Cn_plot(LineWidth, FontSize, LabelSize):
    print("plotting negative as function of lambda and Cn")
    #foldername="../data/Ne1/Feb19_2021"
    foldername="../data/Ne1/Mar3_2021"
    #foldername="../data/Ne1/Mar26_2021"
    q=2.0
    lams=np.arange(4.0,10.1,0.5)
    Cns=np.arange(3.0,10.1,0.5)
    cpar,IK_ave,IK_err,uuc_ave,uuc_err = get_negK_related_data(foldername,q,Cns)
    lamp,Cnp=np.meshgrid(lams,Cns)
    print("cpar-lamp",cpar-lamp)
    negK_ave = -IK_ave/(4*np.pi)
    negK_err = IK_err/(4*np.pi)

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.1))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')

    ax = plt.subplot2grid((1,1), (0, 0))

    axconf = ax.inset_axes([-0.1, 0.9, 1.3, 0.8])
    # plot some configuration
    lamCns=[(4,8.5),(6,6.5),(8,4.5),(10,3.5)]
    filenames,rots=get_negK_conf_data(foldername,q,lamCns)
    nr = len(filenames)
    arrowpos=[(4.5,10.9),(6.5,11),(8.5,11),(10.5,11)]
    for i in range(nr):
        ax_config_xy(axconf,filenames[i],Color="gray",LineWidth=LineWidth*0.6,rotation=rots[i],bead=0,mesh=1,rod=0,xshift=30*i)
        ax_config_xy(axconf,filenames[i],Color="gray",LineWidth=LineWidth*0.6,rotation=rots[i],bead=1,mesh=1,rod=1,d=0.8,xshift=30*i,yshift=20)
        ax.annotate(" ",xy=lamCns[i], xytext=arrowpos[i],arrowprops=dict(arrowstyle="->",lw=LineWidth,color="gray"),annotation_clip=False)
        #axconf.text(i/rn+0.1,0, r"$%.0f$"%lams[i], fontsize=FontSize)
    x1, y1 = -0.035, 0.8
    axconf.text(axconf.get_xlim()[1]*x1+axconf.get_xlim()[0]* (1-x1),axconf.get_ylim()[1]*y1+axconf.get_ylim()[0]* (1-y1), r"(a)", fontsize=FontSize)


    im = ax.pcolormesh(lamp,Cnp,negK_ave,shading="auto",cmap=cm.get_cmap("rainbow"))
    #negK_smooth=gaussian_filter(negK_ave,sigma=1)
    #im = ax.pcolormesh(lamp,Cnp,negK_smooth,shading="auto")


    cbar=fig.colorbar(im, ax=ax)
    #cbar=plt.colorbar(cmap=cm.get_cmap("jet_r"),ax=ax,orientation='horizontal')
    cbar.set_label(r"$-\int K \dd{A}/(4\pi)$",fontsize=FontSize)
    #cbar.set_label(r"negative Gaussian curvature",fontsize=FontSize)

    # dashed line seperating different orders of Enneper's surface
    lam_dash = np.linspace(4,10,100)
    #ax.plot(lam_dash,0.8*lam_dash,"--k")
    #ax.plot(lam_dash,1.2*lam_dash,":k")

    ax.set_xlim(3.75,10.25)
    ax.set_ylim(2.75,10.25)
    ax.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelsize=LabelSize)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_xlabel(r"$\lambda$",fontsize=FontSize)
    #ax.set_xlabel(r"line tension",fontsize=FontSize)
    ax.yaxis.set_major_locator(MultipleLocator(2.0))
    ax.yaxis.set_minor_locator(MultipleLocator(1.0))
    ax.set_ylabel(r"$C=\epsilon_{LL}$",fontsize=FontSize)
    #ax.set_ylabel(r"tilt coupling",fontsize=FontSize)
    x1, y1 = -0.15, 0.9
    ax.text(ax.get_xlim()[1]*x1+ax.get_xlim()[0]* (1-x1),  ax.get_ylim()[1]*y1+ax.get_ylim()[0]* (1-y1), r"(b)", fontsize=FontSize)
    #plot curves, use yaxis as negK instead of color
    '''
    for i in range(len(Cns)):
        if(i%2==0):
            axcurv.errorbar(lamp[i],negK_ave[i],yerr=negK_err[i],label=r"$C=%.0f$"%Cns[i])
    '''
    plt.tight_layout(pad=0.5)
    plt.savefig("figures/negK.pdf")


def negK_lam_q_plot(LineWidth, FontSize, LabelSize):
    print("plotting negative as function of lambda and q")
    #foldername="../data/Ne1/Mar28_2021"
    foldername="../data/Ne1/Apr4_2021"
    q=[1.5,1.9,2.3]
    mcolors = ["red","green","blue"]
    markers = ["^","o","s"]
    msize=4
    malpha=0.7
    Cn,Kd=5.0,5.0
    cpar,IK_ave,IK_err,uuc_ave,uuc_err = get_negK_related_data(foldername,q,Cn,mode="q")
    negK_ave = -IK_ave/(4*np.pi)
    negK_err = IK_err/(4*np.pi)

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.6))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')

    ax = plt.subplot2grid((1,2), (0, 1))
    axconf = plt.subplot2grid((1,2), (0, 0))
    # plot some configuration

    lamCn=(6.0,Cn) # Kd = Cn in this case, just for config data
    filenames,rots=get_negK_conf_data(foldername,q,lamCn,mode="q")
    # assign it's own rot based on observation
    rots=[(0.2,-0.5),(0.5,-1),(0,-1)]
    nr = len(filenames)
    for i in range(nr):
        ax_config_xy(axconf,filenames[i],Color="gray",LineWidth=LineWidth*0.6,rotation=rots[i],bead=1,mesh=1,rod=1,yshift=-22*i)
        ax_config_xy(axconf,filenames[i],Color="gray",LineWidth=LineWidth*0.6,rotation=rots[i],bead=0,mesh=1,rod=0,xshift=30,yshift=-22*i)
        #ax.annotate(" ",xy=lamCns[i], xytext=arrowpos[i],arrowprops=dict(arrowstyle="->",lw=LineWidth,color="gray"),annotation_clip=False)
        #axconf.text(i/rn+0.1,0, r"$%.0f$"%lams[i], fontsize=FontSize)
    x1, y1 = 0.1, 0.97
    axconf.text(axconf.get_xlim()[1]*x1+axconf.get_xlim()[0]* (1-x1),axconf.get_ylim()[1]*y1+axconf.get_ylim()[0]* (1-y1), r"(a)", fontsize=FontSize)

    # negK vs. lam curves
    ns=1
    for i in range(len(q)):
        ax.errorbar(cpar[i][::ns],negK_ave[i][::ns],yerr=negK_err[i][::ns],marker=markers[i],ms=msize,color=mcolors[i],alpha=malpha,linestyle="None",label=r"$%.1f$"%q[i])

    #ax.set_xlim(3.75,10.25)
    #ax.set_ylim(2.75,10.25)
    ax.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelsize=LabelSize)
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_xlabel(r"$\lambda$",fontsize=FontSize)
    #ax.set_xlabel(r"line tension",fontsize=FontSize)
    ax.set_ylim(0.0,1.7)
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.set_ylabel(r"$-\int K \dd{A}/(4\pi)$",fontsize=FontSize)
    #handles, labels = ax.get_legend_handles_labels()
    #ax.legend(handles[::-1], labels[::-1],title=r"$k_c^*=$", ncol=1,handlelength=1.0, columnspacing=1.0,frameon=False,fontsize=FontSize)
    ax.legend(title=r"$k_c^*=$", ncol=1,handlelength=1.0, columnspacing=1.0,frameon=False,fontsize=FontSize)
    x1, y1 = 0.2, 0.9
    ax.text(ax.get_xlim()[1]*x1+ax.get_xlim()[0]* (1-x1),  ax.get_ylim()[1]*y1+ax.get_ylim()[0]* (1-y1), r"(b)", fontsize=FontSize)
    plt.tight_layout(pad=0.5)
    plt.savefig("figures/negK_lam.pdf")
