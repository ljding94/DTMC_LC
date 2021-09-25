import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from config import *
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)


def Oq_data_get(foldername,Kds,Cns,msize,malpha):
    #get ana data from the folder
    # include un2below data

    Olabels,Ocolors,Omarkers=[],[],[]
    data,dataun = [],[]
    lampts = []
    colors = ["red","green","blue","royalblue"]
    markers = ["v", "s", "p", "h","o"]
    legend_elements = []
    # include Kds legend
    for i in range(len(Kds)):
        legend_elements.append(Line2D([0], [0],ls="None", color=colors[i] ,marker="o",ms=msize,alpha=malpha, label=r"$\epsilon_{LL}=%.0f$"%Kds[i]))
    for j in range(len(Cns)):
        legend_elements.append(Line2D([0], [0],ls="None", color="black", mfc="None",marker=markers[j], ms=msize,alpha=malpha,label=r"$C=%.0f$"%Cns[j]))
    for i in range(len(Kds)):
        for j in range(len(Cns)):
            Olabels.append(None)
            if(i==0):
                Olabels[-1]=r"$C$=%.0f"%Cns[j]
            elif(j==0):
                Olabels[-1]=r"$\epsilon_{LL}$=%.0f"%Kds[i]
            Ocolors.append(colors[i])
            Omarkers.append(markers[j])
            #different color Kd
            #different marker Cn
            # usual O data
            filename = foldername +"/O_MC_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd%.1f_qs_Cn%.1f_kard0.0_ana.txt"%(Kds[i],Cns[j])
            data.append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))
            # un2below data
            #filename = foldername +"/un2dis_MC_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd%.1f_qs_Cn%.#1f_kard0.0_ana.txt"%(Kds[i],Cns[j])
            #dataun.append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))
            lampts.append(np.ones(len(data[-1][0]))*np.sqrt(Kds[i]/Cns[j]))

    data = np.transpose(np.array(data), axes=(1, 0, 2))
    print("np.shape(data)",np.shape(data))
    qs,IK_ave,IK_err,uuc_ave,uuc_err,un2_ave,un2_err = data[0],data[16],data[18],data[22],data[24],data[25],data[27]
    #dataun = np.transpose(np.array(dataun), axes=(1, 0, 2))
    return (qs,np.array(lampts),IK_ave,IK_err,uuc_ave,uuc_err,un2_ave,un2_err,legend_elements,Ocolors,Omarkers)

def un2disq_data(foldername,Kd,qs,Cn):
    print("getting u2dis data")
    # Apr22_2021 data uses distribution of theta instead

    un2dis_all = []
    un2pdf = []
    un2x = []
    Olabels,Ocolors,Omarkers=[],[],[]
    colors = ["blue","tomato"]
    linestyles = ["--","-",":"]
    legend_elements = []
    for i in range(len(qs)):
        legend_elements.append(Line2D([0], [0],ls=linestyles[i],color=colors[i],marker="None",linewidth=1, label=r"$%.1f$"%qs[i])) # since here Kd=Cn=5, Kd/Cn=1
        filename = foldername +"/un2dis_MC_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd%.1f_q%.1f_Cn%.1f_kard0.0.txt"%(Kd,qs[i],Cn)
        data = np.loadtxt(filename,skiprows=2,delimiter=",",unpack=True)
        bin_num = len(data)
        un2pdf.append(np.average(data,axis=1)*bin_num)
        un2x.append(np.linspace(0.5/bin_num,1-0.5/bin_num,bin_num))
        Ocolors.append(colors[i])
    return (un2x,un2pdf,linestyles,legend_elements,Ocolors)


def O_q_reduced(LineWidth, FontSize, LabelSize):
    # plot uuc versus q*sqrt(Kd/Cn) see if the curves converge on the x axis
    print("uuc_q_norm() in process")
    #foldername="../data/Ne1/Feb25_2021"
    foldername="../data/Ne1/Apr22_2021"
    # separate data based on qc*lampt
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 3*0.45))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    axconf = plt.subplot2grid((3, 5), (0, 0),colspan=2)
    axdis = plt.subplot2grid((3, 5), (0, 2),colspan=3)
    axuc = plt.subplot2grid((3, 5), (1, 0),colspan=5)
    #axPun = plt.subplot2grid((3, 3), (2, 0),colspan=3,sharex=axuc)
    axIK = plt.subplot2grid((3, 5), (2, 0),colspan=5,sharex=axuc)

    qs = [1.1,1.4] # q to plot for configuration and un2 distribution
    Kd,Cn=3.0,3.0
    # plot membrane config
    axconf.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)
    rot = [(0.2,-1.0),(1.4,-1.0)]
    for i in range(len(qs)):
        filename=foldername+"/State_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd%.1f_q%.1f_Cn%.1f_kard0.0.txt"%(Kd,qs[i],Cn)
        ax_config_xy(axconf,filename,Color="gray",LineWidth=LineWidth*0.65,rotation=rot[i],yshift=-21*i,bead=1,mesh=1,rod=1,d=0.65)
        axconf.text(-20,-21*i+9,r"$k_c^*=%.1f$"%qs[i],fontsize=FontSize)

    #cbar=plt.colorbar(cm.ScalarMappable(norm=Normalize(vmin=0,vmax=0.5*np.pi), cmap=cm.get_cmap("jet_r")),ax=axconf,orientation='horizontal',ticks=[0,0.25*np.pi,0.5*np.pi])
    #cbar.ax.set_xticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"],fontsize=FontSize)
    #cbar.ax.set_xticklabels([r"$0$",r"$0.25$",r"$0.5$",r"$0.75$",r"$1$"],fontsize=FontSize)
    #cbar.ax.tick_params(direction="in",labelsize=LabelSize)
    #cbar.ax.set_title(r"$\arccos{|\vu{u}\cdot\vu{n}|}$",fontsize=FontSize)
    x1, y1 = -0.3, 0.2
    axconf.text(axconf.get_xlim()[1]*x1+axconf.get_xlim()[0]* (1-x1),  axconf.get_ylim()[1]*y1+axconf.get_ylim()[0]* (1-y1), r"(a)", fontsize=FontSize)

    # plot distribution data
    un2x,un2pdf,linestyles,un2dislegend,Ocolors = un2disq_data(foldername,Kd,qs,Cn)
    for i in range(len(un2x)):
        axdis.plot(un2x[i]*0.5*np.pi,un2pdf[i]*2/np.pi,linestyle=linestyles[i],color=Ocolors[i])
        # Apr22 data used code such that undis is actually measuring the distribution of theta
        #axdis.plot(un2pdf[i],un2x[i],linestyle=linestyles[i],color=Ocolors[i])

        #theta=np.arccos(np.sqrt(un2x[i]))
        #pdftheta=un2pdf[i]*2*np.sqrt(un2x[i])*np.sqrt(1-un2x[i])
        #axdis.plot(theta,un2pdf[i]*np.sin(2*theta),linestyle=linestyles[i],color=Ocolors[i])
        #axdis.plot(theta,pdftheta,linestyle=linestyles[i],color=Ocolors[i])
        print("normalized?: int un2pdf", np.sum(un2pdf[i])*(un2x[i][1]-un2x[i][0]))
    axdis.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelsize=LabelSize)
    #axdis.yaxis.set_label_position("right")
    #axdis.yaxis.tick_right()
    #axdis.set_ylabel(r"$p((\vu{u}\cdot\vu{n})^2)$",fontsize=FontSize)
    axdis.set_ylabel(r"$p(\theta)$",fontsize=FontSize)
    axdis.xaxis.set_label_position("top")
    axdis.xaxis.tick_top()
    #axdis.set_xlim(0.0,1.6)
    axdis.set_xticks([0,0.125*np.pi,0.25*np.pi,0.375*np.pi,0.5*np.pi])
    axdis.set_xticklabels([r"$0$",r"$\pi/8$",r"$\pi/4$",r"$3\pi/8$",r"$\pi/2$"])
    axdis.set_xlabel(r"$\theta=\arccos{|\vu{u}\cdot\vu{n}|}$",fontsize=FontSize)
    axdis.legend(title=r"$k_c^*=$",handles=un2dislegend, handlelength=1.0, columnspacing=1.0,frameon=False,fontsize=FontSize)
    axdis.yaxis.set_major_locator(MultipleLocator(1))
    axdis.yaxis.set_minor_locator(MultipleLocator(0.5))
    #axdis.xaxis.set_major_locator(MultipleLocator(0.4))
    axdis.xaxis.set_minor_locator(MultipleLocator(np.pi/16))
    x1, y1 = 0.1, 0.2
    axdis.text(axdis.get_xlim()[1]*x1+axdis.get_xlim()[0]* (1-x1),  axdis.get_ylim()[1]*y1+axdis.get_ylim()[0]* (1-y1), r"(b)", fontsize=FontSize)

    # axuc uc vs. q
    Kds = np.arange(3.0,7.1,2.0)
    Cns = np.arange(3.0,7.1,2.0)
    msize,malpha=4,0.7
    qs,lampts,IK_ave,IK_err,uucs,uuc_errs,un2s,un2_errs,Olegend,Ocolors,Omarkers = Oq_data_get(foldername,Kds,Cns,msize,malpha)
    # pre process data
    qns = qs*lampts

    for i in range(len(qns)):
        axuc.errorbar(qns[i],uucs[i],yerr=uuc_errs[i],linestyle="None",marker=Omarkers[i],mfc="None",ms=msize,color=Ocolors[i],alpha=malpha)
        # highlight cholesteric phase
        cho = uucs[i]<0.1
        axuc.errorbar(qns[i][cho],uucs[i][cho],yerr=uuc_errs[i][cho],linestyle="None",marker=Omarkers[i],ms=msize,color=Ocolors[i],alpha=malpha)
    #axuc.set_ylabel(r"$\left<(\vu{u}_i\cross\vu{u}_j)\cdot\vu{r}_{ij} (\vu{u}_i\cdot \vu{u}_j)\right>_{(i,j)}$",fontsize=LabelSize)
    axuc.set_ylabel(r"$T_s$",fontsize=LabelSize)
    axuc.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False, labelsize=LabelSize)
    axuc.legend(handles=Olegend,ncol=2, handlelength=1.5,handletextpad=0.1, columnspacing=0.2,frameon=False,fontsize=LabelSize)
    #axuc.set_ylim(0.0,0.3)
    axuc.yaxis.set_major_locator(MultipleLocator(0.1))
    axuc.yaxis.set_minor_locator(MultipleLocator(0.05))
    x1, y1 = 0.05, 0.55
    axuc.text(axuc.get_xlim()[1]*x1+axuc.get_xlim()[0]* (1-x1),  axuc.get_ylim()[1]*y1+axuc.get_ylim()[0]* (1-y1), r"(c)", fontsize=FontSize)
    plt.subplots_adjust(hspace=0.5)
    # plot IK versus q*
    for i in range(len(qns)):
        axIK.errorbar(qns[i],-IK_ave[i]/(4*np.pi),yerr=IK_err[i]/(4*np.pi),linestyle="None",marker=Omarkers[i],mfc="None",ms=msize,color=Ocolors[i],alpha=malpha)
        # highlight cholesteric phase
        cho = uucs[i]<0.1
        axIK.errorbar(qns[i][cho],-IK_ave[i][cho]/(4*np.pi),yerr=IK_err[i][cho]/(4*np.pi),linestyle="None",marker=Omarkers[i],ms=msize,color=Ocolors[i],alpha=malpha)
    axIK.set_ylabel(r"$-\int K \dd{A} /(4\pi)$",fontsize=FontSize)
    axIK.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelsize=LabelSize)
    axIK.set_xlabel(r"$k_c^*=k_c\sqrt{\epsilon_{LL}/C}$",fontsize=FontSize)
    axIK.set_ylim(0.0,1.0)
    axIK.yaxis.set_major_locator(MultipleLocator(0.2))
    axIK.yaxis.set_minor_locator(MultipleLocator(0.1))
    axIK.xaxis.set_major_locator(MultipleLocator(1.0))
    axIK.xaxis.set_minor_locator(MultipleLocator(0.2))
    x1, y1 = 0.05, 0.55
    axIK.text(axIK.get_xlim()[1]*x1+axIK.get_xlim()[0]* (1-x1),  axIK.get_ylim()[1]*y1+axIK.get_ylim()[0]* (1-y1), r"(d)", fontsize=FontSize)

    # axPun un2below vs. q
    # replaced with IK for displacing shape formation
    #qs,lampts,un2below,un2below_err,legend_elements,Ocolors,Omarkers = un2below_q_data_get()
    #qns = qs*lampts
    '''
    for i in range(len(qns)):
        axPun.errorbar(qns[i],un2below[i],yerr=un2below_err[i],linestyle="None",marker=Omarkers[i],mfc="None",ms=2,color=Ocolors[i])
    axPun.set_ylabel(r"$P((\hat{u}\cdot\hat{n})^2<1/2)$",fontsize=FontSize)
    axPun.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelsize=LabelSize)
    axPun.set_xlabel(r"$q^*=q(K_d/C_n)^{1/2}$",fontsize=FontSize)
    axPun.yaxis.set_minor_locator(MultipleLocator(0.1))
    axPun.xaxis.set_major_locator(MultipleLocator(1.0))
    axPun.xaxis.set_minor_locator(MultipleLocator(0.2))
    x1, y1 = 0.05, 0.85
    axPun.text(axPun.get_xlim()[1]*x1+axPun.get_xlim()[0]* (1-x1),  axPun.get_ylim()[1]*y1+axPun.get_ylim()[0]* (1-y1), r"(d)", fontsize=FontSize)
    '''


    plt.tight_layout(pad=0.1)
    plt.savefig("figures/O_kc_reduced.pdf",format="pdf")


