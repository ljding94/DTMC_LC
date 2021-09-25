import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from config import *
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


def get_phase_related_data(foldername,Cns,Kds):
    print("getting phase diagram related data")
    #foldername="../data/Ne1/Feb12_2021"
    foldername="../data/Ne1/Feb24_2021"
    #get ana data from the folder
    data = []
    for i in range(len(Cns)):
        for j in range(len(Kds)):
            # usual O data
            filename = foldername +"/O_MC_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd%.1f_qs_Cn%.1f_kard0.0_ana.txt"%(Kds[j],Cns[i])
            data.append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))

    data = np.transpose(np.array(data), axes=(1, 0, 2))
    print("np.shape(data)",np.shape(data))
    qs,p2uu_ave,p2uu_err,uuc_ave,uuc_err,un2_ave,un2_err = data[0],data[19],data[21],data[22],data[24],data[25],data[27]
    return (qs,p2uu_ave,p2uu_err,uuc_ave,uuc_err,un2_ave,un2_err)


def phase_diagram_plot(LineWidth, FontSize, LabelSize):
    #qs = np.array([0.0,0.4,0.8,1.0,1.2,1.4,1.6,2.0,2.4,2.8])
    qs = np.arange(0.0,2.71,0.3)
    Kds =np.arange(0.0,3.61,0.4)
    #Kds = np.array([0.0,0.4,0.8,1.2,1.6,2.0,3.0,4.0,5.0,6.0])
    Cns = np.array([1.0,3.0,5.0,7.0])

    Cnp,Kdp,qp = np.meshgrid(Cns,Kds,qs,indexing="ij")
    qp,Kdp,Cnp=qp.flatten(),Kdp.flatten(),Cnp.flatten()
    # getting observables
    #foldername="../data/Ne1/Feb12_2021"
    foldername="../data/Ne1/Feb24_2021"
    qdata,p2uu_ave,p2uu_err,uuc_ave,uuc_err,un2_ave,un2_err = get_phase_related_data(foldername,Cns,Kds)
    p2uu_ave,uuc_ave = p2uu_ave.flatten(),uuc_ave.flatten()
    #print(np.transpose([Cnp,Kdp,qp])[:200])
    #print("qp[0][1]",qp[0][1])
    #print("Kdp",Kdp)
    #print("Cnp",Cnp)
    iso_bool = p2uu_ave<0.5
    ani_bool = p2uu_ave>=0.5
    q_iso,Kd_iso,Cn_iso=qp[iso_bool],Kdp[iso_bool],Cnp[iso_bool]
    q_ani,Kd_ani,Cn_ani=qp[ani_bool],Kdp[ani_bool],Cnp[ani_bool]
    nem_bool=np.logical_and(ani_bool,uuc_ave<0.1)
    cho_bool=np.logical_and(ani_bool,uuc_ave>=0.1)
    q_nem,Kd_nem,Cn_nem = qp[nem_bool],Kdp[nem_bool],Cnp[nem_bool]
    q_cho,Kd_cho,Cn_cho = qp[cho_bool],Kdp[cho_bool],Cnp[cho_bool]

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1))
    #fig = plt.figure(figsize=plt.figaspect(2.))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    #axconf = plt.subplot2grid((1, 3), (0, 0))
    #ax = plt.subplot2grid((1, 3), (0, 1),colspan=2)
    #ax = plt.subplot2grid((1, 3), (0, 1),colspan=2,projection='3d',proj_type="ortho")
    ax = plt.subplot2grid((1,1), (0, 0),projection='3d')
    axconf = ax.inset_axes([-0.1, 0.0, 0.3, 1])

    #plot configuration
    axconf.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)

    filename,rot = [],[]
    filename.append(foldername+"/State_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd0.0_q0.0_Cn3.0_kard0.0.txt")
    rot.append([0.0,-1]) # isotropic
    filename.append(foldername+"/State_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd3.2_q0.9_Cn3.0_kard0.0.txt")
    rot.append([0.0,-1]) # smectic-A
    filename.append(foldername+"/State_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd3.2_q1.5_Cn3.0_kard0.0.txt")
    rot.append([-0.2,-1.2]) # cholesteric

    for i in range(3):
        ax_config_xy(axconf,filename[i],Color="gray",LineWidth=LineWidth*0.65,rotation=rot[i],xshift=0,yshift=26*i,bead=1,mesh=1,rod=1,d=0.65)
        #axconf.plot([-20],[20*i],marker=markers[i],color=mcolors[i],ms=msize)
        #axconf.text(-15,20*i,mlabels[i],fontsize=FontSize)
    markers=["+","o","s"]
    mcolors=["green","tomato","blue"]
    mlabels=["isotropic","smectic-A","cholesteric"]
    msize=15

    ax.set_box_aspect((1, 1, 1.7))
    #ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1, 1, 1.5, 1]))
    ax.scatter(q_cho,Kd_cho,Cn_cho,marker=markers[2],s=msize,color=mcolors[2],facecolors="None",label=mlabels[2])
    ax.scatter(q_nem,Kd_nem,Cn_nem,marker=markers[1],s=msize,color=mcolors[1],facecolors="None",label=mlabels[1])
    ax.scatter(q_iso,Kd_iso,Cn_iso,marker=markers[0],s=msize,color=mcolors[0],label=mlabels[0])

    ax.tick_params(labelsize=LabelSize,pad=-4)
    ax.xaxis.set_major_locator(MultipleLocator(0.9))
    ax.xaxis.set_minor_locator(MultipleLocator(0.3))
    ax.set_xlabel(r"$k_c$",fontsize=FontSize,labelpad=-8)
    #ax.set_xlabel(r"twist",fontsize=FontSize,labelpad=-8)
    ax.yaxis.set_major_locator(MultipleLocator(1.2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.4))
    ax.set_ylabel(r"$\epsilon_{LL}$",fontsize=FontSize,labelpad=-10)
    #ax.set_ylabel(r"$\epsilon$",fontsize=FontSize,labelpad=-10)
    ax.zaxis.set_major_locator(MultipleLocator(2))
    ax.set_zticks([1,3,5,7])
    ax.set_zlabel(r"$C$",fontsize=FontSize,labelpad=-10)
    #ax.set_zlabel(r"tilt coupling",fontsize=FontSize,labelpad=-10)
    ax.view_init(elev=20., azim=-40)
    ax.legend(loc="center left",bbox_to_anchor=(-0.12, 0.6), frameon=0,labelspacing=4.5,handletextpad=0.1,fontsize=FontSize)
    plt.tight_layout()
    #plt.savefig("test_phase_diagram.pdf")
    plt.savefig("figures/phase_diagram.pdf")
    #plt.show()


