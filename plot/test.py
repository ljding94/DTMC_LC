import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from config import *
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


def is_isotropic(Kd):
    iso =0
    if(Kd<1.5):
        iso = 1
    return iso


def phase_diagram_test_plot(LineWidth, FontSize, LabelSize):
    qs = np.arange(0.0,2.1,0.2)
    Kds = np.arange(0.0,5.1,0.5)
    Cns = np.array([1.0,3.0,5.0,7.0])
    qp,Kdp,Cnp = np.meshgrid(qs,Kds,Cns)

    iso_bool = Kdp<1.5
    ani_bool = Kdp>=1.5
    q_iso,Kd_iso,Cn_iso=qp[iso_bool],Kdp[iso_bool],Cnp[iso_bool]
    q_ani,Kd_ani,Cn_ani=qp[ani_bool],Kdp[ani_bool],Cnp[ani_bool]
    nem_bool = q_ani*np.sqrt(Kd_ani/Cn_ani)<1.2
    q_nem,Kd_nem,Cn_nem = q_ani[nem_bool],Kd_ani[nem_bool],Cn_ani[nem_bool]
    cho_bool =q_ani*np.sqrt(Kd_ani/Cn_ani)>=1.2
    q_cho,Kd_cho,Cn_cho = q_ani[cho_bool],Kd_ani[cho_bool],Cn_ani[cho_bool]

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1))
    #fig = plt.figure(figsize=plt.figaspect(2.))
    plt.rc('text', usetex=True)
    #axconf = plt.subplot2grid((1, 3), (0, 0))
    #ax = plt.subplot2grid((1, 3), (0, 1),colspan=2)
    #ax = plt.subplot2grid((1, 3), (0, 1),colspan=2,projection='3d',proj_type="ortho")
    ax = plt.subplot2grid((1,1), (0, 0),projection='3d')
    fig.subplots_adjust(hspace=0)
    axconf = ax.inset_axes([-0.1, 0.0, 0.3, 1])

    #plot configuration
    axconf.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)

    filename,rot = [],[]
    filename.append("../data/Ne1/Feb5_2021/State_N500_Ne1_L-1_kar100_karg0.0_lam5.0_Kd0.0_q0.0_Cn2.0_kard0.0.txt")
    rot.append([0.0,-1]) # isotropic
    filename.append("../data/Ne1/Feb5_2021/State_N500_Ne1_L-1_kar100_karg0.0_lam5.0_Kd4.0_q1.2_Cn4.0_kard0.0.txt")
    rot.append([0.0,-1]) # nematic
    filename.append("../data/Ne1/Feb5_2021/State_N500_Ne1_L-1_kar100_karg0.0_lam5.0_Kd4.0_q1.6_Cn4.0_kard0.0.txt")
    rot.append([0.5,-1.2]) # cholesteric

    for i in range(3):
        ax_config_xy(axconf,filename[i],Color="gray",LineWidth=LineWidth*0.6,rotation=rot[i],xshift=0,yshift=26*i,bead=1,mesh=1,rod=1)
        #axconf.plot([-20],[20*i],marker=markers[i],color=mcolors[i],ms=msize)
        #axconf.text(-15,20*i,mlabels[i],fontsize=FontSize)
    markers=["+","o","s"]
    mcolors=["green","tomato","blue"]
    mlabels=["isotopic","nematic","cholesteric"]
    msize=15

    ax.set_box_aspect((1, 1, 1.7))
    #ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1, 1, 1.5, 1]))
    ax.scatter(q_cho,Kd_cho,Cn_cho,marker=markers[2],s=msize,color=mcolors[2],facecolors="None",label="cholesteric")
    ax.scatter(q_nem,Kd_nem,Cn_nem,marker=markers[1],s=msize,color=mcolors[1],facecolors="None",label="nematic")
    ax.scatter(q_iso,Kd_iso,Cn_iso,marker=markers[0],s=msize,color=mcolors[0],label="isotropic")

    ax.tick_params(pad=-4)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.set_xlabel(r"$q$",labelpad=-8)
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.set_ylabel(r"$K_d$",labelpad=-10)
    ax.zaxis.set_major_locator(MultipleLocator(2))
    ax.set_zticks([1,3,5,7])
    ax.set_zlabel(r"$C_n$",labelpad=-10)
    ax.view_init(elev=20., azim=-40)
    ax.legend(loc="center left",bbox_to_anchor=(-0.12, 0.6), frameon=0,labelspacing=5,handletextpad=0.1,fontsize=FontSize)
    plt.tight_layout()
    #plt.savefig("test_phase_diagram.pdf")
    plt.savefig("test_phase_diagram.pdf")
    #plt.show()


