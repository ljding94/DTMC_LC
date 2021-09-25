import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
from matplotlib import colors
from matplotlib.colors import Normalize
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import numpy as np
from scipy import special
from matplotlib.colors import ListedColormap
from config import *


def tilt_disk(q0):
    # tilt of disk membrane
    #if(q0==0):
    #    un2=1
    #else:
    un2 = 0.5*(1+special.j1(2*q0)/q0)
    tilt = np.pi*(1-un2)
    # checked by plotting and compare with mathematica
    return tilt

def peri_disk():
    # perimeter of disk, which is 2pi apparently
    return 2*np.pi

def E_Disk(lam,Cn,q):
    tiltdisk=tilt_disk(q)
    peridisk=peri_disk()
    return lam*peridisk+0.5*Cn*tiltdisk


def get_minms(foldername,Cn,ms,qs,dq,lams,ncut=-1):
    #qp,lamp = np.meshgrid(qs,lams[:ncut],indexing="ij")
    mp = [0]
    Emin = []
    for i in range(len(qs)):
        # iterating over all m
        Emin.append([E_Disk(lams[:ncut],Cn,qs[i])])
        for j in range(len(ms)):
            mp.append(ms[j])
            # iterating over all q
            filename=foldername+"/Emin_lams_Cn%.0f_m%d"%(Cn,ms[j])+"_q%.*f.txt"%(dq,qs[i])
            data = np.loadtxt(filename,delimiter=",",skiprows=1,unpack=True)
            Emin[i].append(data[1][:ncut])
    minms = np.argmin(Emin,axis=1)
    return minms
def get_Emin(foldername,Cn,m,qs,dq,lam):
    Emin = []
    nlam=-1
    for i in range(len(qs)):
        # iterating over q
        filename=foldername+"/Emin_lams_Cn%.0f_m1"%Cn+"_q%.*f.txt"%(dq,qs[i])
        data = np.loadtxt(filename,delimiter=",",skiprows=1,unpack=True)
        if(nlam==-1):
            nlam=np.where(data[0]==lam)
        Emin.append(data[1][nlam])
    return np.array(Emin)

def Emin_phase_diagram(LineWidth, FontSize, LabelSize):
    Cn=100
    ms=[1,2]
    # plotting the general phase diagram
    foldername = "../data/pydata/Feb27_2021"
    qs=np.arange(1.00,5.01,0.05)
    lams = np.arange(0.000,5.01,0.05)
    ncut=-1
    qp,lamp = np.meshgrid(qs,lams[:ncut],indexing="ij")
    mp = [0]
    for j in range(len(ms)):
        mp.append(ms[j])
    #qp,lamp = qp.flatten(),lamp.flatten()
    minms=get_minms(foldername,Cn,ms,qs,2,lams,ncut)
    minms_masked = [np.ma.masked_outside(minms,mp[i]-0.5,mp[i]+0.5) for i in range(len(mp))]

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1*1.2))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    axconfENP = plt.subplot2grid((3,2), (0, 0))
    axconfDisk = plt.subplot2grid((3,2), (0, 1))
    ax = plt.subplot2grid((3,2), (1, 0),colspan=2,rowspan=2)

    # plot demo configurations
    # data : q=2.5 lam=3.000000,Emin=104.401054,r1opt=1.201352,qtheta=0.785362
    ax_numeric_config_tilt(axconfENP, m=1, r1=1.201352, q=2.5, qtheta=0.785362)
    ax_numeric_config_tilt(axconfDisk, m=0, r1=1.201352, q=2.5, qtheta=0.785362)

    #axzoom = ax.inset_axes(0.6,0.6,0.3,0.3)
    mcolors=["tomato","limegreen","royalblue"]
    mhatchs=["//","++","xx"]

    cmap =  colors.ListedColormap(["tomato","limegreen","royalblue"])
    #cmap =  colors.ListedColormap(mcolors)
    bounds=np.arange(mp[0]-0.5,mp[-1]+0.51,1.0)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    none_map = ListedColormap(['none'])
    for i in range(len(mp)):
        #ax.pcolor(np.transpose(lamp), np.transpose(qp), np.transpose(minms),cmap=cmap,shading="auto")
        #ax.pcolor(np.transpose(lamp), np.transpose(qp), np.transpose(minms_masked[i]),cmap=cmap,shading="auto")
        ax.pcolor(np.transpose(lamp), np.transpose(qp), np.transpose(minms_masked[i]),cmap=none_map, hatch=mhatchs[i],edgecolor=mcolors[i],lw=0,shading="auto")
        #ax.contour(np.transpose(lamp), np.transpose(qp), np.transpose(minms),2,colors="k",linewidths=1)
    #cbar = plt.colorbar(im,ax=ax, boundaries=bounds, ticks=[0, 1, 2])
    #cbar.set_label(r"$m$")
    ax.set_xlabel(r"$q$")
    ax.set_ylabel(r"$\lambda$")
    plt.tight_layout(pad=0.5)
    plt.savefig("figures/phase_Emin.pdf")
    plt.close()


def Emin_fill_between_plot(LineWidth, FontSize, LabelSize):
    print("plotting q lam phase diagram using fill between method")
    Cn=100
    ms=[1]
    # plotting the general phase diagram
    #foldername = "../data/pydata/Feb27_2021"
    foldername = "../data/pydata/Mar13_2021"
    #qs=np.arange(1.00,5.01,0.05)
    qs=np.arange(1.50,4.70,0.05)
    lams = np.arange(0.000,5.01,0.05)
    ncut=-1
    qp,lamp = np.meshgrid(qs,lams[:ncut],indexing="ij")
    mp = [0]
    for j in range(len(ms)):
        mp.append(ms[j])
    minms=get_minms(foldername,Cn,ms,qs,2,lams,ncut)
    # tracing the critical curves for m=2
    qc,lamc = qs,lams[:ncut]
    lamc /= Cn
    lam_m0l=[] # lower edge of m=0 region
    lam_m2u=[] # upper edge of m=2 region
    for i in range(len(qc)):
        if(minms[i][0]==0):
            lam_m0l.append(lamc[0])
            lam_m2u.append(lamc[0])
            continue
        elif(minms[i][0]==1):
            lam_m2u.append(lamc[0])
        for j in range(len(lamc)-1):
            if(minms[i][j]==2 and minms[i][j+1]==0):
                lam_m0l.append((lamc[j]+lamc[j+1])/2)
                lam_m2u.append((lamc[j]+lamc[j+1])/2)

            if(minms[i][j]==2 and minms[i][j+1]==1):
                lam_m2u.append((lamc[j]+lamc[j+1])/2)
            elif(minms[i][j]==1 and minms[i][j+1]==0):
                lam_m0l.append((lamc[j]+lamc[j+1])/2)

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1*0.8))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    ax = plt.subplot2grid((1,1), (0, 0))
    colors=["tomato","limegreen","royalblue"]
    #ax.plot(qc,lam_m2u,color=colors[1])
    #ax.plot(qc,lam_m0l,color=colors[0])

    ax.fill_betweenx(qc,lamc[-1]*np.ones(len(qc)),lam_m0l,color=colors[0],hatch="///",facecolor="None",label="disk",linewidth=LineWidth)
    ax.plot(lam_m2u,qc,"-k",linewidth=LineWidth)
    ax.plot(lam_m0l,qc,"-k",linewidth=LineWidth)
    ax.fill_betweenx(qc,lam_m2u,lam_m0l,color=colors[1],hatch="xxx",facecolor="None",label=r"$m=1$",linewidth=LineWidth)
    ax.fill_betweenx(qc,np.zeros(len(qc)),lam_m2u,color=colors[2],hatch="+++",facecolor="None",label=r"$m=2$",linewidth=LineWidth)

    ax.tick_params(which="both",direction="in", top="on", right="on", labelsize=LabelSize)
    ax.xaxis.set_major_locator(MultipleLocator(0.01))
    ax.xaxis.set_minor_locator(MultipleLocator(0.002))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_ylim(qc[0],qc[-1])
    ax.set_ylabel(r"$q (A/\pi)^{1/2}$",fontsize=FontSize)
    #ax.set_ylabel(r"$q$",fontsize=FontSize)
    ax.set_xlim(lamc[0],lamc[-1])
    ax.set_xlabel(r"$\lambda(\pi/A)^{1/2} /C$",fontsize=FontSize)
    #ax.set_xlabel(r"$\lambda/C$",fontsize=FontSize)
    ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3,frameon=False,fontsize=LabelSize)
    plt.tight_layout(pad=0.5)
    plt.savefig("figures/numerical_enneper.pdf")
    plt.close()


def Emin_fill_between_plot_m1(LineWidth, FontSize, LabelSize):
    print("plotting q lam phase diagram using fill between method")
    Cn=100
    ms=[1]
    # plotting the general phase diagram
    #foldername = "../data/pydata/Feb27_2021"
    foldername = "../data/pydata/Mar13_2021"
    #qs=np.arange(1.00,5.01,0.05)
    qs=np.arange(1.50,3.16,0.05)
    lams = np.arange(0.000,5.01,0.05)
    lamE=3.5
    ncut=-1
    qp,lamp = np.meshgrid(qs,lams[:ncut],indexing="ij")
    mp = [0]
    for j in range(len(ms)):
        mp.append(ms[j])
    minms=get_minms(foldername,Cn,ms,qs,2,lams,ncut)
    # tracing the critical curves for m=2
    qc,lamc = qs,lams[:ncut]
    lamc /= Cn
    lam_m0l=[] # lower edge of m=0 region
    #lam_m2u=[] # upper edge of m=2 region
    for i in range(len(qc)):
        if(minms[i][0]==0):
            lam_m0l.append(lamc[0])
            #lam_m2u.append(lamc[0])
            continue
        elif(minms[i][0]==1):
            for j in range(len(lamc)-1):
                if(minms[i][j]==1 and minms[i][j+1]==0):
                    lam_m0l.append((lamc[j]+lamc[j+1])/2)

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1*1.2))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    axconfDisk = plt.subplot2grid((3,6), (0, 0),colspan=3,projection='3d')
    axconfENP = plt.subplot2grid((3,6), (0, 3),colspan=3,projection='3d')

    ax = plt.subplot2grid((3,6), (1, 0),colspan=4,rowspan=2)
    axE = plt.subplot2grid((3,6), (1, 4),colspan=2,rowspan=2,sharey=ax)

    # plot demo configurations
    # data : q=2.5 lam=3.000000,Emin=104.401054,r1opt=1.201352,qtheta=0.785362
    ax_numeric_config_tilt(axconfENP, m=1, r1=1.201352, q=2.5, qtheta=0.785362)
    #ax_numeric_config_tilt(axconfENP, m=1, r1=1.304924, q=3.0, qtheta=0.785400)
    axconfENP.view_init(elev=54., azim=66)
    #q=3.0 lam=2.0 1.304924,0.785400
    lims=(-0.64,0.64)
    axconfENP.set_xlim(lims)
    axconfENP.set_ylim(lims)
    axconfENP.set_zlim(lims)
    axconfENP.set_frame_on(False)
    axconfENP.set_axis_off()
    ax_numeric_config_tilt(axconfDisk, m=0, r1=1.201352, q=2.5, qtheta=0.785362)
    #ax_numeric_config_tilt(axconfDisk, m=0, r1=1.304924, q=3.0, qtheta=0.785400)
    axconfDisk.set_xlim(lims)
    axconfDisk.set_ylim(lims)
    axconfDisk.set_zlim(lims)
    axconfDisk.set_frame_on(False)
    axconfDisk.set_axis_off()
    x1,y1=-0.3,0.7
    axconfDisk.text2D(x1,y1, r"(a)", fontsize=FontSize, transform=axconfDisk.transAxes)
    x1,y1=-0.3,0.7
    axconfENP.text2D(x1,y1, r"(b)", fontsize=FontSize, transform=axconfENP.transAxes)

    #ax = plt.subplot2grid((1,3), (0, 0),colspan=2)
    #axE = plt.subplot2grid((1,3), (0, 2),sharey=ax)
    colors=["tomato","royalblue"]
    ax.fill_betweenx(qc,lamc[-1]*np.ones(len(qc)),lam_m0l,color=colors[0],hatch="//",facecolor="None",label="disk",linewidth=LineWidth)
    #ax.plot(lam_m2u,qc,"-k",linewidth=LineWidth)
    ax.plot(lam_m0l,qc,"-k",linewidth=LineWidth)
    ax.fill_betweenx(qc,lamc[0]*np.ones(len(qc)),lam_m0l,color=colors[1],hatch="xx",facecolor="None",label=r"$m=1$",linewidth=LineWidth)
    #ax.fill_betweenx(qc,np.zeros(len(qc)),lam_m2u,color=colors[2],hatch="+++",facecolor="None",label=r"$m=2$",linewidth=LineWidth)

    ax.tick_params(which="both",direction="in", top="on", right="on", labelsize=LabelSize)
    ax.xaxis.set_major_locator(MultipleLocator(0.01))
    ax.xaxis.set_minor_locator(MultipleLocator(0.002))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_ylim(qc[0],qc[-1])
    ax.set_ylabel(r"$q R_0$",fontsize=FontSize)
    #ax.set_ylabel(r"$q$",fontsize=FontSize)
    ax.set_xlim(lamc[0],lamc[-1])
    ax.set_xlabel(r"$\lambda/(R_0 C)$",fontsize=FontSize)
    #ax.set_xlabel(r"$\lambda/C$",fontsize=FontSize)
    ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2,frameon=False,fontsize=LabelSize)
    x1,y1=-0.2,0.9
    ax.text(ax.get_xlim()[1]*x1+ax.get_xlim()[0]* (1-x1),  ax.get_ylim()[1]*y1+ax.get_ylim()[0]* (1-y1), r"(c)", fontsize=FontSize)

    # Emin and E_disk
    Edisk=E_Disk(lamE,Cn,qs)
    axE.plot(Edisk/Cn,qs,color=colors[0],linestyle="-",linewidth=LineWidth,label=r"disk")
    Emin=get_Emin(foldername, Cn, 1, qs, 2, lamE)
    axE.plot(Emin/Cn,qs,color=colors[1],linestyle="--",linewidth=LineWidth,label=r"$m=1$")
    axE.tick_params(which="both",direction="in", top="on", right="on",labelleft=False, labelsize=LabelSize)
    axE.set_xlabel(r"$E'/(CA)$",fontsize=FontSize)
    axE.xaxis.set_major_locator(MultipleLocator(0.1))
    axE.xaxis.set_minor_locator(MultipleLocator(0.05))
    axE.set_title(r"$\lambda/(R_0 C)=%.3f$"%(lamE/Cn),fontsize=LabelSize)
    axE.legend(loc="center left",handlelength=1.2,handletextpad=0.5, columnspacing=0.,frameon=False,fontsize=LabelSize)
    x1,y1=0.1,0.9
    axE.text(axE.get_xlim()[1]*x1+axE.get_xlim()[0]* (1-x1),  axE.get_ylim()[1]*y1+axE.get_ylim()[0]* (1-y1), r"(d)", fontsize=FontSize)

    #plt.show()
    plt.tight_layout(pad=0.5)
    plt.savefig("figures/numerical_enneper.pdf")
    plt.close()
