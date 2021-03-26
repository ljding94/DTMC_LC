import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
from matplotlib import colors
from matplotlib.colors import Normalize
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from cal import *
from os import path
def xyznENP_plot(m,r1,q0,qtheta):
    # plot Enneper's surface and face normal
    fig = plt.figure()
    ax=fig.gca(projection="3d")
    r,phi=np.meshgrid(np.linspace(0,r1,100),np.linspace(0,2*np.pi,100))
    R = R_nmlz(m,r1)
    X,Y,Z=xyzENP(m,R,r,phi)
    X,Y,Z=X.flatten(),Y.flatten(),Z.flatten()
    ax.plot_trisurf(X,Y,Z,linewidth=0,shade=0,color="gray",alpha=0.5)
    #ax.plot_surface(X,Y,Z,linewidth=0,shade=0,color="gray",alpha=0.5)

    # plot face normal
    r,phi=np.meshgrid(np.linspace(0,r1,10),np.linspace(0,2*np.pi,40))
    r,phi=r.flatten(),phi.flatten()
    x,y,z=xyzENP(m,R,r,phi)
    d=0.2
    nx,ny,nz=n(m,R,r,phi)
    #ax.quiver(x-0.5*nx,y-0.5*ny,z-0.5*nz,,nx,ny,nz,length=ds,color="tomato",label=r"$n$")
    x,y,z=x.flatten(),y.flatten(),z.flatten()
    ux,uy,uz=u(q0,qtheta,x,y)
    cmap = cm.get_cmap("jet_r")
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    abs_un = np.abs(nx*ux+ny*uy+nz*uz)
    deg = np.arccos(abs_un)
    for i in range(len(x)):
        #ax.plot3D([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],[z[i]-0.5*d*uz[i],z[i]+0.5*d*uz[i]],"-",linewidth=1,color=cmap(norm(abs_un[i])),label=r"$u$")
        ax.plot3D([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],[z[i]-0.5*d*uz[i],z[i]+0.5*d*uz[i]],"-",linewidth=1,color=cmap(norm(deg[i])),label=r"$u$")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    #sm.set_array([])
    #cbar=plt.colorbar(sm, ticks=[0,0.25,0.5,0.75,1])
    cbar=plt.colorbar(sm, ticks=[0,0.25*np.pi,0.5*np.pi])
    cbar.ax.set_yticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"])
    cbar.ax.set_ylabel(r"$\arccos{|u\cdot n|}$")

    ax.axes.set_xlabel("x")
    ax.axes.set_ylabel("y")
    ax.axes.set_zlabel("z")
    xmin=1.2*X.min()
    xmax=1.2*X.max()
    ax.axes.set_xlim3d(left=xmin, right=xmax)
    ax.axes.set_ylim3d(bottom=xmin, top=xmax)
    ax.axes.set_zlim3d(bottom=xmin, top=xmax)
    plt.title("m=%d,r1=%.1f,q0=%.1f,qtheta=%.1f"%(m,r1,q0,qtheta))
    plt.savefig("config_m%d_r1%.1f_q0%.1f_qtheta%.1f.pdf"%(m,r1,q0,qtheta))



def DelE_r1_plot(lams,Cn,m,q0,qtheta,r1s,bn_r,bn_phi):
    plt.figure()
    #qthetas = np.pi*np.linspace(0,0.5,10)
    for lam in lams:
        Es=[]
        for i in range(len(r1s)):
            print(i,len(r1s))
            Es.append(DelE(lam,Cn,m,r1s[i],q0,qtheta,bn_r,bn_phi))
        plt.plot(r1s,Es,label=r"$\lambda=%.1f,(q_0,\theta)=(%.1f,%.1f)$"%(lam,q0,qtheta))
    plt.plot(r1s,np.zeros(len(r1s)),"k--")
    plt.legend()
    plt.title(r"$m=%.0f,C_n=%.1f$"%(m,Cn))
    plt.xlabel(r"$r_1$")
    plt.ylabel(r"$\Delta(\lambda\oint ds+C_n/2\int (1-(u\cdot n)^2) dA)$")
    plt.savefig("DelE_m%.0f_r1_q.pdf"%m)

def DelE_r1_qtheta_plot(lam,Cn,m,r1s,q0,qthetas,bn_r,bn_phi):
    # note: r1>0
    r1p,qthetap=np.meshgrid(r1s,qthetas)
    Es =np.zeros(np.shape(r1p))
    for i in range(len(Es)):
        for j in range(len(Es[0])):
            #Es[i][j]=E_tilt_ENP(m,r1p[i][j],q0,qthetap[i][j],bn_r,bn_phi)
            Es[i][j]=DelE(lam,Cn,m,r1p[i][j],q0,qthetap[i][j],bn_r,bn_phi)

    plt.figure()
    plt.pcolormesh(r1p,qthetap,Es,shading="auto")
    cbar=plt.colorbar()
    cbar.set_label(r"$E_{ENP}-E_{D}$")
    #plt.scatter(r1p,qthetap,c=Es,marker="s",s=1000)
    plt.title(r"$\lambda=%.1f,C_n=%.1f,m=%d,q_0=%.1f$"%(lam,Cn,m,q0))
    plt.xlabel(r"$r_1$")
    plt.ylabel(r"$\theta$")
    plt.savefig("DelE_mesh_r1_qtheta_lam%.1f_Cn%.1f_m%d_q0%.1f.pdf"%(lam,Cn,m,q0))
    plt.close()

def DelEmin_lam_q0_plot(lams,q0,Cn,ms,bn_r,bn_phi,method):
    DelEmin,r1min,qthetamin = [],[],[]
    for i in range(len(ms)):
        DelEmin.append([])
        r1min.append([])
        qthetamin.append([])
        for j in range(len(lams)):
            opt = opt_r1_qtheta(lams[j],Cn,ms[i],q0,bn_r,bn_phi,method)
            print(opt.x)
            DelEmin[i].append(opt.fun)
            r1min[i].append(opt.x[0])
            qthetamin[i].append(opt.x[1])

    fig = plt.figure(figsize=(5,8))
    axE=plt.subplot2grid((3, 1), (0, 0))
    axr1=plt.subplot2grid((3, 1), (1, 0),sharex=axE)
    axqtheta=plt.subplot2grid((3, 1), (2, 0),sharex=axE)
    for i in range(len(ms)):
        axE.plot(lams,DelEmin[i],"+--",label=r"$m=%.0f$"%ms[i])
        axr1.plot(lams,r1min[i],"+--")
        axqtheta.plot(lams,qthetamin[i],"+--")
    axE.plot(lams,np.zeros(len(lams)),"k:")
    axE.set_ylabel(r"$\Delta E_{min}$")
    axE.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False)
    axr1.set_ylabel(r"$r_1$")
    axr1.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False)
    axqtheta.set_ylabel(r"$\theta$")
    axqtheta.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True)
    axqtheta.set_xlabel(r"$\lambda$")
    axE.legend(frameon=False)
    plt.tight_layout(pad=0.5)
    #plt.savefig("DelEmin_lam_q0_Cn%.1f_m%.0f.pdf"%(Cn,m))
    plt.savefig("DelEmin_lam_q0%.1f_Cn%.1f_ms.pdf"%(q0,Cn))
    plt.close()

def Emin_lam_plot(foldername,Cn,ms,qs):
    lams_all,Emin_all,r1min_all,qthetamin_all = [],[],[],[]
    ncut=120
    for i in range(len(qs)):
        # iterating over all m
        lams_all.append([])
        Emin_all.append([])
        r1min_all.append([])
        qthetamin_all.append([])
        for j in range(len(ms)):
            # iterating over all q
            filename=foldername+"/Emin_lams_Cn%.0f_m%d_q%.2f.txt"%(Cn,ms[j],qs[i])
            data = np.loadtxt(filename,delimiter=",",skiprows=1,unpack=True)
            lams_all[i].append(data[0][:ncut])
            Emin_all[i].append(data[1][:ncut])
            r1min_all[i].append(data[2][:ncut])
            qthetamin_all[i].append(data[3][:ncut])

    # plot style setting
    #mcolors=["red","green","blue","purple"]
    #qmin,qmax=np.min(qs),np.max(qs)
    #qalphas = 0.8*(qs-qmin)/(qmax-qmin)+0.2
    #malphas = 0.4*(ms[-1]-np.array(ms))/(ms[-1]-ms[0])+0.6
    alpha_light=0.2
    mlss = ["-","--","-.",":"]
    mmarkers=["x","v","s","o"]
    cmap = cm.get_cmap("rainbow")
    norm=Normalize(vmin=qs[0],vmax=qs[-1])
    qcolors = cmap(norm(qs))
    mnorm=Normalize(vmin=ms[0],vmax=ms[-1])
    mcolors = cmap(mnorm(ms))
    legend_elements = []
    for i in range(len(qs)):
        pass
        #legend_elements.append(Line2D([0], [0],ls="None", color=qcolors[i] ,marker="o",ms=3,label=r"$q=%.1f$"%qs[i]))
    for j in range(len(ms)):
        legend_elements.append(Line2D([0], [0],ls=mlss[j], color="black", mfc="None",marker="None",label=r"$m=%d$"%ms[j]))
        #legend_elements.append(Line2D([0], [0], color=mcolors[j], mfc="None",marker=mmarkers[j],label=r"$m=%d$"%ms[j]))

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 2*0.8))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    ax = plt.subplot2grid((2,1), (0, 0))
    axr1 = plt.subplot2grid((2,1), (1, 0),sharex=ax)
    #axqtheta = plt.subplot2grid((3,1), (2, 0),sharex=ax)
    # disk related
    E_disks = []
    E_bolds, lam_bolds ,Edisk_bolds= [],[],[]
    mindiffpos = []
    for i in range(len(qs)):
        lams = lams_all[i][0]
        E_disks.append(E_Disk(lams,Cn,qs[i]))
        E_bolds.append([])
        lam_bolds.append([])
        Edisk_bolds.append([])
        Es = Emin_all[i]
        Es.append(E_disks[-1])
        minms = np.argmin(Es,axis=0)
        print(minms)
        print("np.transpose(Es)",np.transpose(Es))
        print("np.transpose(Es)[minms==1]",np.transpose(Es)[minms==4])
        for j in range(len(ms)):
            E_bolds[i].append(np.transpose(np.transpose(Es)[minms==j])[j])
            lam_bolds[i].append(lams[minms==j])
            Edisk_bolds[i].append(E_disks[i][minms==j])
            print("np.shape(E_bolds[i][j])",np.shape(E_bolds[i][j]))
            print("np.shape(lam_bolds[i][j])",np.shape(lam_bolds[i][j]))
        #ax.plot(lams,E_Disk(lams,Cn,qs[j]),color="black",alpha=qalphas[j])
    #plot enneper related
    ax.plot(lams_all[0][0],np.zeros(len(lams_all[0][0])),":k")
    for i in range(len(qs)):
        for j in range(len(ms)):
            ax.plot(lam_bolds[i][j],E_bolds[i][j]-Edisk_bolds[i][j],ls=mlss[j],color=qcolors[i],alpha=1.0)
            #ax.plot(lam_bolds[i][j],qs[i]*np.ones(len(lam_bolds[i][j])),linestyle="None",marker=mmarkers[j],color=mcolors[j],ms=3,alpha=1.0)
            lams = lams_all[i][j]
            if(i==0):
                ax.plot(lams,Emin_all[i][j]-E_disks[i],ls=mlss[j],color=qcolors[i],alpha=alpha_light,label="m=%d"%ms[j])
                #ax.plot(lams,Emin_all[i][j]-E_disks[i],ls=mlss[j],color=qcolors[i],alpha=alpha_light,label="m=%d"%ms[j])
            else:
                axplt=ax.plot(lams,Emin_all[i][j]-E_disks[i],ls=mlss[j],color=qcolors[i],alpha=alpha_light)
            if(j==0):
                #ax.plot([lams_all[i][j][mindiffpos[j]]],[Emin_all[i][j][mindiffpos[j]]-E_disks[j][mindiffpos[j]]],"X",ms=10,color=qcolors[j])
                #ax.plot([lams_all[i][j][mindiffpos[j]],lams_all[i][j][mindiffpos[j]]],[-15,Emin_all[i][j][mindiffpos[j]]-E_disks[j][mindiffpos[j]]],":",color=qcolors[j])
                axr1.plot(lams_all[i][j],r1min_all[i][j],ls=mlss[j],color=qcolors[i],label="q=%.1f"%qs[i])
            else:
                axr1.plot(lams_all[i][j],r1min_all[i][j],ls=mlss[j],color=qcolors[i])
            #axqtheta.plot(lams_all[i][j],qthetamin_all[i][j],ls=mlss[i],color=qcolors[j],alpha=malphas[i])
    ax.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False)
    ax.set_ylabel(r"$\Delta E$")
    ax.legend(handles=legend_elements,ncol=2, handlelength=1.5,handletextpad=0.5, columnspacing=0.2,frameon=False)
    axr1.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True)
    axr1.set_ylabel(r"$r_1$")
    axr1.legend()
    axr1.set_xlabel(r"$\lambda$")
    #axqtheta.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True)
    #axqtheta.set_ylabel(r"$\theta$")
    #axqtheta.set_xlabel(r"$\lambda$")
    plt.tight_layout(pad=0.5)
    plt.savefig(foldername+"/Emin_lams.pdf")
    plt.close()



def get_minms(foldername,Cn,ms,qs,dq,lams,ncut=-1):
    qp,lamp = np.meshgrid(qs,lams[:ncut],indexing="ij")
    mp = [0]
    Emin = []
    for i in range(len(qs)):
        # iterating over all m
        Emin.append([E_Disk(lams[:ncut],Cn,qs[i])])
        for j in range(len(ms)):
            mp.append(ms[j])
            # iterating over all q
            filename=foldername+"/Emin_lams_Cn%.0f_m%d"%(Cn,ms[j])+"_q%.*f.txt"%(dq,qs[i])
            if(path.exists(filename)):
                data = np.loadtxt(filename,delimiter=",",skiprows=1,unpack=True)
                Emin[i].append(data[1][:ncut])
            else:
                Emin[i].append(np.zeros(ncut))
    minms = np.argmin(Emin,axis=1)
    return minms

def phase_of_Emin_pcolormesh(foldername,Cn,ms,qs,qd,lams):
    print("phase diagram (numerical) as function of (lam,q)")
    ncut=100
    qp,lamp = np.meshgrid(qs,lams[:ncut],indexing="ij")
    mp = [0]
    for j in range(len(ms)):
        mp.append(ms[j])
    #qp,lamp = qp.flatten(),lamp.flatten()
    minms=get_minms(foldername,Cn,ms,qs,qd,lams,ncut)

    print("np.shape(minms)",np.shape(minms))
    print("np.shape(qp)",np.shape(qp))
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1*0.8))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    ax = plt.subplot2grid((1,1), (0, 0))

    cmap =  colors.ListedColormap(["tomato","limegreen","royalblue","cyan"])
    #cmap =  colors.ListedColormap(mcolors)
    bounds=np.arange(mp[0]-0.5,mp[-1]+0.51,1.0)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    minm0 = np.ma.masked_not_equal(minms,0)
    im=ax.pcolor(qp,lamp,minms,cmap=cmap,shading="auto")
    #ax.pcolor(qp,lamp, minm0, hatch='///',alpha=0.,shading="auto")
    cbar = plt.colorbar(im,ax=ax, boundaries=bounds, ticks=[0, 1, 2])
    cbar.set_label(r"$m$")
    ax.set_xlabel(r"$q$")
    ax.set_ylabel(r"$\lambda$")
    plt.tight_layout(pad=0.5)
    plt.savefig(foldername+"/phase_Emin.pdf")
    plt.close()

def phase_E_fill_between_plot(foldername,Cn,ms,qs,lams,LineWidth, FontSize, LabelSize):
    print("plotting q lam phase diagram using fill between method")
    ncut=-1
    qp,lamp = np.meshgrid(qs,lams[:ncut],indexing="ij")
    mp = [0]
    for j in range(len(ms)):
        mp.append(ms[j])
    minms=get_minms(foldername,Cn,ms,qs,2,lams,ncut)
    print("np.shape(minms)",np.shape(minms))
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
    print("np.shape(qc)",np.shape(qc))
    print("np.shape(lam_m0l)",np.shape(lam_m0l))
    ax.fill_betweenx(qc,lamc[-1]*np.ones(len(qc)),lam_m0l,color=colors[0],hatch="///",facecolor="None",label="disk",linewidth=LineWidth)
    ax.plot(lam_m2u,qc,"-k",linewidth=LineWidth)
    ax.plot(lam_m0l,qc,"-k",linewidth=LineWidth)
    ax.fill_betweenx(qc,lam_m2u,lam_m0l,color=colors[1],hatch="xxx",facecolor="None",label=r"$m=1$",linewidth=LineWidth)
    ax.fill_betweenx(qc,np.zeros(len(qc)),lam_m2u,color=colors[2],hatch="+++",facecolor="None",label=r"$m=2$",linewidth=LineWidth)

    ax.tick_params(which="both",direction="in", top="on", right="on", labelsize=LabelSize)
    ax.xaxis.set_major_locator(MultipleLocator(0.01))
    ax.xaxis.set_minor_locator(MultipleLocator(0.002))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.2))
    ax.set_ylim(qc[0],qc[-1])
    ax.set_ylabel(r"$q$",fontsize=FontSize)
    ax.set_xlim(lamc[0],lamc[-1])
    ax.set_xlabel(r"$\lambda/C$",fontsize=FontSize)
    ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3,frameon=False,fontsize=LabelSize)
    plt.tight_layout(pad=0.5)
    plt.savefig(foldername+"/numerical_enneper.pdf")
    plt.close()
