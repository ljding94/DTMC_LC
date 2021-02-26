import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from cal import *
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
    for i in range(len(ms)):
        # iterating over all m
        lams_all.append([])
        Emin_all.append([])
        r1min_all.append([])
        qthetamin_all.append([])
        for j in range(len(qs)):
            # iterating over all q
            filename=foldername+"/Emin_lams_Cn%.0f_m%d_q%.1f.txt"%(Cn,ms[i],qs[j])
            data = np.loadtxt(filename,delimiter=",",skiprows=1,unpack=True)

            lams_all[i].append(data[0])
            Emin_all[i].append(data[1])
            r1min_all[i].append(data[2])
            qthetamin_all[i].append(data[3])
    mcolors=["red","green","blue"]
    qmin,qmax=np.min(qs),np.max(qs)
    qalphas = 0.8*(qs-qmin)/(qmax-qmin)+0.2
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 3*0.6))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    ax = plt.subplot2grid((3,1), (0, 0))
    axr1 = plt.subplot2grid((3,1), (1, 0),sharex=ax)
    axqtheta = plt.subplot2grid((3,1), (2, 0),sharex=ax)
    #plot disk related
    #for j in range(len(qs)):
        #lams = lams_all[0][j]
        #ax.plot(lams,E_Disk(lams,Cn,qs[j]),color="black",alpha=qalphas[j])
    #plot enneper related
    ax.plot(lams_all[0][0],np.zeros(len(lams_all[0][0])),"--k")
    for i in range(len(ms)):
        for j in range(len(qs)):
            lams = lams_all[i][j]
            ax.plot(lams,Emin_all[i][j]-E_Disk(lams,Cn,qs[j]),color=mcolors[i],alpha=qalphas[j],label="m=%d,q=%.1f"%(ms[i],qs[j]))
            axr1.plot(lams_all[i][j],r1min_all[i][j],color=mcolors[i],alpha=qalphas[j])
            axqtheta.plot(lams_all[i][j],qthetamin_all[i][j],color=mcolors[i],alpha=qalphas[j])
    ax.set_ylabel(r"$\Delta E$")
    axr1.set_ylabel(r"$r_1$")
    axqtheta.set_ylabel(r"$\theta$")
    axqtheta.set_xlabel(r"$\lambda$")
    plt.tight_layout(pad=0.5)
    plt.savefig("Emin_lam.pdf")
    plt.close()






