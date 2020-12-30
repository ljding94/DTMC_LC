import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import Circle
from matplotlib import cm
from matplotlib.colors import Normalize


def find_nu(ux,uy,uz):
    Q=np.array([[1.5*np.average(ux*ux)-0.5,np.average(ux*uy),1.5*np.average(ux*uz)],[0,1.5*np.average(uy*uy)-0.5,np.average(uy*uz)],[0,0,1.5*np.average(uz*uz)-0.5]])
    Q[1,0]=Q[0,1]
    Q[2,0]=Q[0,2]
    Q[2,1]=Q[1,2]
    w,v=np.linalg.eig(Q)
    w_max=np.max(w)
    for i in range(len(w)):
        if(w[i]==w_max):
            return np.transpose(v)[i]


def config_plot_xyz(filename,mesh=False,rod=True,tag="", Format="pdf",lim=15,fix_index=None):
    print("plotting",filename)
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    x,y,z,sx,sy,sz,dA,I2H,ds,dAK,un2,enum, en0, en1 = data[:14]
    d=1
    #sx,sy,sz=d*sx,d*sy,d*sz
    #x,y,z,sx,sy,sz, enum, en0, en1 = data[5:14]
    x_min, x_max = np.min(x),np.max(x)
    y_min, y_max = np.min(y),np.max(y)
    z_min, z_max = np.min(z),np.max(z)
    alpha_xy = 0.9*(z-z_min+0.1)/(z_max-z_min+0.1)+0.1
    alpha_zx = 0.9*(y-y_min+0.1)/(y_max-y_min+0.1)+0.1
    ns = np.transpose(data[14:])
    #ns = np.transpose(data[14:])
    ens = np.array([en0, en1])
    fig = plt.figure(figsize=(10, 5))
    #ax_xy = fig.add_subplot(111, aspect="equal")
    ax_xy = fig.add_subplot(121, aspect="equal")
    ax_zx = fig.add_subplot(122, aspect="equal")
    # bulk bond

    # track bead ind #
    if(0):
        for i in range(len(x)):
            pass
            ax_xy.annotate(i, (x[i], y[i]), fontsize=5)

    if(mesh):
        bonds = []
        for i in range(len(ns)):
            for j in range(len(ns[0])):
                if ns[i, j] != -1:
                    if(i<ns[i,j]):
                        bonds.append((i,int(ns[i, j])))
                    else:
                        bonds.append((int(ns[i, j]),i))
        bonds = set(bonds)
        for bond in bonds:
            pass
            a,b = bond
            ax_xy.plot([x[a],x[b]], [y[a],y[b]], color="tomato",alpha=alpha_xy[a])
            ax_zx.plot([z[a],z[b]], [x[a],x[b]], color="tomato",alpha=alpha_xy[a])
    # edge bond
    ecolors = ["blue","purple","green"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax_xy.plot([x[j], x[int(ens[i, j])]], [
                    y[j], y[int(ens[i, j])]], "-",linewidth=0.5, color=ecolors[int(enum[j])], alpha=alpha_xy[j])
                ax_zx.plot([z[j], z[int(ens[i, j])]],[x[j], x[int(ens[i, j])]], "-",linewidth=0.5, color=ecolors[int(enum[j])], alpha=alpha_zx[j])
    # spin vector
    nu=find_nu(sx,sy,sz)
    #nu=[0,0,1]
    x_ave,y_ave,z_ave = np.average(x),np.average(y),np.average(z)
    D_ave = 3
    ax_xy.plot([x_ave-D_ave*nu[0],x_ave+D_ave*nu[0]],[y_ave-D_ave*nu[1],y_ave+D_ave*nu[1]],"-",linewidth=3.0,color="k")
    ax_zx.plot([z_ave-D_ave*nu[2],z_ave+D_ave*nu[2]],[x_ave-D_ave*nu[0],x_ave+D_ave*nu[0]],"-",linewidth=3.0,color="k")

    #snu=sx*nu[0]+sy*nu[1]+sz*nu[2]
    #deg = np.arccos(np.absolute(snu))
    deg = np.arccos(np.absolute(un2))
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    cmap = cm.get_cmap("jet_r")
    if(rod):
        for i in range(len(x)):
            ax_xy.plot([x[i]-0.5*d*sx[i],x[i]+0.5*d*sx[i]],[y[i]-0.5*d*sy[i],y[i]+0.5*d*sy[i]],"-",linewidth=1.5,color=cmap(deg[i]))
            ax_zx.plot([z[i]-0.5*d*sz[i],z[i]+0.5*d*sz[i]],[x[i]-0.5*d*sx[i],x[i]+0.5*d*sx[i]],"-",linewidth=1.5,color=cmap(deg[i]))
    #plot fixed bead (if have)
    if fix_index:
        ax_xy.plot([x[fix_index[0]], x[fix_index[1]]], [y[fix_index[0]], y[fix_index[1]]], marker="o",linestyle="None", color="purple")
        ax_zx.plot([z[fix_index[0]], z[fix_index[1]]], [x[fix_index[0]], x[fix_index[1]]], marker="o",linestyle="None", color="purple")
    ax_xy.set_xlim(x_min-2, x_max+2)
    ax_xy.set_ylim(y_min-2, y_max+2)
    ax_xy.set_title("XY  "+tag, fontsize=25)
    ax_zx.set_xlim(z_min-2, z_max+2)
    ax_zx.set_ylim(x_min-2, x_max+2)
    ax_zx.set_title("ZX")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar=plt.colorbar(sm, ticks=[0,0.25*np.pi,0.5*np.pi])
    cbar.ax.set_yticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"])
    ax_xy.legend(title=tag)
    # plt.savefig(filename[:-4] + "_xy.pdf", dpi=300, format="pdf")
    plt.savefig(filename[:-4] + "_xyz."+Format, dpi=100,
                format=Format, bbox_inches='tight',transparent=False)
    plt.close()

def config_plot_xyz_seq(filename,Seq):
    for i in range(Seq):
        config_plot_xyz(filename[:-4]+"_%d.txt"%i,Format="png")


def config_plot3D(filename,mesh=0,rod=0):
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    #x,y,z,sx,sy,sz,enum, en0, en1 = data[5:14]
    x,y,z,sx,sy,sz,dA,I2H,ds,dAK,un2,enum, en0, en1 = data[:14]
    d=2
    #sx,sy,sz=d*sx,d*sy,d*sz
    x,y,z=x-np.average(x),y-np.average(y),z-np.average(z)
    x_min, x_max = np.min(x),np.max(x)
    y_min, y_max = np.min(y),np.max(y)
    z_min, z_max = np.min(z),np.max(z)
    max_range_half = max([x_max-x_min,y_max-y_min,z_max-z_min])*0.5
    alpha_xy = 0.9*(z-z_min+0.1)/(z_max-z_min+0.1)+0.1
    alpha_zx = 0.9*(y-y_min+0.1)/(y_max-y_min+0.1)+0.1
    ns = np.transpose(data[14:])
    ens = np.array([en0, en1])
    fig = plt.figure(figsize=(5, 5))
    ax = plt.axes(projection="3d")

    if(mesh):
        for i in range(len(ns)):
            for j in range(len(ns[0])):
                if ns[i, j] != -1:
                    pass
                    ax.plot3D([x[i], x[int(ns[i, j])]], [y[i], y[int(ns[i, j])]], [z[i], z[int(ns[i, j])]], "-",  color="tomato")

    ecolors = ["blue","purple","green"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax.plot3D([x[j], x[int(ens[i, j])]], [
                    y[j], y[int(ens[i, j])]], [
                    z[j], z[int(ens[i, j])]], "-", color=ecolors[int(enum[j])], alpha=0.7)
    # director
    nu=find_nu(sx,sy,sz)
    print(nu,np.dot(nu,nu))
    x_ave,y_ave,z_ave = np.average(x),np.average(y),np.average(z)
    D_ave = 2
    #ax.plot3D([x_ave-D_ave*nu[0],x_ave+D_ave*nu[0]],[y_ave-D_ave*nu[1],y_ave+D_ave*nu[1]],[z_ave-D_ave*nu[2],z_ave+D_ave*nu[2]],"-",linewidth=3.0,color="k")
    snu=sx*nu[0]+sy*nu[1]+sz*nu[2]
    deg = np.arccos(np.absolute(snu))
    cmap = cm.get_cmap("jet_r")
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    #deg = np.arccos(np.sqrt(un2))
    if(rod):
        for i in range(len(sx)):
            ax.plot3D([x[i]-0.5*d*sx[i],x[i]+0.5*d*sx[i]],[y[i]-0.5*d*sy[i],y[i]+0.5*d*sy[i]],[z[i]-0.5*d*sz[i],z[i]+0.5*d*sz[i]],"-",color=cmap(norm(deg[i])))
    #ax.set_xlim(-t*max_xy, t*max_xy)
    #ax.set_ylim(-t*max_xy, t*max_xy)
    #ax.set_zlim(-t*max_xy, t*max_xy)
    ax.set_xlim(-max_range_half,max_range_half)
    ax.set_ylim(-max_range_half,max_range_half)
    ax.set_zlim(-max_range_half,max_range_half)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar=plt.colorbar(sm, ticks=[0,0.25*np.pi,0.5*np.pi])
    cbar.ax.set_yticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"])

    #ax.scatter3D(x, y, z, s=[1 for i in range(len(x))])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.show()
    #plt.savefig(filename[:-4] + "_3D.png", dpi=300)
    plt.close()

def O_kar_lam_MCstep_plot(Ls, i2Hs, Ens, savefile):
    MCstep = np.arange(len(Ls))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(MCstep, Ls, "d", label=r"$L\sigma_0$")
    ax.plot(MCstep, i2Hs, "s", label=r"$\int (2H)^2 dA$")
    ax.plot(MCstep, Ens, "s", label=r"$E/k_BT$")
    ax.set_ylim(0,500)
    ax.set_xlabel("MC steps")
    plt.legend()
    plt.savefig(savefile,format="pdf",transparent=True)
    plt.close()

def O_kappa_En_dis_lam_plot(Ens,lams,N,savefile):
    dlamN = int(len(lams)/N)
    plt.figure(figsize=(4,2*dlamN))
    for i in range(dlamN):
        plt.subplot(dlamN,1,i+1)
        Ens_p = Ens[i::dlamN]
        lams_p = lams[i::dlamN]
        for j in range(len(Ens_p)):
            plt.hist(Ens_p[j],bins=100,histtype="step",label=r"$\lambda=%.1f$"%lams_p[j])
        plt.legend()
    plt.xlabel(r"$E/k_B T$")
    plt.ylabel("distribution")
    plt.tight_layout()
    plt.savefig(savefile)
    plt.close()


def autocorrelation_plot(rho,tau_int,savefile):
    t = np.linspace(0,1000,1000)
    plt.figure()
    plt.plot(range(1000),rho[:1000],"d")
    plt.plot(t,np.exp(-t/tau_int),"--")
    plt.savefig(savefile,format="pdf",transparent=True)
    plt.close()
