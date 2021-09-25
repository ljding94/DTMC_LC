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

def mean_filter(heat,ns):
    # smooth the heat map by averaging with neighbors bases on ns
    mean_heat = []
    for i in range(len(heat)):
        mean_heat.append(heat[i])
        nei_count=0
        for j in range(len(ns[0])):
            if(ns[i,j]!=-1):
                mean_heat[i]+=heat[int(ns[i,j])]
                nei_count+=1
        mean_heat[i]/=1+nei_count
    return mean_heat


def config_plot_xyz(filename,mesh=0,rod=1,cvt_map="",cmap_smooth=0,tag="", Format="pdf",lim=15,fix_index=None):
    print("plotting",filename)
    ftail = "_xyz"
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    x,y,z,sx,sy,sz,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:14]
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
        #color mapping
    cmap = cm.get_cmap("jet_r")
    if(cvt_map=="Mean"):
        ftail+="_mmap"
        heat=dA*d2H*d2H
        print("mean curvature heat",heat)
        for m in range(cmap_smooth):
            heat = mean_filter(heat,ns)
        norm=Normalize(vmin=0,vmax=0.7)
        ax_xy.scatter(x,y,c=cmap(norm(heat)))
        ax_zx.scatter(z,x,c=cmap(norm(heat)))
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=norm)
        sm.set_array([])
        cbar=plt.colorbar(sm)
        #cbar=plt.colorbar(sm, ticks=[0,0.5,0.7])
        #cbar.ax.set_yticklabels(["0","0.5","0.7"])
    elif(cvt_map=="Gaussian"):
        ftail+="_gmap"
        heat = dAK
        print("Gaussian curvature heat",heat)
        for m in range(cmap_smooth):
            heat = mean_filter(heat,ns)
        norm=Normalize(vmin=-0.1,vmax=0.0)
        ax_xy.scatter(x,y,c=cmap(norm(heat)))
        ax_zx.scatter(z,x,c=cmap(norm(heat)))
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=norm)
        sm.set_array([])
        cbar=plt.colorbar(sm)
        #cbar=plt.colorbar(sm, ticks=[-0.1,-0.05,0])
        #cbar.ax.set_yticklabels(["-0.1","-0.05","0"])
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


    snu=sx*nu[0]+sy*nu[1]+sz*nu[2]
    snu=sz
    #print(np.sqrt(un2)-np.absolute(snu))
    #deg = np.arccos(np.absolute(snu))
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    deg = np.arccos(np.sqrt(un2))
    cmap = cm.get_cmap("jet_r")
    if(rod):
        for i in range(len(x)):
            ax_xy.plot([x[i]-0.5*d*sx[i],x[i]+0.5*d*sx[i]],[y[i]-0.5*d*sy[i],y[i]+0.5*d*sy[i]],"-",linewidth=1.5,color=cmap(norm(deg[i])))

            ax_zx.plot([z[i]-0.5*d*sz[i],z[i]+0.5*d*sz[i]],[x[i]-0.5*d*sx[i],x[i]+0.5*d*sx[i]],"-",linewidth=1.5,color=cmap(norm(deg[i])))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar=plt.colorbar(sm, ticks=[0,0.25*np.pi,0.5*np.pi])
        cbar.ax.set_yticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"])
    #plot fixed bead (if have)
    if fix_index:
        ax_xy.plot([x[fix_index[0]], x[fix_index[1]]], [y[fix_index[0]], y[fix_index[1]]], marker="o",linestyle="None", color="purple")
        ax_zx.plot([z[fix_index[0]], z[fix_index[1]]], [x[fix_index[0]], x[fix_index[1]]], marker="o",linestyle="None", color="purple")
    ax_xy.set_xlim(x_min-2, x_max+2)
    ax_xy.set_ylim(y_min-2, y_max+2)
    ax_xy.set_title("XY  "+tag+" smooth=%d"%cmap_smooth, fontsize=25)
    ax_zx.set_xlim(z_min-2, z_max+2)
    ax_zx.set_ylim(x_min-2, x_max+2)
    ax_zx.set_title("ZX")
    #ax_xy.legend(title=tag)
    # plt.savefig(filename[:-4] + "_xy.pdf", dpi=300, format="pdf")
    plt.savefig(filename[:-4] + ftail+"."+Format, dpi=100,
                format=Format, bbox_inches='tight',transparent=False)
    plt.close()

    #plot the tilt distribution based on only one configuration
    plt.figure()
    plt.hist(deg,bins=50,histtype="step",density=True,label="hist")
    theta=np.linspace(0,0.5*np.pi,100)
    plt.plot(theta,np.sin(theta),label=r"$\sin(\theta)$")
    plt.xlabel(r"$\theta$")
    plt.ylabel("count")
    plt.legend()
    plt.title(filename[-17:-12])
    plt.savefig(filename[:-4]+"tilt_hist."+Format, dpi=100,
                format=Format, bbox_inches='tight',transparent=False)
    plt.close()

def config_plot_xyz_seq(filename,Seq):
    for i in range(Seq):
        config_plot_xyz(filename[:-4]+"_%d.txt"%i,Format="png")

def config_plot3D(filename,mesh=0,rod=0,cvt_map="",cmap_smooth=0):
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    #x,y,z,sx,sy,sz,enum, en0, en1 = data[5:14]
    x,y,z,sx,sy,sz,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:14]
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
    cmap = cm.get_cmap("jet_r")
    if(cvt_map=="Mean"):
        #ftail+="_mmap"
        #norm=Normalize(vmin=0,vmax=0.5*np.pi)
        heat = dA*d2H
        for m in range(cmap_smooth):
            heat = mean_filter(heat,ns)
        ax.scatter3D(x,y,z,c=cmap(heat))
        sm = plt.cm.ScalarMappable(cmap=cmap)
        sm.set_array([])
        cbar=plt.colorbar(sm)

    # director
    nu=find_nu(sx,sy,sz)
    print(nu,np.dot(nu,nu))
    x_ave,y_ave,z_ave = np.average(x),np.average(y),np.average(z)
    D_ave = 2
    #ax.plot3D([x_ave-D_ave*nu[0],x_ave+D_ave*nu[0]],[y_ave-D_ave*nu[1],y_ave+D_ave*nu[1]],[z_ave-D_ave*nu[2],z_ave+D_ave*nu[2]],"-",linewidth=3.0,color="k")
    snu=sx*nu[0]+sy*nu[1]+sz*nu[2]
    #deg = np.arccos(np.absolute(snu))
    #deg = np.arccos(np.sqrt(un2))
    cmap = cm.get_cmap("jet_r")
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    deg = np.arccos(np.sqrt(un2))
    if(rod):
        for i in range(len(sx)):
            ax.plot3D([x[i]-0.5*d*sx[i],x[i]+0.5*d*sx[i]],[y[i]-0.5*d*sy[i],y[i]+0.5*d*sy[i]],[z[i]-0.5*d*sz[i],z[i]+0.5*d*sz[i]],"-",color=cmap(norm(deg[i])))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar=plt.colorbar(sm, ticks=[0,0.25*np.pi,0.5*np.pi])
        cbar.ax.set_yticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"])
    ax.set_xlim(-max_range_half,max_range_half)
    ax.set_ylim(-max_range_half,max_range_half)
    ax.set_zlim(-max_range_half,max_range_half)
    #ax.scatter3D(x, y, z, s=[1 for i in range(len(x))])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.show()
    #plt.savefig(filename[:-4] + "_3D.png", dpi=300)
    plt.close()

def config_nu2_dis(filename,bin_num=20):
    ftail = "_dis"
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    x,y,z,sx,sy,sz,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:14]
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.hist(un2,bins=bin_num,histtype="step",density=True)
    ax1.hist(dA*d2H*d2H,bins=bin_num,histtype="step",density=True)
    ax2.scatter(un2,dA*d2H*d2H)
    plt.show()

def autocorrelation_plot(rho,tau_int,savefile):
    t = np.linspace(0,1000,1000)
    plt.figure()
    plt.plot(range(1000),rho[:1000],"d")
    plt.plot(t,np.exp(-t/tau_int),"--")
    plt.savefig(savefile,format="pdf",transparent=True)
    plt.close()
