import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
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


def MC_update_plot():
    print("plotting MC update figures")

def config_xy_plot(filenames,Color, LineWidth,rotation=None,mesh=0,bead=0,rod=0,xlim=0,ylim=0):

    data = np.loadtxt(filename, skiprows=11, delimiter=",", unpack=True)
    x,y,z,sx,sy,sz,un2, enum, en0, en1 = data[5:15]
    x,y,z=x-np.average(x),y-np.average(y),z-np.average(z)
    if(rotation):
        phi,theta=rotation
        xr=x*np.cos(phi)-y*np.sin(phi)
        yr = np.cos(theta)*(x*np.sin(phi)+y*np.cos(phi))-z*np.sin(theta)
        zr = np.sin(theta)*(x*np.sin(phi)+y*np.cos(phi))+z*np.cos(theta)
        x,y,z=xr,yr,zr
    x_min, x_max = np.min(x),np.max(x)
    y_min, y_max = np.min(y),np.max(y)
    z_min, z_max = np.min(z),np.max(z)
    alpha_xy = 0.9*(z-z_min+0.1)/(z_max-z_min+0.1)+0.1
    ns = np.transpose(data[15:])
    ens = np.array([en0, en1])
    fig = plt.figure(figsize=(5, 5))
    ax_xy = fig.add_subplot(111, aspect="equal")
    # bulk bond
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
            ax_xy.plot([x[a],x[b]], [y[a],y[b]], color="gray",alpha=alpha_xy[a])

    for i in range(len(x)):
            pass
            #ax_xy.annotate(i, (x[i], y[i]), fontsize=5)
    if(bead):
        r = 0.5 * np.ones(len(x))
        for j in range(len(x)):
            xc, yc, rc = x[j], y[j], r[j]
            ax_xy.add_artist(Circle(xy=(xc, yc), radius=rc,color="gray", alpha=alpha_xy[j]))
    # edge bond
    ecolors = ["blue","purple","green"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax_xy.plot([x[j], x[int(ens[i, j])]], [
                    y[j], y[int(ens[i, j])]], "-",linewidth=1.5, color=ecolors[int(enum[j])], alpha=alpha_xy[j])
    #rotate 254
    th=1.0
    #sx[254]=sx[254]*np.cos(th)-sy[254]*np.sin(th)
    #sy[254]=sx[254]*np.sin(th)+sy[254]*np.cos(th)
    if(rod):
        nu=find_nu(sx,sy,sz)
        d=1
        snu=sx*nu[0]+sy*nu[1]+sz*nu[2]
        deg = np.arccos(np.absolute(snu))
        norm=Normalize(vmin=0,vmax=0.5*np.pi)
        cmap = cm.get_cmap("jet_r")
        for i in range(len(x)):
            ax_xy.plot([x[i]-0.5*d*sx[i],x[i]+0.5*d*sx[i]],[y[i]-0.5*d*sy[i],y[i]+0.5*d*sy[i]],"-",linewidth=5.75,color=cmap(norm(deg[i])))
    ax_xy.set_xlim(x_min-2, x_max+2)
    ax_xy.set_ylim(y_min-2, y_max+2)
    ax_xy.set_frame_on(False)
    if(xlim):
        ax_xy.set_frame_on(True)
        ax_xy.set_xlim(xlim)
        ax_xy.set_ylim(ylim)
    ax_xy.set_xticks([])
    ax_xy.set_yticks([])
    #plt.show()
    plt.savefig("xy"+tail+".png", dpi=200,format="png", bbox_inches='tight',transparent=False)
    plt.close()




def config_3d_plot(filename,rotation=None,mesh=0,bead=0,rod=0,xlim=0,ylim=0):
    data = np.loadtxt(filename, skiprows=11, delimiter=",", unpack=True)
    x,y,z,sx,sy,sz,un2, enum, en0, en1 = data[5:15]
    x,y,z=x-np.average(x),y-np.average(y),z-np.average(z)
    if(rotation):
        phi,theta=rotation
        xr=x*np.cos(phi)-y*np.sin(phi)
        yr = np.cos(theta)*(x*np.sin(phi)+y*np.cos(phi))-z*np.sin(theta)
        zr = np.sin(theta)*(x*np.sin(phi)+y*np.cos(phi))+z*np.cos(theta)
        x,y,z=xr,yr,zr
    x_min, x_max = np.min(x),np.max(x)
    y_min, y_max = np.min(y),np.max(y)
    z_min, z_max = np.min(z),np.max(z)
    alpha_xy = 0.9*(z-z_min+0.1)/(z_max-z_min+0.1)+0.1
    ns = np.transpose(data[15:])
    ens = np.array([en0, en1])
    fig = plt.figure(figsize=(5, 5))
    ax_xy = plt.axes(projection="3d")
    # bulk bond
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
            ax_xy.plot3D([x[a],x[b]], [y[a],y[b]],[z[a],z[b]], color="gray",alpha=alpha_xy[a])

    if(bead):
        r = 0.5 * np.ones(len(x))
        for j in range(len(x)):
            xc, yc, rc = x[j], y[j], r[j]
            ax_xy.add_artist(Circle(xy=(xc, yc), radius=rc,color="gray", alpha=alpha_xy[j]))
    # edge bond
    ecolors = ["blue","purple","green"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax_xy.plot([x[j], x[int(ens[i, j])]], [
                    y[j], y[int(ens[i, j])]], "-",linewidth=1.5, color=ecolors[int(enum[j])], alpha=alpha_xy[j])

    if(rod):
        nu=find_nu(sx,sy,sz)
        d=1
        snu=sx*nu[0]+sy*nu[1]+sz*nu[2]
        deg = np.arccos(np.absolute(snu))
        norm=Normalize(vmin=0,vmax=0.5*np.pi)
        cmap = cm.get_cmap("jet_r")
        for i in range(len(x)):
            ax_xy.plot([x[i]-0.5*d*sx[i],x[i]+0.5*d*sx[i]],[y[i]-0.5*d*sy[i],y[i]+0.5*d*sy[i]],"-",linewidth=1.25,color=cmap(norm(deg[i])))
        ax_xy.set_xlim(x_min-2, x_max+2)
        ax_xy.set_ylim(y_min-2, y_max+2)
        ax_xy.set_frame_on(False)

    '''
    for ii in range(1):
        ax_xy.view_init(elev=60., azim=36*ii)
        plt.tight_layout()
        plt.savefig("xyz_ani%d"%ii+".png")
    '''
    ax_xy.view_init(elev=60., azim=0)
    #plt.show()
    #plt.savefig("xy.pdf", format="pdf", bbox_inches='tight',transparent=False)
    plt.close()


def ax_config_xy(ax,filename, Color, LineWidth,rotation=None, xshift=0, yshift=0,bead=0,mesh=0,rod=1,d=1,xlim=0,ylim=0,cvt_map="",cmap_smooth=0,fix_index=None):
    print("adding 2d configuration plot to ax")
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    x,y,z,ux,uy,uz,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:14]
    x,y,z=x-np.average(x),y-np.average(y),z-np.average(z)
    # reorientation
    if(rotation):
        phi,theta=rotation
        # rotate position vector
        xr=x*np.cos(phi)-y*np.sin(phi)
        yr = np.cos(theta)*(x*np.sin(phi)+y*np.cos(phi))-z*np.sin(theta)
        zr = np.sin(theta)*(x*np.sin(phi)+y*np.cos(phi))+z*np.cos(theta)
        x,y,z=xr,yr,zr

        # rotate director
        uxr=ux*np.cos(phi)-uy*np.sin(phi)
        uyr = np.cos(theta)*(ux*np.sin(phi)+uy*np.cos(phi))-uz*np.sin(theta)
        uzr = np.sin(theta)*(ux*np.sin(phi)+uy*np.cos(phi))+uz*np.cos(theta)
        ux,uy,uz=uxr,uyr,uzr

    x_min, x_max = np.min(x),np.max(x)
    y_min, y_max = np.min(y),np.max(y)
    z_min, z_max = np.min(z),np.max(z)
    alpha_xy = 0.9*(z-z_min+0.1)/(z_max-z_min+0.1)+0.1
    if(xlim):
        zselect=z[np.logical_and(xlim[0]<x,x<xlim[1])]
        z_min, z_max = np.min(zselect),np.max(zselect)
        alpha_xy = 0.4*(z-z_min+0.1)/(z_max-z_min+0.1)+0.6

    ns = np.transpose(data[14:])
    ens = np.array([en0, en1])
    x = x + xshift
    y = y + yshift

    if(0):
        for i in range(len(x)):
            pass
            #ax_xy.annotate(i, (x[i], y[i]), fontsize=5)
    if(bead):
        r = 0.5 * np.ones(len(x))
        for j in range(len(x)):
            xc, yc, rc = x[j], y[j], r[j]
            ax.add_artist(Circle(xy=(xc, yc), linewidth=0,radius=rc,edgecolor="None",facecolor=Color, alpha=alpha_xy[j]))

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
            ax.plot([x[a],x[b]], [y[a],y[b]], color=Color,linewidth=LineWidth/4,alpha=alpha_xy[a])
        #color mapping
    cmap = cm.get_cmap("jet_r")
    if(cvt_map=="Mean"):
        ftail+="_mmap"
        heat=dA*d2H
        print("mean curvature heat",heat)
        for m in range(cmap_smooth):
            heat = mean_filter(heat,ns)
        norm=Normalize(vmin=0,vmax=0.7)
        ax_xy.scatter(x,y,c=cmap(norm(heat)))
        ax_zx.scatter(z,x,c=cmap(norm(heat)))
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=norm)
        sm.set_array([])
        cbar=plt.colorbar(sm, ticks=[0,0.5,0.7])
        cbar.ax.set_yticklabels(["0","0.5","0.7"])
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
        cbar=plt.colorbar(sm, ticks=[-0.1,-0.05,0])
        cbar.ax.set_yticklabels(["-0.1","-0.05","0"])
    # edge bond
    ecolors = ["black","purple","green","blue"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax.plot([x[j], x[int(ens[i, j])]], [
                    y[j], y[int(ens[i, j])]], "-",linewidth=LineWidth/2, color=ecolors[int(enum[j])], alpha=alpha_xy[j])
    # director
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    #norm=Normalize(vmin=0,vmax=1)
    deg = np.arccos(np.sqrt(un2))
    sqrtun2=np.sqrt(un2)
    cmap = cm.get_cmap("jet_r")
    if(rod):
        for i in range(len(x)):
            ax.plot([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],"-",linewidth=LineWidth,color=cmap(norm(deg[i])),solid_capstyle='round')
            #ax.plot([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],"-",linewidth=LineWidth,color=cmap(norm(deg[i])),solid_capstyle="butt")
    #plot fixed bead (if have)
    if fix_index:
        ax.plot([x[fix_index[0]], x[fix_index[1]]], [y[fix_index[0]], y[fix_index[1]]], marker="o",linestyle="None", color="purple")

    ax.set_aspect("equal")
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_frame_on(False)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    #ax.set_ylim(y_min-2, y_max+2)

def config_demo_plot(LineWidth, FontSize, LabelSize):
    print("plotting configuration demo plot")
    filename="../data/Ne1/Feb19_2021/State_N500_Ne1_L-1_kar100_karg0.0_lam10.0_Kd5.0_q2.0_Cn5.0_kard0.0.txt"
    #filename="../data/Ne1/Feb14_2021/State_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd9.0_q0.7_Cn5.0_kard0.0.txt" # just for NECF talk

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.4))
    #fig = plt.figure(figsize=plt.figaspect(2.))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    ax = plt.subplot2grid((1,3), (0, 0),colspan=3)
    #axzoom = plt.subplot2grid((1,3), (0, 2))
    #axzoom = ax.inset_axes([0.8, 0.6, 0.5, 0.5])
    rots=[(0.0,0.0),(0.0,-0.8),(0.0,-1.5)]
    rot=(0.0,-1.2)
    #ax_config_xy(ax,filename,Color="gray",LineWidth=LineWidth,rotation=rots[0],xshift=0,yshift=0,bead=1,mesh=1,rod=1)
    ax_config_xy(ax,filename,Color="gray",LineWidth=LineWidth,rotation=rot,xshift=0,yshift=0,bead=0,mesh=1,rod=0)
    ax_config_xy(ax,filename,Color="gray",LineWidth=LineWidth,rotation=rot,xshift=30,yshift=0,bead=1,mesh=1,rod=1,d=0.65)
    cbar=plt.colorbar(cm.ScalarMappable(norm=Normalize(vmin=0,vmax=0.5*np.pi), cmap=cm.get_cmap("jet_r")),ax=ax,ticks=[0,0.25*np.pi,0.5*np.pi])
    cbar.ax.set_yticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"],fontsize=FontSize)
    cbar.ax.tick_params(direction="in",labelsize=LabelSize)
    cbar.ax.set_title(r"$\arccos{|\vu{u}\cdot\vu{n}|}$",fontsize=FontSize)

    x1, y1 = 0.02, 0.9
    ax.text(ax.get_xlim()[1]*x1+ax.get_xlim()[0]* (1-x1),  ax.get_ylim()[1]*y1+ax.get_ylim()[0]* (1-y1), r"(a)", fontsize=FontSize)
    x1, y1 = 0.52, 0.9
    ax.text(ax.get_xlim()[1]*x1+ax.get_xlim()[0]* (1-x1),  ax.get_ylim()[1]*y1+ax.get_ylim()[0]* (1-y1), r"(b)", fontsize=FontSize)
    #ax_config_xy(axzoom,filename,Color="gray",LineWidth=LineWidth,rotation=rot,xshift=0,yshift=0,bead=1,mesh=1,rod=1,d=0.65,xlim=(11,14),ylim=(-6,-3))

    plt.tight_layout(pad=0.5)
    plt.savefig("figures/config_demo.pdf")
    plt.close()


def meron_config_demo(LineWidth, FontSize, LabelSize):
    print("plotting meron lattice configuration demo plot")
    folder = "../data/Ne1/Mar23_2021"
    Kd,Cn=3.0,3.0
    qs=[2.0,3.0,5.0]
    #filename="../data/Ne1/Feb19_2021/State_N500_Ne1_L-1_kar100_karg0.0_lam10.0_Kd5.0_q2.0_Cn5.0_kard0.0.txt"
    #filename="../data/Ne1/Feb14_2021/State_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd9.0_q0.7_Cn5.0_kard0.0.txt" # just for NECF talk

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.5))
    #fig = plt.figure(figsize=plt.figaspect(2.))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    ax = plt.subplot2grid((1,1), (0, 0))
    rots=[(0.0,0.0),(0.0,0.0),(0.0,0.0)]
    for i in range(len(qs)):
        filename=folder+"/State_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd%.1f_q%.1f_Cn%.1f_kard0.0.txt"%(Kd,qs[i],Cn)
        ax_config_xy(ax,filename,Color="gray",LineWidth=LineWidth,rotation=rots[i],xshift=30*i,yshift=5,bead=1,mesh=1,rod=1,d=0.65)
    for i in range(len(qs)):
        xl,yl=i/len(qs)+0.1,-0.05
        ax.text(ax.get_xlim()[1]*xl+ax.get_xlim()[0]* (1-xl),  ax.get_ylim()[1]*yl+ax.get_ylim()[0]*(1-yl), r"$k_c=%.1f$"%qs[i], fontsize=FontSize)
    #cbar=plt.colorbar(cm.ScalarMappable(norm=Normalize(vmin=0,vmax=0.5*np.pi), cmap=cm.get_cmap("jet_r")),ax=ax,ticks=[0,0.25*np.pi,0.5*np.pi])
    #cbar.ax.set_yticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"],fontsize=FontSize)
    #cbar.ax.tick_params(direction="in",labelsize=LabelSize)
    #cbar.ax.set_title(r"$\arccos{|\vu{u}\cdot\vu{n}|}$",fontsize=FontSize)
    #x1, y1 = 0.02, 0.9
    #ax.text(ax.get_xlim()[1]*x1+ax.get_xlim()[0]* (1-x1),  ax.get_ylim()[1]*y1+ax.get_ylim()[0]* (1-y1), r"(a)", fontsize=FontSize)
    plt.tight_layout(pad=0.5)
    plt.savefig("figures/config_meron.pdf")
    plt.close()


def shape_Enneper_plot(LineWidth, FontSize, LabelSize):
    print("plotting examples of Enneper's surfaces")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.35))
    ax1 = fig.add_subplot(1, 3, 1, projection='3d')
    ax2 = fig.add_subplot(1, 3, 2, projection='3d')
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    #ax4 = fig.add_subplot(1, 4, 4, projection='3d')
    ax_Enneper_config(ax1,1,1.0)
    ax1.view_init(elev=30., azim=-60)
    #plt.subplots_adjust(wspace = -0.5)
    ax_Enneper_config(ax2,2,1.0)
    ax2.view_init(elev=30., azim=-40)
    #plt.subplots_adjust(wspace = -0.5)
    ax_Enneper_config(ax3,3,1.0)
    ax3.view_init(elev=30., azim=-70)
    #ax_Enneper_config(ax3,2,0.5)
    #ax_Enneper_config(ax4,2,1.0)
    for ax in [ax1,ax2,ax3]:
        ax.set_frame_on(False)
        ax.set_axis_off()
    plt.tight_layout(pad=-0.5)
    #plt.show()
    plt.savefig("figures/Enneper_shape.pdf")

def xyzENP(m,R,r,phi):
    # cartesian cooredinate of Enneper's surface
    x = r*np.cos(phi)-np.power(r,2*m+1)/(2*m+1)*np.cos((2*m+1)*phi)
    x*=R
    y = -r*np.sin(phi)-np.power(r,2*m+1)/(2*m+1)*np.sin((2*m+1)*phi)
    y*=R
    z = 2*np.power(r,m+1)/(m+1)*np.cos((m+1)*phi)
    z*=R
    return np.array([x,y,z])

def R_nmlz(m,r1):
    # normalized R for A=1
    R =1/r1
    R/=np.sqrt(1+2*np.power(r1,2*m)/(1+m)+np.power(r1,4*m)/(1+2*m))
    return R

def nENP(m,R,r,phi):
    # face normal vector of Enneper's surface
    rm,r2m=np.power(r,m),np.power(r,2*m)
    nx = 2*rm*np.cos(m*phi)/(1+r2m)
    ny = 2*rm*np.sin(m*phi)/(1+r2m)
    nz = (r2m-1)/(1+r2m)
    return np.array([nx,ny,nz])

def u(q0,qtheta,x,y):
    # director field of cholesteric phase
    qx,qy=q0*np.cos(qtheta),q0*np.sin(qtheta)
    rq=qx*x+qy*y
    ux=-np.sin(qtheta)*np.sin(rq)
    uy=np.cos(qtheta)*np.sin(rq)
    uz=np.cos(rq)
    return np.array([ux,uy,uz])

def ax_Enneper_config(ax,m,r1):
    # plot Enneper's surface and face normal
    r,phi=np.meshgrid(np.linspace(0,r1,100),np.linspace(0,2*np.pi,100))
    R = R_nmlz(m,r1)
    X,Y,Z=xyzENP(m,R,r,phi)
    #ax.plot_wireframe(X,Y,Z,color="black",linewidth=1,rstride=5, cstride=20)
    ax.plot_surface(X,Y,Z,color="red",edgecolor="black",linewidth=0.5,rstride=4, cstride=20,alpha=0.7)
    ax.axes.set_xlabel("x")
    ax.axes.set_ylabel("y")
    ax.axes.set_zlabel("z")
    xmin=1*X.min()
    xmax=1*X.max()
    ax.axes.set_xlim3d(left=xmin, right=xmax)
    ax.axes.set_ylim3d(bottom=xmin, top=xmax)
    ax.axes.set_zlim3d(bottom=xmin, top=xmax)
    # plot tilt

    #plt.title("m=%d,r1=%.1f,q0=%.1f,qtheta=%.1f"%(m,r1,q0,qtheta))
    #plt.savefig("config_m%d_r1%.1f_q0%.1f_qtheta%.1f.pdf"%(m,r1,q0,qtheta))

def xyzDisk(r,phi):
    x=r*np.cos(phi)
    y=r*np.sin(phi)
    z=0*np.zeros(np.shape(r))
    return np.array([x,y,z])

def nDisk(r,phi):
    nx=0
    ny=0
    nz=1
    return np.array([nx,ny,nz])

def ax_numeric_config_tilt(ax,m,r1,q,qtheta,d=0.1):

    if(m==0):
        # disk
        r,phi=np.meshgrid(np.linspace(0,1,100),np.linspace(0,2*np.pi,100))
        X,Y,Z=xyzDisk(r, phi)
        r,phi=r.flatten(),phi.flatten()
        nx,ny,nz=nDisk(r,phi)
    elif(m==1):
        r,phi=np.meshgrid(np.linspace(0,r1,100),np.linspace(0,2*np.pi,100))
        R = R_nmlz(m,r1)
        X,Y,Z=xyzENP(m,R,r,phi)
        r,phi=r.flatten(),phi.flatten()
        nx,ny,nz=nENP(m,R,r,phi)

    ax.plot_surface(X,Y,Z,linewidth=0.5,shade=0,color="gray",edgecolor="black",alpha=0.5,rstride=25, cstride=25)

    print("np.shape(X)",np.shape(X))
    print("np.shape(Y)",np.shape(Y))
    print("np.shape(Z)",np.shape(Z))
    x,y,z=X.flatten(),Y.flatten(),Z.flatten()
    ux,uy,uz=u(q,qtheta,x,y)
    deg=np.arccos(np.abs(ux*nx+uy*ny+uz*nz))
    #deg=np.arccos(np.abs(uz))

    cmap = cm.get_cmap("jet_r")
    norm=Normalize(vmin=0,vmax=0.5*np.pi)

    for i in range(len(x))[::5]:
        #ax.plot3D([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],[z[i]-0.5*d*uz[i],z[i]+0.5*d*uz[i]],"-",linewidth=1,color=cmap(norm(abs_un[i])),label=r"$u$")
        ax.plot3D([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],[z[i]-0.5*d*uz[i],z[i]+0.5*d*uz[i]],"-",linewidth=1,color=cmap(norm(deg[i])),label=r"$u$")

    #cbar=plt.colorbar(cm.ScalarMappable(norm=Normalize(vmin=0,vmax=0.5*np.pi), cmap=cm.get_cmap("jet_r")),ax=ax,ticks=[0,0.25*np.pi,0.5*np.pi])
    #cbar.ax.set_yticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"],fontsize=FontSize)
    #cbar.ax.tick_params(direction="in",labelsize=LabelSize)
    #cbar.ax.set_title(r"$\arccos{|\vu{u}\cdot\vu{n}|}$",fontsize=FontSize)


def rippling_config_demo(LineWidth, FontSize, LabelSize):
    print("configurations of rippling with 2 pi walls")
    foldername="../data/Ne1/Mar3_2021"
    #foldername="../data/Ne1/Mar26_2021"
    q=1.5
    lams=[4.0,6.0,8.0,10.0]
    Kd,Cn=5.0,5.0
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.5))
    #fig = plt.figure(figsize=plt.figaspect(2.))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    ax = plt.subplot2grid((1,1), (0, 0))
    rots=[(0.9,-1.0),(0.9,-1.0),(0.6,-1.0),(0.3,-1.0)]
    for i in range(len(lams)):
        lam=lams[i]
        filename = foldername+"/State_N500_Ne1_L-1_kar100_karg0.0_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_kard0.0.txt"%(lam,Kd,q,Cn)
        ax_config_xy(ax,filename,Color="gray",LineWidth=LineWidth,rotation=rots[i],xshift=30*i,yshift=0,bead=1,mesh=1,rod=1,d=0.65)
        ax_config_xy(ax,filename,Color="gray",LineWidth=LineWidth,rotation=rots[i],xshift=30*i,yshift=-20,bead=0,mesh=1,rod=0)
    for i in range(len(lams)):
        xl,yl=i/len(lams)+0.1,-0.05
        ax.text(ax.get_xlim()[1]*xl+ax.get_xlim()[0]* (1-xl),  ax.get_ylim()[1]*yl+ax.get_ylim()[0]*(1-yl), r"$\lambda=%.0f$"%lams[i], fontsize=FontSize)
    plt.tight_layout(pad=0.5)
    plt.savefig("figures/config_rippling.pdf")
    plt.close()

def config3d_demo_plot(LineWidth, FontSize, LabelSize):
    print("plotting configuration demo plot for slide")
    filename="../data/Ne1/Feb19_2021/State_N500_Ne1_L-1_kar100_karg0.0_lam10.0_Kd5.0_q2.0_Cn5.0_kard0.0.txt"
    mesh,rod=1,0
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    x,y,z,ux,uy,uz,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:14]
    x,y,z=x-np.average(x),y-np.average(y),z-np.average(z)
    d=1
    x_min, x_max = np.min(x),np.max(x)
    y_min, y_max = np.min(y),np.max(y)
    z_min, z_max = np.min(z),np.max(z)
    alpha_xy = 0.6*(z-z_min+0.1)/(z_max-z_min+0.1)+0.4

    ns = np.transpose(data[14:])
    ens = np.array([en0, en1])

    ppi = 72
    fig = plt.figure(figsize=(5, 5))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    ax = plt.axes(projection="3d")


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
            ax.plot3D([x[a],x[b]], [y[a],y[b]],[z[a],z[b]], color="gray",linewidth=LineWidth/4,alpha=alpha_xy[a])
        #color mapping
    cmap = cm.get_cmap("jet_r")
    # edge bond
    ecolors = ["black","purple","green","blue"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax.plot3D([x[j], x[int(ens[i, j])]], [
                    y[j], y[int(ens[i, j])]],[
                    z[j], z[int(ens[i, j])]], "-",linewidth=LineWidth/2, color=ecolors[int(enum[j])], alpha=alpha_xy[j])
    # director
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    #norm=Normalize(vmin=0,vmax=1)
    deg = np.arccos(np.sqrt(un2))
    sqrtun2=np.sqrt(un2)
    cmap = cm.get_cmap("jet_r")
    if(rod):
        for i in range(len(x)):
            ax.plot3D([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],[z[i]-0.5*d*uz[i],z[i]+0.5*d*uz[i]],"-",linewidth=LineWidth,color=cmap(norm(deg[i])),solid_capstyle='round')

    ax.set_xlim(x_min-2, x_max+2)
    ax.set_ylim(y_min-2, y_max+2)
    ax.set_zlim(z_min-2, z_max+2)
    ax.set_frame_on(False)
    ax.set_axis_off()
    for ii in range(20):
        ax.view_init(elev=30., azim=36/2*ii)
        plt.tight_layout()
        #plt.savefig("figures/config3d_demo/demo_rod_ani%d"%ii+".png",dpi=300)
        plt.savefig("figures/config3d_demo/demo_ani%d"%ii+".png",dpi=300)
    plt.close()

def director_update_demo(LineWidth, FontSize, LabelSize):
    print("plottin director update demo")

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1))
    #fig = plt.figure(figsize=plt.figaspect(2.))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')
    ax = plt.subplot2grid((1,1), (0, 0),projection="3d")
    #ax.add_artist(Circle(xy=(0, 0), radius=0.5,color="gray"))
    lt=1.5
    # axis
    ax.quiver(0,0,0,1,0,0,length=lt,color="gray",linewidth=LineWidth,alpha=0.7)
    ax.quiver(0,0,0,0,1,0,length=lt,color="gray",linewidth=LineWidth,alpha=0.7)
    #ax.quiver(0,0,0,0,0,1,length=lt,color="gray",linewidth=LineWidth,alpha=0.7)
    ax.quiver(0,0,0,0,0,1,length=lt,color="black",linestyle="--",linewidth=LineWidth,alpha=0.7)

    # director
    d=1
    theta,phi=0.5,0
    phis=np.linspace(0,2*np.pi,100)
    ax.plot(np.sin(theta)*np.cos(phis),np.sin(theta)*np.sin(phis),np.cos(theta),"--",color="gray")
    ux,uy,uz=np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)
    ax.quiver(0,0,0,ux,uy,uz,length=d,color="black",linewidth=2*LineWidth)

    ax.set_xlim(-0.1*lt,0.9*lt)
    ax.set_ylim(-0.1*lt,0.9*lt)
    ax.set_zlim(-0.1*lt,0.9*lt)

    ax.set_frame_on(False)
    ax.set_axis_off()
    #plt.show()
    ax.view_init(elev=15., azim=15)
    plt.savefig("figures/director_demo.pdf",transparent=1)
    plt.close()

