#!/usr/local/bin/python3
import numpy as np
from scipy import special, optimize, integrate


# Enneper's surface related
#TODO: try write in terms of r/r1\in(0,1), instead of r
def xyzENP(m,R,r,phi):
    # cartesian cooredinate of Enneper's surface
    x = r*np.cos(phi)-np.power(r,2*m+1)/(2*m+1)*np.cos((2*m+1)*phi)
    x*=R
    y = -r*np.sin(phi)-np.power(r,2*m+1)/(2*m+1)*np.sin((2*m+1)*phi)
    y*=R
    z = 2*np.power(r,m+1)/(m+1)*np.cos((m+1)*phi)
    z*=R
    return np.array([x,y,z])

def AENP(m,R,r1):
    # Area of Enneper's surface
    A = np.pi*np.power(R,2)*np.power(r1,2)
    A*=(1+2*np.power(r1,2*m)/(1+m)+np.power(r1,4*m)/(1+2*m))
    return A

def jcbENP(m,R,r):
    # discrete area at (r,phi)
    J = r
    J*=np.power((1+np.power(r,2*m)),2)
    J*=R*R
    return J

def dAENP(m,R,r,dr,dphi):
    # discrete area at (r,phi)
    J = r
    J*=np.power((1+np.power(r,2*m)),2)
    J*=R*R
    return J*dr*dphi

def R_nmlz(m,r1):
    # normalized R for A=1
    R =1/r1
    R/=np.sqrt(1+2*np.power(r1,2*m)/(1+m)+np.power(r1,4*m)/(1+2*m))
    return R

def IntKdA(m,r1):
    KdA = -4*m*np.pi
    r2m = np.power(r1,2*m)
    KdA*=r2m/(1+r2m)
    return KdA

def n(m,R,r,phi):
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



def ndoctu_ENP(m,R,r,phi,q0,qtheta):
    # calculate n dot u directly from simplified analytical formula
    # confirmed with preview result
    r2m = np.power(r,2*m)
    cosphipqtheta = np.cos(phi+qtheta)
    cosphimqtheta = np.cos((2*m+1)*phi-qtheta)
    res = (r2m-1)*np.cos(q0*r*R*(cosphipqtheta-r2m*cosphimqtheta/(1+2*m)))
    res+=2*np.power(r,m)*np.sin(m*phi-qtheta)*np.sin(q0*r*R*(cosphipqtheta-r2m*cosphimqtheta/(1+2*m)))
    res/=(r2m+1)
    # difference goes to zero
    return res

def tilt_ENP_2_int(r,phi,m,R,q0,qtheta):
    ndotu=ndoctu_ENP(m,R,r,phi,q0,qtheta)
    tilt = (1-np.power(ndotu,2))*jcbENP(m,R,r)
    return tilt

def tilt_ENP(m,r1,q0,qtheta,bn_r,bn_phi):
    r,phi=np.meshgrid(r1*np.linspace(0.5/bn_r,1-0.5/bn_r,bn_r),2*np.pi*np.linspace(0.5/bn_phi,1-0.5/bn_phi,bn_phi))
    r,phi= r.flatten(),phi.flatten()
    dr=r1/bn_r
    dphi=2*np.pi/bn_phi
    R = R_nmlz(m,r1)
    dA = dAENP(m,R,r,dr,dphi)
    #print("dA.sum()",dA.sum()) 3.141592...
    ndotu=ndoctu_ENP(m,R,r,phi,q0,qtheta)
    tilt = (1-np.power(ndotu,2))*dA
    #print("tilt.sum()",tilt.sum())
    #tilt_integral=integrate.dblquad(tilt_ENP_2_int,0,r1,0,2*np.pi,args=(m,R,q0,qtheta))
    # scipy integral is slow as turtle, and inaccurate
    return tilt.sum()

def peri_ENP(m,r1):
    # perimeter of Enneper surface
    peri = (1+np.power(r1,2*m))*R_nmlz(m,r1)*r1
    peri*=2*np.pi
    return peri

def tilt_disk(q0):
    # tilt of disk membrane
    if(q0==0):
        un2=1
    else:
        un2 = 0.5*(1+special.j1(2*q0)/q0)
    tilt = np.pi*(1-un2)
    # checked by plotting and compare with mathematica
    return tilt

def peri_disk():
    # perimeter of disk, which is 2pi apparently
    return 2*np.pi

def E_ENP(lam,Cn,m,r1,q0,qtheta,bn_r,bn_phi):
    tiltENP=tilt_ENP(m,r1,q0,qtheta,bn_r,bn_phi)
    periENP=peri_ENP(m,r1)
    return lam*periENP+0.5*Cn*tiltENP

def obj_E_ENP(xs,lam,Cn,m,q0,bn_r,bn_phi):
    r1,qtheta=xs
    return E_ENP(lam,Cn,m,r1,q0,qtheta,bn_r,bn_phi)

def E_Disk(lam,Cn,q):
    tiltdisk=tilt_disk(q)
    peridisk=peri_disk()
    return lam*peridisk+0.5*Cn*tiltdisk

def DelE(lam,Cn,m,r1,q0,qtheta,bn_r,bn_phi):
    Del_E = E_ENP(lam,Cn,m,r1,q0,qtheta,bn_r,bn_phi)-E_Disk(lam,Cn,m,r1,q0,qtheta,bn_r,bn_phi)
    return Del_E

def obj_DelE(xs,lam,Cn,m,q0,bn_r,bn_phi):
    r1,qtheta=xs
    #lam,Cn,m,q0,bn_r,bn_phi=paras
    return DelE(lam,Cn,m,r1,q0,qtheta,bn_r,bn_phi)

def opt_r1_qtheta_fun(fun,lam,Cn,m,q0,bn_r,bn_phi,method):
    # find the minimum value of DelE, and corresponding r1 and qtheta
    paras=(lam,Cn,m,q0,bn_r,bn_phi)
    r1qtheta0=(1.0,0.0)
    r1is_m = {1:np.sqrt(3),2:np.sqrt(2),3:np.power(2*np.sqrt(7)-1,1/6),4:1.20393,5:1.16169}
    r1bound=r1is_m[m]-0.1
    if method:
        opt = optimize.minimize(fun,r1qtheta0,args=paras,bounds=((0.8,r1bound),(0,np.pi)),method=method)
    else:
        opt = optimize.minimize(fun,r1qtheta0,args=paras,bounds=((0.8,r1bound),(0,np.pi)))
    return opt

def error_test():
    r1s = np.linspace(0.01,0.11,10)
    m=1
    bn_r=100
    rr1=np.linspace(0.5/bn_r,1-0.5/bn_r,bn_r)
    for i in range(len(r1s)):
        r1=r1s[i]
        r=r1*rr1
        R=R_nmlz(m,r1)
        print("Rr",R*r)
        print("r/r1*R*r1", rr1/np.sqrt(1+2*np.power(r1,2*m)/(1+m)+np.power(r1,4*m)/(1+2*m)))
        print("difference:",R*r-rr1/np.sqrt(1+2*np.power(r1,2*m)/(1+m)+np.power(r1,4*m)/(1+2*m)))