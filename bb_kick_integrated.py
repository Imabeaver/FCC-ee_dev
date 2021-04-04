# beam-beam kick for leptons in FCC-ee

import numpy as np
from scipy.special import wofz
from scipy.constants import physical_constants

#particle properties, could be an input option
_re, _, _ = physical_constants['classical electron radius']
_me, _ ,_ = physical_constants['electron mass energy equivalent in MeV']
_me *= 1e6

# calculate the electric field along x and y
def BassErsk(x,y,sigmax,sigmay):
    sx2 = 2*sigmax*sigmax
    sy2 = 2*sigmay*sigmay
    S= np.sqrt(sx2-sy2)
    factBE=np.sqrt(np.pi)/S
    etaBE=(sigmay/sigmax)*x+1j*(sigmax/sigmay)*y
    zetaBE=x+1j*y 
    #ignore RunTime warning 
    with np.errstate(invalid='ignore'):    
        val=factBE*(wofz(zetaBE/S)-np.exp(-x*x/sx2-y*y/sy2)*wofz(etaBE/S))  
    Ex=val.imag
    Ey=val.real    
    return Ex, Ey

# convert the electric field into deflection angles
def kick_calc(x,y,sigma_x,sigma_y,Ee,Ne,sigma_s,theta_x,theta_y):
    gamma=Ee/_me
    K = -2*_re*Ne/gamma
    sx = np.sqrt(1+(sigma_s/sigma_x*np.tan(theta_x/2))**2)
    sy = np.sqrt(1+(sigma_s/sigma_y*np.tan(theta_y/2))**2)
    dx,dy = BassErsk(x,y,sigma_x*sx,sigma_y*sy)
    return dx*K, dy*K

#compute head-on linear beam-beam parameters
def bb_param(Ne, E0,beta_x,beta_y,sigma_x,sigma_y,sigma_s,theta_x,theta_y):
    sx = np.sqrt(1+(sigma_s/sigma_x*np.tan(theta_x/2))**2)
    sy = np.sqrt(1+(sigma_s/sigma_y*np.tan(theta_y/2))**2)
    bbparamx = Ne*_re*_me*beta_x/(E0*2*np.pi*sigma_x*sx*(sigma_x*sx+sigma_y*sy))
    bbparamy = Ne*_re*_me*beta_y/(E0*2*np.pi*sigma_y*sy*(sigma_x*sx+sigma_y*sy))
    return bbparamx,bbparamy,sx,sy

#apply the kick to a particle array without integration in the lattice
def apply_kick(rin,sigx,sigy,E0,Ne,sigma_s,theta_x,theta_y):
    dx,dy = kick_calc(rin[0],rin[2],sigx,sigy,E0,Ne,sigma_s,theta_x,theta_y)
    #conversion from x' to px/p0
    rin[1] +=  dx*(1+rin[4])
    rin[3] +=  dy*(1+rin[4])

#apply the kick to a particle with pyElement integrated in the lattice
#Warning AT passmethods take vectors as input, not 2D arrays
def passmethod(rin,elem=None):
    sigx, sigy, sigma_s = elem.sigma
    theta_x, theta_y = elem.theta
    E0 = elem.E0
    Ne = elem.Ne
    dx,dy = kick_calc(rin[0::6],rin[2::6],sigx,sigy,E0,Ne,sigma_s,theta_x,theta_y)
    #conversion from x' to px/p0
    rin[1::6] +=  dx*(1+rin[4::6])
    rin[3::6] +=  dy*(1+rin[4::6])
    return rin




