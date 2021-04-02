# 4D beam-beam kick calculation for electrons in FCC-ee

import numpy as np
from numpy import sqrt, pi, exp, sign
from scipy.special import wofz

# Calculates the electric field along x and y
def BassErsk(sigx,sigy,x,y):
    
    # sigmas squared for faster computation
    sigx2 = sigx*sigx
    sigy2 = sigy*sigy
    
    #Constant factor in front of the fields is changed to have it in units of re and gamma (K=2.*rp/gamma)
    eps0 = 1.0

    S = sqrt(2*(sigx2-sigy2))
    factBE = sqrt(pi)*2/(2*eps0*S)
    etaBE = (sigy/sigx)*x + 1j*(sigx/sigy)*y
    zetaBE = x+1j*y
        
    val = factBE*(wofz(zetaBE/S)-exp(-x*x/(2*sigx2)+y*y/(2*sigy2))*wofz(etaBE/S))
           
    Ex=val.imag
    Ey=val.real
     
    return Ex, Ey

# The beam-beam (BB) function, returns deflection angles in x and y
# x and y input in (m), angle deflection output in (m)
# sigx and sigy inputs are beam sizes given in (m)
# Ee is the initial beam energy given in (eV)
# Ne is the number of electrons per bunch
def kick_calc(x, y, Ee, sigx, sigy, Ne):

    # Electron classical radius (m)
    re = 2.82e-15

    # Electron rest energy (eV)
    Eer = 0.511 * 1e6

    # Relativistic gamma factor
    gamma = Ee/Eer

    K = 2*re/gamma
    
    # Scaling for crossing angle
    theta = 30 * 1e-3 # rad
    sigs = 3.15 * 1e-3 # m
    S = 1/np.sqrt(1 + (sigs/sigx * np.tan(theta/2))**2)
    sigx = sigx/S
    
    # Electric fields
    Ex,Ey = BassErsk(sigx,sigy,x,y)

    Dfleix = -K*Ne*Ex # Deflection angle in x (m)
    Dfleiy = -K*Ne*Ey # Deflection angle in y (m)

    return Dfleix,Dfleiy

# end
