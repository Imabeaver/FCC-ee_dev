# beam-beam kick for leptons in FCC-ee

import numpy as np
from numpy import sqrt, pi, exp,sign
from scipy import stats
import sys
from scipy.interpolate import interp1d
from scipy.integrate import quad, dblquad
import matplotlib.pyplot as plt
from scipy.special import wofz

# calculate the electric field along x and y
def BassErsk(Csigx,Csigy,sepx,sepy):

    # Separations are moved from mm to meters
    x=abs(sepx/10**3);
    y=abs(sepy/10**3);

    # Calculating back in meter the sigmax and sigmay of paper CERN-ISR-TH/80-06
    sigmax = np.sqrt((Csigx*(10**-6))**2)
    sigmay = np.sqrt((Csigy*(10**-6))**2)


    #I change the constant factor in front of the fields to have it in units of rp and gamma (K=2.*rp/gamma)
    eps0= 1.  

    S=sqrt(2*(sigmax*sigmax-sigmay*sigmay));
    factBE=sqrt(pi)*2/(2*eps0*S);
    etaBE=(sigmay/sigmax)*x+1j*(sigmax/sigmay)*y;
    zetaBE=x+1j*y;
        
    val=factBE*(wofz(zetaBE/S)-exp(-x*x/(2*sigmax*sigmax)-y*y/(2*sigmay*sigmay))*wofz(etaBE/S));
           
    Ex=abs(val.imag)*sign(sepx);
    Ey=abs(val.real)*sign(sepy);
     
    return Ex, Ey

# the beam-beam (BB) function, returns deflection angles and orbit deflections in x and y
def BB(Csigx,Csigy,sepx,sepy,betax,betay,tunex,tuney,Ne,Ee):

    # Electron mass in atomic units u
    A=0.000549

    # Electron classical radius (m)
    re=2.82e-15

    # Electron rest energy (eV)
    Eer=A*0.511 * 1e6

    # relativistic gamma factor 
    gamma=Ee/Eer

    K = 2*re/gamma

    Ex,Ey = BassErsk(Csigx,Csigy,sepx,sepy)

    Dfleix = -K*Ne*Ex # delta_x
    Dfleiy = -K*Ne*Ey # delta_y 

    Orbx = Dfleix*betax*(1./(2.*np.tan(pi*tunex))) # change in x due to BB, called orbit x
    Orby = Dfleiy*betay*(1./(2.*np.tan(pi*tuney))) # chanfe in y due to BB

    return Dfleix*10**6,Dfleiy*10**6,Orbx*10**6,Orby*10**6

def kick_calc(x,y,Ee, sigx, sigy, betax, betay, tunex, tuney, Ne):
    # x, y must be a single value! If given as array the code doesn't perform as it should.
    
    Csigx = np.sqrt(np.power(sigx,2))
    Csigy = np.sqrt(np.power(sigy,2))
    
    sepx = x * 1e3 # mm
    sepy = y * 1e3 # mm
		
    Dfleix,Dfleiy,Orbx,Orby = BB(Csigx,Csigy,sepx,sepy,betax,betay,tunex,tuney,Ne,Ee)

    return Dfleix, Dfleiy, x*1e6/sigx

# end


