# beam-beam kick for leptons in FCC-ee (RUN with python NOT python3!)

import numpy as np
from numpy import sqrt, pi, exp,sign
from scipy import stats
import sys
from scipy.interpolate import interp1d
from scipy.integrate import quad, dblquad
import matplotlib.pyplot as plt
from errffor import errf

# see if you can use scipy.special.erfcx instead of errf
def wfun(z):
    x=z.real
    y=z.imag
    wx,wy=errf(x,y)

    return wx+1j*wy

# calculate the electric field along x and y
def BassErsk(Csigx,Csigy,sepx,sepy):

		# Separations are moved from mm to meters
    x=abs(sepx/10**3);
    y=abs(sepy/10**3);

		# Calculating back in meter the sigmax and sigmay of paper CERN-ISR-TH/80-06
    sigmax = sqrt((Csigx*(10**-6))**2)
    sigmay = sqrt((Csigy*(10**-6))**2)


		#I change the constant factor in front of the fields to have it in units of rp and gamma (K=2.*rp/gamma)
    eps0= 1.  

    S=sqrt(2*(sigmax*sigmax-sigmay*sigmay));
    factBE=sqrt(pi)*2/(2*eps0*S);
    etaBE=(sigmay/sigmax)*x+1j*(sigmax/sigmay)*y;
    zetaBE=x+1j*y;
        
    val=factBE*(wfun(zetaBE/S)-exp(-x*x/(2*sigmax*sigmax)-y*y/(2*sigmay*sigmay))*wfun(etaBE/S));
           
    Ex=abs(val.imag)*sign(sepx);
    Ey=abs(val.real)*sign(sepy);
     
    return Ex, Ey

# the beam-beam (BB) function, returns deflection angles and orbit deflections in x and y
def BB(Csigx,Csigy,sepx,sepy,betax,betay,tunex,tuney,Np):

    #Mass in atomic units u
    A=1
    # FCC-ee collision energy in eV (there are 4, here shown for the Z)
    Ee=45.6 *1e9

    # Electron classical radius in meter
    re=2.82e-15

    # Electron rest energy
    Eer=A*0.511 * 1e6

    # relativistic gamma factor 
    gamma=Ee/Eer
		
    K = 2*re/gamma

    Ex,Ey = BassErsk(Csigx,Csigy,sepx,sepy)

    Dfleix = K*Np*Ex # delta_x
    Dfleiy = K*Np*Ey # delta_y 

    Orbx = Dfleix*betax*(1./(2.*np.tan(pi*tunex))) # change in x due to BB, called orbit x
    Orby = Dfleiy*betay*(1./(2.*np.tan(pi*tuney))) # chanfe in y due to BB

    return Dfleix*10**6,Dfleiy*10**6,Orbx*10**6,Orby*10**6



# test (WORKS)
Csigx = 6.4 # in micro meters
Csigy = 28 * 1e-3

Scan = Csigx*10**-6

# all at Z
betax = 0.15 # m
betay = 0.8 # mm 
tunex = 269.139
tuney = 269.219

x = np.zeros(12000)
y1 = np.zeros(12000)
y2 = np.zeros(12000)
y3 = np.zeros(12000)
y4 = np.zeros(12000)

Ne = 16640*1.7*1e11 

i = 0

for d in np.arange(-6.*Scan,+6.*Scan,0.001*Scan):

	sepx=d*(10**3)
	sepy=0*Scan*(10**3)

	Dfleix,Dfleiy,Orbx,Orby = BB(Csigx,Csigy,sepx,sepy,betax,betay,tunex,tuney,Ne)

	x[i] = d*10**6;
	y1[i] = Dfleix ;
	y2[i] = Dfleiy ;
	y3[i] = Orbx ;
	y4[i] = Orby ;

	i+=1

plt.figure(dpi = 100)
plt.plot(x, y1, color='r', label = 'Kick_x')
plt.plot(x, y2, color='b', label = 'Kick_y')
plt.xlabel('sep units Capsigma')
plt.ylabel(r'BB Deflection angle [$\mu$ rad]')
plt.legend()

plt.figure(dpi = 100)
plt.plot(x, y3, color='r', label = 'Orbit_x')
plt.plot(x, y4, color='b', label = 'Orbit_y')
plt.xlabel('sep units Capsigma')
plt.ylabel(r'BB Orbit deflection [$\mu$ m]');
plt.legend()
#plt.show()



