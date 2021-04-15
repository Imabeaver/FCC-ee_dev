# beam-beam kick integration in the FCC-ee lattice
import at
import numpy as np
import  bb_kick_integrated as bb
import matplotlib.pyplot as plt
from at.tracking import element_pass, lattice_pass
from at.physics import linopt, fast_ring
from at.lattice import DConstant, get_refpts, RFCavity, Dipole
from at.lattice import PyElement
import time

# load lattices, step sizes adjusted to avoid systematics on lattice calculations
xy_step = 1.0e-9
dp_step = 1.0e-6

ring = at.load_lattice('./fcc_lattice/fcct_norad.mat', mat_key='ring')

#init rings - deepcopy() required to remove all shared references
ring_norad = ring.deepcopy()
ring_norad.radiation_off(copy=True)
ring_rad = ring.deepcopy()

E0 = ring_norad.energy
print('Nominal energy:', E0, 'eV')

ring_rad.radiation_on(quadrupole_pass='auto')
ring_rad.set_cavity_phase()

ring_rad.tapering(niter = 2, XYStep=xy_step, DPStep=dp_step)

# emittance
env = ring_rad.envelope_parameters()
epsilon_x = env.emittances[0]
epsilon_y = env.emittances[1]
sigma_s = env.sigma_l

if (E0 == 120*1e9) is True:
    epsilon_y = 1.3*1e-12 # ZH
    Ne = 1.8*1e11 # H(ZH)
    sigma_s = 5.3e-3 #with BS
elif (E0 == 45.6*1e9) is True:
    epsilon_y = 1.0*1e-12 # Z
    Ne = 1.7*1e11 # Z
    sigma_s = 12.1e-3 #with BS
elif (E0 == 175*1e9) is True:
    epsilon_y = 2.7*1e-12 # TT
    Ne = 2.2*1e11 # TT
    sigma_s = 3.82e-3 #with BS
elif (E0 == 182.5*1e9) is True:
    epsilon_y = 2.9*1e-12 # TT
    Ne = 2.3*1e11 # TT
    sigma_s = 2.54e-3 #with BS
elif (E0 == 80*1e9) is True:
    epsilon_y = 1.7*1e-12 # WW
    Ne = 1.5*1e11 # H(ZH)
    sigma_s = 6.0e-3 #with BS
    
print('Emittance_x:',epsilon_x,'\n','Emittance_y', epsilon_y)

lindata0, _, _, _ = ring_norad.linopt(coupled=False, XYStep=xy_step, DPStep=dp_step)
beta_x, beta_y = lindata0['beta']
print('Beta_x:', beta_x,'\n','Beta_y:',beta_y)

# beam size
sigma_x = np.sqrt(epsilon_x* beta_x)
sigma_y = np.sqrt(epsilon_y*beta_y)
print('Sigma_x:', sigma_x, '\n', 'Sigma_y:', sigma_y, '\n', 'Sigma_s (with BS):', env.sigma_l, sigma_s)
      

#values found online: Boscolo talk
theta_x = 30e-3
theta_y = 0.0
bbparamx, bbparamy,sx,sy = bb.bb_param(Ne, E0, beta_x, beta_y, sigma_x, sigma_y, sigma_s, theta_x, theta_y)
print('theta_x:', theta_x, '\n', 'theta_y:', theta_y)
print('zeta_x:', bbparamx, '\n', 'zeta_y:', bbparamy)

#build beam beam element
params = {'sigma':np.array([sigma_x,sigma_y,sigma_s]),
          'theta':np.array([theta_x,theta_y]),
          'E0':E0,
          'Ne':Ne}         
bb_elem = PyElement('pyelem','bb_kick_integrated','passmethod',**params)
#concatenation makes a copy
ring_bb = ring_norad + [bb_elem]

#Compute tune shift
_, q0, _, _ = ring_norad.linopt(coupled=False, XYStep=xy_step, DPStep=dp_step)
_, qbb, _, _ = ring_bb.linopt(coupled=False, XYStep=xy_step, DPStep=dp_step)
print('Linear BB tune shift: ',qbb-q0)

# define number of turns and number of particles
# and range f in sigma or tracking, equidistant particles
nturns = 100
npart = 10
mina = 0.1
maxa = 5.0

#track with beam-beam
x = np.linspace(mina,maxa, npart)*sigma_x*sx
rin = np.zeros((6,npart))
rin[0,:]=x
t0 = time.time()
#keep_lattice = False required to properly initialize the lattice
rout_bb = np.squeeze(lattice_pass(ring_bb, rin, nturns=nturns, keep_lattice=False))
print('Tracking with bb took: ',time.time()-t0)

#track without beam-beam
x = np.linspace(mina,maxa, npart)*sigma_x*sx
rin = np.zeros((6,npart))
rin[0,:]=x
t0 = time.time()
#keep_lattice = False required to properly initialize the lattice
rout_nobb = np.squeeze(lattice_pass(ring_norad, rin, nturns=nturns, keep_lattice=False))
print('Tracking without bb took: ',time.time()-t0)

# phase space plotting
plt.figure(figsize=(8,6))
plt.scatter(rout_bb[0,:], rout_bb[1,:], color = 'blue', marker='o', s=1)
plt.scatter(rout_nobb[0,:], rout_nobb[1,:], color = 'red', marker='o', s=1)
plt.xlabel('x [m]')
plt.ylabel(r'$x^{prime}$ [rad]')
plt.xlim(np.array([-6.0, 6.0])*sigma_x*sx)
plt.ylim(np.array([-6.0, 6.0])*sigma_x*sx/beta_x)

# kick force plot and its derivative, test that kick is good
# derivative(0) should be equal to the beam beam parameter
npart = 1000
zz = np.zeros(npart)
x = np.linspace(-10,10, npart)*sigma_x*sx
y = np.linspace(-10,10, npart)*sigma_y*sy
xp, _ = bb.kick_calc(x,zz, sigma_x, sigma_y, E0, Ne, sigma_s, theta_x, theta_y)
_, yp = bb.kick_calc(zz,y, sigma_x, sigma_y, E0, Ne, sigma_s, theta_x, theta_y)
 
fig,ax=plt.subplots(figsize=(8,6))
ax2 = ax.twinx()  
ax.scatter(x, xp, color = 'blue', marker='o', s=1, label='Kick X')
ax2.scatter(x, np.gradient(xp,x)*beta_x/(4*np.pi), color = 'green', marker='o', s=1, label='Derivative X')
ax.set_xlabel('x [m]')
ax.set_ylabel(r'dx'' [$\mu$ rad]')
ax2.set_ylabel(r'dx''/dx [$\mu$ rad]')

fig,ax=plt.subplots(figsize=(8,6))
ax2 = ax.twinx()  
ax.scatter(y, yp, color = 'red', marker='o', s=1, label='Kick Y')
ax2.scatter(y, np.gradient(yp,y)*beta_y/(4*np.pi), color = 'black', marker='o', s=1, label='Derivative Y')
ax.set_xlabel('x [m]')
ax.set_ylabel(r'dy'' [$\mu$ rad]')
ax2.set_ylabel(r'dy''/dy [$\mu$ rad]')

plt.show()
