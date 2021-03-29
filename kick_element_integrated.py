# Implementation of beam-beam kick element in the FCC-ee lattice

import at
import numpy as np
from bb_kick_integrated import kick_calc
from at import load_mat, get_refpts
from at.physics import linopt, linopt_rad
import matplotlib.pyplot as plt
import at.plot.specific

def b2_kick(rin, E0, sigx, sigy, Ne):
    kick = kick_calc(rin[0],rin[2],E0, sigx, sigy, Ne)
    
    rin[1] += kick[0]
    rin[3] += kick[1]
    
    return rin

# load lattices
rin = np.array([1.0e-6,0,0,0,0,0])

ring_norad = at.load_mat('./Lattices/fcch_norad.mat', mat_key='ring')
ring_rad = at.load_mat('./Lattices/fcch_rad.mat', mat_key='ring')
ring_tapered = at.load_mat('./Lattices/fcch_rad_tapered.mat', mat_key='ring')
ring_tapered.radiation_on()

# tune, chromaticity and energy
t = ring_norad.get_tune()
q = ring_norad.get_chrom()
E0 = ring_norad.energy

print('Nominal energy:', E0, 'eV')
print('Tune:', t, '\n', 'Chromaticity:', q)

# emittance
env = ring_tapered.envelope_parameters()
epsilon_x = env.emittances[0]
epsilon_y = env.emittances[1]

if (E0 == 120*1e9) is True:
    epsilon_y = 1.3*1e-12 # ZH
elif (E0 == 45.6*1e9) is True:
    epsilon_y = 1.0*1e-12 # Z
elif (E0 == 182.5*1e9) is True:
    epsilon_y = 2.9*1e-12 # TT
elif (E0 == 80*1e9) is True:
    epsilon_y = 1.7*1e-12 # WW
    
print('Emittance_x:',epsilon_x,'\n','Emittance_y', epsilon_y)

# optics
IP_indexes = get_refpts(ring_norad, 'IP*')
lindata0, tune, chrom, lindata = ring_norad.linopt(get_chrom=True, refpts=IP_indexes, coupled=False)
beta_x = lindata['beta'][0,0]
beta_y = lindata['beta'][0,1]
print('Beta_x:', beta_x,'\n','Beta_y:',beta_y)
ring_norad.plot_beta()

# beam size
sigma_x = np.sqrt(epsilon_x* beta_x)  # m
sigma_y = np.sqrt(epsilon_y*beta_y) #  m
print('Sigma_x:', sigma_x, '\n', 'Sigma_y:', sigma_y)

# number of e- per bunch (bunch population)
if (E0 == 120*1e9) or (E0 == 80*1e9) is True:
        Ne = 1.5*1e11 # H(ZH)
elif (E0 == 45.6*1e9) is True:
        Ne = 1.7*1e11 # Z
elif (E0 == 182.5*1e9) is True:
        Ne = 2.3*1e11 # TT

# defining output
nturns = 100
rout_norad = np.zeros((nturns, 6))
rout_rad = np.zeros((nturns, 6))

## ----------------------------------------------------------------------
# phase space ploting with no bb or rad
rin = np.array([1.0e-6,0,0,0,0,0])
rout_no_kick = np.squeeze(at.track.lattice_pass(ring_norad, rin, nturns=100))

plt.figure(figsize=(12,6))
plt.subplot(121)
plt.scatter(rout_no_kick[0], rout_no_kick[1], color = 'blue', marker='o', s=1, label='Radiation off, Kick off')
plt.xlabel('x [m]')
plt.ylabel(r'$x^{prime}$ [rad]')
plt.legend(loc = 1, prop={'size': 8})

plt.subplot(122)
plt.scatter(rout_no_kick[2], rout_no_kick[3], color = 'blue', marker='o', s=1, label='Radiation off, Kick off')
plt.xlabel('y [m]')
plt.ylabel(r'$y^{prime}$ [rad]')
plt.legend(loc = 1, prop={'size': 8})
plt.show()
## ----------------------------------------------------------------------

# radiation off, bb on
rin = np.array([1.0e-6,0,0,0,0,0])
for i in range(nturns):
    rin_norad = np.squeeze(at.track.lattice_pass(ring_norad, rin, nturns=1, keep_lattice=False, omp_num_threads=None))
    rin = b2_kick(rin_norad, E0, sigma_x, sigma_y, Ne)
    rout_norad[i,:] = rin
    #print(rin_norad)

# radiation on
#rin = np.array([1.0e-6,0,0.0e-6,0,0,0])
#for i in range(nturns):
#    rin = np.squeeze(at.track.lattice_pass(ring_rad, rin, nturns=1))
#    rout_rad[i,:] = b2_kick(rin, E0, sigma_x, sigma_y, Ne)

## ----------------------------------------------------------------------
# phase space plotting with bb and rad
plt.figure(figsize=(12,6))
plt.subplot(121)
plt.scatter(rout_norad[:,0], rout_norad[:,1], color = 'blue', marker='o', s=1, label='Radiation off, Kick on')
#plt.scatter(rout_rad[:,0], rout_rad[:,1], color = 'red', marker='o', s=1, label='Radiation on')
plt.xlabel('x [m]')
plt.ylabel(r'$x^{prime}$ [rad]')
plt.legend(loc = 1, prop={'size': 8})
#plt.savefig('phase_space_x.png')

plt.subplot(122)
plt.scatter(rout_norad[:,2], rout_norad[:,3], color = 'blue', marker='o', s=1, label='Radiation off, Kick on')
#plt.scatter(rout_rad[:,0], rout_rad[:,1], color = 'red', marker='o', s=1, label='Radiation on')
plt.xlabel('y [m]')
plt.ylabel(r'$y^{prime}$ [rad]')
plt.legend(loc = 1, prop={'size': 8})
plt.show()
## ----------------------------------------------------------------------

# kick force plot, test that kick is good
x = np.linspace(-10, 10, 1000)*sigma_x
xp = np.zeros(len(x))
for i in range(len(x)):
    kick = kick_calc(x[i],0, E0, sigma_x, sigma_y, Ne)
    xp[i] = kick[0]
    
plt.figure(figsize=(8,6))
plt.scatter(x, xp*1e6, color = 'blue', marker='o', s=1, label='Kick')
plt.scatter(x, np.gradient(xp,x)*beta_x/(4*np.pi), color = 'green', marker='o', s=1, label='Derivative X')
plt.xlabel('x [m]')
plt.ylabel(r'dx [$\mu$ rad]')
plt.legend()
#plt.savefig('Kick_force.png')
plt.show()

## ----------------------------------------------------------------------


# Beam-beam parameter calculation
re = 2.82e-15
Eer = 0.511 * 1e6
gamma = E0/Eer
theta = 30 * 1e-3 # rad
sigma_s = 3.15 * 1e-3 # m
S = 1/np.sqrt(1 + (sigma_s/sigma_x * np.tan(theta/2))**2)

eps_x = Ne*re*beta_x*S**2/(2*np.pi*gamma*sigma_x*(sigma_x+S*sigma_y))
eps_y = Ne*re*beta_y*S/(2*np.pi*gamma*sigma_y*(sigma_x+S*sigma_y))

print('Horisontal beam-beam:', eps_x, '\n', 'Vertical beam-beam:', eps_y)
