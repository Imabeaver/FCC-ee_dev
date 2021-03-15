# beam-beam kick integration in the FCC-ee lattice
import at
import numpy as np
from bb_kick_integrated import kick_calc
from at import load_mat, get_refpts
import matplotlib.pyplot as plt
import at.plot.specific

def b2_kick(rin, E0, sigx, sigy, betax, betay, tunex, tuney, Ne):
    kick = kick_calc(rin[0],rin[2],E0, sigx, sigy, betax, betay, tunex, tuney, Ne)
    rout = rin
    rout[1] = rin[1] + kick[0]*1.0e-6
    rout[3] = rin[3] + kick[1]*1.0e-6
    return rout

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
sigma_x = np.sqrt(epsilon_x* beta_x) *1e6 # micro m
sigma_y = np.sqrt(epsilon_y*beta_y) * 1e6 # micro m
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

# radiation off
for i in range(nturns):
    rin_norad = np.squeeze(at.track.lattice_pass(ring_norad, rin, nturns=1))
    rout_norad[i,:] = b2_kick(rin_norad, E0, sigma_x, sigma_y, beta_x, beta_y, t[0], t[1], Ne)

# radiation on
rin = np.array([1.0e-6,0,0,0,0,0])
for i in range(nturns):
    rin_rad = np.squeeze(at.track.lattice_pass(ring_rad, rin, nturns=1))
    rout_rad[i,:] = b2_kick(rin_rad, E0, sigma_x, sigma_y, beta_x, beta_y, t[0], t[1], Ne)

# phase space plotting
plt.figure(figsize=(8,6))
plt.scatter(rout_norad[:,0], rout_norad[:,1], color = 'blue', marker='o', s=1, label='Radiation off')
plt.scatter(rout_rad[:,0], rout_rad[:,1], color = 'red', marker='o', s=1, label='Radiation on')
plt.xlabel('x [m]')
plt.ylabel(r'$x^{prime}$ [rad]')
plt.legend(loc = 1, prop={'size': 8})
#plt.savefig('phase_space_x.png')
plt.show()

# kick force plot, test that kick is good
x = np.linspace(-6*1e-5,6*1e-5, 1000)
xp = np.zeros(len(x))
for i in range(len(x)):
    rin = [x[i], 0, 0, 0, 0, 0]
    kick = kick_calc(rin[0],rin[2], E0, sigma_x, sigma_y, beta_x, beta_y, t[0], t[1], Ne)
    xp[i] = kick[0]
    
plt.figure(figsize=(8,6))
plt.scatter(x, xp, color = 'blue', marker='o', s=1, label='Kick')
plt.xlabel('x [m]')
plt.ylabel(r'dx [$\mu$ rad]')
#plt.savefig('Kick_force.png')
plt.show()
