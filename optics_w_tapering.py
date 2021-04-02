import at
from at import load_mat
from at.physics import linopt, linopt_rad
import os
import matplotlib.pyplot as plt
import numpy as np
import at.plot.specific

path = './Lattices/'
filename = 'fcch_norad.mat'
key = 'ring'
latticef = path+filename
ring = at.load.matfile.load_mat(latticef, mat_key=key)
ring = at.Lattice(ring)

# optics of the original lattice
ring.radiation_off()
l0, q, qp, l = linopt(ring,dp=0,refpts=range(len(ring)),get_chrom = True, coupled=False)

# optics of tapered lattice
ring.radiation_on(quadrupole_pass='auto')
ring.set_cavity_phase()
ring.tapering(niter = 2, quadrupole=True, sextupole=True)
l0t,qt,qpt,lt = linopt_rad(ring, refpts=range(len(ring)), get_chrom=True, coupled=False)

# s along the track
spos = ring.get_s_pos(range(len(ring)))

# optics comparison
plt.figure(figsize=(12,6), dpi=200)
plt.subplot(121)
plt.plot(spos,(lt.beta[:,0]-l.beta[:,0])/l.beta[:,0], label='x plane')
plt.xlabel('s [m]')
plt.ylabel(r'$\frac{\beta - \beta_0}{\beta_0}$')
plt.legend(loc = 1, prop={'size': 8})

plt.subplot(122)
plt.plot(spos,(lt.beta[:,1]-l.beta[:,1])/l.beta[:,1], label='y plane')
plt.xlabel('s [m]')
plt.ylabel(r'$\frac{\beta - \beta_0}{\beta_0}$')
plt.legend(loc = 1, prop={'size': 8})
plt.show()

plt.figure(figsize=(12,6), dpi=200)
plt.plot(spos, lt.dispersion[:,0] - l.dispersion[:,0], label='x_plane')
plt.xlabel('s [m]')
plt.ylabel(r'$D - D_0$')
plt.legend(loc = 1, prop={'size': 8})
plt.show()
