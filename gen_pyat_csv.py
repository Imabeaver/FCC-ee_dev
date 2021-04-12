import at, os
import matplotlib.pyplot as plt
import numpy as np
from at.physics import linopt_rad
import pandas as pd

# import lattice
path = './'
filename = 'fcch_norad.mat'
key = 'ring'

xy_step = 1.0e-9
dp_step = 1.0e-6

latticef = path+filename
ring = at.load.matfile.load_mat(latticef,key=key)
ring = at.Lattice(ring)
spos = ring.get_s_pos(range(len(ring)))

get_chrom = True

ring.radiation_on(quadrupole_pass='auto')
ring.set_cavity_phase()
ring.tapering(niter = 2, quadrupole=True, sextupole=True, XYStep=xy_step, DPStep=dp_step)
l0rt,qrt,qprt,lrt = linopt_rad(ring,refpts=range(len(ring)),get_chrom=get_chrom,
                               coupled=False, XYStep=xy_step, DPStep=dp_step)

# export pandas file
lattice = 'h' # change accordingly

no = np.arange(len(spos))
export = {'EL. NUMBER': no, 'S': spos.T, 'BETX_PYAT': lrt.beta[:,0].T, 'BETY_PYAT': lrt.beta[:,1].T, 'DX_PYAT': lrt.dispersion[:,0].T, 'DY_PYAT': lrt.dispersion[:,1].T}
df = pd.DataFrame(export, columns=['EL. NUMBER','S','BETX_PYAT','BETY_PYAT','DX_PYAT','DY_PYAT'])
df.to_csv(path + lattice + '/data_tapered_pyat.csv', index = False, header = True)
