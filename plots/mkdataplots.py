"""mkdataplots

Script to make basics triangle plots of data.
"""
import os
import sys
import data_triangle
from matplotlib import pyplot as plt
import numpy as np
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import pm_to_velocities as pmvels

# UCAC Data Load
pmucac = pmvels.PMmeasurements()
pmucac.to_space_velocties()
pmucac.UVW_to_galcen()
# PPXML Data Load
pmppmxl = pmvels.PMmeasurements(pmcatalog='PPMXL')
pmppmxl.to_space_velocties()
pmppmxl.UVW_to_galcen()

# First Triangle Plot
# Contains uncertainties in obs. PMs, R,Z, age
quants = ['PMRA_ERR', 'PMDEC_ERR', 'RC_GALR', 'RC_GALZ', 'cannon_AGE']
tridata = data_triangle.get_data(pmucac, *quants)
labels = [r'$\delta_{\mathrm{PMRA}} [mas/yr]$',
          r'$\delta_{\mathrm{PMDEC}} [mas/yr]$', 'R [kpc]', 'Z [kpc]',
          'age [Gyr]', 'W [km/s]']

plotdata = np.column_stack([tridata, pmucac.get_Ws()])
ranges1 = [(0, 5.0), (0, 5.0), 0.99, 0.99, 0.99, 0.99]
figure = data_triangle.triangle_plot(plotdata, labels, range=ranges1,
                                     quantiles=[.16, .5, .85])
figname = 'pmerr_r_z_age'
figure.savefig('{0}.eps'.format(figname), format='eps')
figure.savefig('{0}.png'.format(figname), format='png')

# U/D, V/D for outer disk
# UCAC First Get outer disk
glon_mask = (pmucac.get_col('GLON') > 165.0) & (pmucac.get_col('GLON') < 195.0)
glat_mask = (pmucac.get_col('GLAT') > -15.0) & (pmucac.get_col('GLAT') < 15.0)
outerdisk_mask = glon_mask & glat_mask
pmra_err = pmucac.get_col('PMRA_ERR')
pmdec_err = pmucac.get_col('PMDEC_ERR')
pmra_mask = pmra_err < 10.0
pmdec_mask = pmdec_err < 10.0
pm_mask = pmra_mask & pmdec_mask

fig = plt.figure()
ax = fig.add_subplot(111)
distances = pmucac.get_col('RC_DIST')
ax.scatter((pmucac.vRvTvZ_g[0]/distances)[outerdisk_mask & pm_mask],
           (pmucac.vRvTvZ_g[1]/distances)[outerdisk_mask & pm_mask])
ax.set_ylim(-50, 260)
ax.set_xlim(-55, 55)
ax.set_xlabel('U/D [km/s/kpc]', fontsize=14)
ax.set_ylabel('V/D [km/s/kpc]', fontsize=14)
fileout = 'U_V_outerdisk_ucac'
plt.savefig('{0}.png'.format(fileout), format='png')
plt.savefig('{0}.eps'.format(fileout), format='eps')

# U/D, V/D for outer disk
# PPXML First Get outer disk
glon_mask = (pmppmxl.get_col('GLON') > 165.0) & (pmppmxl.get_col('GLON') < 195.0)
glat_mask = (pmppmxl.get_col('GLAT') > -15.0) & (pmppmxl.get_col('GLAT') < 15.0)
outerdisk_mask = glon_mask & glat_mask
pmra_err = pmppmxl.get_col('PMRA_ERR')
pmdec_err = pmppmxl.get_col('PMDEC_ERR')
pmra_mask = pmra_err < 10.0
pmdec_mask = pmdec_err < 10.0
pm_mask = pmra_mask & pmdec_mask

fig = plt.figure()
ax = fig.add_subplot(111)
distances = pmppmxl.get_col('RC_DIST')
ax.scatter((pmppmxl.vRvTvZ_g[0]/distances)[outerdisk_mask & pm_mask],
           (pmppmxl.vRvTvZ_g[1]/distances)[outerdisk_mask & pm_mask])
ax.set_ylim(-50, 260)
ax.set_xlim(-55, 55)
ax.set_xlabel('U/D [km/s/kpc]', fontsize=14)
ax.set_ylabel('V/D [km/s/kpc]', fontsize=14)
fileout = 'U_V_outerdisk_ppxml'
plt.savefig('{0}.png'.format(fileout), format='png')
plt.savefig('{0}.eps'.format(fileout), format='eps')

# Positions, Velocities Triangle Plots
plotdata2 = np.column_stack([pmucac.get_col('RC_GALR'),
                             pmucac.get_col('RC_GALZ'), pmucac.vRvTvZ_g[0],
                             pmucac.vRvTvZ_g[1], pmucac.vRvTvZ_g[2],
                             pmucac.spacevel_uncer_var_tensor[:, 0, 0],
                             pmucac.spacevel_uncer_var_tensor[:, 1, 1],
                             pmucac.spacevel_uncer_var_tensor[:, 2, 2],
                             pmucac.get_col('cannon_AGE')])
labels2 = ['R', 'Z', 'U', 'V', 'W', r'$\delta_{\mathrm{U}}$',
           r'$\delta_{\mathrm{V}}$', r'$\delta_{\mathrm{W}}$', 'age']
figure2 = data_triangle.triangle_plot(plotdata2, labels2, range=[.99]*9,
                                      quantiles=[.16, .5, .85])
figname = 'pos_vel_age'
figure2.savefig('{0}.eps'.format(figname), format='eps')
figure2.savefig('{0}.png'.format(figname), format='png')
