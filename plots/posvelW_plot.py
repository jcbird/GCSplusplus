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
# pmppmxl = pmvels.PMmeasurements(pmcatalog='PPMXL')
# pmppmxl.to_space_velocties()
# pmppmxl.UVW_to_galcen()

# Positions, Velocities Triangle Plots
plotdata2 = np.column_stack([pmucac.get_col('RC_GALR'),
                             pmucac.get_col('RC_GALZ'),
                             pmucac.vRvTvZ_g[2],
                             pmucac.spacevel_uncer_var_tensor[:, 2, 2],
                             pmucac.get_col('cannon_AGE')])
labels2 = ['R', 'Z', 'W', r'$\delta_{\mathrm{W}}$', 'age']
ranges = [.97, .97, .97, (0, 300), .97]
figure2 = data_triangle.triangle_plot(plotdata2, labels2, range=ranges,
                                      quantiles=[.16, .5, .85])
figname = 'posvelage_W'
figure2.savefig('{0}.eps'.format(figname), format='eps')
figure2.savefig('{0}.png'.format(figname), format='png')
