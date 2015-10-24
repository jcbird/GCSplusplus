"""pmuncer_limits

Makes a plot of PM uncertainty threshold vs. number of stars with
20% fractional W errors and total number of stars. Used to determine the
best threshold.
"""
import os
import sys
from matplotlib import pyplot as plt
import numpy as np
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import pm_to_velocities as pmvels

# UCAC Data Load
pmucac = pmvels.PMmeasurements(maxheight=0.75)

pmuncers = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]  # mas/yr

allfrac_errW = []

for ii, pmerr in enumerate(pmuncers):
    pmucac.set_maxpmuncer(pmerr)
    pmucac._generatemask()
    pmucac.to_space_velocties()
    pmucac.UVW_to_galcen()
    W = pmucac.get_Ws()
    errW = pmucac.get_sigma2Ws()
    errWfrac = errW/W**2
    allfrac_errW.append(errWfrac)

numltpt2 = [np.sum(efracW < .2) for efracW in allfrac_errW]
ntot = [len(efracW) for efracW in allfrac_errW]

with plt.style.context(['ggplot', 'h277paper']):
    fig = plt.figure()
    ax = plt.subplot2grid((5, 1), (0, 0), rowspan=4)
    ax.plot(pmuncers, numltpt2, label=r'$\sigma_W/W < 0.2$')
    ax.plot(pmuncers, ntot, label='Total stars in sample')
    ax.legend(loc='upper left', frameon=True)
    ax.set_ylabel('Number of Stars')
    fig.subplots_adjust(hspace=0)
    ax.tick_params(axis='x', bottom='off', labelbottom='off')

    ax.grid(b=True, axis='both', c='.7')
    ax2 = plt.subplot2grid((5, 1), (4, 0), sharex=ax)
    ax2.plot(pmuncers,
             np.array(numltpt2, dtype=float)/np.array(ntot, dtype=float),
             c='k')
    ax2.set_ylim(0, 1.0)
    ax2.grid(b=False, axis='y')
    ax2.grid(b=True, c='.7', axis='x')
    ax2.set_ylabel(r'N($\sigma_W/W< 0.2$)/N$_{tot}$', fontsize=10)
    ax2.set_xlabel(r'$\sigma_{\mu_\alpha}, \sigma_{\mu_\beta} < [mas/yr]$',
                   labelpad=10)
    fname = 'pmuncer_limits'
    plt.savefig(fname+'.png', format='png')
    plt.savefig(fname+'.eps', format='eps')
