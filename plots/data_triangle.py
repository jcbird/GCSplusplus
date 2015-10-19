"""Triangle plots in data space

TODO Use PMtovels class to mask out bad data
"""
# Modules
import os
import sys
import numpy as np
import triangle
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pm_to_velocities as pmtovel


def get_data(pm_instance, *quantities):
    """

    # Inputs
    args (list): list of strings, must match fields in combined catalog
    pm_instance (class): class instance of PMmeasurements. Return from
                         load_cat

    # Returns
    array: numpy array for triangle plot. shape = (nsamp, ndim)
           nsamp is number of data pts. ndim is len(quantities)
    """
    data = [np.asarray(pm_instance.get_col(qty)) for qty in quantities]
    return np.column_stack(data)


def load_cat(**PMkwargs):
    """
    Loads up and returns PMmeasurements class from pm_to_velocities
    """
    pmcont = pmtovel.PMmeasurements(**PMkwargs)
    return pmcont


if (__name__ == '__main__') and (__package__ is None):
    # get_data('PMRA_ERR', etc.)
    pms = load_cat(maxpmuncer=20.0)  # uas/yr
    pms.to_space_velocties()
    pms.UVW_to_galcen()
    Wvel = pms.get_Ws()
    quants = ['PMRA_ERR', 'PMDEC_ERR', 'RC_GALR', 'RC_GALZ', 'cannon_AGE']
    tridata = get_data(pms, *quants)
    tridata = np.column_stack([tridata, Wvel])
    figure = triangle.corner(tridata, labels=quants+['W'])
    figure.savefig('data_correlations.png', format='png')
