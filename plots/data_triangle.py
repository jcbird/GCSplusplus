"""Triangle plots in data space

TODO Use PMtovels class to mask out bad data
"""
# Modules
import os
import sys
import numpy as np
import triangle
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils

cat = utils.returncat()


def get_data(*quantities):
    """

    # Inputs
    args (list): list of strings, must match fields in combined catalog

    # Returns
    array: numpy array for triangle plot. shape = (nsamp, ndim)
           nsamp is number of data pts. ndim is len(quantities)
    """
    data = [np.asarray(cat[qty]) for qty in quantities]
    return np.column_stack(data)

if (__name__ == '__main__') and (__package__ is None):
    # get_data('PMRA_ERR', etc.)
    pass
