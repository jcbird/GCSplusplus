#######
# module to combine RC catalog and Ness Catalog with Cannon ages.
######

from __future__ import print_function
import os
##Environment variables for apogee
import sys
os.environ['APOGEE_REDUX']='v603'
os.environ['APOGEE_DATA']='/Users/jquark/obs_data/apogee/'
import apogee.tools.read as apread
import numpy as np
#import galpy.util.bovy_coords as bovy_coords
#import matplotlib as mpl
#from matplotlib import pyplot as plt
#import numpy as np
import numpy.lib.recfunctions as rfn


def returncat():
    homedir = os.path.expanduser('~')
    apogee_data_dir = os.path.join(homedir, 'obs_data/apogee')
    ness_age_cat_file = os.path.join(apogee_data_dir, 'HWR_redclump_sample.txt')

    RCcat = apread.rcsample()
    #Now need to read in ascii Ness Sample
    fields = ['APOGEE_ID', 'RC_DIST', 'RC_GALR', 'RC_GALPHI', 'RC_GALZ', 'TEFF', 'LOGG', 'FE_H', 'ALPHAFE', 'RC_GALR']
    #Using the Rccatalog to provide data types for making record array with ness catalog
    # The last 'RC_GALR' is just to steal the datatype to make the age array
    ness_dtypes=[RCcat[field].dtype for field in fields]
    ness_fields = ['APOGEE_ID', 'RC_DIST', 'RC_GALR', 'RC_GALPHI', 'RC_GALZ', 'TEFF', 'LOGG', 'FE_H', 'ALPHAFE', 'AGE']
    ndtypes = [(i,j) for i,j in zip(ness_fields, ness_dtypes)]
    ness_cat = np.genfromtxt(ness_age_cat_file, dtype=ndtypes)

    ###can write out record array here and now

    ###
    ### After this, change the columns for the merging
    print((RCcat['APOGEE_ID']==ness_cat['APOGEE_ID']).all())
    print('RC catalog and Cannon catalog have matching APOGEE_IDs')
    print([(np.allclose(RCcat[i],ness_cat[i])) for i in ness_fields[1:5]])
    print('RC catalog and Cannon catalog have matching {}'.format(','.join(ness_fields[1:5])))
    newnames = ["cannon_"+i for i in ness_fields]
    ness_cat.dtype.names=newnames

    RC_cannon_cat = rfn.merge_arrays([RCcat, ness_cat], asrecarray=True, flatten=True)
    return RC_cannon_cat




