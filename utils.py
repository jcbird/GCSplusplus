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

_cannon_cat_file = "redclump_sample_A.txt"
_homedir = os.path.expanduser('~')
_apogee_data_dir = os.path.join(_homedir, 'obs_data/apogee')

def read_ness_cat():
    """
    Just APOGEEID and Ages
    """
    ness_age_cat_file = os.path.join(_apogee_data_dir, _cannon_cat_file)
    ness_cat = np.genfromtxt(ness_age_cat_file, dtype="|S18,f4", usecols=(0,-2), names=['ID', 'AGE']) # ID, Gyr
    return ness_cat

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
    ###ness_cat = np.genfromtxt(ness_age_cat_file, dtype=ndtypes)
    ness_cat = read_ness_cat()

    #TEST NESS CAT
    if (ness_cat['ID']==RCcat['APOGEE_ID']).all():
        print("Ness Catalog and RC Catalog successfully merged.")  #proper ordering
        print('RC catalog and Cannon catalog have matching APOGEE_IDs')
    else:
        print("Ness Catalog merge unsuccessful. You lose!")
        sys.exit()

    newnames = ["cannon_"+i for i in ness_cat.dtype.names]
    ness_cat.dtype.names=newnames

    ### Just add in delta corrections from PPMXL
    pm_corr_file = os.path.join(apogee_data_dir, 'qso_gal_pm_corrections')
    corr_fields = ['corr_RA', 'corr_DEC', 'dgaluRA', 'dgaluDEC', 'dqsouRA', 'dqsouDEC']
    corr_dtypes=ness_dtypes[1:7]
    ndtypes_corr = [(i,j) for i,j in zip(corr_fields, corr_dtypes)]
    pmcorr_cat = np.genfromtxt(pm_corr_file, dtype=ndtypes_corr)
    pmcorr_cat.sort(order= 'corr_RA')
    #### Need to sort all against RA, then join

    print ('Sum difference of RA in correction file and RC cat: {0}'.format(np.sum(pmcorr_cat['corr_RA'] - RCcat['RA'])))
    print ('Sum difference of DEC in correction file and RC cat: {0}'.format(np.sum(pmcorr_cat['corr_DEC'] - RCcat['DEC'])))

    RC_cannon_cat = rfn.merge_arrays([RCcat, ness_cat, pmcorr_cat], asrecarray=True, flatten=True)
    return RC_cannon_cat




