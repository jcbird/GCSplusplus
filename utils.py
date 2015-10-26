#######
# module to combine RC catalog and Ness Catalog with Cannon ages.
######

from __future__ import print_function
import os
import sys
import numpy.lib.recfunctions as rfn
import numpy as np
import astropy.io.fits as pyfits

# Data directories and files
_homedir = os.path.expanduser('~')
_apogee_data_dir = os.path.join(_homedir, 'obs_data/apogee')
_cannon_cat_file = "redclump_sample_A.txt"
_ucac_apogee_file = "allStar-v603_UCAC4.fits"
if os.path.exists(_apogee_data_dir):
    pass
else:  # Some other computer
    from datapath import datadir as _apogee_data_dir

# Environment variables for apogee module
os.environ['APOGEE_REDUX'] = 'v603'
os.environ['APOGEE_DATA'] = _apogee_data_dir
import apogee.tools.read as apread


def read_anders():
    ucac_apogee_file = os.path.join(_apogee_data_dir, _ucac_apogee_file)
    ucac_fits = pyfits.open(ucac_apogee_file)
    # change first H header name, it was duplicated
    ucac_fits[1].columns['H'].name = 'H_APO'
    ucac_data = ucac_fits[1].data
    return ucac_data


def read_ness_cat():
    """
    Just APOGEEID and Ages
    """
    ness_age_cat_file = os.path.join(_apogee_data_dir, _cannon_cat_file)
    ness_cat = np.genfromtxt(ness_age_cat_file, dtype="|S18,f4",
                             usecols=(0, -2), names=['ID', 'AGE'])  # ID, Gyr
    return ness_cat


def returncat():
    RCcat = apread.rcsample()
    # Now need to read in ascii Ness Sample
    fields = ['APOGEE_ID', 'RC_DIST', 'RC_GALR', 'RC_GALPHI', 'RC_GALZ',
              'TEFF', 'LOGG', 'FE_H', 'ALPHAFE', 'RC_GALR']
    # Using the Rccatalog to provide data types for making record array
    # with ness catalog
    # The last 'RC_GALR' is just to steal the datatype to make the age array
    ness_dtypes = [RCcat[field].dtype for field in fields]
    ness_cat = read_ness_cat()

    # TEST NESS CAT
    if (ness_cat['ID'] == RCcat['APOGEE_ID']).all():
        print("Ness Catalog and RC Catalog successfully merged.")
        print('RC catalog and Cannon catalog have matching APOGEE_IDs')
    else:
        print("Ness Catalog merge unsuccessful. You lose!")
        sys.exit()

    newnames = ["cannon_"+i for i in ness_cat.dtype.names]
    ness_cat.dtype.names = newnames

    # Just add in delta corrections from PPMXL
    pm_corr_file = os.path.join(_apogee_data_dir, 'qso_gal_pm_corrections')
    corr_fields = ['corr_RA', 'corr_DEC', 'dgaluRA', 'dgaluDEC',
                   'dqsouRA', 'dqsouDEC']
    corr_dtypes = ness_dtypes[1:7]
    ndtypes_corr = [(i, j) for i, j in zip(corr_fields, corr_dtypes)]
    pmcorr_cat = np.genfromtxt(pm_corr_file, dtype=ndtypes_corr)
    pmcorr_cat.sort(order='corr_RA')
    # Need to sort all against RA, then join

    print ('Sum difference of RA in correction file and RC cat: {0}'.format(
           np.sum(pmcorr_cat['corr_RA'] - RCcat['RA'])))
    print ('Sum difference of DEC in correction file and RC cat: {0}'.format(
           np.sum(pmcorr_cat['corr_DEC'] - RCcat['DEC'])))

    RC_cannon_cat = rfn.merge_arrays([RCcat, ness_cat, pmcorr_cat],
                                     asrecarray=True, flatten=True)
    return RC_cannon_cat
