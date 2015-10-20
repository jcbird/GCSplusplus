import numpy as np
from galpy.util import bovy_coords as bcoords
import utils
from collections import namedtuple

catalog = utils.returncat()


class PMmeasurements(object):  # New style class
    """
    Class to contain PMs and errors. This will make it each to switch between
    UCAC and PMXXL later.

    # Inputs
    - biascorect  (string) 'dqsou', 'dgalu', 'ppmxl'
                  'dqsou' and 'dgalu' refer to interpolated corrections
                  from Bertrand (see HWR email from 2015-08-04) using either
                  quasars or galaxies as a reference objects.
                  'ppmxl' uses mean velocity of entire ppmxl catalog as a
                  function of position
                  PPMXL NOT YET IMPLIMENTED
    # Bugs
    - galpy:   rewrites the input l,b; need to fork and fix
        + Hack: change degree keyword to False for conversions (NASTY)

    # Decisions, todos

    - REWRITE data section to use getters
    - MAJOR Propagate uncertainty tensor from UVW to galactocentric frame
    - 2015-08-03 Assuming 5% distance errors to calc spacevels.
    - 2015-08-03 Fork Galpy and rewrite conv to galcencycl coords.
    """
    colmod_lookup = {'PPMXL': '_PPMXL', 'UCAC': ''}  # Dict lookup
    # for column names in RCcatalog

    def __init__(self, pmcatalog='UCAC', RCcatalog=catalog, biascorrect=None,
                 degree=True, maxheight=None, maxpmuncer=None):
        self.pmcatalog = pmcatalog
        self.catalog = catalog    # RC catalog w/ages
        self._pmcolname_mod = PMmeasurements.colmod_lookup.get(pmcatalog, '')
        # default: UCAC
        self.set_maxheight(maxheight)  # [kpc] to cut sample close to the plane
        self.set_maxpmuncer(maxpmuncer)  # [mas/yr] cut sample by PM uncer
        self.degree = degree  # RA,DEC in degrees, set to false otherwise
        self.biascorrect = biascorrect  # ignoring this, not using PPMXL
        self._generatemask()  # create mask over data to be used

    # Write now explicity writing getter/setter, should do property!
    def set_maxpmuncer(self, maxpmuncer):
        """
        Cut the catalog in PM uncertainty. Take everything by default
        Auto-updates mask.
        """
        if maxpmuncer is None:
            maxpmuncer = np.max((self.catalog[self._pmcolname('PMRA_ERR')],
                                 self.catalog[self._pmcolname('PMDEC_ERR')]))
        self.maxpmuncer = maxpmuncer
        # if getattr(self, 'maxheight', None) is not None:
        #     self._generatemask()

    def set_maxheight(self, maxheight):
        """
        Cut the catalog by Z height above plane.
        Auto-updates mask if maxpmuncer exists
        """
        if maxheight is None:
            maxheight = np.abs(self.catalog['RC_GALZ']).max()
        self.maxheight = maxheight
        # if getattr(self, 'maxpmuncer', None) is not None:
        #     self._generatemask()

    def _pmcolname(self, basename):
        """
        Puts PM catalog ID in string for getters
        """
        return '{0}{1}'.format(basename, self._pmcolname_mod)

    def _generatemask(self, addmask=None, init=False):
        """
        Uses 'match' arrays to build boolean mask that finds all stars with PM
        matches corresponding to the catalog used. Mask applied to catalog
        before conversion.

        Also makes sure those PMs are reasonable with > -100 mas/yr cutoff
        # Inputs
        - init: [False]  if True will reset mask and start anew

        # Checked quantities
        - pm uncertainties: always checked, by default 100 mas/yr cutoff
        - maxheight: optional, height of stars from midplane, uses RCcat
        - matches: Demands match from single PM catalog (for now)

        # If AGE mask needed, put here
        """
        matches = (self.catalog['PMMATCH{}'.format(self._pmcolname_mod)] == 1)
        # setting above 0 below makes sure of no -9999 vals
        goodpmRAerr = np.logical_and(self.catalog[self._pmcolname('PMRA_ERR')] >=
                                     0.,
                                     self.catalog[self._pmcolname('PMRA_ERR')] <=
                                     self.maxpmuncer)  # mas/yr
        goodpmDECerr = np.logical_and(self.catalog[self._pmcolname('PMDEC_ERR')] >=
                                      0.,
                                      self.catalog[self._pmcolname('PMDEC_ERR')] <=
                                      self.maxpmuncer)  # mas/yr
        goodHeight = (np.abs(self.catalog['RC_GALZ']) < self.maxheight)
        if addmask is not None:
            self.mask = np.logical_and.reduce((goodpmRAerr, goodpmDECerr,
                                               matches, goodHeight, addmask))
        else:
            self.mask = np.logical_and.reduce((goodpmRAerr, goodpmDECerr,
                                               matches, goodHeight))

    def get_pm_corr(self):
        """
        Getter for proper motion corrections (subtracting off mean PM of PPMXL)
        # Inputs
        - source: (string) 'dqsou', 'dgalu', 'ppmxl'
                  'dqsou' and 'dgalu' refer to interpolated corrections from Bertrand (see HWR email
                  from 2015-08-04) using either quasars or galaxies as a reference objects.
                  'ppmxl' uses mean velocity of entire ppmxl catalog as a function of position

        # Outputs
        - dpmra: (array)  delta offset in RA [mas/yr]
        - dpmdec: (array) delta offset in DEC [mas/yr]
        self.mask is applied to each array

        # TODO
        'ppmxl' NOT IMPLIMENTED YET
        """
        if self.biascorrect is not 'ppmxl':
            dpmra = self.get_col(self.biascorrect+'RA')
            dpmdec = self.get_col(self.biascorrect+'DEC')
            dpmra[np.isnan(dpmra)] = 0.0  # mas/yr , if NO correction available, set to 0.0
            dpmdec[np.isnan(dpmdec)] = 0.0  # mas/yr , if NO correction available, set to 0.0
            print('Bias-correction for PMs\nSource: {0}'.format(self.biascorrect))
            return(dpmra,dpmdec)
        else:
            print('ppmxl mean PM correction not available yet')

    def calc_covar_pmrapmdec(self):
        """
        Calculates covariance matrix of PM uncertainties at the input
        measurement level (RA,DEC)

        # Notes
        - Locally renamed variables to be clear about uncertainties vs. errors
        """
        pmra_uncer = self.get_col('PMRA_ERR')
        pmdec_uncer = self.get_col('PMDEC_ERR')
        pmrapmdec_corrcoefs = np.corrcoef(pmra_uncer, pmdec_uncer)
        # step 2: propogate corresponding errors
        covpm = np.array([[pmra_uncer**2,
                           pmrapmdec_corrcoefs[0, 1]*pmra_uncer*pmdec_uncer],
                          [pmrapmdec_corrcoefs[1, 0]*pmra_uncer*pmdec_uncer,
                           pmdec_uncer**2]]).T  # transpose for shape for galpy
        self.covar_pmradec = covpm  # For style reasons

    def conv_pmrapmdec_to_pmllpmbb(self):
        # step 1: convert PM radec to PM l,b
        pmra = self.get_col('PMRA')
        pmdec = self.get_col('PMDEC')
        if self.biascorrect is not None:
            (dpmra, dpmdec) = self.get_pm_corr()
            pmra += dpmra
            pmdec += dpmdec
        self.pmll_pmbb = bcoords.pmrapmdec_to_pmllpmbb(pmra, pmdec,
                                                       self.get_col('RA'),
                                                       self.get_col('DEC'),
                                                       degree=self.degree,
                                                       epoch=2000.0)

    def calc_covar_pmllpmbb(self):
        covpmllbb = bcoords.cov_pmrapmdec_to_pmllpmbb(self.covar_pmradec,
                                                      self.get_col('RA'),
                                                      self.get_col('DEC'),
                                                      degree=self.degree,
                                                      epoch=2000.0)
        self.covar_pmllpmbb = covpmllbb

    def calc_spacevel(self):
        """
        Calculates U,V,W
        """
        uvw = bcoords.vrpmllpmbb_to_vxvyvz(self.get_col('VHELIO_AVG'),
                                           self.pmll_pmbb[:, 0],
                                           self.pmll_pmbb[:, 1],
                                           self.get_col('GLON'),
                                           self.get_col('GLAT'),
                                           self.get_col('RC_DIST'),
                                           degree=self.degree)

        self.spacevels = uvw

    def calc_spacevel_uncer_var_tensor(self):
        """
        Right now assumes ZERO error in RV (km/s) and 5% distance errors
        """
        dist_uncer = 0.05 * self.get_col('RC_DIST')
        uncer_tensor = bcoords.cov_dvrpmllbb_to_vxyz(self.get_col('RC_DIST'),
                               dist_uncer,
                               np.zeros_like(self.get_col('VHELIO_AVG')),
                               self.pmll_pmbb[:, 0], self.pmll_pmbb[:, 1],
                               self.covar_pmllpmbb, self.get_col('GLON'),
                               self.get_col('GLAT'), degree=self.degree)
        self.spacevel_uncer_var_tensor = uncer_tensor

    def to_space_velocties(self):
        """
        Wrapper around galpy, and wrapper to go through steps to calc UVW
        """
        self.conv_pmrapmdec_to_pmllpmbb()
        self.calc_covar_pmrapmdec()
        self.calc_covar_pmllpmbb()
        self.XYZ = np.array(bcoords.lbd_to_XYZ(self.get_col('GLON'),
                                               self.get_col('GLAT'),
                                               self.get_col('RC_DIST'),
                                               degree=True))
        self.calc_spacevel()
        self.calc_spacevel_uncer_var_tensor()

    def UVW_to_galcen(self):
        """
        Converts UVW space vels to galactocentric frame
        Assumes R0=8 kpc, z0 = 0.025 kpc, vsun=[-11.1,30.24*8.,7.2] km/s
        Uses R, phi, z in kpc from RC catalog
        This matches JoBo's RC catalog assumptions

        # BUGS

        - galpy maipulates input GLON, GLAT
        - HACK: puts in back in degrees
        - EVENTUAL FIX: PR galpy with not mutating attribute
        # Output

        vRg, vTg, vZg  # km/s

        """
        self.vRvTvZ_g = bcoords.vxvyvz_to_galcencyl(self.spacevels[:, 0],
                                                    self.spacevels[:, 1],
                                                    self.spacevels[:, 2],
                                                    self.get_col('RC_GALR'),
                                                    self.get_col('RC_GALPHI'),
                                                    self.get_col('RC_GALZ'),
                                                    vsun=[-11.1, 30.24*8., 7.2],
                                                    galcen=True)
        # BEWARE, HACK TIME

    def get_col(self, colname):
        """
        Important function to grab column from catalog while applying
        the PMMATCH mask.

        #TODO
        - propagate use of this function everywhere internally
        """
        if ('PM' in colname) or ('GALV' in colname):
            # If PM measurement, decide if UCAC or PPXML
            return self.catalog[self._pmcolname(colname)][self.mask]
        else:
            return self.catalog[colname][self.mask]

    def get_ages(self):
        """
        Convenient wrapper, nothing more
        """
        return self.get_col('cannon_AGE')

    # THESE Getters should be properties ! TODO
    def get_Ws(self):
        return self.vRvTvZ_g[2]  #GALVZ equivalent #incorporates bias corrections!1
        #return self.get_col('GALVZ') #km/s from Jo

    def get_radii(self):
        return self.get_col('RC_GALR') #km/s from Jo

    def get_sigma2Ws(self):
        return self.spacevel_uncer_var_tensor[:,2,2] #(km.s)**2

    def height_cut(self, maxheight=0.5, heightcol='RC_GALZ'):
        """
        Function that updates data mask and cuts on height
        Eventually could look at different calculations of height

        # Input
        - maxheight: number of kpc (default = 0.5)

        # Notes
        - reference to self.catalog is correct. When creating mask, always
          do so on full catalog
        """
        self.set_maxheight(maxheight)

    def get_tau_radii_vZg_sigma2Ws_container(self,
                                             max_uncer_variance=10000.):
        """
        WEIRD NAME: reminds me that need to propagate sigma2Ws to sigma2vZg
        Propagates velocity uncertainty cut to data
        returns shape(4,N)
        (ages, radii, Ws, sigma2Ws)
        """
        self.sigma2W_uncer_cut = max_uncer_variance
        sigma2Ws_mask = self.get_sigma2Ws() < max_uncer_variance
        age_cut = np.logical_and(self.get_ages() > 0.1, self.get_ages() < 14)
        combined_mask = sigma2Ws_mask & age_cut
        self.combined_mask = combined_mask
        print("N:{}".format(np.sum(self.mask)))
        data_container = {'ages': self.get_ages(),
                          'radii': self.get_radii(),
                          'Ws': self.get_Ws(),
                          'sigma2Ws': self.get_sigma2Ws()}
        return data_container
