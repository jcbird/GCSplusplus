import numpy as np
import galpy
from galpy.util import bovy_coords as bcoords
import utils
from collections import namedtuple

catalog = utils.returncat()


class PMmeasurements:
    """
    Class to contain PMs and errors. This will make it each to switch between UCAC and PMXXL later.

    # Inputs
    - biascorect  (string) 'dqsou', 'dgalu', 'ppmxl'
                  'dqsou' and 'dgalu' refer to interpolated corrections from Bertrand (see HWR email
                  from 2015-08-04) using either quasars or galaxies as a reference objects.
                  'ppmxl' uses mean velocity of entire ppmxl catalog as a function of position
                  PPMXL NOT IMPLIMENTED
    # Bugs
    - galpy:   rewrites the input l,b; need to fork and fix
        + Hack: change degree keyword to False for conversions (NASTY)

    # Decisions, todos

    - REWRITE data section to use getters
    - MAJOR Propagate uncertainty tensor from UVW to galactocentric frame
    - 2015-08-03 Assuming 5% distance errors to calc spacevels.
    - 2015-08-03 Fork Galpy and rewrite conv to galcencycl coords.
    """
    def __init__(self,pmcatalog='PPMXL', RCcatalog=catalog, biascorrect=None):
        self.pmcatalog = pmcatalog
        self.catalog = catalog
        self._colname_mod = '_PPMXL' if pmcatalog=='PPMXL' else ''
        self.degree = True #RA,DEC in degrees
        self.biascorrect = biascorrect
        self._generatemask()
        self.Data = namedtuple('Data', ['ages', 'radii', 'Ws', 'sigma2Ws'])
        self.grabdata()

    def _generatemask(self):
        """
        Uses 'match' arrays to build boolean mask that finds all stars with PM matches corresponding to the catalog used. Mask applied to catalog before conversion
        """
        self.mask = (self.catalog['PMMATCH{}'.format(self._colname_mod)]==1)

    def _colname(self,basename):
        return '{0}{1}'.format(basename,self._colname_mod)

    def grabdata(self, degree=True):
        """
        assign reasonable variable names to PM columns of RC catalog depending upon which survey is being used. Then pass to Galpy.
        """
        self.data = {}
        self.data['GALVZ'] = self.catalog[self._colname('GALVZ')][self.mask]
        self.data['pmra'] = self.catalog[self._colname('PMRA')][self.mask]
        self.data['pmra_uncer'] = self.catalog[self._colname('PMRA_ERR')][self.mask]
        self.data['pmdec'] = self.catalog[self._colname('PMDEC')][self.mask]
        self.data['pmdec_uncer'] = self.catalog[self._colname('PMDEC_ERR')][self.mask]
        self.data['distance'] = self.catalog['RC_DIST'][self.mask]
        self.data['vlos'] = self.catalog['VHELIO_AVG'][self.mask]
        self.data['RA'] = self.catalog['RA'][self.mask]
        self.data['DEC'] = self.catalog['DEC'][self.mask]
        self.data['GLAT'] = self.catalog['GLAT'][self.mask]
        self.data['GLON'] = self.catalog['GLON'][self.mask]
        self.data['RC_GALR'] = self.catalog['RC_GALR'][self.mask]
        self.data['RC_GALZ'] = self.catalog['RC_GALZ'][self.mask]
        self.data['RC_GALPHI'] = self.catalog['RC_GALPHI'][self.mask]

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
        pmrapmdec_corrcoefs = np.corrcoef(self.data['pmra_uncer'], self.data['pmdec_uncer'])
        pmra_uncer = self.data['pmra_uncer']
        pmdec_uncer = self.data['pmdec_uncer']
        ###step 2: propogate corresponding errors
        self.covar_pmradec = np.array([[pmra_uncer**2, pmrapmdec_corrcoefs[0,1]*pmra_uncer*pmdec_uncer], [pmrapmdec_corrcoefs[1,0]*pmra_uncer*pmdec_uncer, pmdec_uncer**2]]).T  #transpose for shape for galpy

    def conv_pmrapmdec_to_pmllpmbb(self):
        ###step 1: convert PM radec to PM l,b
        if self.biascorrect is not None:
            (dpmra,dpmdec) = self.get_pm_corr()#self.biascorrect)
            pmra = self.get_col('PMRA') + dpmra
            pmdec = self.get_col('PMDEC') + dpmdec
            self.pmll_pmbb = bcoords.pmrapmdec_to_pmllpmbb(pmra, pmdec, self.get_col('RA'),self.get_col('DEC'), degree=self.degree, epoch=2000.0)
        else:
            self.pmll_pmbb = bcoords.pmrapmdec_to_pmllpmbb(self.data['pmra'], self.data['pmdec'], self.data['RA'],self.data['DEC'], degree=self.degree, epoch=2000.0)

    def calc_covar_pmllpmbb(self):
        self.covar_pmllpmbb = bcoords.cov_pmrapmdec_to_pmllpmbb(self.covar_pmradec,self.data['RA'], self.data['DEC'], degree=self.degree,epoch=2000.0)

    def calc_spacevel(self):
        """
        Calculates U,V,W
        """
        self.spacevels = bcoords.vrpmllpmbb_to_vxvyvz(self.data['vlos'],self.pmll_pmbb[:,0]*np.cos(np.radians(self.data['GLAT'])), self.pmll_pmbb[:,1], self.data['GLON'], self.data['GLAT'], self.data['distance'], degree=self.degree)

    def calc_spacevel_uncer_var_tensor(self):
        """
        Right now assumes ZERO error in RV (km/s) and 5% distance errors
        """
        dist_uncer = 0.05 * self.data['distance']
        self.spacevel_uncer_var_tensor = bcoords.cov_dvrpmllbb_to_vxyz(self.data['distance'], dist_uncer, np.zeros_like(self.data['vlos']), self.pmll_pmbb[:,0], self.pmll_pmbb[:,1], self.covar_pmllpmbb, self.data['GLON'], self.data['GLAT'], degree=self.degree)

    def to_space_velocties(self):
        """
        Wrapper around galpy, and wrapper to go through steps to calc UVW
        """
        self.conv_pmrapmdec_to_pmllpmbb()
        self.calc_covar_pmrapmdec()
        self.calc_covar_pmllpmbb()
        self.XYZ = np.array(bcoords.lbd_to_XYZ(self.data['GLON'], self.data['GLAT'], self.data['distance'], degree=True))
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
        self.vRvTvZ_g = bcoords.vxvyvz_to_galcencyl(self.spacevels[:,0], self.spacevels[:,1], self.spacevels[:,2], self.data['RC_GALR'], self.data['RC_GALPHI'], self.data['RC_GALZ'], vsun=[-11.1,30.24*8.,7.2], galcen=True)
        ###BEWARE, HACK TIME
        self.data['GLON'] = np.degrees(self.data['GLON'])
        self.data['GLAT'] = np.degrees(self.data['GLAT'])

    def get_col(self, colname):
        """
        Important function to grab column from catalog while applying
        the PMMATCH mask.

        #TODO
        - propagate use of this function everywhere internally
        """
        return self.catalog[colname][self.mask]

    def get_ages(self):
        """
        HACK!!!! Taking absolute ages until ln(age) available
        """
        return np.abs(self.catalog['cannon_AGE'][self.mask])

    def get_Ws(self):
        return self.data['GALVZ'] #km/s from Jo

    def get_radii(self):
        return self.data['RC_GALR'] #km/s from Jo

    def get_sigma2Ws(self):
        return self.spacevel_uncer_var_tensor[:,2,2] #(km.s)**2

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
        age_cut = np.logical_and(self.get_ages()>0.1,self.get_ages()<12)
        combined_mask = sigma2Ws_mask & age_cut
        print(np.sum(combined_mask))
        data_container = self.Data(ages=self.get_ages()[combined_mask],
                radii=self.get_radii()[combined_mask],
                Ws=self.get_Ws()[combined_mask],
                sigma2Ws=self.get_sigma2Ws()[combined_mask])
        return data_container

###Abandoned, left for dead until further notice    
#print (self.data['GLON'][:10], "8")
#self.rectgal = bcoords.sphergal_to_rectgal(self.data['GLON'], self.data['GLAT'], self.data['distance'], self.data['vlos'], self.pmll_pmbb[:,0], self.pmll_pmbb[:,1], degree=False)
#




#need to convert PMs and errors into velocities and errors

###use pmrapmdec_to_pmllpmbb(pmra,pmdec,ra,dec,degree=False,epoch=2000.0)
###cov_pmrapmdec_to_pmllpmbb(cov_pmradec,ra,dec,degree=False,epoch=2000.0)
###step 3: convert l,b,d,rv to space velocities
###vrpmllpmbb_to_vxvyvz(vr,pmll,pmbb,l,b,d,XYZ=False,degree=False)
###step 4: propagate errors as welll
###cov_dvrpmllbb_to_vxyz(d,e_d,e_vr,pmll,pmbb,cov_pmllbb,l,b,
###                                  plx=False,degree=False)
###


# def vrpmllpmbb_to_vxvyvz(vr,pmll,pmbb,l,b,d,XYZ=False,degree=False):
