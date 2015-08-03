
import numpy as np
import galpy
from galpy.util import bovy_coords as bcoords
import combinecats

catalog = combinecats.returncat()

class PMmeasurements:
    """
    Class to contain PMs and errors. This will make it each to switch between UCAC and PMXXL later.
    """
    def __init__(self,pmcatalog='PPMXL', RCcatalog=catalog):
        self.pmcatalog = pmcatalog
        self.catalog = catalog
        self._colname_mod = '_PPMXL' if pmcatalog=='PPMXL' else ''
        self.degree = True #RA,DEC in degrees
        self._generatemask()
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
        self.data['pmra'] = self.catalog[self._colname('PMRA')][self.mask]
        self.data['pmra_err'] = self.catalog[self._colname('PMRA_ERR')][self.mask]
        self.data['pmdec'] = self.catalog[self._colname('PMDEC')][self.mask]
        self.data['pmdec_err'] = self.catalog[self._colname('PMDEC_ERR')][self.mask]
        self.data['distance'] = self.catalog['RC_DIST'][self.mask]
        self.data['vlos'] = self.catalog['VHELIO_AVG'][self.mask]
        self.data['RA'] = self.catalog['RA'][self.mask]
        self.data['DEC'] = self.catalog['DEC'][self.mask]
        self.data['GLAT'] = self.catalog['GLAT'][self.mask]
        self.data['GLON'] = self.catalog['GLON'][self.mask]
        self.data['RC_GALR'] = self.catalog['RC_GALR'][self.mask]
        self.data['RC_GALZ'] = self.catalog['RC_GALZ'][self.mask]
        self.data['RC_GALPHI'] = self.catalog['RC_GALPHI'][self.mask]

    def covar_pmrapmdec(self):
        pmrapmdec_corrcoefs = np.corrcoef(self.data['pmra_err'], self.data['pmdec_err'])
        pmra_err = self.data['pmra_err']
        pmdec_err = self.data['pmdec_err']
        ###step 2: propogate corresponding errors
        self.covar_pmradec = np.array([[pmra_err**2, pmrapmdec_corrcoefs[0,1]*pmra_err*pmdec_err], [pmrapmdec_corrcoefs[1,0]*pmra_err*pmdec_err, pmdec_err**2]]).T  #transpose for shape for galpy

    def pmrapmdec_to_pmllpmbb(self):
        ###step 1: convert PM radec to PM l,b
        self.pmll_pmbb = bcoords.pmrapmdec_to_pmllpmbb(self.data['pmra'], self.data['pmdec'], self.data['RA'],self.data['DEC'], degree=self.degree, epoch=2000.0)


    def to_space_velocties(self):
        """
        Wrapper around galpy, and wrapper to go through steps
        """
        #Build covar_pmrapmdec
        print (self.data['GLON'][:10], "Two")
        pmrapmdec_corrcoefs = np.corrcoef(self.data['pmra_err'], self.data['pmdec_err'])
        #easier to read
        print (self.data['GLON'][:10], "Three")
        pmra_err = self.data['pmra_err']
        pmdec_err = self.data['pmdec_err']
        ###step 2: propogate corresponding errors
        covar_pmradec = np.array([[pmra_err**2, pmrapmdec_corrcoefs[0,1]*pmra_err*pmdec_err], [pmrapmdec_corrcoefs[1,0]*pmra_err*pmdec_err, pmdec_err**2]]).T  #transpose for shape for galpy
        print (self.data['GLON'][:10], "4")
        self.covar_pmllpmbb = bcoords.cov_pmrapmdec_to_pmllpmbb(covar_pmradec,self.data['RA'], self.data['DEC'], degree=self.degree,epoch=2000.0)
        #### Make galactocentric X,Y,Z for spacevel calc. Include Zsun
        print (self.data['GLON'][:10], "5")

        #galXYZ = bcoords.cyl_to_rect(self.data['RC_GALR'], self.data['RC_GALPHI'], self.data['RC_GALZ'])
        self.XYZ = np.array(bcoords.lbd_to_XYZ(self.data['GLON'], self.data['GLAT'], self.data['distance'], degree=True))
        ###step 3: convert l,b,d,rv to space velocities
        ### Replacing l,b,d with galactocentric coords
        print (self.data['GLON'][:10], "6")
        self.spacevels = bcoords.vrpmllpmbb_to_vxvyvz(self.data['vlos'],self.pmll_pmbb[:,0]*np.cos(np.radians(self.data['GLAT'])), self.pmll_pmbb[:,1], self.data['GLON'], self.data['GLAT'], self.data['distance'], degree=self.degree)
        ###step 4: propagate errors as welll
        ## assume 5% distance errors
        print (self.data['GLON'][:10], "6")
        dist_err = 0.05 * self.data['distance']
        self.spacevel_err = bcoords.cov_dvrpmllbb_to_vxyz(self.data['distance'], dist_err, np.zeros_like(self.data['vlos']), self.pmll_pmbb[:,0], self.pmll_pmbb[:,1], self.covar_pmllpmbb, self.data['GLON'], self.data['GLAT'], degree=self.degree)
        print (self.data['GLON'][:10], "7")

        self.spacevel_galcen = bcoords.vxvyvz_to_galcencyl(self.spacevels[:,0], self.spacevels[:,1], self.spacevels[:,2], self.data['RC_GALR'], self.data['RC_GALPHI'], self.data['RC_GALZ'], vsun=[-11.1,30.24*8.,7.2], galcen=True)
        print (self.data['GLON'][:10], "8")
        self.rectgal = bcoords.sphergal_to_rectgal(self.data['GLON'], self.data['GLAT'], self.data['distance'], self.data['vlos'], self.pmll_pmbb[:,0], self.pmll_pmbb[:,1], degree=False)





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
