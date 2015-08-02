import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pm_to_velocities as pmtospacevel
import pandas as pd

pmstars = pmtospacevel.PrMoMeasuremnts(RCcatalog=pmtospacevel.catalog)
pmstars.to_space_velocties()

galr = np.array(pmtospacevel.catalog['RC_GALR'][pmstars.mask]).byteswap().newbyteorder()
galz = np.array(pmtospacevel.catalog['RC_GALZ'][pmstars.mask]).byteswap().newbyteorder()
age = np.array(pmtospacevel.catalog['cannon_AGE'][pmstars.mask]).byteswap().newbyteorder()

bovy_vZ = np.array(pmtospacevel.catalog['GALVZ_PPMXL'][pmstars.mask]).byteswap().newbyteorder()
#Build spatial masks
#central locations to take measurement

cenlocs = [(i,0.0) for i in np.arange(7.0,9.1,0.5)]
cen_radii= np.arange(7.0,9.1,0.5)
distfromcen = [np.sqrt((galr-cenR)**2 + (galz-0.)**2) for cenR in cen_radii]

central_masks = [dists<1.0 for dists in distfromcen]
vW = pmstars.spacevels[:,2] #actually Z velocity, but same (Z)
vW_err = pmstars.spacevel_err[:,2,2]

#all endianness needs to change based of Ness stuff
#end_vW = vW.byteswap().newbyteorder()
#end_vW_err = vW_err.byteswap().newbyteorder()
#end_age = age.byteswap().newbyteorder()

data_dict = {'vW':vW, 'vW_err':vW_err, 'age':age, 'bovy_vz':bovy_vZ, 'galr':galr, 'galz':galz}
dataframe = pd.DataFrame(data_dict)

### Try just doing GCS first
velerr_cut = vW_err<100.
df = dataframe[central_masks[1] & velerr_cut]
dfbyage = df.groupby(pd.qcut(df.age,10, labels=False))
plt.figure()
plt.scatter(dfbyage.age.mean(), dfbyage.vW.std())
plt.show()

