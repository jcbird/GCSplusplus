import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pm_to_velocities as pmtospacevel
import pandas as pd

pmstars = pmtospacevel.PMmeasurements(RCcatalog=pmtospacevel.catalog, pmcatalog='UCAC')
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
velerr_cut = vW_err<75.
labels = ['{}'.format(cenloc) for cenloc in cen_radii]

fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(len(cen_radii)):
    df = dataframe[central_masks[i] & velerr_cut]
    equallogbins = 10**np.arange(0,1.14,0.125)  #half the estimated error
    equallogbins2 = 10**np.arange(0.0625,1.24,0.125)  #half the estimated error
    dfbyage = df.groupby(pd.cut(df.age,equallogbins, labels=False))
    dfbyage2 = df.groupby(pd.cut(df.age,equallogbins2, labels=False))
    sortage = np.argsort(np.hstack((dfbyage.age.mean(),dfbyage2.age.mean())))
    ax.plot(np.hstack((dfbyage.age.mean(),dfbyage2.age.mean()))[sortage], np.hstack((dfbyage.vW.std(),dfbyage2.vW.std()))[sortage], label=labels[i])
    #plt.scatter(dfbyage.age.mean(), dfbyage.vW.std(), label=labels[i])
    #plt.scatter(dfbyage2.age.mean(), dfbyage2.vW.std(), label=labels[i])
plt.legend(loc='upper left')
plt.show()
