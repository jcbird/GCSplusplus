import numpy as np
import synthetic_data as synth
import pm_to_velocities as pmvels

# sample with flat age and radii ranges
# sdata_flat_b35_sw0100_iRd_07
N = 5000
rand1 = np.random.random(N)
rand2 = np.random.random(N)
minage = 0.1
maxage = 13.0
dage = maxage - minage
ages = minage + rand1*dage
minrad = 4.0
maxrad = 12.0
drad = maxrad - minrad
radii = minrad + rand2*drad

# Eventually do ranges of params
betas = [.2, .35, .5]
sw0s = [100., 200., 400.0]
inv_Rds = [0.0, .15, .3]

# for now
betas = betas[1:]
sw0s = sw0s[1]
inv_Rds = inv_Rds[1]

sigs = [0.0, 49.0, 100.0, 400.0, 900.0]

for beta in betas:
    fakedata = synth.fidsynth(ages, radii, beta=beta, SW0=sw0s, inv_Rd=inv_Rds)
    fakedata.trueWs()
    fakedata.gen_sigma2Ws(scale=5.0, sig2Wmin=1000.0)
    # First, 1 m/s errors in velocity
    for sigma in sigs:
        fakedata.sigma2Ws = sigma * np.random.rand(N)
        fakedata.addnoise()
        fakedata.run_emcee()
        summary = 'flat ages; flat radii; flat sigma lim={0}'.format(sigma)
        fakedata.savetoascii(description=summary, asciifile='datatests.stats5')

# Now with ages, radii pulled from data
pmucac = pmvels.PMmeasurements(maxheight=0.7, maxpmuncer=2.5)
pmucac.to_space_velocties()
pmucac.UVW_to_galcen()
pmdata = pmucac.get_tau_radii_vZg_sigma2Ws_container()

for beta in betas:
    fakedata = synth.fidsynth(pmdata.ages, pmdata.radii, beta=beta, SW0=sw0s,
                              inv_Rd=inv_Rds)
    fakedata.trueWs()
    fakedata.gen_sigma2Ws(scale=5.0, sig2Wmin=1000.0)
    # First, 1 m/s errors in velocity
    for sigma in sigs:
        fakedata.sigma2Ws = sigma * np.random.rand(fakedata.N)
        fakedata.addnoise()
        fakedata.run_emcee()
        summary1 = 'ages; radii from pmucac 2.5mas/yr uncer;'
        summary2 = ' flat sigma lim={0}'.format(sigma)
        summary = summary1 + summary2
        fakedata.savetoascii(description=summary, asciifile='datatests.stats5')
