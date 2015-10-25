import numpy as np
import synthetic_data as synth
import model
import pm_to_velocities as pmvels

# sample with flat age and radii ranges
# sdata_flat_b35_sw0100_iRd_07
N = 2500
rand1 = np.random.random(N)
rand2 = np.random.random(N)
minage = 0.1
maxage = 13.0
dage = maxage - minage
ages = minage + rand1*dage
minrad = 6.0
maxrad = 10.0
drad = maxrad - minrad
radii = minrad + rand2*drad

def mkplotname(datacont, source='flat'):
    name = 'sdata_{}_b{}_sw0{}_iRd{}_sig{}'.format(source, datacont.beta,
                                                   datacont.SW0,
                                                   datacont.inv_Rd,
                                                   datacont.sig2Wmax)
    return name

def mkflatplots():
    fakedata = synth.fidsynth(ages, radii, beta=.35, inv_Rd=.07, SW0=100.0)
    fakedata.genWs(sig2Wmax=350.0)
    fmod = model.fidmodel(plotname='sdata_flat_b35_sw0100_iRd_07_sig350')
    fmod.init_emcee()
    hparams = model.HyperParams(fakedata.sigma2Ws)
    fmod.run_emcee(fakedata, hparams)
    fmod.plot_emcee_params(truths=[100.0, .35, .07])

    # same age, radius sample, and params, smaller errors
    fakedata = synth.fidsynth(ages, radii, beta=.35, inv_Rd=.07, SW0=100.0)
    fakedata.genWs(sig2Wmax=150.0)
    fmod = model.fidmodel(plotname='sdata_flat_b35_sw0100_iRd_07_sig150')
    fmod.init_emcee()
    hparams = model.HyperParams(fakedata.sigma2Ws)
    fmod.run_emcee(fakedata, hparams)
    fmod.plot_emcee_params(truths=[100.0, .35, .07])

    # same age, radius sample, and params, larger errors
    fakedata = synth.fidsynth(ages, radii, beta=.35, inv_Rd=.07, SW0=100.0)
    fakedata.genWs(sig2Wmax=600.0)
    fmod = model.fidmodel(plotname='sdata_flat_b35_sw0100_iRd_07_sig600')
    fmod.init_emcee()
    hparams = model.HyperParams(fakedata.sigma2Ws)
    fmod.run_emcee(fakedata, hparams)
    fmod.plot_emcee_params(truths=[100.0, .35, .07])

    # same age, radius sample different params
    # sdata_flat_b2_sw0250_iRd_2

    fakedata = synth.fidsynth(ages, radii, beta=.2, inv_Rd=.2, SW0=250.0)
    fakedata.genWs(sig2Wmax=350.0)
    fmod = model.fidmodel(plotname='sdata_flat_b2_sw0250_iRd_2_sig350')
    fmod.init_emcee()
    hparams = model.HyperParams(fakedata.sigma2Ws)
    fmod.run_emcee(fakedata, hparams)
    fmod.plot_emcee_params(truths=[250.0, .2, .2])

# Now using pmucac values for age, radii
pmucac = pmvels.PMmeasurements(maxheight=0.7, maxpmuncer=2.5)
pmucac.to_space_velocties()
pmucac.UVW_to_galcen()
pmdata = pmucac.get_tau_radii_vZg_sigma2Ws_container()

def mkfromdata_plots():
    # Typical errors
    fakedata = synth.fidsynth(pmdata.ages, pmdata.radii, beta=.2,
                            inv_Rd=.2, SW0=250.0)
    fakedata.genWs(sig2Wmax=350.0)
    fmod = model.fidmodel(plotname='sdata_frompm_b2_sw0250_iRd_2_sig350')
    fmod.init_emcee()
    hparams = model.HyperParams(fakedata.sigma2Ws)
    fmod.run_emcee(fakedata, hparams)
    fmod.plot_emcee_params(truths=[250.0, .2, .2])

    # smaller errors
    fakedata = synth.fidsynth(pmdata.ages, pmdata.radii, beta=.2,
                            inv_Rd=.2, SW0=250.0)
    fakedata.genWs(sig2Wmax=75.0)
    fmod = model.fidmodel(plotname='sdata_frompm_b2_sw0250_iRd_2_sig75')
    fmod.init_emcee()
    hparams = model.HyperParams(fakedata.sigma2Ws)
    fmod.run_emcee(fakedata, hparams)
    fmod.plot_emcee_params(truths=[250.0, .2, .2])

    # larger errors
    fakedata = synth.fidsynth(pmdata.ages, pmdata.radii, beta=.2,
                            inv_Rd=.2, SW0=250.0)
    fakedata.genWs(sig2Wmax=600.0)
    fmod = model.fidmodel(plotname='sdata_frompm_b2_sw0250_iRd_2_sig600')
    fmod.init_emcee()
    hparams = model.HyperParams(fakedata.sigma2Ws)
    fmod.run_emcee(fakedata, hparams)
    fmod.plot_emcee_params(truths=[250.0, .2, .2])


def mklowerr_plots():
    fakedata = synth.fidsynth(ages, radii, beta=.2, inv_Rd=.2, SW0=250.0)
    fakedata.genWs(sig2Wmax=50.0)
    fname = mkplotname(fakedata, source='flat')
    fmod = model.fidmodel(plotname=fname)
    fmod.init_emcee()
    hparams = model.HyperParams(fakedata.sigma2Ws)
    fmod.run_emcee(fakedata, hparams)
    fmod.plot_emcee_params(truths=[250.0, .2, .2])

if __name__ == '__main__':
    mklowerr_plots()
