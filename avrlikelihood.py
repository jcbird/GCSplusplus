from __future__ import print_function
import numpy as np
import pm_to_velocities as pm_to_vels
import emcee

class HyperParams():
    """
    """
    def __init__(self, sigma2Ws):
        self.age0 = 1.0 #Gyr
        self.sigma2Ws = sigma2Ws
        self.set_sigma2Ws(sigma2Ws)

    def set_sigma2Ws(self, sigma2Ws):
        self.sigmaWs = np.sqrt(sigma2Ws)
        self.sigma2Ws = self.sigma2Ws

    def get_age0(self):
        return self.age0

    def get_sigma2Ws(self):
        return self.sigma2Ws

def run_emcee(ndim, nwalkers, lnprob_func, lnprob_args):
    if (len(lnprob_args)==2):
        ivar = 1. / np.random.rand(ndim)
        p0 = [np.random.rand(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_func,
                args=lnprob_args)
        sampler.run_mcmc(p0,500)
        return sampler

def lnlh(data, params, hyperparams):
    """
    # Inputs
    - data: ages, galacto-centric radius, vertical velocity
            shape (3,N)
    - params: variance_W0, beta_W, scalelen_W
              shape (3)
    - hyperpars: Various getters,read source

    # Bugs
    - This is really brittle, it is only there to serve a single customer

    """
    ages = data.get_ages()
    radii = data.get_radii()
    Ws = data.get_Ws()
    variance_W0, beta_W, inv_scalelen_W = params
    print ('betaW: {}'.format(beta_W))
    variances_W = (variance_W0 * np.power((ages/hyperparams.get_age0()), beta_W) * np.exp(-2.*radii*inv_scalelen_W) + hyperparams.get_sigma2Ws())
    return -0.5*np.sum(Ws**2/variances_W) -0.5*np.sum(np.log(variances_W))

def lnprior(params, hyperparams):
    """
    Just enforcing constraints (for now)
    """
    variance_W0, beta_W, inv_scalelen_W = params #you're pathetic
    if variance_W0 < 0.:
        return -np.Inf
    if variance_W0 > 10000.:  # (km/s)**2
        return -np.Inf
    if beta_W < 0.15:   # need to justify prior (repeat everywhere)
        return -np.Inf
    if beta_W > 1.:
        return -np.Inf
    if inv_scalelen_W < -.1 :# kpc**-1
        return -np.Inf
    if inv_scalelen_W > .1 : # kpc**-1
        return -np.Inf
    return 0.

def lnprob(params, data, hyperparams):
    lnp = lnprior(params, hyperparams)
    if not np.isfinite(lnp):
        return -np.Inf
    return lnp + lnlh(data, params, hyperparams)


if __name__ == "__main__":
    data = pm_to_vels.PMmeasurements(RCcatalog = pm_to_vels.catalog)
    data.to_space_velocties()
    data.UVW_to_galcen()
    params = (10.,0.4,0.05) #km/s,beta,kpc**-1 
    #(variance_W0, beta_W, inv_scalelen_W) = params
    hyperparams = HyperParams(data.get_sigma2Ws())
    emcee_sampler = run_emcee(ndim=3, nwalkers=10, lnprob_func=lnprob,
            lnprob_args=[data,hyperparams])
    #result = lnprob(data,params,hyperparams)
    #print ('lnprob: {}'.format(result))
    #pass
