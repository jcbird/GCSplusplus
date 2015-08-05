from __future__ import print_function
import numpy as np
import pm_to_velocities as pm_to_vels
import emcee
import triangle
from collections import namedtuple
import sys

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

def init_emcee(init_guess, nwalkers):
    """
    - init_guess: shape(ndim), initial guesses of params
    """
    randperturb_di = [np.random.normal(0,0.015,nwalkers)
            for i in range(ndim)]
    p0 = [(init_guess[0]*(1.+i),init_guess[1]*(1.+j),init_guess[2]*(1.+k)) for i,j,k in zip(*randperturb_di)] # (km/s,beta,kpc**-1)
    return p0

def run_emcee(ndim, nwalkers, p0, lnprob_func, lnprob_args):
    if (len(lnprob_args)==2):
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_func,
                args=lnprob_args)
        sampler.run_mcmc(p0,45000)
        return sampler

def plot_and_run_emcee(nexec,ndim, nwalkers, p0, lnprob_func, lnprob_args):
    """
    - nexec: Number of times to continue looping emcee and making plot
    """
    sampler = run_emcee(ndim, nwalkers, p0, lnprob_func, lnprob_args)
    mk_triangle_plot(sampler)
    if nexec>1:
        for i in xrange(nexec):
            sampler = run_emc(ndim, nwalkers, p0, lnprob_func, lnprob_args)
            mk_triangle_plot(sampler)

def lnlh(data, params, hyperparams):
    """
    # Inputs
    - data: ages, galacto-centric radius, vertical velocity
            shape (3,N)
    - params: variance_W0, beta_W, inv_scalelen_W
              shape (3)
    - hyperpars: Various getters,read source

    # Bugs
    - This is really brittle, it is only there to serve a single customer

    """
    ages = data.ages
    radii = data.radii
    Ws = data.Ws
    variance_W0, beta_W, inv_scalelen_W = params
    #print ('betaW: {}'.format(beta_W))
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
    if beta_W < 0.05:   # need to justify prior (repeat everywhere)
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

def mk_triangle_plot(sampler, nstart=500):
    samples = sampler.chain[:,nstart:,:].reshape((-1,3))
    fig = triangle.corner(samples, labels=[r"$S_{W_0}$ [km$^2$/s$^2$]",
        r"$\beta$", r"$R_{W_0}^{-1}$ [kpc]"])
    fig.savefig("triangle0.png")


if __name__ == "__main__":
    pmdata = pm_to_vels.PMmeasurements(RCcatalog = pm_to_vels.catalog, biascorrect='dqsou')
    pmdata.to_space_velocties()
    pmdata.UVW_to_galcen()
    data = pmdata.get_tau_radii_vZg_sigma2Ws_container(max_uncer_variance=200.)
    #data = pmdata.get_tau_radii_vZg_sigma2Ws_container()#max_uncer_variance=200.)
    hyperparams = HyperParams(data.sigma2Ws)
    ndim = 3
    nwalkers = 20
    init_guess = (225.0, 0.2, 0.06) #km/s,beta,kpc**-1 param guesses
    p0 = init_emcee(init_guess,nwalkers)
    plot_and_run_emcee(nexec=1, ndim=ndim, nwalkers=nwalkers, p0=p0, lnprob_func=lnprob, lnprob_args=[data,hyperparams])
    #emcee_sampler = run_emcee(ndim=ndim, nwalkers=nwalkers, p0=p0, lnprob_func=lnprob, lnprob_args=[data,hyperparams])
    #mk_triangle_plot(emcee_sampler)


##ndim = 3   #HACK, need to fix!!!
####samples = emcee_sampler.chain[:,1000:,:].reshape((-1,ndim))

#result = lnprob(data,params,hyperparams)
#print ('lnprob: {}'.format(result))
#pass
