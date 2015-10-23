from __future__ import print_function
import numpy as np
import pm_to_velocities as pm_to_vels
import emcee
import corner


class HyperParams():
    """
    """
    def __init__(self, sigma2Ws):
        self.age0 = 1.0  # Gyr
        self.sigma2Ws = sigma2Ws
        self.set_sigma2Ws(sigma2Ws)

    def set_sigma2Ws(self, sigma2Ws):
        self.sigmaWs = np.sqrt(sigma2Ws)
        self.sigma2Ws = self.sigma2Ws

    def get_age0(self):
        return self.age0

    def get_sigma2Ws(self):
        return self.sigma2Ws


class fidmodel(object):
    """
    fidmodel

    Container for fiducial model and emcee setup. Emcee wrappers will
    call on this class.

    Model parameters are: SW0, Beta, R^-1

    SW0 = Normalization of variance of W vels
    Beta = Power law index of growth of variance as func. of age
    R^-1 = Inverse scale length of variance; assume variance is an exp.
           function of age.

    Model
    ---------

    SW(R, tau) = SW0 * (tau/tau0)^(2*Beta) * e^(-2. * (R - R0) * R^-1) + sigW^2

    R:      Galactocentric Radius
    tau:    Stellar age
    tau0:   Min. age or normalization of age (tau0=1 Gyr)
    R0:     Solar radius (R0=8 kpc)
    sigW:   Uncertainty in W velocity
    """

    def __init__(self, data, ndim=3, nwalkers=20, guess=(200.0, 0.35, 0.1),
                 plotname='triangle0'):
        """
        Initialize fidmodel instance.

        Parameters
        ----------------

        data : namedtuple
            Object containing radii, W vels, ages, and sigma_W^2
            attrs of object are 'radii','Ws','ages','sigma2Ws', respectively.
        ndim : int
            Number of dimensions (# of parameters)
        nwalkers : int
            Number of walkers for emcee
        guess : iterable
            Initial guess of parameters SW0, Beta, and R^-1.
            Units are (km/s)^2, unitless, and kpc, respectively

        Notes
        ----------------

        Hyperparams automatically assigned sigma2Ws, this can be made more
        versitile.

        """
        self.data = data
        self.ndim = ndim
        self.nwalkers = nwalkers
        self.guess = guess
        self.plotname = plotname
        self.hparams = HyperParams(data.sigma2Ws)

    def init_emcee(self):
        """
        - init_guess: shape(ndim), initial guesses of params
        """
        randperturb_di = [np.random.normal(0, 0.015, self.nwalkers)
                          for i in range(self.ndim)]
        p0 = [(self.guess[0]*(1.+i), self.guess[1]*(1.+j),
               self.guess[2]*(1.+k)) for i, j, k in zip(*randperturb_di)]
        # Units (km/s,beta,kpc**-1)
        self.p0 = p0

    def lnlh(self, params):
        """
        # Inputs
        - params: variance_W0, beta_W, inv_scalelen_W
                shape (3)
        # Bugs
        - Rewrite to class, this could break!!!
        - This is really brittle, it is only there to serve a single customer
        """
        ages = self.data.ages
        radii = self.data.radii
        Ws = self.data.Ws
        variance_W0, beta_W, inv_scalelen_W = params
        variances_W = (variance_W0 * np.power((ages / self.hparams.get_age0()),
                                              2. * beta_W) *
                       np.exp(-2. * (radii - 8.0) * inv_scalelen_W) +
                       self.hparams.get_sigma2Ws())
        return (-0.5 * np.sum(Ws * Ws / variances_W) -
                0.5 * np.sum(np.log(variances_W)))

    @staticmethod
    def lnprior(params):
        """
        Just enforcing constraints (for now)

        Can access hyperparams through self.hparams
        Current priors are flat priors in "acceptable" regions.
        Need to justify priors.

        0 < variance_W0 < 10000.  (km/s)^2
        .05 < beta_W < 1
        -.1 < inv_scalelen_W < .4  (kpc^-1)
        """

        variance_W0, beta_W, inv_scalelen_W = params
        if variance_W0 < 0:
            return -np.inf
        if variance_W0 > 10000.:
            return -np.inf
        if beta_W < 0.05:
            return -np.inf
        if beta_W > 1.0:
            return -np.inf
        if inv_scalelen_W < -.1:
            return -np.inf
        if inv_scalelen_W > .4:
            return -np.inf
        return 0.0

    def lnprob(self, params):
        lnp = self.lnprior(params)
        if not np.isfinite(lnp):
            return -np.Inf
        return lnp + self.lnlh(params)

    def run_emcee(self, threads=1):
        lnprob_args = [self.data, self.hparams]
        lnprob_args = []
        lnprob_func = self.lnprob

        if (len(lnprob_args) == 0):
            sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim,
                                            lnprob_func, args=lnprob_args,
                                            threads=threads)
            sampler.run_mcmc(self.p0, 20000)
            self.sampler = sampler

    def plot_emcee_params(self, nstart=500, **corner_kwargs):
        samples = self.sampler.chain[:, nstart:, :].reshape((-1, 3))
        fig = corner.corner(samples, labels=[r"$S_{W_0}$ [km$^2$/s$^2$]",
                                             r"$\beta$",
                                             r"$R_{W_0}^{-1}$ [kpc]"],
                            quantiles=[.16, .5, .84], show_titles=True,
                            title_kwargs={'fontsize': 12}, **corner_kwargs)
        fig.savefig(self.plotname+".png", format='png')
        fig.savefig(self.plotname+".eps", format='eps')


if __name__ == "__main__":
    pmdata = pm_to_vels.PMmeasurements(RCcatalog=pm_to_vels.catalog,
                                       pmcatalog='UCAC',
                                       maxpmuncer=2.0, maxheight=0.7)
    pmdata.to_space_velocties()
    pmdata.UVW_to_galcen()
    data = pmdata.get_tau_radii_vZg_sigma2Ws_container()
    model = fidmodel(data, plotname='current_tri0')
    model.init_emcee()
    model.run_emcee()
    model.plot_emcee_params()

#    def plot_emcee(self, nexec=1):
#    """
#    - nexec: Number of times to continue looping emcee and making plot
#             Not used yet, hook
#    """
#    try:
#        sampler = self.sampler
#    except AttributeError:
#        run_emcee(self)
#    mk_triangle_plot(sampler)
#    if nexec > 1:
#        for i in xrange(nexec):
#            sampler = run_emcee(ndim, nwalkers, p0, lnprob_func, lnprob_args)
#            mk_triangle_plot(sampler)
# #ndim = 3   #HACK, need to fix!!!
# ###samples = emcee_sampler.chain[:,1000:,:].reshape((-1,ndim))

# result = lnprob(data,params,hyperparams)
# print ('lnprob: {}'.format(result))
# pass
