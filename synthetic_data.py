from __future__ import print_function
import numpy as np
import random


class fidsynth(object):
    """
    fidsynth

    Class to create synethic data according to fiducal model. Creates a data
    container that will replace the pm_to_velocities data container when
    checking synthetic data.

    Model parameters are: SW0, Beta, R^-1

    SW0 = Normalization of variance of W vels
    Beta = Power law index of growth of variance as func. of age
    R^-1 = Inverse scale length of variance; assume variance is an exp.
           function of age.

    Model
    ---------

    SW(R, tau) = SW0 * (tau/tau0)^(2*Beta) * e^(-2. * (R - R0) * R^-1)

    R:      Galactocentric Radius
    tau:    Stellar age
    tau0:   Min. age or normalization of age (tau0=1 Gyr)
    R0:     Solar radius (R0=8 kpc)
    """

    def __init__(self, ages, radii, beta=.3, inv_Rd=.1, SW0=10.0):
        """
        Initialize class

        Parameters
        ----------------

        ages : array [Gyr]
            Array of input ages, can be taken from data or created by user
        radii : array [kpc]
            Array of input radii, can be taken from data or created by user
        beta : float
            "Truth" of beta for synthetic data
        inv_Rd : tuple [kpc]
            Inverse scale length of velocity dispersion

        Returns
        ----------------
        fidsynth instance with attributes Ws, sigma2Ws, ages, radii.
        Can be fed directly into model class.

        """
        self.ages = ages
        self.radii = radii
        self.beta = beta
        self.inv_Rd = inv_Rd
        self.N = len(ages)
        self.SW0 = SW0
        # assume you just want to make the data
        print('Generating synthetic data using:')
        print('SW0: {0}, Beta: {1}, Rd^-1: {2}'.format(self.SW0, self.beta,
                                                       self.inv_Rd))

    def trueWs(self, age0=1.):
        """
        For each data point, samples a  gaussian whose shape is determined
        by the model and the 'truth'. Generates W velocities.

        age0 : float
            normalizing age [Gyr]
        """
        variances_W = (self.SW0 *
                       np.power((self.ages / age0), 2. * self.beta) *
                       np.exp(-2. * (self.radii - 8.0) * self.inv_Rd))
        std_W = np.sqrt(variances_W)
        Ws = np.array([random.normalvariate(mu=0.0, sigma=stdev) for stdev
                       in std_W])
        self.trueWs = Ws

    def sigma2Ws(self, sig2Wmax=50.0):
        """
        Given Ws, generate sigma_W**2 using the relationship between W
        and sigma2W found in data.

        Model:
        dW^2 / W^2 = N / |W| * RV_beta(.65,6.0)

        N : normalization of upper envelope of error. This is sig2Wmax.
        W : W velocity
        RV_beta(.65, 6.0): random variate from beta distribution with shape
        parameters alpha = 0.6, beta = 6.0.

        Parameters
        ----------------

        sig2Wmax : float (50.0)
            Normalization of the upper envelope of fractional error in W^2.

        Notes
        ----------------
        For now, the shape parameters of the beta distribution are decided by
        eye. Plots show this is a very good description.

        Returns
        ----------------
        None, but self.sigma2Ws is created.

        """
        self.sig2Wmax = sig2Wmax

        envelope = sig2Wmax / np.abs(self.trueWs)
        betaRVs = np.random.beta(.65, 6.0, len(envelope))
        frac_err = envelope * betaRVs
        self.sigma2Ws = frac_err * self.trueWs**2.

    def addnoise(self):
        """
        Given sigma2Ws, perturb the trueWs to yield measured Wvels.
        Assume gaussian noise.
        """
        sigmaWs = np.sqrt(self.sigma2Ws)  # sigW*sigW = sigma2W
        newWs = np.array([random.normalvariate(mu=W, sigma=sigW) for
                          W, sigW in zip(self.trueWs, sigmaWs)])
        self.Ws = newWs
