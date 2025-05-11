import numpy as np
from scipy.stats import halfnorm
import warnings

class Prior:
# Define class of prior distribution objects containing type of distribution all important statistical attributes
    def __init__(self,name, distribution, param_1, param_2, weight_LHS = True, inp_cov = True):
        # Initialise attributes of prior based upon the type of distribution
        self.name = name
        self.distribution = distribution
        self.weight_LHS = weight_LHS
        self.inp_cov = inp_cov
        if self.distribution == 'Gaussian':
            # If Gaussian, the prior is described by a mean and Coefficient of Variation
            self.mu = param_1
            self.COV = param_2
            # Cannot define via Coefficient of Variation if the mean is zero
            if self.mu == 0:
                self.inp_cov = False
        elif self.distribution == 'Lognormal':
            self.mu = param_1
            self.s = param_2
        elif self.distribution == 'Halfnormal':
            self.mu = param_1
            self.scale = param_2
        elif self.distribution == 'Uniform' or self.distribution == 'Loguniform':
            # If Uniform or loguniform, the prior is described by a lower and upper bound
            self.lb = param_1
            self.ub = param_2
        elif self.distribution == 'Loggamma':
            # Loggamma described by shape and rate
            self.shape = param_1
            self.rate = param_2

    @property
    def sigma(self):
        if self.distribution == 'Gaussian':
            if self.inp_cov:
                return self.mu*self.COV/100.0
            else:
                return self.COV
        elif self.distribution == 'Halfnormal':
            # Not the standard deviation, but a parameter I refer to as "sigma"
            return self.scale
        else:
            return None
            # Shouldn't be used, but if necessary could actually specify this
        
    @property
    def min(self):
        if self.distribution == 'Gaussian':
            return self.mu - 3.5*self.sigma
        elif self.distribution == 'Lognormal':
            # Defined by +/- 3 standard deviations of log(x)
            return np.exp(self.mu - 3.5*self.s)
        elif self.distribution == 'Halfnormal':
            return self.mu
        elif self.distribution == 'Uniform' or self.distribution == 'Loguniform':
            return self.lb
        else:
            warnings.warn("Min and max have not been implemented for this distribution yet")
            return None

    @property
    def max(self):
        if self.distribution == 'Gaussian':
            return self.mu + 3.5*self.sigma
        elif self.distribution == 'Lognormal':
            return np.exp(self.mu + 3.5*self.s)
        elif self.distribution == 'Halfnormal':
            return halfnorm.ppf(0.999)*self.sigma + self.mu
        elif self.distribution == 'Uniform' or self.distribution == 'Loguniform':
            return self.ub
        else:
            warnings.warn("Min and max have not been implemented for this distribution yet")
            return None