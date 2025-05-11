import numpy as np
import pandas as pd
from scipy.stats import norm, uniform, loguniform, lognorm, halfnorm
from Prior import Prior

def sample_prior(input_list, N, inp_cov = True):
    
# Generate a set of N random samples for a set of prior parameters, distributed
# with parameters stored in input list. Each entry of this list contains a list
# of parameters for each random variable, with entries:
# 0 : Name
# 1 : Distribution
# 2 : Parameter 1
# 3 : Parameter 2
# 4 : Boolean indicating whether or not to weight according to the RV PDF (optional)
# Where parameters 1 and 2 depend upon the distribution
# Alternatively, input_list can be a pandas dataframe wherein the index is the 
# variable name, and the following columns are as above

    # Create a Prior object for each item in input_list
    if isinstance(input_list, pd.DataFrame):
        d = input_list.shape[0]
        priors = []
        for index, item in input_list.iterrows():
            priors.append(Prior(index, item["distribution"], item["param_1"], item["param_2"], inp_cov = inp_cov))

    elif isinstance(input_list, list):
        d = len(input_list) # Number of inputs
        priors = []
        for item in input_list:
            priors.append(Prior(item[0],item[1],item[2],item[3], inp_cov = inp_cov))
    else:
        raise Exception("Please input either a Pandas DataFrame or list")
    
    x = np.empty([N, d])
    for i, inp in enumerate(priors):
        # Would be better if this were a method of the Priors object...
        # It might be possible to do this via the rv_continuous class of scipy
        # without the conditional
        if inp.distribution == 'Gaussian':
            x[:,i] = norm(loc = inp.mu, scale = inp.sigma).rvs(size = N)
        elif inp.distribution == 'Lognormal':
            # scipy respresents lognormal distribution in an unusual way
            x[:,i] = lognorm(inp.s, loc = 0.0, scale = np.exp(inp.mu)).rvs(size = N)
        elif inp.distribution == 'Halfnormal':
            x[:,i] = halfnorm().rvs(size = N)*inp.sigma + inp.mu
        elif inp.distribution == 'Uniform':
            # Note, the scipy uniform distribution is bounded by [loc, loc + scale]
            x[:,i] = uniform(loc = inp.lb, scale = (inp.ub - inp.lb)).rvs(size = N)
        elif inp.distribution == 'Loguniform':
            x[:,i] = loguniform.rvs(inp.lb, inp.ub, size = N) 
            
    return(x)