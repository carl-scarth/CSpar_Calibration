# Monte Carlo of Priors

import sys
#import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

header_dir = "C:\\Users\\cs2361\\Documents\\Bayesian_Model_Calibration\\source" # directory of sampling headers
sys.path.append(header_dir) # add header directory to path
from sample_prior import *  # Import Latin Hypercube module
from utils import set_plot_params

def combine_priors(name, params, N, plot_pdfs = True, units = "GPa"):
    inputs = [[name, 'Gaussian', item[0], item[1], True] for item in params]
    samples = sample_prior(inputs, N)
    x_all = pd.DataFrame(samples.flatten(),columns=[name])
    if plot_pdfs:
        x = pd.DataFrame(samples, columns = ["_".join((name,str(i))) for i in range(len(inputs))])
        #fig, axes = plt.subplots(1,x.shape[1]+1)
        fig, ax = plt.subplots()
        #for ax, column in zip(axes, x.columns.values):
        #    x.hist(column=column, bins=25, ax=ax, density=True)
        #    ax.set_xlabel(" ".join((column, units)))
        #x_all.hist(bins=25, ax=axes[-1], density = True)
        x_all.hist(bins=25, ax=ax, density = True)
    desc_df = x_all.describe()
    desc_df = pd.concat((desc_df, pd.DataFrame((x_all.std()/x_all.mean()*100.0).to_numpy(),columns=[name], index=["CoV%"])), axis=0)
    return(desc_df)

if __name__ == "__main__":
    properties = {
        "E1C"   :   [[150.0, 2.525],[140.9, 1.55],[154.489, 3.5]],
        "E1T"   :   [[164.0, 1.83],	[162.1, 2.27],[171.42, 1.39]],
        "E2T"   :   [[8.96, 3.37], [12, 2.2], [9.08, 1.03]],
        "G12"   :   [[4.69, 3.27], [5.29, 2.53]],
        "XT"    :   [[2.538, 4.72], [2.558, 4.1], [2.724, 4.40], [2.625, 4.56], [2.3262, 5.8]],
        "XC"    :   [[1.69, 4.75],	[1.731, 3.22],	[1.194, 5.2], [1.2001,12.1], [1.017, 5.2]],
        "YT"    :   [[0.064, 9.47], [0.0602, 9.432], [0.0623, 8.5]],
        "YC"    :   [[0.2857, 4.5],	[0.4006, 4.148],	[0.1998, 10.2]],
        "S"     :   [[0.12, 0.876], [0.0923, 0.7], [0.091, 1.6], [0.1127, 0.933]]
        }
    
    print(properties)
    N = 10000 # Number of samples per distribution
    for key, value in properties.items():
        desc_df = combine_priors(key, value, N)
        print(desc_df)
    plt.show()