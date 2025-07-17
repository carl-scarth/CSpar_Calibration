# Collection of useful functions used across all examples
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from matplotlib import rcParams

def estimate_mode(x, n_points = 10000):
    # Estimate mode of a numpy array of random samples, belonging to some continuous
    # distribution, using kernel density estimate
    x_vals = np.linspace(np.min(x), np.max(x), num=n_points)
    y_vals = gaussian_kde(x).evaluate(x_vals)
    return x_vals[np.argmax(y_vals)]

def estimate_mode_df(df, **kwargs):
    # Estimates mode for a Pandas DataFrame, the columns of which contain samples from some 
    # underlying continuous random variable
    modes = pd.Series([estimate_mode(df[col], **kwargs) for col in df], index = [col for col in df])
    return modes

def set_plot_params():
    # Set commonly used plot parameters to desired values
    # Plotting parameters
    rcParams.update({'figure.figsize' : (8,6),
                    'font.size': 16,
                    'figure.titlesize' : 18,
                    'axes.labelsize': 18,
                    'xtick.labelsize': 15,
                    'ytick.labelsize': 15,
                    'legend.fontsize': 12})    
    
def set_plot_params_elsevier():
    # Set commonly used plot parameters to desired values
    # Plotting parameters
    rcParams.update({
                    'figure.figsize'   : (3.54331,2.6575), # 90 mm width in inches, 4:3 aspect ratio
                    'font.serif'        : 'Times New Roman',
                    'font.family'       : 'serif',
                    'font.size': 10,
                    'figure.titlesize' : 12,
                    'axes.labelsize': 11,
                    'xtick.labelsize': 10,
                    'ytick.labelsize': 10,
                    'legend.fontsize': 10})    