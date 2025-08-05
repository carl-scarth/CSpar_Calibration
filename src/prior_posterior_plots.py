import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import halfnorm, gamma, invgamma
from math import floor, ceil, log2

def ff_ls_hist(ls, p_eta = None, inp_labels=None, prior_shape = 4.0, prior_rate = 4.0):
    # Produce plots of posterior and prior distributions of correlation lengths
    # (ls) for full-field emulator. Arrange with a row for each principal 
    # and a column for each input.
    # inp_labels is an optional input with list of input labels for plot titles placed at head of column
    # Option to pass shape and rate parameters of the inverse gamma prior 
    # with default values as recommended in:
    #https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#priors-for-gaussian-processes
    
    # Get number of inputs and principal components
    ls_names = ls.columns.values
    if p_eta is None:
        p_eta = int(max([name.split("_")[1] for name in ls_names])) + 1
        q = int(max([name.split("_")[2] for name in ls_names])) + 1
    else:
        q = ls.shape[1]//p_eta # Use if passing p_eta as an input

    ls_plot = np.linspace(0.0,4.0,100) # Range of values for prior plot
    # Set-up subplots.
    fig, axes = plt.subplots(p_eta,q)# figsize=[10,6])
    fig.tight_layout(pad = 1.01)
    # Loop over each row of the plot
    for i, ax_row in enumerate(axes):
        # Loop over the columns
        for j, ax in enumerate(ax_row):
            colname = "_".join(("ls",str(i),str(j)))
            # Plot histogram of posterior distribution
            ax.hist(ls[colname], bins = 25, density = True, color = "tomato", edgecolor="black")
            # Overlay plot of prior distribution
            # Note: I think the scipy documentation is incorrect about the inverse gamma pdf
            # I've compared with my own implementation with the pymc docs and it matches, 
            # so I'm confident the below correctly matches the prior
            ax.plot(ls_plot, invgamma.pdf(ls_plot, prior_shape, loc = 0.0, scale = prior_rate), lw = 2, color = "blue")
            
            ax.set_xlabel(colname)
            ax.set_ylabel("PDF")
            if inp_labels is not None and i == 0:
                # Only include title on first row
                title = inp_labels[j]
                ax.set_title(title)
  

def ff_sigma_hist(sigma, p = None, prior_scale = 1.0):
    # Produce plots of posterior and prior distributions of GP scale hyperparameters 
    # (sigma) for full-field output. Each plot corresponds to a different basis function
    # Option to pass scale parameter of halfnormal prior
    
    # Get p if not passed as input
    if p is None:
        p = sigma.shape[1]

    # Split over multiple rows if more than 5 basis functions
    if p <= 5:
        fig, axes = plt.subplots(1,p)
        axes = [axes]
    else:
        fig, axes = plt.subplots(p//5,5)
    
    fig.tight_layout()
    sigma_plot = np.linspace(0.0,10.0,200) # Range of values for prior plot
    for i, sigma_i in enumerate(sigma):
        # plot posterior histogram
        axes[i//5][i%5].hist(sigma[sigma_i], bins = 25, density = True, color="tomato", edgecolor="black")
        # Overlay plot of prior
        axes[i//5][i%5].plot(sigma_plot, halfnorm.pdf(sigma_plot, scale = prior_scale), lw = 2, color = "blue")
        # Just use column label
        axes[i//5][i%5].set_xlabel(sigma_i)
        axes[i//5][i%5].set_ylabel("PDF")
        
    fig.suptitle("Emulator standard deviation")
    
def lambda_hist(lambda_sam, prior_shape = 5.0, prior_rate = 5.0, label = "lambda", adj_prior_shape = None, adj_prior_rate = None):
    # Plots precision hyper-parameter lambda for scalar-valued quantities
    # Optional inputs prior_shape and prior_rate give the properties of the gamma
    # prior. Default values are from Higdon et al. 
    # If adj_prior_shape and adj_prior_rate are passed a second prior plot is 
    # overlaid with these parameter values.
    
    n_lambda = lambda_sam.shape[1]
    fig, axes = plt.subplots(1,n_lambda)
    if n_lambda == 1:
        axes = [axes]

    # Get x-limits for plots
    xlim_prior = [gamma.ppf(0.01, prior_shape, scale = 1.0/prior_rate), 
                  gamma.ppf(0.99, prior_shape, scale = 1.0/prior_rate)]
    
    max_x = xlim_prior[1]
    if adj_prior_shape is not None:
        xlim_adj_prior = [gamma.ppf(0.01, adj_prior_shape, scale = 1.0/adj_prior_rate),
                          gamma.ppf(0.99, adj_prior_shape, scale = 1.0/adj_prior_rate)]
        max_x = max([max_x, xlim_adj_prior[1]])
          
    lambda_plot = np.linspace(xlim_prior[0], xlim_prior[1], 1000) # set of x values for prior plot
    for i, lambda_i in enumerate(lambda_sam):
        if n_lambda > 1:
            header = "_".join(label, str(i))
        else:
            header = label

        # Plot posterior histogram
        axes[i].hist(lambda_sam[lambda_i], bins = 25, density = True, color="tomato", edgecolor="black")
        # Overlay with prior distribution
        axes[i].plot(lambda_plot, gamma.pdf(lambda_plot, prior_shape, scale = 1.0/prior_rate), lw = 2, color = "blue")

        axes[i].set_title(header)
        axes[i].set_xlabel(label)
        axes[i].set_ylabel("PDF")
        axes[i].set_xlim([0, floor(2**log2(max_x))])
        if adj_prior_shape is not None:
            # Produce additional prior plot if needed
            adj_lambda_plot = np.linspace(xlim_adj_prior[0], xlim_adj_prior[1], 1000)
            axes[i].plot(adj_lambda_plot, gamma.pdf(adj_lambda_plot, adj_prior_shape, scale = 1.0/adj_prior_rate))
    