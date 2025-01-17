import os, sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Add the src directory to the pythonpath for loading shared modules
src_path = "source/"
sys.path.insert(0, src_path)
from sample_prior import *

# Write the Pearson correlation coefficient to the centre of the axes for each plot
def label_input(x, **kwargs):
    ax = plt.gca()
    ax.annotate(x.name, xy = (0.5, 0.5), size = 12, xycoords = ax.transAxes, horizontalalignment = "center", verticalalignment = "center")#, weight="bold")
    # Can't seem to control axis or tick appearance from within function
    # May need to do this when defining the PairGrid
    # Not a disaster by any means as can just delete them in inkscape, but rather wouldn't
    # Quickly google but don't waste time
    # Remember - lower triangular might make the most sense
    #ax.grid(visible=True)
    #ax.tick_params(color = "w")
    # Make axes invisible (somehow???)
    # CANT SEEM TO DO THIS FOR THE DIAGONAL TERMS
    # GOOGLE?
    # COULD JUST DELETE - FEWER PLOTS HERE IN INKSCAPE

def format_axes(x,y,**kwargs):
    ax = plt.gca()
    ax.axes.tick_params(direction = "in")

def label_coorcoeff(x, y, **kwargs):
    # Calculate the value
    coeff = np.corrcoef(x, y)[0][1]
    
    # Make the label
    plot_text = r'$\rho$ = ' + str(round(coeff, 2))
    
    # Add the label to the plot
    ax = plt.gca()

    if kwargs["label"] == "Prior":
        #ax.annotate(plot_text, xy = (0.5, 0.525), size = 10, xycoords = ax.transAxes, color = kwargs["color"], horizontalalignment = "center", verticalalignment = "bottom")
        #sns.kdeplot(x=x,y=y,color=kwargs["color"], fill = True, ax=ax, bw_adjust=1.5, cut = 1, levels = [0.05, 0.25, 0.5, 0.75, 0.95, 1.0]) # higher bw_adjust = higher smoothing
        pass
    else:
        # Could also specify quantiles/percentiles
        #sns.kdeplot(x=x,y=y,color=kwargs["color"], fill = True, ax=ax, bw_adjust=1.5, cut = 1, levels = [0.05, 0.25, 0.5, 0.75, 0.95, 1.0]) # higher bw_adjust = higher smoothing
        # ax.annotate(plot_text, xy = (0.5, 0.475), size = 10, xycoords = ax.transAxes, color = kwargs["color"], horizontalalignment = "center", verticalalignment = "top")
        ax.annotate(plot_text, xy = (0.5, 0.05), size = 12, xycoords = ax.transAxes, color = "k", horizontalalignment = "center", verticalalignment = "bottom")

    #ax.set_axis_off()
    #ax.tick_params(color = "w")

# wd = "E:\\Calibration_outputs_for_paper"

plt.rcParams.update({'font.serif'        : 'Times New Roman',
                    'font.family' : 'serif',
                    'font.size': 12,
                    'figure.titlesize' : 12,
                    'axes.labelsize': 12,
                    'xtick.labelsize': 12,
                    'ytick.labelsize': 12,
                    'legend.fontsize': 12})    

#plt.rcParams["axes.labelsize"] = 16

wd = os.getcwd()
in_str = "LHSDesign100x8"

postsam_file = os.path.join(wd, in_str + "_posterior_samples.csv")
priorparam_file = os.path.join(wd, in_str + "_tf_param_mod.csv")

# Load in posterior samples
post_samples = pd.read_csv(postsam_file)
N = post_samples.shape[0]
inp_labels = post_samples.columns.values

# Load in prior distribution parameters
prior_param = pd.read_csv(priorparam_file, index_col = 0)
prior_samples = sample_prior(prior_param, N, inp_cov = False)
prior_samples = pd.DataFrame(prior_samples, columns=inp_labels)

# Plot histograms and print summary statistics for the two datasets
print(prior_samples.describe())
print(post_samples.describe())

# Combine the two datasets into a single dataframe with categorical 
# variable indicating which set each point is from
prior_samples["Category"] = "Prior"
post_samples["Category"] = "Posterior"
samples = pd.concat((prior_samples,post_samples),axis = 0)

# sns.pairplot(post_samples, vars = inp_labels, kind = "reg")#, plot_kws=dict(s = 50)) # reg adds a regression line
sns.pairplot(post_samples, vars = inp_labels, kind = "hist")#, plot_kws=dict(s = 50))
#g = sns.pairplot(samples, vars = inp_labels, kind = "hist", diag_kind = None, hue = "Category",
#             palette = {"Prior" : sns.color_palette("bright")[0], "Posterior" : sns.color_palette("bright")[3]}, corner=True, height = 1.25, aspect = 1)#plot_kws=dict(s = 50)
#print(dir(g))

plt.show()
adsad
# Directly implement own formatting
grid = sns.PairGrid(data = samples, vars = inp_labels, hue = "Category", height = 1.25, aspect = 1.0, layout_pad=0.1, despine=True,
                    palette = {"Prior" : sns.color_palette("bright")[0], "Posterior" : sns.color_palette("bright")[3]}, diag_sharey=False)
#grid = grid.map_lower(sns.histplot, bins=50)#, shrink = 1.5)#, pthresh = 0.0025) # Bins = 50 might be a slightly smaller filesize... Likewise pthresh (minimum probability threshold)
grid = grid.map_upper(sns.histplot, bins=50)#, shrink= 1.5)#, pthresh = 0.0025)
grid = grid.map_lower(sns.histplot, bins=50)
#grid = grid.map_lower(label_coorcoeff)
#grid = grid.map_diag(label_input)
#grid.add_legend()

limits = [[120, 190], [6.9, 11.1], [3.75, 5.85], [0.105, 0.15], [-2.35, 2.35], [5.9, 12.5], [4.0, 12.0], [0.0, 1.75]]
ticks = [[125, 150, 175], [7, 8, 9, 10], [4, 4.75, 5.5], [0.11, 0.125, 0.14], [-2, 0, 2], [6, 8, 10, 12], [5, 7.5, 10], [0, 0.5, 1, 1.5]]
ticklabel = [[str(label) for label in plot_ticks] for plot_ticks in ticks] # This produces the simplest string for each tick label

# How to hide the diagonal axes - use the same format 
# to do other things. E.g. rotate text on last column
# for rad_thin to hug closer to lower diagonal
# So you can just index the axes as if an array
# So the below code is a little silly - 
# just use single variable to loop over the variables
# and use both indices ie grid.axes[i,i]
for i, y_var in enumerate(grid.y_vars):
    for j, x_var in enumerate(grid.x_vars):
        if j == i:
            grid.axes[i, j].set_visible(False)
            #grid.axes[i, j].annotate(samples.columns.values[i], xy = (0.5, 0.5), size = 12, xycoords = grid.axes[i,j].transAxes, horizontalalignment = "center", verticalalignment = "center")#, weight="bold")
        else:
            #grid.axes[i, j].grid(visible = True)
            grid.axes[i,j].set_xlim(limits[j][0],limits[j][1])
            grid.axes[i,j].set_ylim(limits[i][0],limits[i][1]) 
            grid.axes[i,j].set_xticks(ticks[j], labels = ticklabel[j])
            grid.axes[i,j].set_yticks(ticks[i], labels = ticklabel[i])
            grid.axes[i,j].axes.tick_params(direction = "out", length = 5)#, labelsize=12)
            
            # If over-tracing
            # Reduce the spine width
            plt.setp(grid.axes[i,j].spines.values(), linewidth=1.0)
            # Also the ticks
            grid.axes[i,j].xaxis.set_tick_params(width=1.0)
            grid.axes[i,j].yaxis.set_tick_params(width=1.0)

for i in range(len(grid.y_vars)):
    # Remove axis labels
    grid.axes[i,0].axes.tick_params(labelleft = False)
    grid.axes[-1,i].axes.tick_params(labelbottom = False)
    grid.axes[i,0].yaxis.set_label_text(samples.columns.values[0],visible=False)
    grid.axes[-1,i].xaxis.set_label_text(samples.columns.values[0],visible=False)
# Add missing labels due to removing the diagonals
grid.axes[0,1].axes.tick_params(labelleft = False)
grid.axes[-2,-1].axes.tick_params(labelbottom = False)
grid.axes[0,1].yaxis.set_label_text(samples.columns.values[0],visible=False)
grid.axes[-2,-1].xaxis.set_label_text(samples.columns.values[-2],visible=False)
#plt.savefig("prior_posterior_pairs.png", transparent = True, dpi = 1200, pad_inches = 0)
plt.show()

# Otherwise have KDE on lower triangular? with correlation coefficient overlaid
# remember there's a shade option with kdeplot
# lots of other options to look at for individual seaborn functions as well as the overall pairgrid settings
# KDE STILL HAS ADVANTAGE OF REDUCING FILE SIZE...
# As would maybe some of the settings of histplot

#plt.cm.coolwarm(np.linspace(0,1,model_disp.shape[1])) # get colourmap
#sns.pairplot(samples, vars = inp_labels, diag_kind=None, plot_kws=dict(s = 50), hue="Category", 
#             palette = {"Prior" : sns.color_palette("bright")[0], "Posterior" : sns.color_palette("bright")[3]})
