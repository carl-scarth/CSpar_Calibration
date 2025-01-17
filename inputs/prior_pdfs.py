import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
header_dir = "C:\\Users\\cs2361\\Documents\\Bayesian_Model_Calibration\\source" # directory of headers
sys.path.append(header_dir) # add header directory to path
from utils import set_plot_params

set_plot_params()
infile = "LHSDesign5000x3"
x = pd.read_csv(infile + ".csv")
fig, axes = plt.subplots(1,x.shape[1])
units = ["(GPa)", "(mm)", ""]
for i, column in enumerate(x.columns.values):
    x.hist(column=column, bins=50, ax=axes[i], density=True)
    axes[i].set_xlabel(" ".join((column, units[i])))

axes[0].set_ylabel("PDF")

# Add spring study (degree of clamping) 
K = np.loadtxt('spring_study.csv', delimiter = ',', skiprows=1)
K = np.log(K[:,2])
N = K.shape[0]

uvw = np.loadtxt('spring_study_displacements.csv', delimiter = ',', skiprows = 1)
uvw = uvw[1,:].reshape([N,3])
norm_disp = (uvw[:,2] - np.min(uvw[:,2]))/(np.max(uvw[:,2])-np.min(uvw[:,2]))
ax2 = axes[2].twinx()
ax2.plot(K,norm_disp,'-r')
ax2.set_ylabel("Degree of clamping")
plt.show()