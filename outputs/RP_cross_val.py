import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({'figure.figsize' : (8,6),
                'font.size' : 14,
                'font.family' : 'serif',
                'figure.titlesize' : 16,
                'axes.labelsize': 16,
                'xtick.labelsize': 14,
                'ytick.labelsize': 14,
                'legend.fontsize': 14})

infile = "LHSDesign40x4_1"
# Load training data
# DoE = np.loadtxt(infile + "_RP_displacements.csv", delimiter=",", skiprows=0)
DoE = np.loadtxt(infile + "_u_max_displacements.csv", delimiter=",", skiprows=0)
# Load gp mean and standard deviation (used a different name for training gp - fix later)
gp_mean = np.loadtxt("LHSDesign40x4_RP_eta_mu_mu.csv", delimiter=",", skiprows=0)
gp_sd = np.loadtxt("LHSDesign40x4_RP_eta_sigma_mu.csv", delimiter=",", skiprows=0)

# Create plot
fig, ax = plt.subplots()
# Cherry pick from output data
# ind = [28]
# ind = range(0,10)
# ind = [2]
ind = range(31,40)
# ind = [30]
DoE = DoE[ind,:]
gp_mean = gp_mean[ind,:]
gp_sd = gp_sd[ind,:]
  
force = 200.0 # Maximum applied load
# Plot Training data
y_DoE = np.linspace(0.0, force, DoE.shape[1])
for displacement in DoE:
    ax.plot(-displacement, y_DoE, "k")

y_GP = np.linspace(0.0, force, gp_mean.shape[1])
# Repeat for gp mean
for displacement in gp_mean:
    ax.plot(-displacement, y_GP, "r")
# Plot gp mean +/- two standard deviations
for i, sd in enumerate(gp_sd):
    ax.plot(-gp_mean[i,:]-2*sd, y_GP, "b")
    ax.plot(-gp_mean[i,:]+2*sd, y_GP, "b")
        
ax.set_ylabel("Force (kN)")
ax.set_xlabel("Displacement (mm)")
plt.show()