import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Used to produce force displacement plots for validating emulator

rcParams.update({'figure.figsize' : (12,9),
                'font.size' : 14,
                'font.family' : 'serif',
                'figure.titlesize' : 16,
                'axes.labelsize': 16,
                'xtick.labelsize': 14,
                'ytick.labelsize': 14,
                'legend.fontsize': 14})

gpfile = "LHSDesign50x3_2"
infile = "LHSDesign50x3_3"
disp_str = disp_str = ["u", "v", "w"]
# Load training data
RP_DoE = np.loadtxt(infile + "_RP_displacements.csv", delimiter=",", skiprows=0)
max_DoE = {key : np.loadtxt("_".join((infile,key,"max_displacements.csv")),delimiter=",", skiprows=0) for key in disp_str}

# DoE = np.loadtxt(infile + "_u_max_displacements.csv", delimiter=",", skiprows=0)
# Load gp mean and standard deviation (used a different name for training gp - fix later)
gp_mean_RP = np.loadtxt(gpfile + "_RP_eta_mu_mu.csv", delimiter=",", skiprows=0)
gp_sd_RP = np.loadtxt(gpfile + "_RP_eta_sigma_mu.csv", delimiter=",", skiprows=0)
gp_mean_max = {key : np.loadtxt(gpfile+"_max_eta_mu_mu_"+key+".csv", delimiter=",", skiprows=0) for key in disp_str}
gp_sd_max = {key : np.loadtxt(gpfile+"_max_eta_sigma_mu_"+key+".csv", delimiter=",", skiprows=0) for key in disp_str}
# Create plot
fig, ax = plt.subplots()
# Cherry pick from output data

# ind = [2]
ind = range(20,25)
# ind = [5]
# ind = [10]
RP_DoE = RP_DoE[ind,:]
max_DoE = {key : max_DoE[key][ind,:] for key in disp_str}

gp_mean_RP = gp_mean_RP[ind,:]
gp_sd_RP = gp_sd_RP[ind,:]
gp_mean_max = {key : gp_mean_max[key][ind,:] for key in disp_str}
gp_sd_max = {key : gp_sd_max[key][ind,:] for key in disp_str}

force = 200.0 # Maximum applied load
# Plot Training data
y_DoE = np.linspace(0.0, force, RP_DoE.shape[1])
for displacement in RP_DoE:
    ax.plot(-displacement, y_DoE, "k")

y_GP = np.linspace(0.0, force, gp_mean_RP.shape[1])
# Repeat for gp mean
for displacement in gp_mean_RP:
    ax.plot(-displacement, y_GP, "r")
# Plot gp mean +/- two standard deviations
for i, sd in enumerate(gp_sd_RP):
    ax.plot(-gp_mean_RP[i,:]-2*sd, y_GP, "b")
    ax.plot(-gp_mean_RP[i,:]+2*sd, y_GP, "b")
        
ax.set_ylabel("Force (kN)")
ax.set_xlabel("Displacement (mm)")
ax.set_title("Longitudinal dispacement at the RP")

fig2, axs2 = plt.subplots(1,len(disp_str))
if len(disp_str) == 1:
    axs2 = [axs2]
for ax2, key in zip(axs2, disp_str):
    for displacement in max_DoE[key]:
        ax2.plot(-displacement, y_DoE, "k")
    for displacement in gp_mean_max[key]:
        ax2.plot(-displacement, y_GP, "r")
        # Plot gp mean +/- two standard deviations
    for i, sd in enumerate(gp_sd_max[key]):
        ax2.plot(-gp_mean_max[key][i,:]-2*sd, y_GP, "b")
        ax2.plot(-gp_mean_max[key][i,:]+2*sd, y_GP, "b")
        
    ax2.set_ylabel("Force (kN)")
    ax2.set_xlabel("Displacement (mm)")
    ax2.set_title(key)
fig2.suptitle("Predictions at location of maximum mean displacement across training data")
plt.show()