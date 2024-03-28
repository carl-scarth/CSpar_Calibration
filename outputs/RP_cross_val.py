import matplotlib.pyplot as plt
import numpy as np
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

gpfile = "LHSDesign60x6_4"
infile = "LHSDesign60x6_5"
disp_str = ["w"]

apply_force = True # Is force applied, or displacement?

# Load training data
RP_DoE = np.loadtxt(infile + "_RP_displacements.csv", delimiter=",", skiprows=0)
max_DoE = {key : np.loadtxt("_".join((infile,key,"max_displacements.csv")),delimiter=",", skiprows=0) for key in disp_str}

if apply_force:
    force = 200.0 # Maximum applied load
    RP_Force = np.tile(np.linspace(0.0, force, RP_DoE.shape[1]), (RP_DoE.shape[0],1))
    
else:
    RP_Force = np.loadtxt(infile + "_RP_Force.csv", delimiter=",", skiprows=0)
    RP_Force = np.tile(RP_Force, (RP_DoE.shape[0],1))

# at a push could do separate emulator - not using in the calibration?
# but could...

# Need to replace y_DoE with RP force and loop over array below
# DoE = np.loadtxt(infile + "_u_max_displacements.csv", delimiter=",", skiprows=0)
# Load gp mean and standard deviation (used a different name for training gp - fix later)
gp_mean_RP = np.loadtxt(gpfile + "_RP_eta_mu_mu.csv", delimiter=",", skiprows=0)
gp_sd_RP = np.loadtxt(gpfile + "_RP_eta_sigma_mu.csv", delimiter=",", skiprows=0)
gp_mean_max = {key : np.loadtxt(gpfile+"_max_eta_mu_mu_"+key+".csv", delimiter=",", skiprows=0) for key in disp_str}
gp_sd_max = {key : np.loadtxt(gpfile+"_max_eta_sigma_mu_"+key+".csv", delimiter=",", skiprows=0) for key in disp_str}

if apply_force:
    gp_force = np.tile(np.linspace(0.0, force, gp_mean_RP.shape[1]), (gp_mean_RP.shape[0],1))
else:
    gp_force = RP_Force # Cheating to use these as we already know the output - to do this properly we'd also need to predict the reaction force
# Create plot
fig, ax = plt.subplots()
# Cherry pick from output data

# ind = [4]
ind = range(10,20)
# ind = range(5,10)
# ind = [0,1,4]
# ind = [5]0
# ind = [10]
RP_DoE = RP_DoE[ind,:]
max_DoE = {key : max_DoE[key][ind,:] for key in disp_str}

gp_mean_RP = gp_mean_RP[ind,:]
gp_mean_RP_mean = np.mean(gp_mean_RP,axis=0)
gp_sd_RP = gp_sd_RP[ind,:]
gp_sd_RP_mean = np.mean(gp_sd_RP,axis=0)
RP_DoE_mean = np.mean(RP_DoE,axis=0)
gp_mean_max = {key : gp_mean_max[key][ind,:] for key in disp_str}
gp_sd_max = {key : gp_sd_max[key][ind,:] for key in disp_str}

print(RP_DoE)
print(RP_Force)
for i, (displacement, force) in enumerate(zip(RP_DoE,RP_Force)):
    if i == 0:
        ax.plot(-displacement, force, "k", label="Model")
    else:
        ax.plot(-displacement, force, "k")

# Repeat for gp mean
for i, (displacement, sd, force) in enumerate(zip(gp_mean_RP, gp_sd_RP, gp_force)):
    if i == 0:
        ax.plot(-displacement, force, "r",label="Emulator Mean")
        ax.plot(-displacement-2*sd, force, "b", label = "+/- 2 Emulator Standard Deviations")
    else:
        ax.plot(-displacement, force, "r")
        ax.plot(-displacement-2*sd, force, "b")
    ax.plot(-displacement+2*sd, force, "b")

# Plot gp mean +/- two standard deviatins
# for i, forsd in enumerate(gp_sd_RP):
#    if i == 0:
#        ax.plot(-gp_mean_RP[i,:]-2*sd, force, "b", label = "+/- 2 Emulator Standard Deviations")
#    else:
#        ax.plot(-gp_mean_RP[i,:]-2*sd, force, "b")
#    ax.plot(-gp_mean_RP[i,:]+2*sd, force, "b")

ax.set_ylabel("Force (kN)")
ax.set_xlabel("Displacement (mm)")
ax.set_title("Longitudinal dispacement at the RP")
ax.legend()

fig2, axs2 = plt.subplots(1,len(disp_str))
if len(disp_str) == 1:
    axs2 = [axs2]
for j, (ax2, key) in enumerate(zip(axs2, disp_str)):
    for displacement, force in zip(max_DoE[key],RP_Force):
        ax2.plot(-displacement, force, "k")
    for displacement, sd, force in zip(gp_mean_max[key], gp_sd_max[key], gp_force):
        ax2.plot(-displacement, force, "r")
        # Plot gp mean +/- two standard deviations
        ax2.plot(-displacement-2*sd, force, "b")
        ax2.plot(-displacement+2*sd, force, "b")
        
    if j == 0:
        ax2.set_ylabel("Force (kN)")
    ax2.set_xlabel("Displacement (mm)")
    ax2.set_title(key)
fig2.suptitle("Predictions at location of maximum mean displacement across training data")

#fig3, ax3 = plt.subplots()
#ax3.plot(-gp_mean_RP_mean, y_GP, "r")
#ax3.plot(-gp_mean_RP_mean + 2*gp_sd_RP_mean, y_GP, "b")
#ax3.plot(-gp_mean_RP_mean - 2*gp_sd_RP_mean, y_GP, "b")
#ax3.plot(-RP_DoE_mean, y_GP, "k")
#ax3.set_ylabel("Force (kN)")
#ax3.set_xlabel("Displacement (mm)")
# plot average


plt.show()

