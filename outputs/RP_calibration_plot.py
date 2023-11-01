import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

exp_file = "Crosshead_Force_Disp.csv" # Force-displacement from test machine
GP_file = "LHSDesign40x4_RP_pred.csv"

exp_data = pd.read_csv(exp_file,sep=";")
exp_force = exp_data["Force [kV]"] # Subtract first entry??
exp_disp = exp_data["Displacement [mm]"]


# Load gp mean and standard deviation (used a different name for training gp - fix later)
gp_mean = np.loadtxt("LHSDesign40x4_RP_eta_mu_mu.csv", delimiter=",", skiprows=0)
gp_sd = np.loadtxt("LHSDesign40x4_RP_eta_sigma_mu.csv", delimiter=",", skiprows=0)

# should probably subtract first entry for force
# also look at taking fixed points from DIC data - think about this...
sdfsdfds
# Create plot
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)
# Cherry pick from output data
# ind = [28]
ind = range(1,10)
# ind = [1]
# ind = range(15,30)

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
        
label_font = {'family': 'serif', 'size': 16,}
ax.set_ylabel("Force (kN)", fontdict = label_font)
ax.set_xlabel("Displacement (mm)", fontdict = label_font)

plt.show()
