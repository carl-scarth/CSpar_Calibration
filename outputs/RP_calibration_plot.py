import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

# Produces force-displacement plots comparing DIC data with calibrated model

# Thoughts - could be better to load in actual force disp rather than a large time steps
# Bit of a mess. Commented unused code. Delete later if not needed
rcParams.update({'figure.figsize' : (12,9),
                'font.size' : 14,
                'font.family' : 'serif',
                'figure.titlesize' : 16,
                'axes.labelsize': 16,
                'xtick.labelsize': 14,
                'ytick.labelsize': 14,
                'legend.fontsize': 14})

file_str = "LHSDesign40x4"
exp_file = "Crosshead_Force_Disp.csv" # Force-displacement from test machine
# GP_file = "LHSDesign40x4_RP_pred.csv"
# DIC_file = "Interpolated_DIC.csv"
# For if I want to load in actual data
# Folder Containing DIC data
DIC_folder = "..\\..\\Geir_Olafsson\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2\\Data_rad_trimmed"
# Folder Containg GP predictions, interpolated to DIC coordinates
GP_DIC_folder = "interp_gp_output" # folder containng interpolated DIC
# n_incs = 17 # number of increments

exp_data = pd.read_csv(exp_file,sep=";")
exp_force = exp_data["Force [kV]"] # Subtract first entry??
exp_disp = exp_data["Displacement [mm]"]
zeroed_force = exp_force - exp_force[0]
zeroed_disp = exp_disp - exp_disp[0]
max_force = 200

# Load gp mean and standard deviation (used a different name for training gp - fix later)
gp_mean = np.loadtxt(file_str+ "_RP_eta_sam_mu.csv", delimiter=",", skiprows=0)
gp_sd = np.loadtxt(file_str + "_RP_eta_sam_sigma.csv", delimiter=",", skiprows=0)
gp_sam = np.loadtxt(file_str + "_RP_eta_sam.csv", delimiter=",", skiprows=0)

# gp_mean = np.loadtxt("LHSDesign40x4_RP_eta_sam_mu_cal.csv", delimiter=",", skiprows=0)
# gp_sd = np.loadtxt("LHSDesign40x4_RP_eta_sam_sigma_cal.csv", delimiter=",", skiprows=0)
# gp_sam = np.loadtxt("LHSDesign40x4_RP_eta_sam_cal.csv", delimiter=",", skiprows=0)

# Load in (interpolated) DIC data
# DIC_data = pd.read_csv(DIC_file)
# DIC_data.columns.values[0] = "point_ind"

# Alternative source of DIC data - contains entire test dataset rather than interpolated data used to train model
for i, file in enumerate(os.listdir(DIC_folder)):
    data = pd.read_csv(DIC_folder + "\\" + file)
    data["Force"] = zeroed_force[i]
    data["point_ind"] = data.index
    if i == 0:
        DIC_all = data
    else:
        DIC_all = pd.concat((DIC_all, data),axis=0)

# Load in GP_predictions for DIC data points
for file in os.listdir(GP_DIC_folder):
    data = pd.read_csv(GP_DIC_folder + "\\" + file)
    data["Increment"] = int(file.strip("Frame_.csv"))
    #data["Increment"] = int(file.strip("Image_Inc_.csv"))
    try:
        GP_DIC = pd.concat((GP_DIC, data),axis=0)
    except:
        GP_DIC = data
        
n_incs = GP_DIC["Increment"].max()
# Create plot
fig, ax = plt.subplots()
#ax = fig.add_subplot(1, 1, 1)

force = 200.0 # Maximum applied load
# Plot emulator predictions
y_gp = np.linspace(0.0, force, gp_mean.shape[0])
ax.set_title("Calibrated model against machine force-displacement")
ax.plot(-gp_sam, y_gp, "c", linewidth=0.25, label="Sample")
ax.plot(-gp_mean, y_gp, "r", linewidth=1.5, label = "Mean")
# Plot gp mean +/- two standard deviations
ax.plot(-gp_mean-2*gp_sd, y_gp, "b", linewidth=1.5, label="95% Interval")
ax.plot(-gp_mean+2*gp_sd, y_gp, "b", linewidth=1.5)
# Plot experimental force-displacement
ax.plot(-zeroed_disp.to_numpy(),-zeroed_force.to_numpy(),"k", linewidth=1.5)     
label_font = {'family': 'serif', 'size': 16}
ax.set_ylabel("Force (kN)", fontdict = label_font)
ax.set_xlabel("Displacement (mm)", fontdict = label_font)
# ax.legend()

# Also plot force (longitudinal) displacement for selected DIC data points
point_subset = [41306, 45855, 43614, 46631, 48792, 49622]
fig2 = plt.figure()
# coord_list = []
ax2 = []
for i, point in enumerate(point_subset):
    # Extract DIC data (across all increments) at current point
    DIC_point = DIC_all[DIC_all["point_ind"] == point]
    # Also extract values from training data
    gp_point = GP_DIC[GP_DIC["point_ind"] == point].sort_values(["Increment"])
    gp_force = gp_point["Increment"].to_numpy()/n_incs*max_force
    # Extract samples, mean and standard deviation from gp predictions
    gp_point_sam = gp_point[[column for column in GP_DIC.columns.values if "eta_sam_w" in column and column != "eta_sam_mu_w" and column != "eta_sam_sigma_w"]]
    gp_point_mu = gp_point[["eta_sam_mu_w"]]
    gp_point_sigma = gp_point[["eta_sam_sigma_w"]]

    # Extract force and displacemnt data from DIC 
    force_disp = DIC_point[["w_rot", "Force"]].to_numpy()
    #DIC_point = DIC_data[DIC_data["point_ind"] == point]
    #force_disp = DIC_point[["w_rot", "Compressive Force"]].to_numpy()
    #coord_list.append(DIC_point[["x_proj","y_proj","z_proj"]].iloc[0].to_list())
    ax2.append(fig2.add_subplot(2, 3, i+1))
    ax2[-1].plot(-gp_point_sam.to_numpy(), gp_force, "c", linewidth = 0.25, label = "sample")
    ax2[-1].plot(-force_disp[:,0], -force_disp[:,1], "k", linewidth=1.5, label="DIC")
    ax2[-1].plot(-gp_point_mu.to_numpy(), gp_force, "r", linewidth = 1.5, label = "mean")
    ax2[-1].plot(-gp_point_mu.to_numpy()-2*gp_point_sigma.to_numpy(), gp_force, "b", linewidth = 1.0, label = "-2 standard deviations")
    ax2[-1].plot(-gp_point_mu.to_numpy()+2*gp_point_sigma.to_numpy(), gp_force, "b", linewidth = 1.0, label = "+2 standard deviations")
    ax2[-1].set_ylabel("Force (kN)", fontdict = label_font)
    ax2[-1].set_xlabel("Displacement (mm)", fontdict = label_font)

#print(coord_list)
fig2.suptitle("Force (longitudinal) displacement plots from DIC near spar end")
plt.show()

# Works up to here.
# Next plot summary u and v.
# Decide on a few points from each
# Points from middle (for u) are [11710, 10738, 10699, 21171, 10650]
# Try and get some for the other bits (nodes)
