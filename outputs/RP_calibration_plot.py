import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

# Produces force-displacement plots comparing DIC data with calibrated model
# Thoughts - could be better to load in actual force disp rather than a large time steps
# Bit of a mess. Commented unused code. Delete later if not needed

# First attempt at coding a general purpose plotting function - place in separate header if useful
def gp_plot(fig = [], ax = [], n_row = 1, n_col = 1, gp_force = [], gp_mean = [], gp_sd = [], gp_sam = [], xval_force = [], xval_y = [], add_legend = False):
    # Individual plot of Gaussian process mean and standard deviation, and samples
    # ax = list of axes on the figure
    if not fig:
        fig = plt.figure()
    ax.append(fig.add_subplot(n_row, n_col, len(ax)+1))
    if len(gp_sam) > 0:
        ax[-1].plot(gp_sam, get_force(gp_sam, force=gp_force), "c", linewidth = 0.25, label = "sample")
    if len(gp_mean) > 0:
        ax[-1].plot(gp_mean, get_force(gp_mean, force=gp_force), "r", linewidth = 1.5, label = "mean")
    if len(gp_sd) > 0 and len(gp_mean) > 0:
        ax[-1].plot(gp_mean-2*gp_sd, get_force(gp_mean, force=gp_force), "b", linewidth = 1.0, label = "-2 standard deviations")
        ax[-1].plot(gp_mean+2*gp_sd, get_force(gp_mean, force=gp_force), "b", linewidth = 1.0, label = "+2 standard deviations")
    print(xval_y)
    if len(xval_y) > 0:
        ax[-1].plot(xval_y, get_force(xval_y, force = xval_force), "k", linewidth=1.5, label="DIC")
    ax[-1].set_ylabel("Force (kN)")
    ax[-1].set_xlabel("Displacement (mm)")
    if add_legend:
        ax.legend()
    return(ax)

def get_force(y, force = []):
    # Check if force is given, if not return range length of y QoI
    if len(force) < 1:
        force = np.arange(y.shape[0])
    return(force)

def subplots_loop(point_inds, DIC, GP_DIC, force_inc, disp_str, n_row, n_col, minus = False):
    # Wrapper function
    # coord_list = []
    fig = plt.figure()
    ax = []
    # for i, point in enumerate(point_inds):
    for point in point_inds:
        # Extract DIC data (across all increments) at current point
        DIC_point = DIC[DIC["point_ind"] == point]
        # Also extract values from training data
        gp_point = GP_DIC[GP_DIC["point_ind"] == point].sort_values(["Increment"])
        gp_force = gp_point["Increment"].to_numpy()*force_inc
        # Extract samples, mean and standard deviation from gp predictions
        gp_point_sam = gp_point[[column for column in GP_DIC.columns.values if "_".join(("eta_sam",disp_str)) in column and column != "_".join(("eta_sam_mu",disp_str)) and column != "_".join(("eta_sam_sigma",disp_str))]]
        gp_point_mu = gp_point[["_".join(("eta_sam_mu", disp_str))]]
        gp_point_sigma = gp_point[["_".join(("eta_sam_sigma", disp_str))]]
        if minus:
            gp_point_sam = -gp_point_sam
            gp_point_mu = -gp_point_mu
        # Extract force and displacemnt data from DIC 
        force_disp = DIC_point[["Force", "_".join((disp_str,"rot"))]].to_numpy()
        force_disp[:,0] = - force_disp[:,0]
        if minus:
            force_disp[:,1] = -force_disp[:,1]


        #coord_list.append(DIC_point[["x_proj","y_proj","z_proj"]].iloc[0].to_list())
        # print(coord_list)
        gp_plot(fig = fig, ax = ax, n_row = n_row, n_col = n_col, gp_force = gp_force, gp_mean = gp_point_mu.to_numpy(), gp_sd = gp_point_sigma.to_numpy(), gp_sam = gp_point_sam.to_numpy(), xval_force = force_disp[:,0], xval_y = force_disp[:,1], add_legend = False)


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
    print(file)
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
ax.set_ylabel("Force (kN)")
ax.set_xlabel("Displacement (mm)")
# ax.legend()

# Also plot force (longitudinal) displacement for selected DIC data points
point_subset_w = [41306, 45855, 43614, 46631, 48792, 49622]
subplots_loop(point_subset_w, DIC_all, GP_DIC, max_force/n_incs, "w", 2, 3, minus = True)
# Plot minmum vertical displacement
point_subset_umin = [11710, 10738, 10699, 21171, 10650]
subplots_loop(point_subset_umin, DIC_all, GP_DIC, max_force/n_incs, "u", 1, 5)
# Plot maximum vertical displacement at both ends
point_subset_umax = [669, 554, 398, 22793, 28281, 38707, 388, 38834, 45144, 45568]
subplots_loop(point_subset_umax, DIC_all, GP_DIC, max_force/n_incs, "u", 2, 5)
# Plot maximum transverse diplacement at both end of the flange tips
point_subset_vmax = [51980, 51513, 50460, 51045, 51126, 51318]
subplots_loop(point_subset_vmax, DIC_all, GP_DIC, max_force/n_incs, "v", 2, 3)
# Decide on a few points from each
# Points from middle (for u) are [11710, 10738, 10699, 21171, 10650]
# Try and get some for the other bits (nodes)
plt.show()