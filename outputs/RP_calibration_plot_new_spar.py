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
        ax[-1].plot(xval_y, get_force(xval_y, force = xval_force), "kx", linewidth=1.5, label="DIC")
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
        DIC_point = DIC[DIC["point_ind"] == point].sort_values(["Increment"])
        # Also extract values from training data
        gp_point = GP_DIC[GP_DIC["point_ind"] == point].sort_values(["Increment"])
        #print(force_inc)
        #gp_force = gp_point["Increment"].to_numpy()*force_inc
        gp_force = gp_point["Compressive Force"].to_numpy()
        #print(gp_force)
        # Extract samples, mean and standard deviation from gp predictions
        gp_point_sam = gp_point[[column for column in GP_DIC.columns.values if "_".join(("eta_sam",disp_str)) in column and column != "_".join(("eta_sam_mu",disp_str)) and column != "_".join(("eta_sam_sigma",disp_str))]]
        gp_point_mu = gp_point[["_".join(("eta_sam_mu", disp_str))]]
        gp_point_sigma = gp_point[["_".join(("eta_sam_sigma", disp_str))]]
        if minus:
            gp_point_sam = -gp_point_sam
            gp_point_mu = -gp_point_mu
        # Extract force and displacemnt data from DIC 
        force_disp = DIC_point[["Compressive Force", "_".join((disp_str,"rot"))]].to_numpy()
        # force_disp[:,0] = - force_disp[:,0]
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

file_str = "LHSDesign100x8"

# GP_file = "LHSDesign40x4_RP_pred.csv"
# DIC_file = "Interpolated_DIC.csv"
# For if I want to load in actual data
# Folder Containing DIC data
#DIC_folders = {"LC": "..\\..\\Geir_Olafsson\\CS02P\\DIC\\Left_Camera_Pair\\Interpolated_Data_150kN", 
#              "RC": "..\\..\\Geir_Olafsson\\CS02P\\DIC\\Right_Camera_Pair\\Interpolated_Data_150kN"}
DIC_folders = {"LC": "E:Calibration_outputs_for_paper\\DIC\\Left_Camera_Pair\\Interpolated_Data_150kN", 
              "RC": "E:Calibration_outputs_for_paper\\DIC\\Right_Camera_Pair\\Interpolated_Data_150kN"}
# Folder Containg GP predictions, interpolated to DIC coordinates
GP_DIC_folder = "E:Calibration_outputs_for_paper\\interp_gp_output" # folder containng interpolated DIC
n_incs = 16 # number of increments
# DIC_force_disp_file = "..\\inputs\\CS02P_mean_disp_412.5mm.csv"
DIC_force_disp_file = "..\\inputs\\CS02P_mean_w_tip_LC.csv"
DIC_force_disp_file2 = "..\\inputs\\CS02P_mean_w_tip_RC.csv"

DIC_force_disp_file_u = "..\\inputs\\CS02P_mean_u_mid_LC.csv"
DIC_force_disp_file_u2 = "..\\inputs\\CS02P_mean_u_mid_RC.csv"

# Load gp mean and standard deviation (used a different name for training gp - fix later)
gp_mean = np.loadtxt("E:Calibration_outputs_for_paper\\Predictions_obserr5e-2\\"+file_str+ "_max_eta_sam_mu_w.csv", delimiter=",", skiprows=0)
gp_sd = np.loadtxt("E:Calibration_outputs_for_paper\\Predictions_obserr5e-2\\"+file_str + "_max_eta_sam_sigma_w.csv", delimiter=",", skiprows=0)
gp_sam = np.loadtxt("E:Calibration_outputs_for_paper\\Predictions_obserr5e-2\\"+file_str + "_max_eta_sam_w.csv", delimiter=",", skiprows=0)

gp_mean_u = np.loadtxt("E:Calibration_outputs_for_paper\\Predictions_obserr5e-2\\"+file_str+ "_max_eta_sam_mu_u.csv", delimiter=",", skiprows=0)
gp_sd_u = np.loadtxt("E:Calibration_outputs_for_paper\\Predictions_obserr5e-2\\"+file_str + "_max_eta_sam_sigma_u.csv", delimiter=",", skiprows=0)
gp_sam_u = np.loadtxt("E:Calibration_outputs_for_paper\\Predictions_obserr5e-2\\"+file_str + "_max_eta_sam_u.csv", delimiter=",", skiprows=0)

# Load in experimental force displacement from loaded end of web
DIC_force_disp = pd.read_csv(DIC_force_disp_file)
DIC_force_disp2 = pd.read_csv(DIC_force_disp_file2)

DIC_force_disp_u = pd.read_csv(DIC_force_disp_file_u)
DIC_force_disp_u2 = pd.read_csv(DIC_force_disp_file_u2)

# Alternative source of DIC data - contains entire test dataset rather than interpolated data used to train model
DIC_all = {}
for key, DIC_folder in DIC_folders.items():
    for i, file in enumerate(os.listdir(DIC_folder)):
        data = pd.read_csv(DIC_folder + "\\" + file)
        # data["point_ind"] = data.index
        #print(data)
        #print(data["index"])
        #data["point_ind"] = data["index"]
        data.columns.values[0] = "point_ind"
        if i == 0:
            DIC_all_i = data
        else:
            DIC_all_i = pd.concat((DIC_all_i, data),axis=0)
    DIC_all[key] = DIC_all_i

# Load in GP_predictions for DIC data points
for file in os.listdir(GP_DIC_folder):
    data = pd.read_csv(GP_DIC_folder + "\\" + file)
    data["Increment"] = int(file.strip("Frame_.csv"))
    # data["Increment"] = int(file.strip("Image_Inc_.csv"))
    try:
        GP_DIC = pd.concat((GP_DIC, data),axis=0)
    except:
        GP_DIC = data

#GP_DIC["Force"] = GP_DIC["Increment"]*10.0
#GP_DIC.loc[GP_DIC["Increment"]==0,"Force"] = 5.39
#print(GP_DIC)

##n_incs = GP_DIC["Increment"].max()
# Create plot
fig, axes = plt.subplots(1,2)
force = 150.0 # Maximum applied load

#<=3.5)
DIC_force_disp = DIC_force_disp[(DIC_force_disp["Load"].abs()<=(force*1.025)) & (DIC_force_disp["Mean DIC Displacement"].abs()<=3.5)]
DIC_force_disp2 = DIC_force_disp2[(DIC_force_disp2["Load"].abs()<=(force*1.025)) & (DIC_force_disp2["Mean DIC Displacement"].abs()<=3.5)]
DIC_force_disp_u = DIC_force_disp_u[(DIC_force_disp_u["Load"].abs()<=(force*1.025)) & (DIC_force_disp_u["Mean DIC Displacement"].abs()<=3.5)]
DIC_force_disp_u2 = DIC_force_disp_u2[(DIC_force_disp_u2["Load"].abs()<=(force*1.025)) & (DIC_force_disp_u2["Mean DIC Displacement"].abs()<=3.5)]

# Plot emulator predictions
y_gp = np.concatenate((np.array([5.39]), np.linspace(10.0, 150.0, 15)))
axes[0].set_title("Calibrated model longitudinal displacement against DIC from near the loaded end")
axes[0].plot(-gp_sam, y_gp, "c", linewidth=0.25, label="Sample")
axes[0].plot(-gp_mean, y_gp, "r", linewidth=1.5, label = "Mean")
# Plot gp mean +/- two standard deviations
axes[0].plot(-gp_mean-2*gp_sd, y_gp, "b", linewidth=1.5, label="95% Interval")
axes[0].plot(-gp_mean+2*gp_sd, y_gp, "b", linewidth=1.5)
axes[0].plot(-DIC_force_disp["Mean DIC Displacement"].to_numpy(),-DIC_force_disp["Load"].to_numpy(),"mx", label = "Experiment (DIC) LC")
axes[0].plot(-DIC_force_disp2["Mean DIC Displacement"].to_numpy(),-DIC_force_disp2["Load"].to_numpy(),"gx", label = "Experiment (DIC) RC")
axes[0].set_ylabel("Force (kN)")
axes[0].set_xlabel("Displacement (mm)")

axes[1].set_title("Calibrated model out-of-plane displacement against DIC from mid-point")
axes[1].plot(-gp_sam_u, y_gp, "c", linewidth=0.25, label="Sample")
axes[1].plot(-gp_mean_u, y_gp, "r", linewidth=1.5, label = "Mean")
# Plot gp mean +/- two standard deviations
axes[1].plot(-gp_mean_u-2*gp_sd_u, y_gp, "b", linewidth=1.5, label="95% Interval")
axes[1].plot(-gp_mean_u+2*gp_sd_u, y_gp, "b", linewidth=1.5)
axes[1].plot(-DIC_force_disp_u["Mean DIC Displacement"].to_numpy(),-DIC_force_disp_u["Load"].to_numpy(),"mx", label = "Experiment (DIC) LC")
axes[1].plot(-DIC_force_disp_u2["Mean DIC Displacement"].to_numpy(),-DIC_force_disp_u2["Load"].to_numpy(),"gx", label = "Experiment (DIC) RC")
axes[1].set_ylabel("Force (kN)")
axes[1].set_xlabel("Displacement (mm)")
# ax.legend()

# Also plot force (longitudinal) displacement for selected DIC data points
point_subset_w = {"LC" : [26805, 34982, 32996], "RC" : [12734 ,33496, 18927]}
point_subset_u = {"LC" : [22841, 23197, 23275, 23294, 23717], "RC" : [8513, 8529, 8566, 15030, 28731]}

# THIS IS GOING TO BE TOO MUCH OF A PAIN ALSO. LOOK AT FINDING THE CLOSEST NODE TO WHERE THE FD IS, EVEN PICK ALL FOUR AND TAKE AN AVERAGE
# extract data tomorrow before leaving. Plot using calibration_res code. Then produce figure, plus observation error figure next week

# print(GP_DIC.columns.values)
for key, points in point_subset_w.items():
    subplots_loop(points, DIC_all[key], GP_DIC[GP_DIC["Camera_pair"]==key], force/n_incs, "w", 1, 3, minus = True)

# Plot minmum vertical displacement
for key, points in point_subset_u.items():
#point_subset_umin = [23275, 23903, 24285] # (first one is maximum of the model, second two are either side of experimental points)
    subplots_loop(points, DIC_all[key], GP_DIC[GP_DIC["Camera_pair"]==key], force/n_incs, "u", 1, 5)
plt.show()
adsadd
# Plot maximum transverse diplacement at both end of the flange tips
point_subset_vmax = [14313, 15408, 22474]
subplots_loop(point_subset_vmax, DIC_all, GP_DIC, force/n_incs, "v", 1, 3)
# Decide on a few points from each
# Points from middle (for u) are [11710, 10738, 10699, 21171, 10650]
# Try and get some for the other bits (nodes)
plt.show()