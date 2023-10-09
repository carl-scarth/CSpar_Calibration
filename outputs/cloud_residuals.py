import pandas as pd
import numpy as np
import os

# Read in csv files containing point cloud of Gaussian process predictions, and produce and output cloud including residuals with experiental data
label = "LHSDesign40x4"
# Booleans indicating which outputs are to be read in, and combined into the relevent csv files
# training_data_mean = False # Plot mean of training data (This is done automatically in R so maybe don't need this)
GP_Mean = False # Average across posterior samples of Gaussian process predictions
GP_Std = False # Standard deviation of posterior predictive samples
GP_Post_Sam = True # Gaussian process predictions across a set of posterior samples
GP_Sam = True # Samples from Gaussian process
file_series = True # Output as a file series for animating in paraview 

cloud_file = "Image_0074.csv"
# Column headings for coordinates and displacements
coord_cols = ["x_0_rot","y_0_rot","z_0_rot"]
disp_col = "w_rot"

file_list = []
if GP_Mean:
    file_list.extend(["eta_mu_y_mu","eta_sigma_y_mu"])
    if GP_Sam:
        file_list.append("eta_sam_y_mu")
if GP_Std:
    file_list.append("eta_mu_y_sigma")
    if GP_Sam:
        file_list.append("eta_sam_y_sigma")
if GP_Post_Sam:
    file_list.extend(["eta_mu_y", "eta_sigma_y"])
    if GP_Sam:
        file_list.append("eta_sam_y")

# Read in point cloud data file
in_cloud = pd.read_csv(cloud_file, sep=",")
out_cloud = in_cloud.filter(items=coord_cols)
uvw = in_cloud[disp_col].values.reshape(-1,1)
out_cloud[disp_col] = uvw
if file_series:
    frame_list = [] # For storing dataframes when
# Read in GP output quantites and append to the data
for i, filestr in enumerate(file_list):
    # load in output from current filename
    output = pd.read_csv(filestr + "_" + label + ".csv")
    res = output.values - uvw   # Residual
    abs_res = np.absolute(res)  # Absolute value of residual
    rel_err = abs_res/uvw*100.0 # Percentage error
    if file_series:
        if "output_csv" not in os.listdir(os.getcwd()):
            os.mkdir("output_csv")
        for j, col in enumerate(output.columns):
            if i == 0:
                frame_list.append(out_cloud.copy())
                
            frame_list[j][filestr] = output[col]
            if "sigma" not in filestr:
                frame_list[j][filestr + "_res"] = res[:,j]
                frame_list[j][filestr + "_abs_res"] = abs_res[:,j]
                frame_list[j][filestr + "_rel_err"] = rel_err[:,j]
    else:
        if output.shape[1] == 1:
            out_cloud[filestr] = output
        else:
            out_cloud = pd.concat((out_cloud, output), axis=1)
        # Calculate residual and append to data if quantity is not a standard deviation
        if "sigma" not in filestr:
            if output.shape[1] == 1:
                out_cloud[filestr+"_res"] = res
                out_cloud[filestr+"_abs_res"] = abs_res
                out_cloud[filestr+"_rel_err"] = rel_err
            else:
                out_cloud = pd.concat((out_cloud, pd.DataFrame(res,columns=[col + "_res" for col in output.columns.values])), axis=1)
                out_cloud = pd.concat((out_cloud, pd.DataFrame(abs_res,columns=[col + "abs_res" for col in output.columns.values])), axis=1)
                out_cloud = pd.concat((out_cloud, pd.DataFrame(rel_err,columns=[col + "_rel_err" for col in output.columns.values])), axis=1)
            
# Write output to csv
if file_series:
    for i, frame in enumerate(frame_list):
        frame.to_csv("output_csv\\Sample_"+str(i)+".csv",sep=",",index=False)
else:
    out_cloud.to_csv("output_csv.csv",sep=",",index=False)