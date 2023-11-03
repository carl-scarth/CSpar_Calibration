import pandas as pd
import numpy as np
import os

split_frames = True # Is output comprised of multiple frames (e.g. load steps)

file_str = "LHSDesign40x4"
# Combines output with point cloud coordinates for visualising in paraview
cloud_file = "..\\inputs\\Interpolated_DIC_downsam_12.csv"
coord_str = ["x_proj","y_proj","z_proj"]
disp_str = ["w_rot"]
# data_file = "training_data_mean_LHSDesign40x4"
# data_files = ["eta_mu_LHSDesign40x4.csv","eta_sigma_LHSDesign40x4.csv","eta_sam_LHSDesign40x4.csv","eta_mu_mu_LHSDesign40x4.csv","eta_sigma_mu_LHSDesign40x4.csv"]
# data_files = ["eta_mu_mu_LHSDesign40x4.csv","eta_sigma_mu_LHSDesign40x4.csv"]
data_files = ["eta_sam_y", 
              "eta_sigma_y",
              "eta_mu_y", 
              "eta_mu_y_mu", 
              "eta_mu_y_sigma", 
              "eta_sam_y_mu",
              "eta_sam_y_sigma",
              "eta_sigma_y_mu"]

res_label = ["eta_sam_y","eta_sam_y_mu"] # label of files for which residual calculations are required, if any

cloud_data = pd.read_csv(cloud_file)

# initialise output frame using cloud coordinates
data = cloud_data[coord_str]
displacement = cloud_data[disp_str].to_numpy()

for file in data_files:
    file_data = pd.read_csv(file+"_"+file_str+".csv")
    data = pd.concat((data, file_data), axis=1)
    if file in res_label:
        print(file)
        res = file_data.to_numpy() - displacement # Residual
        abs_res = np.absolute(res)  # Absolute value of residual
        per_err = abs_res/np.absoluate(displacement)*100.0 # Percentage error
        res = pd.DataFrame(res, columns=[col_label + "_res" for col_label in file_data.columns.values])
        abs_res = pd.DataFrame(abs_res, columns=[col_label + "_absres" for col_label in file_data.columns.values])
        per_err = pd.DataFrame(per_err, columns=[col_label + "_pererr" for col_label in file_data.columns.values])
        data = pd.concat((data, res, abs_res, per_err), axis=1)

# If data is split across several frames, output to multiple csvs.
if split_frames:
    inc_ind = cloud_data.Increment
    frames = inc_ind.unique()
    if "cloud_output" not in os.listdir(os.getcwd()):
        os.mkdir("cloud_output")
    for frame_ind in frames:
        frame_data = data[inc_ind==frame_ind]
        frame_data.to_csv("cloud_output\\frame_"+str(frame_ind)+".csv",sep=",",index=False)
else:
    data.to_csv("cloud_output.csv",sep=",",index=False)
