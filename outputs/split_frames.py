# Reads in a csv with point cloud data for different frames, and splits it into separate csvs (one for each frame)
# for reading into paraview

import pandas as pd
import os

file_str = "LHSDesign40x4"
# quantity = "mean_error"
quantity = "interpolated_basis"

cloud_file = "_".join((file_str, quantity, "nonlinear.csv"))
#cloud_file = "Interpolated_DIC_inc2_downsam4.csv"
out_folder = "_".join((quantity, file_str))
if out_folder not in os.listdir(os.getcwd()):
    os.makedirs(out_folder)

cloud_data = pd.read_csv(cloud_file, sep = ",")
n_inc = cloud_data["Increment"].max() # Number of increments
for i in range(n_inc):
    print(i)
    data_i = cloud_data[cloud_data["Increment"]==(i+1)]
    data_i.to_csv(os.path.join(out_folder,"Increment_"+str(i)+".csv"), sep=",", index=False)
