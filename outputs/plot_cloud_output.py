import pandas as pd

# Combines output with point cloud coordinates for visualising in paraview
cloud_file = "..\\inputs\\Interpolated_DIC_200kN.csv"
coord_str = ["x_proj","y_proj","z_proj"]
# data_file = "training_data_mean_LHSDesign40x4"
data_files = ["eta_mu_LHSDesign40x4.csv","eta_sigma_LHSDesign40x4.csv","eta_sam_LHSDesign40x4.csv","eta_mu_mu_LHSDesign40x4.csv","eta_sigma_mu_LHSDesign40x4.csv"]

cloud_data = pd.read_csv(cloud_file)
# initialise output frame using cloud coordinates
data = cloud_data[coord_str]
for file in data_files:
    file_data = pd.read_csv(file)
    data = pd.concat((data, file_data), axis=1)

data.to_csv("cloud_output.csv",sep=",",index=False)
