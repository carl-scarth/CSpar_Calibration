import pandas as pd
import numpy as np
import os

split_frames = True # Is output comprised of multiple frames (e.g. load steps)
split_cam_pairs = True # Separate output for different camera pairs
#file_str = "LHSDesign100x8"
file_str = "nominal_inputs_translator"
# Combines output with point cloud coordinates for visualising in paraview
# Modify for multiple input point clouds
cloud_dir = "..\\inputs\\"
#cloud_file = "..\\inputs\\Interpolated_DIC_downsam_12.csv"
# Need to keep in same order as calibration code
cloud_files = ["E:Calibration_outputs_for_paper\\Interpolated_DIC_multistep_150kN_LC.csv", "E:Calibration_outputs_for_paper\\Interpolated_DIC_multistep_150kN_RC.csv"]
coord_str = ["x_proj","y_proj","z_proj"]
# Need to account for multiple, if not already
disp_str = ["u","w"]
disp_label = "rot"

#data_files = ["eta_sam_y_mu"]
data_files = ["ABAQUS"]


#res_label = ["eta_sam_y_mu"] # label of files for which residual calculations are required, if any
res_label = ["ABAQUS"]
cloud_data = [pd.read_csv(os.path.join(cloud_dir, cloud_file)) for cloud_file in cloud_files]

# initialise output frame using cloud coordinates
data = pd.DataFrame(columns=cloud_data[0].columns.values)
n_y = []
for frame in cloud_data:
    data = pd.concat((data,frame),axis=0)
    n_y.append(frame.shape[0])


#displacement = cloud_data[disp_str].to_numpy()
displacement = data[[disp+"_"+disp_label for disp in disp_str]].to_numpy()
print("here")
for file in data_files:
    print(file)
    file_data = pd.read_csv("E:Calibration_outputs_for_paper\\"+file+"_"+file_str+".csv")
    print(file_data)
    print(file_data.columns.values)
    file_data_cols = file_data.columns.values
    new_cols = np.empty((data.shape[0],len(disp_str)*file_data.shape[1]))
    for i in range(len(cloud_files)):
        for j, comp in enumerate(disp_str):
            file_data_ij = file_data.iloc[(j*n_y[i]+2*i*sum([n_y[i] for i in range(i)])):(j+1)*n_y[i]+2*i*sum([n_y[i] for i in range(i)])]
            file_data_ij.columns = ['_'.join((col_name, comp)) for col_name in file_data_cols]
            new_cols[i*sum([n_y[i] for i in range(i)]):sum(n_y[i] for i in range(i+1)),j*file_data.shape[1]:(j+1)*file_data.shape[1]] = file_data_ij
            
    new_cols = pd.DataFrame(new_cols,columns= ['_'.join((col_name,comp)) for comp in disp_str for col_name in file_data_cols])
    data = pd.concat((data, new_cols.set_index(data.index)), axis=1)
    
    if file in res_label:
        for i, comp in enumerate(disp_str):
            print(i)
            print(comp)
            # Calculate residual
            res = new_cols[['_'.join((col_name,comp)) for col_name in file_data_cols]].to_numpy() - displacement[:,i].reshape((-1,1))
            abs_res = np.absolute(res)  # Absolute value of residual
            per_err = abs_res/np.absolute(displacement[:,i].reshape((-1,1)))*100.0 # Percentage error
            res = pd.DataFrame(res, columns=["_".join((col_name, comp, "res")) for col_name in file_data_cols])
            abs_res = pd.DataFrame(abs_res, columns=["_".join((col_name, comp, "absres")) for col_name in file_data_cols])
            per_err = pd.DataFrame(per_err, columns=["_".join((col_name, comp, "pererr")) for col_name in file_data_cols])
            data = pd.concat((data, res.set_index(data.index), abs_res.set_index(data.index), per_err.set_index(data.index)), axis=1)

# If data is split across several frames, output to multiple csvs.
if split_cam_pairs:
    cloud_data = [data.iloc[sum([n_y[i] for i in range(i)]):sum([n_y[i] for i in range(i+1)]),:] for i in range(len(n_y))]
else:
    cloud_data = [data]

for i, dataset in enumerate(cloud_data):
    if split_frames:
        inc_ind = dataset.Increment
        frames = inc_ind.unique()
        if len(cloud_data) == 1:
            out_dir = "cloud_output"
        else:
            out_dir = "cloud_output_dataset_" + str(i)
        #if out_dir not in os.listdir(os.getcwd()):
        if out_dir not in os.listdir("E:Calibration_outputs_for_paper"):
            os.mkdir("E:Calibration_outputs_for_paper\\"+out_dir)
        for frame_ind in frames:
            frame_data = dataset[inc_ind==frame_ind]
            frame_data.to_csv(os.path.join("E:Calibration_outputs_for_paper\\"+out_dir,"frame_"+str(frame_ind)+".csv"),sep=",",index=False)
    else:
        if len(cloud_data) == 1:
            dataset.to_csv("cloud_output.csv",sep=",",index=False)
        else:
            dataset.to_csv("cloud_output_dataset_"+str(i)+".csv",sep=",",index=False)