import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json

file_str = "nominal_inputs_new_spar"
max_load = -300.0
max_inc = 0.05
displacements = np.loadtxt("inputs\\" + file_str + "_displacements_load=" + str(max_load) + "_max_inc=" + str(max_inc) +".csv",delimiter=",",skiprows=1)
# displacements = pd.read_csv("inputs\\LHSDesign60x6_displacements_load=-250.0_max_inc=0.1.csv",sep=",")
RFs = np.loadtxt("inputs\\" + file_str + "_RFs_load=" + str(max_load) + "_max_inc=" + str(max_inc) + ".csv", delimiter=',', skiprows=1)
with open ("inputs\\" + file_str + "_incs_load=" + str(max_load) + "_max_inc=" + str(max_inc) + ".txt",'r') as f:
    increments = [[float(increment) for increment in line.strip().split(',')] for line in f.readlines()]

# Play around with different data structures: json or dataframe
displacements_struct = {"Sample" : []}
# Creates an empty dataframe with all the correct column names. Not sure how to add to these in the loop
displacements_df = pd.DataFrame(columns=["Sample","Frame","Increment","Node","RF1","RF2","RF3","U","V","W"])
frame_cnt = 0
# Dataframe seems slower...
for i, sample in enumerate(increments):
    displacements_struct["Sample"].append({"Frame" : []})
    # displacements_df["Sample"].append(i) # Doesn't work. Can append a series though. Look at other ways of constructing pandas dataframes on the fly
    for j, frame in enumerate(sample):
        displacements_struct["Sample"][i]["Frame"].append({"Increment" : frame,
                                                           "RFs" : RFs[3*frame_cnt:3*(frame_cnt+1)].tolist(),
                                                           "Displacements" : displacements[:,3*frame_cnt:3*(frame_cnt+1)]})
        # Altermative format using Pandas dataframe
        uvw_df = pd.DataFrame(displacements[:,3*frame_cnt:3*(frame_cnt+1)], columns=["U","V","W"])
        uvw_df["Sample"] = i
        uvw_df["Frame"] = j
        uvw_df["Increment"] = frame
        uvw_df = pd.concat([uvw_df, pd.DataFrame(np.broadcast_to(RFs[3*frame_cnt:3*(frame_cnt+1)],(uvw_df.shape[0],3)), columns=["RF1","RF2","RF3"])], axis=1)
        uvw_df["Node"] = range(uvw_df.shape[0])
        uvw_df[displacements_df.columns]
        displacements_df = pd.concat([displacements_df,uvw_df])#[displacements_df.columns]])
        
        frame_cnt += 1

# Number of frames in each sample
n_frames = [len(sample["Frame"]) for sample in displacements_struct["Sample"]]
# Final increment for each sample (to check which runs didn't complete)
end_incs = [sample["Frame"][-1]["Increment"] for sample in displacements_struct["Sample"]]

# Alterative formatting using pandas
n_samples = displacements_df["Sample"].max() + 1
# This way is quite slow
n_frames_df = [displacements_df.loc[displacements_df["Sample"]==i]["Frame"].max() + 1 for i in range(n_samples)]
end_incs_df = [displacements_df.loc[displacements_df["Sample"]==i]["Increment"].max() for i in range(n_samples)]

displacements_sam_0 = displacements_df.loc[displacements_df["Sample"]==0]
# Continue the below with Pandas formatting - is it better?

# Create a different structure with increments of 0.8 or lower
# First create a copy of the dictionary to prevent modifying the original
# This is slightly awkward as it's also necessary to change the entries in the dictionary
displacements_struct_subset = displacements_struct.copy()
displacements_struct_subset["Sample"] = displacements_struct_subset["Sample"].copy()
displacements_struct_subset["Sample"] = [sample.copy() for sample in displacements_struct_subset["Sample"]]
for sample in displacements_struct_subset["Sample"]:
    # Pick 0.81 as the increment of interest is very slightly over 0.80
    # 
    sample["Frame"] = [frame for frame in sample["Frame"] if frame["Increment"]<=1.01]
    
# Look at the number of frames and increment at the final frame to check that 
# the following would be suitable for  an emulator trained using fixed increments
n_frames_subset = [len(sample["Frame"]) for sample in displacements_struct_subset["Sample"]]
end_incs_subset = [sample["Frame"][-1]["Increment"] for sample in displacements_struct_subset["Sample"]]

# Plot force-displacement of reference point at the spar tip
# Extract displacement (w) at reference point for all increments in each sample
ref_w = [[-frame["Displacements"][1,2] for frame in sample["Frame"]]for sample in displacements_struct["Sample"]]
RF_z = [[frame["RFs"][2] for frame in sample["Frame"]]for sample in displacements_struct["Sample"]]

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)
for i, sample in enumerate(ref_w):
    ax.plot(sample, RF_z[i])

label_font = {'family': 'serif', 'size': 16,}
ax.set_ylabel("Force (kN)", fontdict = label_font)
ax.set_xlabel("End axial displacement (mm)", fontdict = label_font)

max_u = [[frame["Displacements"][:,0].max() for frame in sample["Frame"]]for sample in displacements_struct["Sample"]]
min_u = [[frame["Displacements"][:,0].min() for frame in sample["Frame"]]for sample in displacements_struct["Sample"]]

fig2 = plt.figure(figsize=(10,8))
ax2 = fig2.add_subplot(1, 1, 1)
for i, sample in enumerate(max_u):
    ax2.plot(sample, RF_z[i])

ax2.set_ylabel("Force (kN)", fontdict = label_font)
ax2.set_xlabel("Maximum vertical displacement (mm)", fontdict = label_font)

fig3 = plt.figure(figsize=(10,8))
ax3 = fig3.add_subplot(1, 1, 1)
for i, sample in enumerate(min_u):
    ax3.plot(sample, RF_z[i])

ax3.set_ylabel("Force (kN)", fontdict = label_font)
ax3.set_xlabel("Minimum vertical displacement (mm)", fontdict = label_font)

plt.show()
for sample in displacements_struct_subset["Sample"]:
    for frame in sample["Frame"]:
        frame["Displacements"] = frame["Displacements"].tolist()

with open("inputs\\" + file_str + "_output_struct.json",'w') as f:
    f.write(json.dumps(displacements_struct_subset))
    # f.write(json.dumps(displacements_struct))
