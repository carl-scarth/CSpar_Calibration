import json
import matplotlib.pyplot as plt
import pandas as pd

infile = "nominal_inputs_new_spar"
output_file = infile + "_output_struct.json"
RP_node = 1 # node index for reference point where load is applied
labels = ["NIAR database properties", "Meng Yi Properties"]

# Open json
with open(output_file, "r") as f:
    # Load in string from file
    in_dict = json.loads(f.readline())  

RFs, w_RP = ([], [])
for sample in in_dict["Sample"]:
    RFs.append([frame["RFs"][2] for frame in sample["Frame"]])
    w_RP.append([-frame["Displacements"][RP_node][2] for frame in sample["Frame"]])

# print(w_RP)

# Load in Meng Yi Force displacement from the fileshare
mengyi_fd = pd.read_csv("FD_MengYi_plus4.csv",header=None).to_numpy()
exp_fd = pd.read_csv("CS02P_Failure.csv").to_numpy()
print(exp_fd)

# Plot reference point displacement
fig, ax = plt.subplots(figsize=(10,4))
max_load = 300.0
for force, disp, label in zip(RFs, w_RP, labels):
    ax.plot(disp, force, linewidth=2.0, label=label)
ax.plot(mengyi_fd[:,0], mengyi_fd[:,1], linewidth=2.0, label = "MengYi Model")
ax.plot(exp_fd[:,0], -exp_fd[:,1], "bx", label = "Experiment") 

ax.set_ylabel("Force (kN)")
ax.set_xlabel("RP Displacement (mm)")
ax.legend()
plt.show()