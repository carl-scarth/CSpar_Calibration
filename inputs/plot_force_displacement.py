import json
import matplotlib.pyplot as plt
import pandas as pd

infile = "nominal_inputs_new_spar"
# output_file = infile + "_output_struct.json" 
output_file_2 = infile + "_output_struct_disp.json"
prior_file = "LHSDesign50x3_2_output_struct_disp.json"
RP_node = 1 # node index for reference point where load is applied
labels = ["E_11,C = 140.9 kN", "E_11,C = 150 kN"]

# Open json
#with open(output_file, "r") as f:
    # Load in string from file
#    in_dict = json.loads(f.readline())  
with open(output_file_2, "r") as f:
    # Load in string from file
    in_dict_2 = json.loads(f.readline())  

with open(prior_file, "r") as f:
    # Load in string from file
    in_dict_prior = json.loads(f.readline())  

#RFs, w_RP = ([], [])
#for sample in in_dict["Sample"]:
#    RFs.append([frame["RFs"][2] for frame in sample["Frame"]])
#    w_RP.append([-frame["Displacements"][RP_node][2] for frame in sample["Frame"]])

RFs2, w_RP2 = ([], [])
for sample in in_dict_2["Sample"]:
    RFs2.append([frame["RFs"][2] for frame in sample["Frame"]])
    w_RP2.append([-frame["Displacements"][RP_node][2] for frame in sample["Frame"]])

RFs_prior, w_RP_prior = ([], [])
for sample in in_dict_prior["Sample"]:
    RFs_prior.append([frame["RFs"][2] for frame in sample["Frame"]])
    w_RP_prior.append([-frame["Displacements"][RP_node][2] for frame in sample["Frame"]])

# print(w_RP)

# Load in Meng Yi Force displacement from the fileshare
mengyi_fd = pd.read_csv("FD_MengYi_plus4.csv",header=None).to_numpy()
exp_fd = pd.read_csv("CS02P_Failure.csv").to_numpy()
print(exp_fd)

# Plot reference point displacement
fig, ax = plt.subplots(figsize=(10,4))
#for force, disp, label in zip(RFs, w_RP, labels):
#    ax.plot(disp, force, linewidth=2.0, label=label)

for force, disp in zip(RFs_prior, w_RP_prior):
    ax.plot(disp, force, linewidth=1.0)

for force, disp, label in zip(RFs2, w_RP2, labels):
    ax.plot(disp, force, linewidth=4.0, label=label)
    
# ax.plot(mengyi_fd[:,0], mengyi_fd[:,1], linewidth=2.0, label = "Meng Yi's Model")
ax.plot(exp_fd[:,0], -exp_fd[:,1], "bx", label = "Experiment (DIC)") 

ax.set_ylabel("Force (kN)")
ax.set_xlabel("RP Displacement (mm)")
ax.legend()
plt.show()