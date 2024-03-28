import json
import matplotlib.pyplot as plt

in_file = "inputs\\LHSDesign60x6_4_output_struct"
#in_file = "outputs\\gp_predictions_nonlinear_LHSDesign40x4.json"
with open (in_file+".json",'r') as f:
    sample_dict = json.loads(f.readline())

# Reduce frame size when convergence issues
for sample in sample_dict["Sample"]:
    # print(sample["Frame"][-1]["Increment"])
    print(len(sample["Frame"]))
    if len(sample["Frame"]) > 16:
        sample["Frame"] = sample["Frame"][0:17]
    print(sample["Frame"][-1]["RFs"][2])
    
# Reference point displacement for old approach
# ref_w = [[-frame["Displacements"][1][2] for frame in sample["Frame"]]for sample in sample_dict["Sample"]]
# Refernce point displacement for new approach
ref_w_pivot = [[-frame["Displacements"][-1][2] for frame in sample["Frame"]] for sample in sample_dict["Sample"]]
ref_w_end = [[-frame["Displacements"][-3][2] for frame in sample["Frame"]] for sample in sample_dict["Sample"]]
RF_z = [[frame["RFs"][2] for frame in sample["Frame"]] for sample in sample_dict["Sample"]]

# =============================================================================
# fig = plt.figure(figsize=(10,8))
# ax = fig.add_subplot(1, 1, 1)
# for i, sample in enumerate(ref_w_pivot):
#     if min(RF_z[i]) < -1:
#         print(i)
#     if max(RF_z[i]) <= 300:
#         ax.plot(sample, RF_z[i])
#         
# =============================================================================

fig, axs = plt.subplots(1,2)
for ax, ref_w in zip(axs, (ref_w_pivot, ref_w_end)):
    for i, sample in enumerate(ref_w):
        if min(RF_z[i]) < -1:
            print(i)
        if max(RF_z[i]) <= 300:
            ax.plot(sample, RF_z[i])



# Only one sample is below 225. Try checking other 60 example
# Mode switch occurs well above 225, maybe I can just ignore?
cczxcxzczxc

with open(in_file+"_200kN.json",'w') as f:
    f.write(json.dumps(sample_dict))