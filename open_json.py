import json
import matplotlib.pyplot as plt
import numpy as np

def lin_interp(before, after, fraction):
    interp = before + (after - before)*fraction
    return(interp)

def frame_corrections(sample, split_ind, missing_incs):
    correct_frames = sample["Frame"][0:split_ind]
    remaining_frames = sample["Frame"][split_ind:]
    for inc in missing_incs:
        before_frame = [frame for frame in remaining_frames if frame["Increment"] <= inc][-1]
        after_frame = [frame for frame in remaining_frames if frame["Increment"] > inc][0]
        fraction = (inc - before_frame["Increment"])/(after_frame["Increment"] - before_frame["Increment"])
        intp_frame = before_frame.copy()
        intp_frame["Increment"] = lin_interp(before_frame["Increment"], after_frame["Increment"], fraction)
        intp_frame["RFs"] = [lin_interp(RF1, RF2, fraction) for RF1, RF2 in zip(before_frame["RFs"],after_frame["RFs"])]
        intp_frame["Displacements"] = [[lin_interp(d1, d2, fraction) for d1, d2 in zip(d1_row, d2_row)] for d1_row, d2_row in zip(before_frame["Displacements"],after_frame["Displacements"])]
        correct_frames.append(intp_frame)
    sample["Frame"] = correct_frames
    
    
# in_file = "inputs\\LHSDesign25x1_1_output_struct_disp"
in_file = "E:Working_folder\\LHSDesign60x6_7_output_struct_200kN"
#in_file = "outputs\\gp_predictions_nonlinear_LHSDesign40x4.json"
with open (in_file+".json",'r') as f:
    sample_dict = json.loads(f.readline())

# Reduce frame size when convergence issues
#for i, sample in enumerate(sample_dict["Sample"]):
#    print([frame["RFs"][2] for frame in sample["Frame"][0:17]])
#    # Interpolation code used for LHSDesign60x6_7
#    if i == 15:
#        missing_incs = [0.75, 0.8]
#        split_ind = 15
#        frame_corrections(sample, split_ind, missing_incs)
#    elif i == 58:
#        missing_incs = [0.8]
#        split_ind = 16
#        frame_corrections(sample, split_ind, missing_incs)
#    else:
#        if len(sample["Frame"]) > 16:
#            sample["Frame"] = sample["Frame"][0:17]
#    print(sample["Frame"][-1]["Increment"])

# Reference point displacement for old approach
# ref_w = [[-frame["Displacements"][1][2] for frame in sample["Frame"]]for sample in sample_dict["Sample"]]
# Refernce point displacement for new approach
ref_w_pivot = [[-frame["Displacements"][-1][2] for frame in sample["Frame"]] for sample in sample_dict["Sample"]]
ref_w_end = [[-frame["Displacements"][-3][2] for frame in sample["Frame"]] for sample in sample_dict["Sample"]]
RF_z = [[frame["RFs"][2] for frame in sample["Frame"]] for sample in sample_dict["Sample"]]

fig, axs = plt.subplots(1,2)
for ax, ref_w in zip(axs, (ref_w_pivot, ref_w_end)):
    for i, sample in enumerate(ref_w):
        if min(RF_z[i]) < -1:
            print(i)
        if max(RF_z[i]) <= 300:
            ax.plot(sample, RF_z[i], marker = "x")

plt.show()
asdsadsad

#init_grad = [sample[1] for sample in RF_z]

with open(in_file+"_200kN.json",'w') as f:
    f.write(json.dumps(sample_dict))

