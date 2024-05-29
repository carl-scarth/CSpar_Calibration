import json
import matplotlib.pyplot as plt
import numpy as np
import os

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
    
    
wd = os.getcwd()
# in_file = "inputs\\nominal_inputs_translator_output_struct"
in_file = "inputs\\nominal_inputs_translator_output_struct"
with open (in_file+".json",'r') as f:
    sample_dict = json.loads(f.readline())

# Reduce frame size when convergence issues
last_frame = -1
last_load = 187.5
for i, sample in enumerate(sample_dict["Sample"]):
    # Automated detection and interpolation code
    print(sample["Frame"][last_frame]["RFs"][2])
    if round(sample["Frame"][last_frame]["RFs"][2],1) != last_load:
#        missing_incs = [0.05*i for i in range(last_frame+1)]
#        split_ind = [j+1 for j,frame in enumerate(sample["Frame"][1:]) if round(frame["RFs"][2] - sample["Frame"][j]["RFs"][2],1) != 12.5 ][0]
#        missing_incs = missing_incs[split_ind:]
        print(str(i))
#        print(missing_incs)
#        print(split_ind)
#        print(sample["Frame"][-1]["RFs"])
#        frame_corrections(sample, split_ind, missing_incs)
    else:
        print("")
        #if len(sample["Frame"]) > last_frame:
        #    sample["Frame"] = sample["Frame"][0:last_frame+1]

for i, sample in enumerate(sample_dict["Sample"]):
    # print(len(sample["Frame"]))
    print(i)
    print(sample["Frame"][-1]["RFs"][2])

# Code for subtracting displacements
# to undo and make new predictions just add value back on
for i, sample in enumerate(sample_dict["Sample"]):
    zero_disp = np.array(sample["Frame"][1]["Displacements"])
    for j, frame in enumerate(sample["Frame"]):
        frame["Unzeroed_disp"] = frame["Displacements"]
        frame["Displacements"] = (np.array(frame["Displacements"]) - zero_disp).tolist()
    # Delete first frame
    sample["Frame"] = sample["Frame"][1:]

        
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

#init_grad = [sample[1] for sample in RF_z]
with open(in_file+"_zeroed.json",'w') as f:
    f.write(json.dumps(sample_dict))
