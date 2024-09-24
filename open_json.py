import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def lin_interp(before, after, fraction):
    interp = before + (after - before)*fraction
    return(interp)

def frame_corrections(sample, split_ind, missing_incs):
    correct_frames = sample["Frame"][0:split_ind]
    remaining_frames = sample["Frame"][split_ind:]
    for inc in missing_incs:
        before_frame = [frame for frame in remaining_frames if frame["Increment"] <= inc][-1]
        print(before_frame)
        print([frame["Increment"] for frame in remaining_frames])
        print(inc)
        after_frame = [frame for frame in remaining_frames if frame["Increment"] > inc][0]
        print(after_frame)
        # This is as far as I've got with the new version of the code - I think it's working up to here, but I haven't had a proper chance to check
        asdsadasd
        fraction = (inc - before_frame["Increment"])/(after_frame["Increment"] - before_frame["Increment"])
        intp_frame = before_frame.copy()
        intp_frame["Increment"] = lin_interp(before_frame["Increment"], after_frame["Increment"], fraction)
        intp_frame["RFs"] = [lin_interp(RF1, RF2, fraction) for RF1, RF2 in zip(before_frame["RFs"],after_frame["RFs"])]
        intp_frame["Displacements"] = [[lin_interp(d1, d2, fraction) for d1, d2 in zip(d1_row, d2_row)] for d1_row, d2_row in zip(before_frame["Displacements"],after_frame["Displacements"])]
        correct_frames.append(intp_frame)
    sample["Frame"] = correct_frames
    
    
wd = os.getcwd()
# in_file = "inputs\\nominal_inputs_translator_output_struct"
# in_file = "Cspar_sam_output"
in_file = "inputs\\LHSDesign100x10_1"
with open (in_file+"_output_struct.json",'r') as f:
    sample_dict = json.loads(f.readline())
asdsads


# Code for deleting samples if required
inputs = pd.read_csv(in_file+".csv")

# Delete multiple samples
keep_ind = [i for i in range(len(sample_dict["Sample"]))]
del_ind = [13,33,48,57,65]
[keep_ind.remove(i) for i in del_ind]
sample_dict["Sample"] = [sample_dict["Sample"][ind] for ind in keep_ind]
[inputs.drop(index = ind, inplace = True) for ind in del_ind]

# Delete only one sample
sample_dict["Sample"].pop(13)
inputs = inputs.drop(index = 13)
inputs.to_csv(in_file + "_subset200kN.csv",index=False,sep=",")


print(sample_dict["Sample"][46]["Frame"][-1]["RFs"])

# Reduce frame size when convergence issues
last_frame = 17
last_load = 160
# n_frames_min = 23
for i, sample in enumerate(sample_dict["Sample"]):
    #print(i)
    # Automated detection and interpolation code
    if round(sample["Frame"][last_frame]["RFs"][2],0) != last_load:# or len(sample["Frame"])> n_frames_min:
        print(i)
        print(sample["Frame"][last_frame]["RFs"][2])
        bvbvcb
        # missing_incs = [0.05*i for i in range(last_frame+1)]
        # split_ind = [j+1 for j,frame in enumerate(sample["Frame"][1:]) if round(frame["RFs"][2] - sample["Frame"][j]["RFs"][2],1) != 12.5 ][0]
        # missing_incs = missing_incs[split_ind:]
        # frame_corrections(sample, split_ind, missing_incs)
        missing_incs = [0.05*i for i in range(last_frame-1)] # for version with first two increments
        split_ind = [j+1 for j,frame in enumerate(sample["Frame"][3:]) if round(frame["RFs"][2] - sample["Frame"][j+2]["RFs"][2],0) != 10 ][0]
    
        # The below is correct for both versions of the model

        missing_incs = missing_incs[split_ind:]

        #print(len(sample["Frame"]))
#        print(missing_incs)
#        print(split_ind)
#        print(sample["Frame"][-1]["RFs"])
        # Uncomment to use with multistep
        #frame_corrections(sample, split_ind+2, missing_incs)
        
        #if i == 13:
        #    # Fudge
        #    print("Fudge!")
        #    correct_frames = sample["Frame"][1:split_ind+2]
        #    final_frame = sample["Frame"][-1]
        #    final_frame["Increment"] = round(final_frame["Increment"],1)
        #    final_frame["RFs"][2] = round(final_frame["RFs"][2],-1)
        #    correct_frames.append(final_frame)
        #    sample["Frame"] = correct_frames
        
    else:
        print("")
        if len(sample["Frame"]) > last_frame:
            sample["Frame"] = sample["Frame"][0:last_frame+1]
            

for i, sample in enumerate(sample_dict["Sample"]):
    print(i)
    print(len(sample["Frame"]))
    print(sample["Frame"][-1]["RFs"][2])
    #if round(sample["Frame"][-1]["RFs"][2]) < 200:
    #    print(i)
    #    print(sample["Frame"][-1]["RFs"][2])
last_frame = [sample["Frame"][-1]["RFs"][2] for sample in sample_dict["Sample"]]
print(min(last_frame))

# Code for subtracting displacements
# to undo and make new predictions just add value back on
for i, sample in enumerate(sample_dict["Sample"]):
    zero_disp = np.array(sample["Frame"][1]["Displacements"])
    for j, frame in enumerate(sample["Frame"]):
        #frame["Unzeroed_disp"] = frame["Displacements"]
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
with open(in_file+"_output_struct_160kN_zeroed.json",'w') as f:
    f.write(json.dumps(sample_dict))
