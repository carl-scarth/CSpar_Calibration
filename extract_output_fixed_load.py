# Extracts output at a fixed load value from a json file
import json
import pandas as pd

specify_load = True
if not specify_load:
    specify_inc = True
else:
    specify_inc = False

load = 200 # Desired applied load Value
inc = 2

file_str = "LHSDesign40x4_1"
with open("inputs\\" + file_str + "_output_struct.json",'r') as f:
    output_dict = json.loads(f.readline())

if specify_load:
    out_frame = pd.DataFrame() # Create empty dataframe, the fill by looping over samples
    for i, sample in enumerate(output_dict["Sample"]):
        disp_i = pd.DataFrame([frame for frame in sample["Frame"] if round(frame["RFs"][2]) == load][0]["Displacements"],columns=["u_"+str(i+1),"v_"+str(i+1),"w_"+str(i+1)])
        out_frame = pd.concat((out_frame, disp_i), axis=1)
        print([frame for frame in sample["Frame"] if round(frame["RFs"][2]) == load][0]["RFs"])
    # Write to csv
    out_frame.to_csv("inputs\\" + file_str +"_fixed_" + str(load) + "kN.csv", index=False)
elif specify_inc:
    out_dict = {"Sample" : []}
    for i, sample in enumerate(output_dict["Sample"]):
        # [print(frame["RFs"]) for j, frame in enumerate(sample["Frame"]) if j % inc == 0]
        out_dict["Sample"].append({"Frame": [frame for j, frame in enumerate(sample["Frame"]) if j % inc == 0]})
    # Write to json
    with open("inputs\\" + file_str + "_downsam_" + str(inc) + ".json",'w') as f:
        f.write(json.dumps(out_dict))        
