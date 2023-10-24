# Extracts output at a fixed load value from a json file
import json
import pandas as pd

file_str = "LHSDesign40x4"
with open("inputs\\" + file_str + "_output_struct.json",'r') as f:
    output_dict = json.loads(f.readline())
    
load = 200 # Desired applied load Value
out_frame = pd.DataFrame() # Create empty dataframe, the fill by looping over samples
for i, sample in enumerate(output_dict["Sample"]):
    disp_i = pd.DataFrame([frame for frame in sample["Frame"] if round(frame["RFs"][2]) == load][0]["Displacements"],columns=["u_"+str(i+1),"v_"+str(i+1),"w_"+str(i+1)])
    print([frame for frame in sample["Frame"] if round(frame["RFs"][2]) == load][0]["RFs"])
    out_frame = pd.concat((out_frame, disp_i), axis=1)

out_frame.to_csv("inputs\\" + file_str +"_fixed_" + str(load) + "kN.csv", index=False)