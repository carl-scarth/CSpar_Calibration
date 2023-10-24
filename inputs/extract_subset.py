import json

# Extract a subset of frames from a json
in_file = "LHSDesign40x4"
with open(in_file + "_output_struct.json","r") as f:
    in_data = json.loads(f.readline()) 

load = 200
for sample in in_data["Sample"]:
    # sample["Frame"] = [sample["Frame"][-1]]
    sample["Frame"] = [frame for frame in sample["Frame"] if round(frame["RFs"][2]) == load]

with open(in_file + "_output_struct_sub.json",'w') as f:
    f.write(json.dumps(in_data))