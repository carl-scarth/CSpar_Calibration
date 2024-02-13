import json

in_file = "inputs\\LHSDesign50x3_3_output_struct.json"
#in_file = "outputs\\gp_predictions_nonlinear_LHSDesign40x4.json"
with open (in_file,'r') as f:
    sample_dict = json.loads(f.readline())

# Reduce frame size when convergence issues
for sample in sample_dict["Sample"]:
    # print(sample["Frame"][-1]["Increment"])
    print(len(sample["Frame"]))
    if len(sample["Frame"]) > 16:
        sample["Frame"] = sample["Frame"][0:16]

with open(in_file,'w') as f:
    f.write(json.dumps(sample_dict))