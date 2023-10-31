import json

# in_file = "inputs\\LHSDesign40x4_downsam_2.json"
in_file = "outputs\\gp_predictions_nonlinear_LHSDesign40x4.json"
with open (in_file,'r') as f:
    sample_dict = json.loads(f.readline())
