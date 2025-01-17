import pandas as pd
import os

model_file = "E:Calibration_outputs_for_paper\\nominal_inputs_translator_interp_150kN"
out_dir = "E:Calibration_outputs_for_paper\\"

#prior_sample = pd.read_csv(model_file+".csv")
# Convert to abaqus notation on increments for consistency with experimental data
camera_pairs = ["LC","RC"]

for pair in camera_pairs:
    if "prior_residuals_" + pair not in os.listdir(out_dir):
        os.mkdir(os.path.join(out_dir,"prior_residuals_" + pair))

prior_sample = {pair : pd.read_csv("_".join((model_file,pair))+".csv") for pair in camera_pairs}
for sample in prior_sample.values():
    sample["Increment"] = sample["Increment"] + 1

experimental_file = "E:Calibration_outputs_for_paper\\Interpolated_DIC_multistep_150kN"
disp_str = ["u","v","w"]
suffix_model = "interp"
suffix_data = "rot"
experimental_data = {pair : pd.read_csv("_".join((experimental_file,pair))+".csv") for pair in camera_pairs}
increments = prior_sample[camera_pairs[0]]["Increment"].unique()

for data in experimental_data.values():
    data.columns.values[0] = "point_ind"
#for data in prior_sample.values():
    #data.columns.values[0] = "point_ind"

for i, increment in enumerate(increments):
    print(i)
    for pair, data in experimental_data.items():
        prior_pair = prior_sample[pair]
        model_i = prior_pair[prior_pair["Increment"] == increment].set_index("point_ind")
        data_i = data[data["Increment"] == increment].set_index("point_ind")
        for comp in disp_str:
            res = model_i["_".join((comp, suffix_model))] - data_i["_".join((comp, suffix_data))]
            model_i[comp+"_res"] = res
            model_i[comp+"_absres"] = res.abs()
            model_i[comp+"_pererr"] = model_i[comp+"_absres"]/data_i["_".join((comp, suffix_data))].abs()*100.0
        model_i.to_csv(os.path.join(out_dir,"prior_residuals_" + pair,"frame_"+str(i)+".csv"), sep=",", index=False)
