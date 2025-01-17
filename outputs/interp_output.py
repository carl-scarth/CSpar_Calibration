# Interpolates model outputs for Design of Experiments to DIC coordinates
import numpy as np
import pandas as pd
import json
import sys
import os
sys.path.append("..\\source")
from interpolate_data import intp_nodes_to_cloud

if __name__ == "__main__":

    # Set up for output across multiple frames, input via json
    infile = "LHSDesign100x8"
    # output_file = "gp_predictions_nonlinear_" + infile + ".json"
    output_file = "E:Calibration_outputs_for_paper\\gp_predictions_nonlinear_" + infile + "obs_err5e-2.json"
    # conn_file = "nominal_shell_mesh_outer_surface_elements" # Connectivity
    # Load inputs from file
    conn_file = "new_spar_mesh_outer_surface_elements"
    multiple_datasets = True
    exp_data_file = "..\\inputs\\Interpolated_DIC_multistep_150kN_combined" # csv containing experimental data, including increment index and point index
    
    # Read in prediction data
    with open(output_file, "r") as f:
        # Load in string from file
        in_dict = json.loads(f.readline())
    # Read in connectivity
    conn = pd.read_csv(conn_file+".csv").to_numpy(dtype=int)
    conn = conn - 1 # Convert to Python indexing from Abaqus
    # Read in experimental data
    exp_data = pd.read_csv(exp_data_file+".csv")
    exp_data.columns.values[0] = "point_ind"
    print(exp_data)

    # Convert increment to python indexing
    exp_data.Increment = exp_data.Increment - 1

    # Get number of predictions and number of frames.
    n_frames = len(in_dict["Frame"])
    # Get labels of outputs and number of nodes
    if "Posterior_Sample" in list(in_dict["Frame"][0].values())[0]:
        out_str = list(list(in_dict["Frame"][0].values())[0]["Posterior_Sample"][0].keys())
        n_nodes = len(list(list(in_dict["Frame"][0].values())[0]["Posterior_Sample"][0].values())[0])
    else:
        out_str = list(list(in_dict["Frame"][0].values())[0].keys())
        n_nodes = len(list(list(in_dict["Frame"][0].values())[0].values())[0])
    d_y = len(out_str) # Number of y components
    
    # Loop over the entries of in_dict and interpolate the data to a point cloud
    #if "interp_gp_output" not in os.listdir(os.getcwd()):
    #    os.mkdir("interp_gp_output")
    if "interp_gp_output" not in os.listdir("E:Calibration_outputs_for_paper"):
        os.mkdir("E:Calibration_outputs_for_paper\\interp_gp_output")

    for i, frame in enumerate(in_dict["Frame"]):
        print(i)
        # Extract experimental data for the current frame
        exp_data_i = exp_data[exp_data.Increment == i]

        # Extract element index and natural coordinates
        el_ind = exp_data_i["Element"].to_numpy(dtype=int)
    
        # gh = exp_data_i[["h","r"]].to_numpy()
        gh = exp_data_i[["g","h"]].to_numpy()
        out_frame = exp_data_i[["point_ind","x_proj","y_proj","z_proj","Camera_pair","Load","Compressive Force"]]
        for QoI, predictions in frame.items():
            print(QoI)
            # Does output have multiple posterior samples?
            if "Posterior_Sample" in predictions:
                n_post_sam = len(predictions["Posterior_Sample"])
                # output = np.empty((len(predictions["Posterior_Sample"][0]),n_post_sam*d_y), dtype=float)
                output = np.empty((n_nodes,0), dtype = float)
                output_names = []
                for j, post_sam in enumerate(predictions["Posterior_Sample"]):
                    # Convert from list to numpy
                    output = np.concatenate((output, np.array([component for component in post_sam.values()],ndmin=2, dtype="float").T), axis=1)
                    output_names.extend(["_".join((QoI, coord, str(j))) for coord in post_sam.keys()])

                interp_output = intp_nodes_to_cloud(el_ind, gh, output, conn, GH = [], skip_nodes=0)
                interp_output = pd.DataFrame(interp_output, columns=output_names)
            else:
                # Otherwise the prediction is an average across the posterior
                output = np.array([component for component in predictions.values()],ndmin=2).T
                interp_output = intp_nodes_to_cloud(el_ind, gh, output, conn, GH = [], skip_nodes=0)
                interp_output = pd.DataFrame(interp_output, columns=["_".join((QoI,coord)) for coord in predictions.keys()])
                interp_output.reset_index(drop=True, inplace=False)

            # Problem with index here
            out_frame.reset_index(drop=True, inplace=True)
            out_frame = pd.concat((out_frame,interp_output), axis=1)

        # Write to csv
        #out_frame.to_csv("interp_gp_output\\Frame_"+str(i)+".csv",sep=",",index=False)
        out_frame.to_csv("E:Calibration_outputs_for_paper\\interp_gp_output\\Frame_"+str(i)+".csv",sep=",",index=False)
        # Add residuals if need be, though leave till next time - doesn't matter for this first run and
        # may not downsample in same way in future