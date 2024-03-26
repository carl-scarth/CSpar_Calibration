# Interpolates model outputs for Design of Experiments to DIC coordinates

import pandas as pd
import sys
sys.path.append("..\\source")

from interpolate_data import intp_nodes_to_cloud

if __name__ == "__main__":
    # Load inputs from file
    conn_file = "nominal_shell_mesh_outer_surface_elements" # Connectivity
    exp_data_file = "Selected_Points_CS02P"                # Experimental data
    DoE_file = "nominal_inputs_new_spar"                  # Model output
    conn = pd.read_csv(conn_file+".csv").to_numpy(dtype=int)
    conn = conn - 1 # Convert to Python indexing from Abaqus
    
    model_out = pd.read_csv(DoE_file+".csv")
    exp_data = pd.read_csv(exp_data_file+".csv")

    # Extract element index and natural coordinates
    el_ind = exp_data["Element"].to_numpy(dtype=int)
    gh = exp_data[["h","r"]].to_numpy()

    interp_outputs = intp_nodes_to_cloud(el_ind, gh, model_out.to_numpy(), conn, GH = [], skip_nodes = 2)
    interp_outputs = pd.DataFrame(interp_outputs, columns=model_out.columns.values)
    interp_outputs.to_csv(DoE_file+"_interp.csv", sep=",", index=False)

    # concatenate with coordinates for plotting
    interp_outputs = pd.concat((exp_data[["x_0_rot","y_0_rot","z_0_rot"]],interp_outputs),axis=1)
    interp_outputs.to_csv(DoE_file+"_interp_plot.csv", sep=",", index=False)