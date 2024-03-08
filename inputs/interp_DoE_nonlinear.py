# Interpolates model outputs for Design of Experiments to DIC coordinates

import pandas as pd
import numpy as np
import json
import sys
sys.path.append("..\\source")

from interpolate_data import intp_nodes_to_cloud

if __name__ == "__main__":
    # Load inputs from file
    in_file = "LHSDesign50x5_1"
    conn_file = "nominal_shell_mesh_outer_surface_elements" # Connectivity
    # exp_data_file = "Selected_Points_CS02P"                 # Experimental data
    exp_data_file = "mid_point"                 # Experimental data
    DoE_file = in_file + "_output_struct_disp.json"         # Model output
    conn = pd.read_csv(conn_file+".csv").to_numpy(dtype=int)
    conn = conn - 1 # Convert to Python indexing from Abaqus
    output_mean = True # Only output mean
    print(conn)
    with open (DoE_file,'r') as f:
        sample_dict = json.loads(f.readline())
    
    exp_data = pd.read_csv(exp_data_file+".csv")

    # Extract element index and natural coordinates
    el_ind = exp_data["Element"].to_numpy(dtype=int)
    gh = exp_data[["g","h"]].to_numpy()

    if output_mean:
        interp_outputs = pd.DataFrame()
    else:
        interp_outputs = pd.DataFrame(columns = ["index", 'x_proj','y_proj','z_proj','x_rot','y_rot','z_rot',"u_interp", "v_interp", "w_interp", "Sample", "Increment", "Force"])
    for i,sample in enumerate(sample_dict["Sample"]):
        if output_mean:
            newcol = np.empty([len(sample["Frame"]),1])
            force = np.empty([len(sample["Frame"]),1])
            crosshead = np.empty([len(sample["Frame"]),1])
        for j,frame in enumerate(sample["Frame"]):
            interp_outputs_i = intp_nodes_to_cloud(el_ind, gh, np.array(frame["Displacements"]), conn, skip_nodes = 2)
            interp_outputs_i = pd.DataFrame(interp_outputs_i, columns = ["u_interp", "v_interp", "w_interp"])
            if output_mean:
                newcol[j,0] = interp_outputs_i["w_interp"].mean()
                force[j,0] = frame["RFs"][2]
                crosshead[j,0] = frame["Displacements"][1][2]
            else:
                interp_outputs_i = pd.concat([exp_data[['index','x_proj','y_proj','z_proj','x_rot','y_rot','z_rot']],interp_outputs_i],axis=1)
                interp_outputs_i["Sample"] = i
                interp_outputs_i["Increment"] = j
                interp_outputs_i["Force"] = frame["RFs"][2]
                interp_outputs = pd.concat([interp_outputs, interp_outputs_i],axis=0,ignore_index=True)

        if output_mean:
            interp_outputs = pd.concat([interp_outputs, pd.DataFrame(newcol, columns = ["disp_" + str(i)]), pd.DataFrame(force, columns = ["force_"+str(i)]), pd.DataFrame(crosshead, columns = ["crosshead_"+str(i)])], axis=1)
    
    interp_outputs.to_csv(in_file+"_interp_mid.csv", sep=",", index=False)