# Plots outputs from analysis as a vtk file

import numpy as np
import pandas as pd
import meshio
import os

shell_mesh = True
label = "LHSDesign40x4"
plot_basis = True # Plot of basis functions
plot_training_data_mean = False # Plot mean of training data
plot_GP_Mean = True # Plot average across posterior samples of Gaussian process predictions
plot_GP_Std = False # Plot standard deviation of posterior predictive samples
plot_GP_Post_Sam = False # Plot Gaussian process predictions across a set of posterior samples
plot_GP_Sam = False # Plot samples from Gaussian process?
file_series = False # Output as a file series for animating in paraview 
out_str = ["u", "v", "w"] # Strings for labelling output components. Used to determine number of quantities in output vector
 
file_list = []
if plot_basis:
    file_list.append("basis")
if plot_training_data_mean:
    file_list.append("training_data_mean")
if plot_GP_Mean:
    file_list.extend(["eta_mu_mu","eta_sigma_mu"])
    if plot_GP_Sam:
        file_list.append("eta_sam_mu")
if plot_GP_Std:
    file_list.append("eta_mu_sigma")
    if plot_GP_Sam:
        file_list.append("eta_sam_sigma")
if plot_GP_Post_Sam:
    file_list.extend(["eta_mu", "eta_sigma"])
    if plot_GP_Sam:
        file_list.append("eta_sam")

# Open the output files
if shell_mesh:
    node_file = "CSpar_sam_shell_mesh_nodes.csv" # Nodes of nominal input (ignores geometric uncertainty)
    element_file = "CSpar_sam_shell_mesh_elements.csv" # Element connectivity
else:
    node_file = "CSpar_sam_mesh_nodes.csv" # Nodes of nominal input (ignores geometric uncertainty)
    element_file = "CSpar_sam_mesh_elements.csv" # Element connectivity
# Add option for model mean

# Read in the element and node definitions
elements = np.loadtxt(element_file, dtype = int, delimiter = ',')
nodes = np.loadtxt(node_file, delimiter = ',')

# Remove the first column of nodes and elements which is just an index
nodes = nodes[:,1:]
elements = elements[:,1:]

# Discard the first two points which are reference points defined on the assembly, and aren't
# referenced in the connectivity file
nodes = nodes[2:,:]

# Define list of indice of the nodes which define each face of the brick (i.e. which column of the connectivity)
# Look at Jean's code - figure out how to write as solid rather than faces
face_nodes = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]]
faces = np.empty([0,4], dtype = int)
for face in face_nodes:
    elements_face = elements[:,face]
    faces = np.concatenate((faces,elements_face))

# Convert connectivity from abaqus to python indexing
faces = faces - 1

dict_list = []
output_dict = {}
for i, filestr in enumerate(file_list):
    # load in output from current filename
    output = pd.read_csv(filestr + "_" + label + ".csv")
    # If the csv contains output for multiple displacement components, separate these out first
    if len(out_str) > 0:
        labels = output.columns.values
        output = output.to_numpy()
        n_out = output.shape[0]//len(out_str)
        # Reshape output and update column labels
        output = pd.DataFrame(output.reshape(n_out,-1, order = 'F'), columns= ["_".join((label, comp)) for label in labels for comp in out_str])
        
    # Discard the first two points from the output as these are reference points
    output = output.loc[2:]
    for j, col in enumerate(output.columns):
        if  file_series:
            # If outputting a file series, create individual vtks for each sample
            if i == 0:
                dict_list.append({filestr:output[col]})
            else:
                dict_list[j][filestr] = output[col]
        else:
            # Otherwise write to a single vtk
            for col in output.columns:
                output_dict[col] = output[col]

# Get data in the correct format for meshio and write to vtk
if file_series:
    if "output_vtk" not in os.listdir(os.getcwd()):
        os.makedirs("output_vtk")
    for  i, output_dict in enumerate(dict_list):
        meshio.Mesh(points = nodes, cells = [("quad",faces)], point_data = output_dict).write(os.path.join("output_vtk","Sample_"+str(i)+".vtk"), file_format="vtk")
else:
    meshio.Mesh(points = nodes, cells = [("quad",faces)], point_data = output_dict).write("output_vtk.vtk", file_format="vtk")

