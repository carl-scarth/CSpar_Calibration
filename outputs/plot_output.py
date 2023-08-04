# Plots outputs from analysis as a vtk file

import numpy as np
import pandas as pd
import meshio

shell_mesh = True
label = "LHSDesign40x4"
plot_basis = False # Plot of basis functions
plot_training_data_mean = False # Plot mean of training data
plot_GP_Mean = True # Plot average across posterior samples of Gaussian process predictions
plot_GP_Post_Sam = False # Plot Gaussian process predictions across a set of posterior samples of hyperparameters
plot_GP_Sam = True # Plot samples from Gaussian process?
 
file_list = []
if plot_basis:
    file_list.append("basis")
if plot_training_data_mean:
    file_list.append("training_data_mean")
if plot_GP_Mean:
    file_list.extend(["eta_mu_mu","eta_sigma_mu"])
    if plot_GP_Sam:
        file_list.append("eta_sam_mu")
if plot_GP_Post_Sam:
    file_list.extend(["eta_mu", "eta_sigma"])
    if plot_GP_Sam:
        file_list.append("eta_sam")

# Files containing output samples
# sample_file_list = ['eta_mu','delta_mu','delta_sam','eta_sam','delta_sigma','eta_sigma']

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

# Loop over different csv files containing output data, and create a dictionary
# No longer used. Delete if I find it definitely isn't necessary
#for file_name in sample_file_list:
#    output_dict = {}
#    output = np.loadtxt(file_name + '.csv', delimiter = ',', skiprows = 1)
#    # Loop over samples of the current output and add to dictionary
#    #for i in range(output.shape[1]):
#    for i in range(10):
#        print(i)
#        output_dict[file_name + str(i+1)] = output[2:,i]
#    meshio.Mesh(points = nodes, cells = [("quad",elements_face)], point_data = output_dict).write(file_name + "_samples.vtk", file_format="vtk")

# Dictionary containing name of output field, and the field values
output_dict = {}
# Loop over all files for which output is required and populate the output dictionary
for filestr in file_list:
    # load in output from current filename
    # output = np.loadtxt(filestr + "_" + label + ".csv", delimiter = ",", skiprows = 1) # legacy code
    # output = output[2:,:]
    # for i in range(output.shape[1]):
    #    output_dict[filestr + "_" + str(i+1)] = output[:,i]

    output = pd.read_csv(filestr + "_" + label + ".csv")

    # Discard the first two points from the output as these are reference points
    output = output.loc[2:]
    # Loop over entries of the output file and store in the output dictionary
    for col in output.columns:
        output_dict[col] = output[col]

# Get data in the correct format for meshio
meshio.Mesh(points = nodes, cells = [("quad",faces)], point_data = output_dict).write("output_vtk.vtk", file_format="vtk")