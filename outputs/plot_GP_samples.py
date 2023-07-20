# Plots outputs from analysis as a vtk file
# Try to set up for any output



import numpy as np
import meshio

shell_mesh = True
label = "LHSDesign40x4"
plot_basis = False
plot_GP_Mean = True

file_list = []
if plot_basis:
    file_list.append("basis")
if plot_GP_Mean:
    # Sort format after done with R
    file_list.extend(["eta_mu_mu","eta_sigma_mu"])

# GOOD PRACTICE WOULD BE TO USE PANDAS, AND WRITE HEADER. THEN I COULD USE FOR INPUT AND OUTPUT
# SEE PLOT DOE IN INPUTS FOR EXAMPLE

# Files containing output samples
# sample_file_list = ['eta_mu','delta_mu','delta_sam','eta_sam','delta_sigma','eta_sigma']
# sample_file_list = ['eta_mu','eta_sam','eta_sigma']
# sample_file_list = []
# Files containing average output
# mean_file_list = ['eta_mu_mu', 'eta_sigma_mu','eta_sam_mu']

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
    output = np.loadtxt(filestr + "_" + label + ".csv", delimiter = ",", skiprows = 1)
    # Discard the first two points from the output as these are reference points
    output = output[2:,:]
    # Loop over entries of the output file and store in the output dictionary
    for i in range(output.shape[1]):
        output_dict[filestr + "_" + str(i+1)] = output[:,i]

# Get data in the correct format for meshio
meshio.Mesh(points = nodes, cells = [("quad",faces)], point_data = output_dict).write("output_vtk.vtk", file_format="vtk")