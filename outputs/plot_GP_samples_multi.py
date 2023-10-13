# Packages outputs from the calibration code, stored in different .csv files, into a vtk file

import numpy as np
import meshio

# Open the file in question
node_file = 'CSpar_sam_mesh_nodes.csv' # consider using the nodal positions of the nominal properties
element_file = 'CSpar_sam_mesh_elements.csv'

suffix = '_50x3_multi' # Identifier which is appended to all input/output strings

# Files containing output samples
# sample_file_list = ['eta_mu','delta_mu','delta_sam','eta_sam','delta_sigma','eta_sigma']
sample_file_list = ['eta_mu','eta_sam','eta_sigma']
sample_file_list = []
sample_file_list = [file_name + suffix for file_name in sample_file_list]

# Files containing average output
mean_file_list = ['eta_sam_mu','eta_mu_mu', 'eta_sigma_mu']
mean_file_list = [file_name + suffix for file_name in mean_file_list] # add suffix

# Read in the element and node definitions
elements = np.loadtxt(element_file, dtype = int, delimiter = ',')
nodes = np.loadtxt(node_file, delimiter = ',')
# first column is just an index, which I don't need here as everything is (I think) in ascending order
# with no gaps
nodes = nodes[:,1:4]
nodes = nodes[2:,:]

# Define list of indice of the nodes which define each face of the brick (i.e. which column of the connectivity)
# Note that first column of elements is just the element number, hence indexing starts at 1
face_nodes = [[1,2,3,4],[5,6,7,8],[1,2,6,5],[2,3,7,6],[3,4,8,7],[4,1,5,8]]
faces = np.empty([0,4], dtype = int)
for face in face_nodes:
    elements_face = elements[:,face]
    faces = np.concatenate((faces,elements_face))
    

# Re-writing from abaqus to python indexing for meshio
elements_face = elements_face - 1

# Loop over different csv files containing output data, and create a dictionary
#for file_name in sample_file_list:
#    output_dict = {}
#    output = np.loadtxt(file_name + '.csv', delimiter = ',', skiprows = 1)
#    # Loop over samples of the current output and add to dictionary
#    for i in range(output.shape[1]):
#        output_dict[file_name + str(i+1)] = output[2:,]
#    meshio.Mesh(points = nodes, cells = [("quad",elements_face)], point_data = output_dict).write(file_name + "_samples.vtk", file_format="vtk")

# For files containing mean output there is no need to append the field name with a sample number
output_dict = {}
for file_name in mean_file_list:
    output = np.loadtxt(file_name + '.csv', delimiter = ',', skiprows = 1)
    output_w = output[0:(output.shape[0]//2),:]
    output_u = output[output.shape[0]//2:,:]
    for i in range(output.shape[1]):
        output_dict[file_name + "_u" + str(i+1)] = output_u[2:,i]
        output_dict[file_name + "_w" + str(i+1)] = output_w[2:,i]

# Get data in the correct format for meshio
meshio.Mesh(points = nodes, cells = [("quad",elements_face)], point_data = output_dict).write("GP_mean" + suffix + ".vtk", file_format="vtk")