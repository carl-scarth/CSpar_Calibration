# Extracts the elements from the surface ply based upon the element-set given in the mesh .inp file
import numpy as np
import meshio

# Open the file in question
node_file = 'CSpar_sam_mesh_nodes.csv' # consider using the nodal positions of the nominal properties
element_file = 'CSpar_sam_mesh_elements.csv'
output_file = 'model_mean_50x3_multi.csv'
output_str = output_file.strip(".csv")


# Read in the element and node definitions
elements = np.loadtxt(element_file, dtype = int, delimiter = ',')
nodes = np.loadtxt(node_file, delimiter = ',')
# first column is just an index, which I don't need here as everything is (I think) in ascending order
# with no gaps
nodes = nodes[:,1:4]
nodes = nodes[2:,]

# Read in basis data
output = np.loadtxt(output_file, delimiter = ',', skiprows = 1, ndmin = 2)

# Define list of indice of the nodes which define each face of the brick (i.e. which column of the connectivity)
# Note that first column of elements is just the element number, hence indexing starts at 1
face_nodes = [[1,2,3,4],[5,6,7,8],[1,2,6,5],[2,3,7,6],[3,4,8,7],[4,1,5,8]]
faces = np.empty([0,4], dtype = int)
for face in face_nodes:
    elements_face = elements[:,face]
    faces = np.concatenate((faces,elements_face))

# Re-writing from abaqus to python indexing for meshio
elements_face = elements_face - 1

output_w = output[0:(output.shape[0]//2),:]
output_u = output[output.shape[0]//2:,:]
# Create dictionary using output data
output_dict = {output_str + "_u" : output_u[2:,]}
output_dict[output_str + "_w"] = output_w[2:,]

# Get data in the correct format for meshio
meshio.Mesh(points = nodes, cells = [("quad",elements_face)], point_data = output_dict).write(output_str + ".vtk", file_format="vtk")