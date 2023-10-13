# Quick sanity check of interpolation

import numpy as np
import meshio

elements = np.array([115304, 115387, 115388, 115305])
elements = elements[np.newaxis,:]

# Read node definitions
node_file = 'CSpar_sam_mesh_nodes.csv' # consider using the nodal positions of the nominal properties
nodes = np.loadtxt(node_file, delimiter = ',')
# first column is just an index, which I don't need here as everything is (I think) in ascending order
# with no gaps
nodes = nodes[:,1:4]

output_file = 'model_mean.csv'
# Read in basis data
output = np.loadtxt(output_file, delimiter = ',', skiprows = 1)

# Re-writing from abaqus to python indexing for meshio
elements = elements - 1

# Create dictionary using basis data
output_dict = {"model_mean": output}

# Get data in the correct format for meshio
meshio.Mesh(points = nodes, cells = [("quad",elements)], point_data = output_dict).write("test.vtk", file_format="vtk")