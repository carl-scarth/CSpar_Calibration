# Plots the basis functions of finite element model output
# Write a general function for plotting all outputs and place in parent directory??

import numpy as np
import meshio
import json
import matplotlib.pyplot as plt

shell_mesh = True # Is the mesh comprised of continuum shells?
infile = "LHSDesign40x4"

# Open the output files
if shell_mesh:
    node_file = "CSpar_sam_shell_mesh_nodes.csv" # Nodes of nominal input (ignores geometric uncertainty)
    element_file = "CSpar_sam_shell_mesh_elements.csv" # Element connectivity
else:
    node_file = "CSpar_sam_mesh_nodes.csv" # Nodes of nominal input (ignores geometric uncertainty)
    element_file = "CSpar_sam_mesh_elements.csv" # Element connectivity

basis_file = "basis_nonlinear_" + infile + ".json"

# Read in the element and node definitions
elements = np.loadtxt(element_file, dtype = int, delimiter = ',')
nodes = np.loadtxt(node_file, delimiter = ',')
# Read in basis data
with open(basis_file, "r") as f:
    # Load in string from file
    in_dict = json.loads(f.readline())
    
# Convert Basis data to numpy
for frame in in_dict["Frame"]:
    frame["Bases"] = np.array(frame["Bases"])

# Extract number of frames and number of bases
n_frames = len(in_dict["Frame"])
n_bases = in_dict["Frame"][0]["Bases"].shape[1]

# Remove the first column of nodes and elements which is just an index
nodes = nodes[:,1:]
elements = elements[:,1:]

# Discard the first two points which are reference points defined on the assembly, and aren't
# referenced in the connectivity file
nodes = nodes[2:,:]
# Extract the second row from each basis, which gives corresponds to 
# displacement at the reference point
basis_RP = np.empty((n_frames,n_bases))
for i, frame in enumerate(in_dict["Frame"]):
    basis_RP[i,:] = frame["Bases"][1,:]
    # Now delete the reference points
    frame["Bases"] = frame["Bases"][2:,:]

# Define list of indice of the nodes which define each face of the brick (i.e. which column of the connectivity)
# Look at Jean's code - figure out how to write as solid rather than faces
face_nodes = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]]
faces = np.empty([0,4], dtype = int)
for face in face_nodes:
    elements_face = elements[:,face]
    faces = np.concatenate((faces,elements_face))

# Convert connectivity from abaqus to python indexing
faces = faces - 1

# Loop over each frame, then write a separate frame which paraview can 
# interpret as a file series
for i, frame in enumerate(in_dict["Frame"]):
    # Create a dictionary using output data for each basis
    basis_dict = {}
    for j in range(n_bases):
        basis_dict["basis_" + str(j+1)] = frame["Bases"][:,j]
        # Get data in the correct format for meshio
        meshio.Mesh(points = nodes, cells = [("quad",faces)], point_data = basis_dict).write("basis_nonlinear" + infile + "_" + str(i) + ".vtk", file_format="vtk")

# Finally, plot the basis function at the reference point
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)
for base in basis_RP.T:
    ax.plot(range(n_frames), base)

label_font = {'family': 'serif', 'size': 16,}
ax.set_ylabel("Increment", fontdict = label_font)
ax.set_xlabel("K_i", fontdict = label_font)
plt.show()
