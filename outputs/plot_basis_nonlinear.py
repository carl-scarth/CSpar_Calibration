# Plots the basis functions of finite element model output
# Write a general function for plotting all outputs and place in parent directory??

import numpy as np
import meshio
import json
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

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
    frame["Training_Data_Mean"] = np.array(frame["Training_Data_Mean"], dtype="float")

# Find the row with maximum mean value (across training data) for extracting relevant output
max_ind = np.argmax(np.abs(in_dict["Frame"][-1]["Training_Data_Mean"]))
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
# displacement at the reference point. Also store output for the row
# with the maximum absolute value in the final frame
basis_RP = np.empty((n_frames,n_bases))
basis_max = np.empty((n_frames,n_bases))
mean_RP = np.empty(n_frames)
mean_max = np.empty(n_frames)
for i, frame in enumerate(in_dict["Frame"]):
    basis_RP[i,:] = frame["Bases"][1,:]
    mean_RP[i] = frame["Training_Data_Mean"][1]
    basis_max[i,:] = frame["Bases"][max_ind,:]
    mean_max[i] = frame["Training_Data_Mean"][max_ind]
    # Now delete the reference point
    frame["Bases"] = frame["Bases"][2:,:]
    frame["Training_Data_Mean"] = frame["Training_Data_Mean"][2:]

# Define list of indice of the nodes which define each face of the brick (i.e. which column of the connectivity)
# Look at Jean's code - figure out how to write as solid rather than faces
face_nodes = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]]
faces = np.empty([0,4], dtype = int)
for face in face_nodes:
    elements_face = elements[:,face]
    faces = np.concatenate((faces,elements_face))

# Convert connectivity from abaqus to python indexing
faces = faces - 1

# Create a new directory for the vtk files, if one does not exist already
if not(os.path.isdir("basis_nonlinear_" + infile)):
    os.mkdir("basis_nonlinear_" + infile)
    
# Loop over each frame, then write a separate frame which paraview can 
# interpret as a file series
for i, frame in enumerate(in_dict["Frame"]):
    # Create a dictionary using output data for each basis
    basis_dict = {}
    for j in range(n_bases):
        basis_dict["basis_" + str(j+1)] = frame["Bases"][:,j]
        # Get data in the correct format for meshio
    
    basis_dict["Training_Data_Mean"] = frame["Training_Data_Mean"]
    meshio.Mesh(points = nodes, cells = [("quad",faces)], point_data = basis_dict).write("basis_nonlinear_" + infile + "\\frame_" + str(i) + ".vtk", file_format="vtk")

# Finally, plot the basis function at the reference point
rcParams.update({'figure.figsize' : (8,6),
                'font.size' : 14,
                'font.family' : 'serif',
                'figure.titlesize' : 16,
                'axes.labelsize': 16,
                'xtick.labelsize': 14,
                'ytick.labelsize': 14,
                'legend.fontsize': 14})

fig, ax = plt.subplots()
for base in basis_RP.T:
    ax.plot(range(n_frames), base)

ax.set_xlabel("Increment")
ax.set_ylabel("K_i")
ax.set_title("Basis values at reference point")

fig2, ax2 = plt.subplots()
ax2.plot(range(n_frames), mean_RP)
ax2.set_xlabel("Increment")
ax2.set_ylabel("Training Data Mean")
ax2.set_title("Training data mean at reference point")

# Plot output at location of maximum absolute mean displacement (across training data)
# in final frame (useful for vertical displacement)
fig3, ax3 = plt.subplots()
for base in basis_max.T:
    ax3.plot(range(n_frames),base)
ax3.set_xlabel("Increment")
ax3.set_ylabel("K_i")
ax3.set_title("Basis at location of maximum absolute displacement across training data")

fig4, ax4 = plt.subplots()
ax4.plot(range(n_frames), mean_max)
ax4.set_xlabel("Increment")
ax4.set_ylabel("Training Data Mean")
ax4.set_title("Mean displacement at location of maximum absolute value across training data")
plt.show()