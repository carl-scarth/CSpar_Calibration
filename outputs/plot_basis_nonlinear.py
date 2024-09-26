# Plots the basis functions of finite element model output
# Write a general function for plotting all outputs and place in parent directory??

import numpy as np
import meshio
import json
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

shell_mesh = True # Is the mesh comprised of continuum shells?
infile = "LHSDesign100x8_6"
new_spar = True # Are we considering the new spar geometry?
flex_support = True # Are we using the model with flexible supports?
flex_ground = False # Are we using the model with flexible ground?

# Open the output files
if shell_mesh:
    if new_spar:
        node_file = "new_spar_mesh_nodes.csv"
        element_file = "new_spar_mesh_elements.csv"
    else:
        node_file = "CSpar_sam_shell_mesh_nodes_new.csv" # Nodes of nominal input (ignores geometric uncertainty)
        element_file = "CSpar_sam_shell_mesh_elements.csv" # Element connectivity
else:
    node_file = "CSpar_sam_mesh_nodes.csv" # Nodes of nominal input (ignores geometric uncertainty)
    element_file = "CSpar_sam_mesh_elements.csv" # Element connectivity

basis_file = "basis_nonlinear_" + infile + "_uw150z.json"

# Read in the element and node definitions
elements = np.loadtxt(element_file, dtype = int, delimiter = ',')
nodes = np.loadtxt(node_file, delimiter = ',')
# Read in basis data
with open(basis_file, "r") as f:
    # Load in string from file
    in_dict = json.loads(f.readline())
  
# Convert Basis and mean data to numpy
for frame in in_dict["Frame"]:
    frame["Bases"] = {label : np.array(value, dtype="float") for label, value in frame["Bases"].items()}
    frame["Training_Data_Mean"] = {label : np.array(value, dtype="float") for label, value in frame["Training_Data_Mean"].items()}

# Get maximum mean displacement for plotting/exporting


max_ind = {label : np.argmax(np.abs(value)) for label, value in in_dict["Frame"][-1]["Training_Data_Mean"].items()}
max_ind = {'u': 4164, 'v': 182, 'w': 19} # New spars (simply supported) (comes out as same value as above)
if flex_support or flex_ground:
    max_ind = {key : value - 2 for key, value in max_ind.items()} # (don't have the reference points in the new model)

# Extract number of frames, number of bases, and number of outputs
n_frames = len(in_dict["Frame"])
out_str = list(in_dict["Frame"][0]["Bases"].keys())
n_bases = in_dict["Frame"][0]["Bases"][out_str[0]].shape[1]
n_out = len(out_str)

# Remove the first column of nodes and elements which is just an index
nodes = nodes[:,1:]
elements = elements[:,1:]
# Discard the first two points which are reference points defined on the assembly, and aren't
# referenced in the connectivity file
nodes = nodes[2:,:]

# If longittudinal displacement w is in the output, extract the second row from
# each basis, which gives corresponds to displacement at the reference point. 
# This is zero for other components
# Also store output for the row with the maximum absolute value in the final frame
if "w" in out_str:
    basis_RP = np.empty((n_frames,n_bases)) # Only one value as displacement components u and v are constrained to zero
    mean_RP = np.empty(n_frames)
basis_max = {label : np.empty((n_frames,n_bases)) for label in out_str}
mean_max = {label : np.empty(n_frames) for label in out_str}
for i, frame in enumerate(in_dict["Frame"]):
    if "w" in out_str:
        if flex_support or flex_ground:
            basis_RP[i,:] = frame["Bases"]["w"][-1,:]
            mean_RP[i] = frame["Training_Data_Mean"]["w"][-1]
        else:
            basis_RP[i,:] = frame["Bases"]["w"][1,:]
            mean_RP[i] = frame["Training_Data_Mean"]["w"][1]
    
    for label in out_str:
        basis_max[label][i,:] = frame["Bases"][label][max_ind[label],:]
        mean_max[label][i] = frame["Training_Data_Mean"][label][max_ind[label]]
        # Now delete the reference point
        if flex_ground:
            frame["Bases"][label] = frame["Bases"][label][:-5,:]
            frame["Training_Data_Mean"][label] = frame["Training_Data_Mean"][label][:-5]    
        elif flex_support:
            frame["Bases"][label] = frame["Bases"][label][:-4,:]
            frame["Training_Data_Mean"][label] = frame["Training_Data_Mean"][label][:-4]
        else:
            frame["Bases"][label] = frame["Bases"][label][2:,:]
            frame["Training_Data_Mean"][label] = frame["Training_Data_Mean"][label][2:]
        
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
    for label in out_str:
        for j in range(n_bases):
            basis_dict["_".join(("basis", label, str(j+1)))] = frame["Bases"][label][:,j]
        # Get data in the correct format for meshio
        basis_dict["Training_Data_Mean_" + label] = frame["Training_Data_Mean"][label]
    
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

if "w" in out_str:
    fig, ax = plt.subplots()
    for base in basis_RP.T:
        ax.plot(range(n_frames), base)
    ax.set_xlabel("Increment")
    ax.set_ylabel("K_i")
    ax.set_title("Basis values in w at reference point")

    fig2, ax2 = plt.subplots()
    ax2.plot(range(n_frames), mean_RP)
    ax2.set_xlabel("Increment")
    ax2.set_ylabel("Training Data Mean w")
    ax2.set_title("Training data mean w at reference point")

# Plot output at location of maximum absolute mean displacement (across training data)
# in final frame (useful for vertical displacement)
fig3, axes3 = plt.subplots(1,n_out)#, sharey=True)
if n_out == 1:
    axes3 = [axes3]
for i, ax3 in enumerate(axes3):
    for j, base in enumerate(basis_max[out_str[i]].T):
        ax3.plot(range(n_frames),base, label="K_"+str(j+1))
        ax3.set_title(out_str[i])
        ax3.set_xlabel("Increment")
        if i == 0:
            ax3.set_ylabel("K_i")
        ax3.legend()
fig3.suptitle("Basis at location of maximum absolute displacement across training data")

fig4, ax4 = plt.subplots()
for label, value in mean_max.items():
    ax4.plot(range(n_frames), value, label = label)
ax4.set_xlabel("Increment")
ax4.set_ylabel("Training Data Mean")
ax4.set_title("Mean displacement at location of maximum absolute value across training data")
ax4.legend()
plt.show()