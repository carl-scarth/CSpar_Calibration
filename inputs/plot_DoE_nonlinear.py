# Plots the basis functions of finite element model output
# Write a general function for plotting all outputs and place in parent directory??

import numpy as np
import meshio
import json
import matplotlib.pyplot as plt
import os

shell_mesh = True # Is the mesh comprised of continuum shells?
infile = "LHSDesign40x4_1"

# Open the output files
if shell_mesh:
    node_file = "CSpar_sam_shell_mesh_nodes.csv" # Nodes of nominal input (ignores geometric uncertainty)
    element_file = "CSpar_sam_shell_mesh_elements.csv" # Element connectivity
else:
    node_file = "CSpar_sam_mesh_nodes.csv" # Nodes of nominal input (ignores geometric uncertainty)
    element_file = "CSpar_sam_mesh_elements.csv" # Element connectivity

output_file = infile + "_output_struct.json"

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

# Read in basis data
with open(output_file, "r") as f:
    # Load in string from file
    in_dict = json.loads(f.readline())  

# Get number of predictions and number of frames
n_sam = len(in_dict["Sample"]) # Number of predictions
n_frames = len(in_dict["Sample"][0]["Frame"]) # Number of frames
increments = [frame["Increment"] for frame in in_dict["Sample"][0]["Frame"]]
keep_ind = np.linspace(0,n_frames-1, 17, dtype="int")
RFs = [frame["RFs"][2] for frame in in_dict["Sample"][0]["Frame"]]
max_ind = {"u" : 4164, "v" : 207, "w" : 1299} # Index of node containing maximum absolute value across training data
# col_ind = {"u" : 0, "v" : 1, "w" : 2} # Column index for each displacement component
output_RP = np.empty((n_sam, n_frames)) # Reference point displacement
output_max = {key : np.empty((n_sam, n_frames)) for key in max_ind.keys()} # Maximum displacement
frame_dict = [{} for i in range(n_frames)] # List of dictionaries containing output for each frame
for i, sample in enumerate(in_dict["Sample"]):
    for j, frame in enumerate(sample["Frame"]):
        # Extract displacements and convert to numpy
        displacements = np.array(frame["Displacements"], dtype = "float")
        # Store RP displacement
        output_RP[i,j] = displacements[1,2]
        # for key, value in output_max.items():
          #  value[i,j] = displacements[max_ind[key], col_ind[key]]
            #output_max[i,j] = displacements[max_ind,0]
        for k, QoI in enumerate(["u","v","w"]):
            output_max[QoI][i,j] = displacements[max_ind[QoI],k]
            # Dictionary containing output
            if j in keep_ind:
                frame_dict[j][QoI+"_"+str(i)] = displacements[2:,k]
                
frame_dict = [frame_dict[i] for i in keep_ind] # Only retain those at specified increments
# Create a new directory for the vtk files, if one does not exist already
if not(os.path.isdir("displacements_" + infile)):
    os.mkdir("displacements_" + infile)
    
# Loop over each frame, then write a separate frame which paraview can 
# interpret as a file series
for i, frame in enumerate(frame_dict):
    meshio.Mesh(points = nodes, cells = [("quad",faces)], point_data = frame).write("displacements_" + infile + "\\frame_" + str(i) + ".vtk", file_format="vtk")
    
    
# Finally, plot the displacements at the reference point

fig, ax = plt.subplots(figsize=(10,4))
for sample in output_RP:
    ax.plot(-sample, [250*inc for inc in increments],linewidth=2.0)

ax.set_ylabel("Force (kN)")
ax.set_xlabel("RP Displacement (mm)")


fig2, axs2 = plt.subplots(1, 3, figsize=(10,4))
for ax2, key in zip(axs2, ["u", "v", "w"]):
    for sample in output_max[key]:
        ax2.plot(-sample, [250*inc for inc in increments],linewidth=2.0)
    ax2.set_ylabel("Force (kN)")
    ax2.set_xlabel("Displacement (mm)")
    ax2.set_title(key+" at node "+str(max_ind[key]))
fig2.suptitle("Displacements at location of maximum mean value across training data")
plt.tight_layout(pad = 1.0)
plt.show()

# Write to a csv file
np.savetxt(infile+"_RP_displacements.csv", output_RP, delimiter=',')
for key in ["u","v","w"]:
    np.savetxt("_".join((infile,key,"max_displacements.csv")), output_max[key], delimiter=',')
