# Plots the basis functions of finite element model output
# Write a general function for plotting all outputs and place in parent directory??

import numpy as np
import meshio
import json
import matplotlib.pyplot as plt
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

output_file = "gp_predictions_nonlinear_" + infile + ".json"

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

# Get number of predictions and number of frames. If there is no prediction key then there is only one prediction
try:
    n_pred = len(in_dict["Prediction"]) # Number of predictions
    n_frames = len(in_dict["Prediction"][0]["Frame"]) # Number of frames
    sam_iter = enumerate(in_dict["Prediction"])
except:
    n_pred = 1
    n_frames = len(in_dict["Frame"])
    sam_iter = [(0, in_dict)]

output_RP = {} # Dictionary for storing reference point info
frame_dict = [{} for i in range(n_frames)] # List of dictionaries containing output for each frame
for i, sample in sam_iter:
    for j, frame in enumerate(sample["Frame"]):
        for QoI, predictions in frame.items():
            if "Posterior_Sample" in predictions:
                for k, post_sam in enumerate(predictions["Posterior_Sample"]):
                    # Convert from list to numpy
                    post_sam = np.array(post_sam, dtype="float")
                    # Initialise entry in the reference point dictionary if it doesn't exist already
                    # (not necessary if doing RP cross-val separately)
                    if QoI+"_sam_"+str(k) not in output_RP:
                        output_RP[QoI+"_sam_"+str(k)] = np.empty((n_pred, n_frames))
                    
                    output_RP[QoI+"_sam_"+str(k)][i, j] = post_sam[1]
                    # Dictionary containing output
                    if n_pred > 1:
                        frame_dict[j][QoI+"_"+str(i)+"_sam_"+str(k)] = post_sam[2:]
                    else:
                        frame_dict[j][QoI+"_sam_"+str(k)] = post_sam[2:]
        
            else:
                # Otherwise the prediction is an average across the posterior
                if QoI not in output_RP:
                    output_RP[QoI] = np.empty((n_pred, n_frames))
                
                output_RP[QoI][i,j] = predictions[1]
                if n_pred > 1:
                    frame_dict[j][QoI+"_"+str(i)] = predictions[2:]
                else:
                    frame_dict[j][QoI] = predictions[2:]
                
# Create a new directory for the vtk files, if one does not exist already
if not(os.path.isdir("gp_predictions_nonlinear_" + infile)):
    os.mkdir("gp_predictions_nonlinear_" + infile)
    
# Loop over each frame, then write a separate frame which paraview can 
# interpret as a file series
for i, frame in enumerate(frame_dict):
    meshio.Mesh(points = nodes, cells = [("quad",faces)], point_data = frame).write("gp_predictions_nonlinear_" + infile + "\\frame_" + str(i) + ".vtk", file_format="vtk")
    

# Finally, plot the GP predictionsat the reference point
# Plots for samples
force = 200
if len([key for key in output_RP.keys() if "eta_mu_sam" in key]) > 0:
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1, 1, 1)
    #ind = range(18)
    ind = [0]
    for key, value in output_RP.items():
        value = value[ind,:]
        if "eta_mu_sam" in key:
            for displacement in value:
                ax.plot(-displacement, np.linspace(0,force,n_frames), "r")
        if "eta_sam_sam" in key:
            for displacement in value:
                ax.plot(-displacement, np.linspace(0,force,n_frames), "c", linewidth=0.25)
        if "eta_sigma_sam" in key:
            mu_key = "eta_mu_sam_" + key.strip("eta_sigma_sam_")
            mu = output_RP[mu_key][ind,:]
            for i, sd in enumerate(value):
                ax.plot(-mu[i,:]-2*sd, np.linspace(0,force,n_frames), "b")
                ax.plot(-mu[i,:]+2*sd, np.linspace(0,force,n_frames), "b")
            
    label_font = {'family': 'serif', 'size': 16,}
    ax.set_ylabel("Force (kN)", fontdict = label_font)
    ax.set_xlabel("Displacement (mm)", fontdict = label_font)

if len([key for key in output_RP.keys() if "eta_mu_mu" in key]) > 0:
    fig2 = plt.figure(figsize=(10,8))
    ax2 = fig2.add_subplot(1, 1, 1)
    # ind = range(10)
    ind = [1]
    for key, value in output_RP.items():
        value = value[ind,:]
        if "eta_mu_mu" in key:
            for displacement in value:
                ax2.plot(-displacement, np.linspace(0,force,n_frames), "r")
        if "eta_sam_mu" in key:
            for displacement in value:
                ax2.plot(-displacement, np.linspace(0,force,n_frames), "g")
        if "eta_sigma_mu" in key:
            mu_key = "eta_mu_mu" + key.strip("eta_sigma_mu_")
            mu = output_RP[mu_key][ind,:]
            for i, sd in enumerate(value):
                ax2.plot(-mu[i,:]-2*sd, np.linspace(0,force,n_frames), "b")
                ax2.plot(-mu[i,:]+2*sd, np.linspace(0,force,n_frames), "b")
            
    label_font = {'family': 'serif', 'size': 16,}
    ax2.set_ylabel("Force (kN)", fontdict = label_font)
    ax2.set_xlabel("Displacement (mm)", fontdict = label_font)
    
# Write data to csv files for plotting separately
for key, value in output_RP.items():
    np.savetxt(infile+"_RP_"+key+".csv", value, delimiter=',')

plt.show()
