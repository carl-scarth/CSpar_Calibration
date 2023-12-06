# Plots the basis functions of finite element model output
# Write a general function for plotting all outputs and place in parent directory??

import numpy as np
import meshio
import json
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import pandas as pd

# This code is getting messy - it would be good to package up some aspects to tidy

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

# Read in prediction data
with open(output_file, "r") as f:
    # Load in string from file
    in_dict = json.loads(f.readline())

# Get number of predictions and number of frames. If there is no prediction key then there is only one prediction
try:
    n_pred = len(in_dict["Prediction"]) # Number of predictions
    n_frames = len(in_dict["Prediction"][0]["Frame"]) # Number of frames
    sam_iter = enumerate(in_dict["Prediction"])
    if "Posterior_Sample" in list(in_dict["Prediction"][0]["Frame"][0].values())[0]:
        out_str = list(list(in_dict["Prediction"][0]["Frame"][0].values())[0]["Posterior_Sample"][0].keys())
    else:
        out_str = list(list(in_dict["Prediction"][0]["Frame"][0].values())[0].keys())
except:
    n_pred = 1
    n_frames = len(in_dict["Frame"])
    sam_iter = [(0, in_dict)]
    if "Posterior_Sample" in list(in_dict["Frame"][0].values())[0]:
        out_str = list(list(in_dict["Frame"][0].values())[0]["Posterior_Sample"][0].keys())
    else:
        out_str = list(list(in_dict["Frame"][0].values())[0].keys())

n_out = len(out_str)
if "w" in out_str:
    output_RP = {} # Dictionary for storing reference point info

output_max = {} # Dictionary for storing outputs at location of maximum displacement
max_ind = {'u': 4164, 'w': 1299} # Index of node containing maximum absolute value across training data
# Ideally I'd caluclate this here but it's a little too messy.
frame_dict = [{} for i in range(n_frames)] # List of dictionaries containing output for each frame
for i, sample in sam_iter:
    for j, frame in enumerate(sample["Frame"]):
        for QoI, predictions in frame.items():
            if "Posterior_Sample" in predictions:
                n_post_sam = len(predictions["Posterior_Sample"])
                for k, post_sam in enumerate(predictions["Posterior_Sample"]):
                    # Convert from list to numpy
                    for coord in out_str:
                        post_sam = np.array(post_sam[coord], dtype="float")
                        # Initialise entry in the reference point dictionary if it doesn't exist already
                        # (not necessary if doing RP cross-val separately)
                        # Consider trying to condense the below by squeezing an array?/re-naming below
                        if n_pred > 1:
                            if QoI+"_sam_"+str(k) not in output_max:
                                output_max[QoI+"_sam_"+str(k)] = {key : np.empty((n_pred, n_frames)) for key in out_str}
                                if "w" in out_str:
                                    output_RP[QoI+"_sam_"+str(k)] = np.empty((n_pred, n_frames))
                            
                            output_max[QoI+"_sam_"+str(k)][coord][i, j] = post_sam[max_ind[coord]]
                            if "w" in out_str:
                                output_RP[QoI+"_sam_"+str(k)][i, j] = post_sam[1]
                        
                        else:
                            if QoI not in output_max:
                                output_max[QoI] = {key : np.empty((n_frames, n_post_sam)) for key in out_str}
                                if coord == "w":
                                    output_RP[QoI] = np.empty((n_frames, n_post_sam))
                            
                            output_max[QoI][coord][j,k] = post_sam[max_ind[coord]]
                            if coord == "w":
                                output_RP[QoI][coord][j,k] = post_sam[1]
                    
                    # Dictionary containing output
                    if n_pred > 1:
                        frame_dict[j]["_".join((QoI,coord,str(i),"sam",str(k)))] = post_sam[2:]
                    else:
                        frame_dict[j]["_".join((QoI,coord,"sam",str(k)))] = post_sam[2:]
        
            else:
                # Otherwise the prediction is an average across the posterior
                for coord in out_str:
                    if QoI not in output_max:
                        if n_pred > 1:
                            output_max[QoI] = {key : np.empty((n_pred, n_frames)) for key in out_str}
                            if "w" in out_str:
                                output_RP[QoI] = np.empty((n_pred, n_frames))

                        else:
                            output_max[QoI] = {key : np.empty(n_frames) for key in out_str}
                            if "w" in out_str:
                                output_RP[QoI] = np.empty(n_frames)

                    if n_pred > 1:
                        output_max[QoI][coord][i,j] = predictions[coord][max_ind[coord]]
                        if coord == "w":
                            output_RP[QoI][i,j] = predictions[coord][1]
                        frame_dict[j]["_".join((QoI,coord,str(i)))] = predictions[coord][2:]
                    else:                        
                        output_max[QoI][coord][j] = predictions[coord][max_ind[coord]]
                        if coord == "w":
                            output_RP[QoI][j] = predictions[1]
                        frame_dict[j]["_".join((QoI, coord))] = predictions[2:]

# Create a new directory for the vtk files, if one does not exist already
if not(os.path.isdir("gp_predictions_nonlinear_" + infile)):
    os.mkdir("gp_predictions_nonlinear_" + infile)
    
# Loop over each frame, then write a separate frame which paraview can 
# interpret as a file series
for i, frame in enumerate(frame_dict):
    meshio.Mesh(points = nodes, cells = [("quad",faces)], point_data = frame).write("gp_predictions_nonlinear_" + infile + "\\frame_" + str(i) + ".vtk", file_format="vtk")

# Finally, plot the GP predictionsat the reference point/max point
# Plots for samples
force = 200
rcParams.update({'figure.figsize' : (8,6),
                'font.size' : 14,
                'font.family' : 'serif',
                'figure.titlesize' : 16,
                'axes.labelsize': 16,
                'xtick.labelsize': 14,
                'ytick.labelsize': 14,
                'legend.fontsize': 14})

# The below code is getting very messy - consider packaging. Possible example would be 
# package code for plotting output at specified node
if n_pred > 1:
    # ind = [0]
    if "w" in out_str:
        if len([key for key in output_RP.keys() if "eta_mu_sam" in key]) > 0:
            fig, ax = plt.subplots()
            for key, value in output_RP.items():
                # value = value[ind,:]
                if "eta_mu_sam" in key:
                    for displacement in value:
                        ax.plot(-displacement, np.linspace(0,force,n_frames), "r")
                        if "eta_sam_sam" in key:
                            for displacement in value:
                                ax.plot(-displacement, np.linspace(0,force,n_frames), "c", linewidth=0.25)
                if "eta_sigma_sam" in key:
                    mu_key = "eta_mu_sam_" + key.strip("eta_sigma_sam_")
                    mu = output_RP[mu_key]#[ind,:]
                    for i, sd in enumerate(value):
                        ax.plot(-mu[i,:]-2*sd, np.linspace(0,force,n_frames), "b")
                        ax.plot(-mu[i,:]+2*sd, np.linspace(0,force,n_frames), "b")   
        
            ax.set_ylabel("Force (kN)")
            ax.set_xlabel("Displacement (mm)")
            ax.set_title("Posterior sample longitudinal displacement at reference point")
    
    if len([key for key in output_max.keys() if "eta_mu_sam" in key]) > 0:
        fig2, axs2 = plt.subplots(1, n_out)
        for key, disp_dict in output_max.items():
            for i, (disp_key, value) in enumerate(disp_dict.items()):
                # value = value[ind,:]
                if "eta_mu_sam" in key:
                    for displacement in value:
                        axs2[i].plot(-displacement, np.linspace(0,force,n_frames), "r")
                if "eta_sam_sam" in key:
                    for displacement in value:
                        axs2[i].plot(-displacement, np.linspace(0,force,n_frames), "c", linewidth=0.25)
                if "eta_sigma_sam" in key:
                    mu_key = "eta_mu_sam_" + key.strip("eta_sigma_sam_")
                    mu = output_max[mu_key][disp_key]#[ind,:]
                    for j, sd in enumerate(value):
                        axs2[i].plot(-mu[j,:]-2*sd, np.linspace(0,force,n_frames), "b")
                        axs2[i].plot(-mu[j,:]+2*sd, np.linspace(0,force,n_frames), "b")   
                    
        for ax2 in axs2:
            ax2.set_ylabel("Force (kN)")
            ax2.set_xlabel("Displacement (mm)")
            ax2.set_title(disp_key + "at node " + str(max_ind))
        fig2.suptitle("Prediction at location of maximum displacement across training data")
adsadsadsadsad
# Output names depends on if the calibration or emulator code is used
elif "eta_mu" in output_RP.keys():
    fig3, ax3 = plt.subplots()
    for key, value in output_RP.items():
        if key == "eta_mu":
            ax3.plot(-value, np.linspace(0,force,n_frames), "r")
        if key == "eta_sam":
            ax3.plot(-value, np.linspace(0,force,n_frames), "c", linewidth=0.25)
        if key == "eta_sigma":
            mu = output_RP["eta_mu"]
            ax3.plot(-mu-2*value, np.linspace(0,force,n_frames), "b")
            ax3.plot(-mu+2*value, np.linspace(0,force,n_frames), "b")   
            
    ax3.set_ylabel("Force (kN)")
    ax3.set_xlabel("Displacement (mm)")
    ax3.set_title("Posterior sample displacement at reference point")

    fig4, ax4 = plt.subplots()
    for key, value in output_max.items():
        if key == "eta_mu":
            ax4.plot(-value, np.linspace(0,force,n_frames), "r")
        if key == "eta_sam":
            ax4.plot(-value, np.linspace(0,force,n_frames), "c", linewidth=0.25)
        if key == "eta_sigma":
            mu = output_max["eta_mu"]
            ax3.plot(-mu-2*value, np.linspace(0,force,n_frames), "b")
            ax3.plot(-mu+2*value, np.linspace(0,force,n_frames), "b")   
            
    ax4.set_ylabel("Force (kN)")
    ax4.set_xlabel("Displacement (mm)")
    ax4.set_title("Posterior sample displacement at node " + str(max_ind))

# Plot mean values
if len([key for key in output_RP.keys() if "eta_mu_mu" in key]) > 0:
    fig5, ax5 = plt.subplots()
    # ind = range(10)
    # ind = [0]
    for key, value in output_RP.items():
        # value = value[ind,:]
        if "eta_mu_mu" in key:
            if value.ndim > 1:
                for displacement in value:
                    ax5.plot(-displacement, np.linspace(0,force,n_frames), "r")
            else:
                ax5.plot(-value, np.linspace(0,force,n_frames), "r")
        if "eta_sam_mu" in key:
            if value.ndim > 1:
                for displacement in value:
                    ax5.plot(-displacement, np.linspace(0,force,n_frames), "g")
            else:
                ax5.plot(-value, np.linspace(0,force,n_frames), "g")
        if "eta_sigma_mu" in key:
            mu_key = "eta_mu_mu" + key.strip("eta_sigma_mu_")
            mu = output_RP[mu_key]#[ind,:]
            if value.ndim > 1: 
                for i, sd in enumerate(value):
                    ax5.plot(-mu[i,:]-2*sd, np.linspace(0,force,n_frames), "b")
                    ax5.plot(-mu[i,:]+2*sd, np.linspace(0,force,n_frames), "b")
            else:
                ax5.plot(-mu-2*value, np.linspace(0,force,n_frames), "b")
                ax5.plot(-mu+2*value, np.linspace(0,force,n_frames), "b")
        if "eta_mu_sigma" in key:
            mu_key = "eta_mu_mu" + key.strip("eta_sigma_mu_")
            mu = output_RP[mu_key]
            if value.ndim > 1: 
                for i, sd in enumerate(value):
                    ax5.plot(-mu[i,:]-2*sd, np.linspace(0,force,n_frames), "c")
                    ax5.plot(-mu[i,:]+2*sd, np.linspace(0,force,n_frames), "c")
            else:
                ax5.plot(-mu-2*value, np.linspace(0,force,n_frames), "c")
                ax5.plot(-mu+2*value, np.linspace(0,force,n_frames), "c")
        
        if "eta_sam_sigma" in key:
            mu_key = "eta_sam_mu" + key.strip("eta_sigma_mu_")
            mu = output_RP[mu_key]
            if value.ndim > 1: 
                for i, sd in enumerate(value):
                    ax5.plot(-mu[i,:]-2*sd, np.linspace(0,force,n_frames), "m")
                    ax5.plot(-mu[i,:]+2*sd, np.linspace(0,force,n_frames), "m")
            else:
                ax5.plot(-mu-2*value, np.linspace(0,force,n_frames), "m")
                ax5.plot(-mu+2*value, np.linspace(0,force,n_frames), "m")
            
    ax5.set_ylabel("Force (kN)")
    ax5.set_xlabel("Displacement (mm)")
    ax5.set_title("Summary statistics across posterior for reference point displacement")

if len([key for key in output_max.keys() if "eta_mu_mu" in key]) > 0:
    fig6, ax6 = plt.subplots()
    # ind = range(10)
    # ind = [0]
    for key, value in output_max.items():
        # value = value[ind,:]
        if "eta_mu_mu" in key:
            if value.ndim > 1:
                for displacement in value:
                    ax6.plot(-displacement, np.linspace(0,force,n_frames), "r")
            else:
                ax6.plot(-value, np.linspace(0,force,n_frames), "r")
        if "eta_sam_mu" in key:
            if value.ndim > 1:
                for displacement in value:
                    ax6.plot(-displacement, np.linspace(0,force,n_frames), "g")
            else:
                ax6.plot(-value, np.linspace(0,force,n_frames), "g")
        if "eta_sigma_mu" in key:
            mu_key = "eta_mu_mu" + key.strip("eta_sigma_mu_")
            mu = output_max[mu_key]#[ind,:]
            if value.ndim > 1: 
                for i, sd in enumerate(value):
                    ax6.plot(-mu[i,:]-2*sd, np.linspace(0,force,n_frames), "b")
                    ax6.plot(-mu[i,:]+2*sd, np.linspace(0,force,n_frames), "b")
            else:
                ax6.plot(-mu-2*value, np.linspace(0,force,n_frames), "b")
                ax6.plot(-mu+2*value, np.linspace(0,force,n_frames), "b")
        if "eta_mu_sigma" in key:
            mu_key = "eta_mu_mu" + key.strip("eta_sigma_mu_")
            mu = output_max[mu_key]
            if value.ndim > 1: 
                for i, sd in enumerate(value):
                    ax6.plot(-mu[i,:]-2*sd, np.linspace(0,force,n_frames), "c")
                    ax6.plot(-mu[i,:]+2*sd, np.linspace(0,force,n_frames), "c")
            else:
                ax6.plot(-mu-2*value, np.linspace(0,force,n_frames), "c")
                ax6.plot(-mu+2*value, np.linspace(0,force,n_frames), "c")
        
        if "eta_sam_sigma" in key:
            mu_key = "eta_sam_mu" + key.strip("eta_sigma_mu_")
            mu = output_max[mu_key]
            if value.ndim > 1: 
                for i, sd in enumerate(value):
                    ax6.plot(-mu[i,:]-2*sd, np.linspace(0,force,n_frames), "m")
                    ax6.plot(-mu[i,:]+2*sd, np.linspace(0,force,n_frames), "m")
            else:
                ax6.plot(-mu-2*value, np.linspace(0,force,n_frames), "m")
                ax6.plot(-mu+2*value, np.linspace(0,force,n_frames), "m")
            
    ax6.set_ylabel("Force (kN)")
    ax6.set_xlabel("Displacement (mm)")
    ax6.set_title("Summary statistics across posterior for displacement at node " + str(max_ind))
    
#if n_pred == 1:
#    # output_RP = {key : value.flatten() for key, value in output_RP.items()}
#    RP_frame = pd.DataFrame(output_RP)
#    RP_frame.to_csv(infile+"_RP_pred.csv",sep=",",index=False)
#else:
#    # Write data to csv files for plotting separately
for key, value in output_RP.items():
    np.savetxt(infile+"_RP_"+key+".csv", value, delimiter=',')
for key, value in output_max.items():
    np.savetxt(infile+"_max_"+key+".csv", value, delimiter=',')

plt.show()
