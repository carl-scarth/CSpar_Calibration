# Create vtk files for 3D plots of calibration outputs

import numpy as np
import meshio
import json
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
from utils import set_plot_params

def plot_basis_nonlinear_vtk(in_dict = {}, file_str = "", out_folder = "", mesh_str = "", from_json = False, skip_nodes = []):
    # Produce vtks of svd basis and training data mean
    # in_dict = dictionary containing structured output data. See abaqus_json.py for details of structure
    # file_str = identifier of .json containing structured output, if used, and for the folder where vtk output will be stored
    # out_folder = folder where output data is to be read from/written to
    # mesh_str = string used to identify .csvs containing node and element definitions
    # from_json: Boolean indicating to read output from a json
    # skip_nodes: Slice of nodes to skip for removing output at reference points which don't appear in the nodes or connectivity
    
    # Load in and process mesh data
    node_file = mesh_str + "_nodes.csv" # File 
    element_file = mesh_str + "_elements.csv"

    # Read in the element and node definitions
    elements = np.loadtxt(element_file, dtype = int, delimiter = ',')
    nodes = np.loadtxt(node_file, delimiter = ',')
    # Remove the first column of nodes and elements which is just an index
    nodes = nodes[:,1:]
    elements = elements[:,1:]

    # Define list of indice of the nodes which define each face of the brick
    face_nodes = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]]
    # Convert Brick element connectivity into 6 shell elements for the faces
    faces = np.empty([0,4], dtype = int)
    for face in face_nodes:
        elements_face = elements[:,face]
        faces = np.concatenate((faces,elements_face))

    # Convert connectivity from abaqus to python indexing
    faces = faces - 1

    # Load in json data if not passed directly
    # Raise Exception if the wrong combination of inputs has been passed
    if not from_json and not in_dict:
        raise Exception("Input either in_dict = dictionary of output data, or set from_json = True and pass file_str to identify .json file containing output data")
    if from_json:
        in_dict = load_from_json(file_str = file_str, output_type = "basis_nonlinear", out_folder = out_folder)
    
    # Convert Basis and mean data to numpy
    for frame in in_dict["Frame"]:
        frame["Bases"] = {label : np.array(value, dtype="float") for label, value in frame["Bases"].items()}
        frame["Training_Data_Mean"] = {label : np.array(value, dtype="float") for label, value in frame["Training_Data_Mean"].items()}
    
    # Extract number of bases, and output labels
    out_str = list(in_dict["Frame"][0]["Bases"].keys())
    n_bases = in_dict["Frame"][0]["Bases"][out_str[0]].shape[1]

    # Delete Reference point nodes from output (meshio raises error otherwise)
    for frame in in_dict["Frame"]:
        for label in out_str:
            frame["Bases"][label] = np.delete(frame["Bases"][label], skip_nodes, axis=0)
            frame["Training_Data_Mean"][label] = np.delete(frame["Training_Data_Mean"][label], skip_nodes, axis=0)

    # Create a new directory for the vtk files, if one does not exist already
    out_subfolder = os.path.join(out_folder, "basis_nonlinear_" + file_str)
    if not(os.path.isdir(out_subfolder)):
        os.mkdir(out_subfolder)

    # Loop over each frame, then write a separate vtk
    for i, frame in enumerate(in_dict["Frame"]):
        # Create a dictionary using output data for each basis
        basis_dict = {}
        for label in out_str:
            for j in range(n_bases):
                basis_dict["_".join(("basis", label, str(j+1)))] = frame["Bases"][label][:,j]
            # Get data in the correct format for meshio
            basis_dict["Training_Data_Mean_" + label] = frame["Training_Data_Mean"][label]
        
        meshio.Mesh(points = nodes, cells = [("quad",faces)], point_data = basis_dict).write(os.path.join(out_subfolder, "frame_" + str(i) + ".vtk"), file_format="vtk")

def plot_basis_nonlinear_point(point_ind, in_dict = {}, file_str = "", out_folder = "", from_json = False, location_title = ""):
    # Produce plots of basis and training data mean against frame number at a fixed point
    # point_ind = dictionary with key matching displacement components, and value giving the corresponding node index for plotting output
    # in_dict = dictionary containing structured output data. See abaqus_json.py for details of structure
    # file_str = identifier of .json containing structured output, if used, and for the folder where vtk output will be stored
    # out_folder = folder where output data is to be read from/written to
    # from_json: Boolean indicating to read output from a json
    # location_title: Label to be placed in figure title

    # Update plot parameters
    set_plot_params()

    # Raise Exception if the wrong combination of inputs has been passed
    if not from_json and not in_dict:
        raise Exception("Input either in_dict = dictionary of output data, or set from_json = True and pass file_str to identify .json file containing output data")
    # Load from json if needed
    if from_json:
        in_dict = load_from_json(file_str = file_str, output_type = "basis_nonlinear", out_folder = out_folder)
    
    # Extract number of frames, number of bases, output labels
    n_frames = len(in_dict["Frame"])
    # out_str = [item for item in list(in_dict["Frame"][0]["Bases"].keys()) if item in point_ind.keys()]
    out_str = list(point_ind.keys())
    n_bases = in_dict["Frame"][0]["Bases"][out_str[0]].shape[1]
    
    # Initialise dictionaries for storing output
    basis_point = {label : np.empty((n_frames,n_bases)) for label in out_str}
    # Loop over each output frame and populate point output arrays
    for i, frame in enumerate(in_dict["Frame"]):
        for label in out_str:
            basis_point[label][i,:] = frame["Bases"][label][point_ind[label],:]

    # Produce plots
    fig, axes = plt.subplots(len(out_str),1, sharex=True)
    if len(out_str) == 1:
        axes = [axes]
    for ax, out_str_ax in zip(axes, out_str):
        for i, base in enumerate(basis_point[out_str_ax].T):
            ax.plot([j for j in range(n_frames)], base, label="K_"+str(i+1))
            ax.set_title(out_str_ax)
            ax.set_xlim(0,n_frames-1)
            ax.set_ylabel("K_i")
            ax.legend()
    axes[-1].set_xlabel("Frame")
    if location_title:
        fig.suptitle("Basis at " + location_title)

def load_from_json(file_str = "", output_type = "", out_folder = ""):
    # Load in output data from a json file
    if not file_str:
        raise Exception("If passing output data via json, also set file_str with the name of the .json file")

    # Get name of input .json if using
    filename = os.path.join(out_folder,"_".join((output_type,file_str)) + ".json")
    # Read in basis data from json
    with open(filename, "r") as f:
        # Load in string from file
        in_dict = json.loads(f.readline())

    return(in_dict)
