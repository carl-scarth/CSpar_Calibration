# Loads in a Design of Experiments and runs the C_spar model for this DoE

import numpy as np
import pandas as pd
import sys
import os

model_dir = "C:\\Users\\cs2361\\Documents\\Jean_Benezech\\CompositesFEMesh\\Csection_mesh" # directory of C Spar model
gridMod_dir = "C:\\Users\\cs2361\\Documents\\Jean_Benezech\\CompositesFEMesh\\C++\\gridModification\\build" # directory of grid modification code
# Add the model directory to the Python path
sys.path.append(model_dir)
# Import the required methods from the model directory
from write_parameters import write_parameters
from write_mesh import *
from write_inp import *
from write_shell_inp import *
from write_shell_parameters import *

infile = "inputs\\LHSDesign75x7"
# infile = "inputs\\LHSDesign30x3" # file in which DoE is stored
# infile = "inputs\\LHSDesign60x6"
# infile = "Problem_run"
# infile = "inputs\\LHSDesign100x7_1" # file in which DoE is stored
# infile = "inputs\\LHSDesign50x4" # file in which DoE is stored
# infile = "inputs\\nominal_inputs" # file in which DoE is stored

change_inc = True # Do I want to play with the increment size?
write_buffer = False # Do I want to write a temporary file to store displacements as I go?
restart = False # Am I restarting a previous analysis?
shell_mesh = True # Is the mesh comprised of shells?

max_inc = 0.1  # maximum increment
init_inc = max_inc # initial increment. Set equal to maximum increment in the hope that this keeps the output regular
min_inc = 1.0e-5 # minimum increment
load = -250.0 # Applied load

# Numpy version
# Load in the Design of Experiments
# x_DoE = np.loadtxt(infile + ".csv", delimiter = ",", skiprows = 1, ndmin = 2)
#N, d = x_DoE.shape # Number of samples and number of inputs
#print(x_DoE)

# Pandas version
x_DoE = pd.read_csv(infile + ".csv", sep = ",")
N, d = x_DoE.shape # Number of samples and number of inputs

if shell_mesh:
    # I think this is correct as the 5 reference nodes from Jean's mesh get omitted from the output - check this!
    n_Nodes = 9353
else:
    # n_Nodes = 116876 # Number of nodes per simulation - Old version before adding new reference nodes
    n_Nodes = 116877 # Number of nodes per simulation (this increased with the additional BCs) 3D Mesh

# Define inputs which are held constant
model_name = 'CSpar_sam' # Name of the Abaqus model
Zlength = 420.0 # effective length of the spar (between end blocks).
height = 55.0 # distance from tip of the flanges to the outer surface of the web

# If we are restarting an aborted set of runs, then reload existing data, otherwise initialise entities used to store data
if restart:
    # I think this might require some correction. Not sure though
    # Note: displacements were transposed before writing to the buffer, and so need to be transposed again for 
    # compatibility with the loop below
    displacements = np.loadtxt('displacement_buffer.csv', delimiter = ',').T
    print(displacements.shape)
    with open ('increment_buffer.csv','r') as f:
        incs = [[float(increment) for increment in line.strip().split()] for line in f.readlines()]
    # Reaction forces are flattened into a single-row numpy array to match format from loop
    RFs = np.empty([1,0])
    with open ('RF_buffer.csv','r') as f:
        for line in f.readlines():
            RFs = np.concatenate((RFs, np.array([float(RF) for RF in line.strip().split()], ndmin = 2)), axis=1)
    N_complete = len(incs) # number of runs completed before being aborted
    # Create iterable object for for loop
    # iterable = enumerate(x_DoE[N_complete:,:], start=N_complete)
    iterable = x_DoE[N_complete:,:].iterrows()
    # It isn't possible to start from a different index with iterrows. Perhaps it doesn't matter
    # as the index is a property of the dataframe anyway. My hope is that this will fix itself.
    # Check if needed....

else:
    # Do we want to vary the increment size? If not, stick with default settings?
    if change_inc:
        # Don't know in advance how many frames there will be
        displacements = np.empty([n_Nodes,0])
        incs = []
        RFs = np.empty([1,0])
    else:
        init_inc = 1.0
        min_inc = 1.0e-5
        max_inc = 1.0
        displacements = np.empty([n_Nodes,3*N])
    # iterable = enumerate(x_DoE)
    iterable = x_DoE.iterrows()

# nodes = np.empty([n_Nodes,3*N])
# Loop over entries of the DoE and run the Abaqus model
for i, x_i in iterable:
    print(i)
    print(x_i)
        
    
    # t_ply = x_i[5] # Extract ply thickness from the DoE
    # Use conditional for now to extract parameter values - this will be neater with Pandas
    # if d == 3 or d == 4:
    #    t_ply = x_i[1] 
    #elif d >= 6 and d <= 8:
    #    t_ply = x_i[5]
    if "t_ply" in x_i:
        t_ply = x_i["t_ply"]
    else:
        t_ply = 0.196

    print(t_ply)
    dsadsad    
    # Write parameters to a text file which will be used to generate mesh
    write_parameters(t_ply = t_ply, Zlength = Zlength, height = height, model_name = model_name)
    sdsadsa
    # Generate the mesh using Gmsh
    write_mesh()
    # Run the gridModification code to create the ramp in the spar
    command = gridMod_dir + "\\gridMod"
    os.system(command)
    sadasdasd
    # Extract other material properties from the DoE
    # Sort with Pandas later    
    if d == 3:
       E11 = x_i[0]
       K = x_i[2]
       # Write input file for Abaqus
       write_inp(E11 = E11, K=K, load=load, init_inc=init_inc, min_inc=min_inc, max_inc=max_inc)
    elif d == 4:
        E11 = x_i[0]
        K = x_i[2]
        x_spring = x_i[3]
        write_inp(E11 = E11, K=K, x_spring=x_spring, load=load, init_inc=init_inc, min_inc=min_inc, max_inc=max_inc)
    elif d == 6:
        E11 = x_i[0]
        E22 = x_i[1]
        nu12 = x_i[2]
        nu23 = x_i[3]
        G12 = x_i[4]
        write_inp(E11 = E11, E22 = E22, nu12 = nu12, nu23 = nu23, G12 = G12, load=load, init_inc=init_inc, min_inc=min_inc, max_inc=max_inc)
    elif d == 7:
        E11 = x_i[0]
        E22 = x_i[1]
        nu12 = x_i[2]
        nu23 = x_i[3]
        G12 = x_i[4]
        K = x_i[6]
        write_inp(E11 = E11, E22 = E22, nu12 = nu12, nu23 = nu23, G12 = G12, K=K, load=load, init_inc=init_inc, min_inc=min_inc, max_inc=max_inc)
    elif d == 8:
        E11 = x_i[0]
        E22 = x_i[1]
        nu12 = x_i[2]
        nu23 = x_i[3]
        G12 = x_i[4]
        K = x_i[6]
        x_spring = x_i[7]
        write_inp(E11 = E11, E22 = E22, nu12 = nu12, nu23 = nu23, G12 = G12, K=K, x_spring=x_spring, load=load, init_inc=init_inc, min_inc=min_inc, max_inc=max_inc)
    
    # Run Abaqus from the command line
    command = "Abaqus Job=" + model_name + " input=\"Abaqus\\" + model_name + ".inp\" interactive ask_delete=OFF cpus=4"
    os.system(command)
    
    # Process the abaqus output to extract displacements and nodal coordinates
    if change_inc:
        # If multiple frames are requested, run the code which processes all frames
        command = "abaqus viewer noGUI=process_outputs_all_frames.py -- " + model_name
    else:
        # If only one increment is required just run the script which processes the final frame
        command = "abaqus viewer noGUI=process_outputs.py -- " + model_name
    
    os.system(command)

    # sense. The increment index outputted by the postprocessing file should help keep track of things.
    # At some point consider looking at numpy options for 3d arrays, and json files. This is going to be an 
    # issue for the experimental data as well so worth getting a decent structure set downl
    disp_i = np.loadtxt(model_name + "_displacement.csv", delimiter = ",", skiprows = 0)
    # nodes_i = np.loadtxt(model_name + "_nodes.csv", delimiter = ",", skiprows = 0)
    if change_inc:
        displacements = np.concatenate((displacements,disp_i), axis=1)
        # Also load in increment indices and reaction forces
        incs_i = np.loadtxt(model_name + "_increments.csv", delimiter = ",", skiprows = 0)
        incs.append(incs_i.tolist())
        RFs_i = np.loadtxt(model_name + "_RFs.csv", delimiter = ",", skiprows = 0)
        RFs = np.concatenate((RFs, RFs_i.reshape((1, -1))), axis=1)
    else:
        displacements[:,3*i:3*(i+1)] = disp_i

    # nodes[:,3*i:3*(i+1)] = nodes_i

    if write_buffer:
        with open('displacement_buffer.csv', 'a') as f:
            np.savetxt(f, disp_i.T, delimiter = ',')
        with open('increment_buffer.csv', 'a') as f:
            np.savetxt(f, incs_i.reshape((1,-1)))
        with open('RF_buffer.csv', 'a') as f:
            np.savetxt(f, RFs_i.reshape(1,-1))
        
# Consider using json format, stick with csv for now. This should be relatively easy if working with python dictionaries

# Write model outputs from all simulations to csv files   
# Create header for csv
if change_inc:      
    head_str = ', '.join(['u_' + str(i+1) + '_' + str(j+1) + ', v_' + str(i+1) + '_' + str(j+1) + ', w_' + str(i+1) + '_' + str(j+1) for i in range(N) for j in range(len(incs[i]))])
else:
    head_str = ', '.join(['u_' + str(i+1)+', v_'+str(i+1)+', w_'+str(i+1) for i in range(N)])

np.savetxt(infile + "_displacements_load=" + str(load) + "_max_inc=" + str(max_inc) + ".csv", displacements, delimiter=",", header = head_str, comments = "")
# head_str = ', '.join(['x_' + str(i+1)+', y_'+str(i+1)+', z_'+str(i+1) for i in range(N)])
# np.savetxt(infile + "_nodes_load=" + str(load) + "_max_inc=" + str(max_inc) + ".csv", nodes, delimiter=",", header = head_str, comments = "")

if change_inc:
    with open(infile + "_incs_load=" + str(load) + "_max_inc=" + str(max_inc) + ".txt", 'w') as f:
        [f.write(', '.join([str(inc) for inc in sam])+"\n") for sam in incs]
    head_str = ', '.join(['Rx_' + str(i+1) + '_' + str(j+1) + ', Ry_' + str(i+1) + '_' + str(j+1) + ', Rz_' + str(i+1) + '_' + str(j+1) for i in range(N) for j in range(len(incs[i]))])
    np.savetxt(infile + "_RFs_load=" + str(load) + "_max_inc=" + str(max_inc) + ".csv", RFs, delimiter=",", header = head_str, comments = "")

# Delete temporary files
if write_buffer:
    for filename in ['displacement_buffer.csv','increment_buffer.csv','RF_buffer.csv']:
        if os.path.isfile(filename): # Check if file exists first to avoid errors
            os.remove(filename)
