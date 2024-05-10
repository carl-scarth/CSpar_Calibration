    # Loads in a Design of Experiments and runs the C_spar model for this DoE

import numpy as np
import pandas as pd
import sys
import os
import json

   
model_dir = "C:\\Users\\cs2361\\Documents\\Jean_Benezech\\CompositesFEMesh\\Csection_mesh" # directory of C Spar model
gridMod_dir = "C:\\Users\\cs2361\\Documents\\Jean_Benezech\\CompositesFEMesh\\C++\\gridModification\\build" # directory of grid modification code
# Add the model directory to the Python path
sys.path.append(model_dir)
# Import the required methods from the model directory
from write_parameters import write_parameters
from write_mesh import *
from write_inp import *
from write_shell_inp_beam import write_inp as write_shell_inp_beam
from write_shell_inp_translator import write_inp as write_shell_inp
from write_shell_parameters import write_parameters as write_shell_parameters

def set_input(x_series, x_name, default_val):
    # Checks if Pandas series "x_series" contains an entry indexed by string "x_name". If True sets the input x to this value,
    # if not uses a default value
    if x_name in x_series:
        x = x_series[x_name]
    else:
        x = default_val
    return x

infile = "inputs\\translator_study"

change_inc = True # Do I want to play with the increment size?
shell_mesh = True # Is the mesh comprised of continuum shells?
store_all_sam = False
rotate_flanges = False
thin_corners = False
apply_load = True
fix_to_ground = True

max_inc = 0.05  # maximum increment
init_inc = max_inc # initial increment. Set equal to maximum increment in the hope that this keeps the output regular
min_inc = 1.0e-5 # minimum increment
load = -250.0 # Applied load
applied_disp = -4.0 # Applied displacement (specified as alternative to load)
sym = True # Representing only half of a symmetric layup
layup = 3*[45.0,45.0,-45.0,-45.0, 90.0, 90.0, 0.0, 0.0] #, 45.0, -45.0, 90.0, 0.0, 45.0, -45.0, 90.0, 0.0] # Layup of new spar
# layup = [45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0] # Half layup for old CSpar
if sym:
    n_plies = 2*len(layup)
else:
    n_plies = len(layup) # Number of plies
# Define inputs which are held constant
model_name = 'CSpar_sam' # Name of the Abaqus model
Zlength = 420.0 # effective length of the spar (between end blocks).
height = 55.0 # distance from tip of the flanges to the outer surface of the web.
rotation_offset = 160.0 # Offset from bearing centre to spar ends

# Read in DoE
x_DoE = pd.read_csv(infile + ".csv", sep = ",")
print(x_DoE)

N, d = x_DoE.shape # Number of samples and number of inputs

if shell_mesh:
    if fix_to_ground:
        n_Nodes = 9354
    else:
        n_Nodes = 9355
else:
    # n_Nodes = 116877 # Number of nodes per simulation (this increased with the additional BCs) 3D Mesh
    n_Nodes = (n_plies + 1)*4675 + 2 # Assumes mesh density is fixed
print(n_Nodes)

# Do we want to vary the increment size? If not, stick with default settings?
if change_inc:
    # Don't know in advance how many frames there will be
    displacements = np.empty([n_Nodes,0])
    incs = []
    RFs = np.empty([1,0])
    out_dict = {"Sample" : []}
else:
    # Otherwise just look at the final frame
    init_inc = 1.0
    min_inc = 1.0e-5
    max_inc = 1.0
    displacements = np.empty([n_Nodes,3*N])

# nodes = np.empty([n_Nodes,3*N]) # Uncomment if storing nodes
# Loop over entries of the DoE and run the Abaqus model
for i, x_i in x_DoE.iterrows():
    print(i)
    print(x_i)
    
    # Extract inputs which govern the geometry, or othewise set to their default values
    t_ply = set_input(x_i, "t_ply", 0.125)
    t_ply_rad = set_input(x_i, "t_ply_rad", 0.125)
    if "Flange_theta" in x_i:
        LFlange_theta = set_input(x_i, "Flange_theta", 0.0)
        RFlange_theta = set_input(x_i, "Flange_theta", 0.0)
    else:
        LFlange_theta = set_input(x_i, "LFlange_theta", 0.0)
        RFlange_theta = set_input(x_i, "RFlange_theta", 0.0)
    rad_thin = set_input(x_i, "rad_thin", 0.0)

    # Give unique filenames to each run?
    if store_all_sam:
        file_str = "_".join((model_name, str(i)))
    else:
        file_str = model_name

    # Pick different input file generating script depending on if the mesh is comprised of shells or not
    # ideally I'd just use the same file taking a Boolean input, but keep for now in case of changes
    if shell_mesh:
        write_shell_parameters(t_ply=t_ply, rad_thin = rad_thin, Zlength=Zlength, height=height+n_plies*t_ply, LFlange_theta=LFlange_theta, RFlange_theta=RFlange_theta, n_plies=n_plies, rotate_flanges = rotate_flanges, thin_corners = thin_corners, model_name=file_str)
    else:
        # Write parameters to a text file which will be used to generate mesh
        # write_parameters(t_ply = t_ply, Zlength = Zlength, height = height, model_name = model_name)
        write_parameters(t_ply = t_ply, Zlength = Zlength, height = height, StackSeq = layup, rotate_flanges = rotate_flanges, model_name = file_str)
  
    # Generate the mesh using Gmsh
    write_mesh()
    # Run the gridModification code to create the ramp in the spar
    command = gridMod_dir + "\\gridMod"
    os.system(command)
    
    # Extract other material properties from the DoE
    E11 = set_input(x_i, "E11", 140.9)
    E22 = set_input(x_i, "E22", 8.96)
    nu12 = set_input(x_i, "nu12", 0.32)
    nu23 = set_input(x_i, "nu23", 0.43)
    G12 = set_input(x_i, "G12", 4.69)
    x_spring = set_input(x_i, "x_spring", 46.27)
    x_spring_error = set_input(x_i, "x_spring_error", 0.0)
    x_misalign_slope = set_input(x_i, "x_misalign_slope", 0.0)
    pivot_offset_error = set_input(x_i, "pivot_offset_error", 0.0)
    # Calculate the eccentricity of the support pivot at either end
    x_spring_fix = x_spring + x_spring_error - x_misalign_slope*(Zlength/2+rotation_offset)
    x_spring_load = x_spring + x_spring_error + x_misalign_slope*(Zlength/2+rotation_offset)
    if x_misalign_slope != 0.0:
        print(x_misalign_slope*(Zlength/2+rotation_offset))
    if "K" in x_i:
        K = x_i["K"]
    elif "log_K" in x_i:
        K = np.exp(x_i["log_K"])
    else:
        K = 10.0
    if "log_K_rig" in x_i:
        K_rig = np.exp(x_i["log_K_rig"])
    else:
        K_rig = 1.0e7
        print(K_trans)
    if "log_K_ground" in x_i:
        K_ground = np.exp(x_i["log_K_ground"])
        print(K_ground)
    else:
        K_ground = []
    
    # write Abaqus input file using the appropriate function, depending on whether or not the mesh is 
    # comprised of shells
    if shell_mesh:
        write_shell_inp(E11=E11, E22=E22, nu12=nu12, nu23=nu23, G12=G12, t_ply=t_ply, K=K, K_rig=K_rig, K_ground=K_ground, x_spring_fix=x_spring_fix, x_spring_load=x_spring_load, load=load, displacement = applied_disp, rotation_offset = rotation_offset+pivot_offset_error, StackSeq = layup, init_inc=init_inc, min_inc=min_inc, max_inc=max_inc, apply_load = apply_load, fix_to_ground = fix_to_ground)
    else:
        write_inp(E11=E11, E22 = E22, nu12 = nu12, nu23 = nu23, G12 = G12, K=K, x_spring=x_spring, load=load, rotation_offset = rotation_offset+pivot_offset_error, init_inc=init_inc, min_inc=min_inc, max_inc=max_inc)
    # Run Abaqus from the command line
    command = "Abaqus Job=" + file_str + " input=\"Abaqus\\" + file_str + ".inp\" interactive ask_delete=OFF cpus=2"
    print(command)
    os.system(command)
    # Process the abaqus output to extract displacements and nodal coordinates
    if change_inc:
        # If multiple frames are requested, run the code which processes all frames
        # command = "abaqus viewer noGUI=process_outputs_all_frames.py -- " + file_str
        command = "abaqus viewer noGUI=process_outputs_all_frames_new.py -- " + str(fix_to_ground) + " " + file_str
    else:
        # If only one increment is required just run the script which processes the final frame
        command = "abaqus viewer noGUI=process_outputs.py -- " + file_str
    
    os.system(command)

    # The increment index outputted by the postprocessing file should help keep track of things.
    # At some point consider looking at numpy options for 3d arrays, and json files. This is going to be an 
    # issue for the experimental data as well so worth getting a decent structure set downl
    # disp_i = np.loadtxt(model_name + "_displacement.csv", delimiter = ",", skiprows = 0)
    disp_i = np.loadtxt(file_str + "_displacement.csv", delimiter = ",", skiprows = 0)
    # nodes_i = np.loadtxt(model_name + "_nodes.csv", delimiter = ",", skiprows = 0)
    if change_inc:
        displacements = np.concatenate((displacements,disp_i), axis=1)
        # Also load in increment indices and reaction forces
        #incs_i = np.loadtxt(model_name + "_increments.csv", delimiter = ",", skiprows = 0, ndmin = 1)
        incs_i = np.loadtxt(file_str + "_increments.csv", delimiter = ",", skiprows = 0, ndmin = 1)
        incs.append(incs_i.tolist())
        # RFs_i = np.loadtxt(model_name + "_RFs.csv", delimiter = ",", skiprows = 0)
        RFs_i = np.loadtxt(file_str + "_RFs.csv", delimiter = ",", skiprows = 0)
        RFs = np.concatenate((RFs, RFs_i.reshape((1, -1))), axis=1)
        # Load in json file if using
        # with open(model_name + "_output.json", "r") as f:
        with open(file_str + "_output.json", "r") as f:
            # Load in string from file
            sample_dict = json.loads(f.readline())
        
        # Should follow the structure of the dictionary defined in the process outputs file
        # If the below code works fine then it should be possible just to add this to another
        # dictionary defined at this level. Then output this file. Could even output at each increment
        # then overwrite to get rid of the complicated buffer code
        out_dict["Sample"].append(sample_dict)
    else:
        displacements[:,3*i:3*(i+1)] = disp_i

    # Store nodes???
    # nodes[:,3*i:3*(i+1)] = nodes_i
    
# Consider using json format, stick with csv for now. This should be relatively easy if working with python dictionaries

# Write model outputs from all simulations to csv files   
# Create header for csv
if change_inc:
    # issue if some runs haven't started as then the stored increment list is a single float and not a list. I've added ndmin to the code which adds to this list
    # in the hope that it fixes this issue
    head_str = ', '.join(['u_' + str(i+1) + '_' + str(j+1) + ', v_' + str(i+1) + '_' + str(j+1) + ', w_' + str(i+1) + '_' + str(j+1) for i in range(N) for j in range(len(incs[i]))])
else:
    head_str = ', '.join(['u_' + str(i+1)+', v_'+str(i+1)+', w_'+str(i+1) for i in range(N)])

if apply_load:
    np.savetxt(infile + "_displacements_load=" + str(load) + "_max_inc=" + str(max_inc) + ".csv", displacements, delimiter=",", header = head_str, comments = "")
else:
    np.savetxt(infile + "_displacements_disp=" + str(applied_disp) + "_max_inc=" + str(max_inc) + ".csv", displacements, delimiter=",", header = head_str, comments = "")
# head_str = ', '.join(['x_' + str(i+1)+', y_'+str(i+1)+', z_'+str(i+1) for i in range(N)])
# np.savetxt(infile + "_nodes_load=" + str(load) + "_max_inc=" + str(max_inc) + ".csv", nodes, delimiter=",", header = head_str, comments = "")

if change_inc:
    if apply_load:
        with open(infile + "_incs_load=" + str(load) + "_max_inc=" + str(max_inc) + ".txt", 'w') as f:
            [f.write(', '.join([str(inc) for inc in sam])+"\n") for sam in incs]
    else:
        with open(infile + "_incs_disp=" + str(applied_disp) + "_max_inc=" + str(max_inc) + ".txt", 'w') as f:
            [f.write(', '.join([str(inc) for inc in sam])+"\n") for sam in incs]
    head_str = ', '.join(['Rx_' + str(i+1) + '_' + str(j+1) + ', Ry_' + str(i+1) + '_' + str(j+1) + ', Rz_' + str(i+1) + '_' + str(j+1) for i in range(N) for j in range(len(incs[i]))])
    if apply_load:
        np.savetxt(infile + "_RFs_load=" + str(load) + "_max_inc=" + str(max_inc) + ".csv", RFs, delimiter=",", header = head_str, comments = "")
    else:
        np.savetxt(infile + "_RFs_disp=" + str(applied_disp) + "_max_inc=" + str(max_inc) + ".csv", RFs, delimiter=",", header = head_str, comments = "")

# Write json file
if apply_load:
    with open(infile + "_output_struct.json",'w') as f:
        f.write(json.dumps(out_dict))
else:
    with open(infile + "_output_struct.json",'w') as f:
        f.write(json.dumps(out_dict))