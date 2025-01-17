import os
import pandas as pd
import numpy as np
import json
import sys
import matplotlib.pyplot as plt

src_path = "source"
sys.path.insert(0, src_path)
sys.path.insert(0, src_path)

from abaqus_json import extract_const_frame, basis_mean_to_json
from transform_input_output import normalise_inputs, standardise_vector_output, rescale_output, rescale_inputs
from dimension_reduction import svd_basis

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation

in_file = "LHSDesign100x8" # File identifier string for input and output files
pred_file = "LHSDesign10x8" # File identifier for input values where predictions are required
in_suffix = "150kN_zeroed" # Other text added to the file identifier

p_eta = 20 # Number of basis functions retained for the emulator from SVD
exp_tol = 1e-4 # Tolerance variance fraction used to assess SVD convergence
disp_str = ["u", "w"] # String which identifies the displacement component of interest (u,v, or w) or list of multiple
#disp_str = ["u"] # String which identifies the displacement component of interest (u,v, or w) or list of multiple
q_y = len(disp_str)
print_svd_output = True # Print diagnostic output of svd to the terminal?
export_modes = True # Calculate modes and means of emulator hyperparameters and write to file?
output_json = True # Output structured data via json
output_vtk = True # Output plottable data via vtk

# Define parameters of the gamma prior on the error associated with truncating
# the series expansion for the model output
a_eta = 1.0     # Shape parameter for the lambda_eta prior
b_eta = 0.0001  # Rate parameter for the lambda_eta prior
#iter = 5000 # Number of samples per chain
#chains = 3 # Number of chains for simulation

#-------------------------------------------------------------------------------

# Set up simulation data

# Load in emulator training data input values from Design of Experiments. 
XT_sim = pd.read_csv(os.path.join("inputs",in_file+".csv"))

# Load in test points at which predictions are required
XT_pred = pd.read_csv(os.path.join("inputs",pred_file+".csv"))
n_pred = XT_pred.shape[0] # number of predictions

# Determine useful quantities from model inputs and outputs. Variable names 
# match the notation of Higdon et al. 2008
m, q = XT_sim.shape # sample size of simulation data, and number of calibration inputs

# Load in training data output displacement values from Abaqus.
# Outputs are structured in a json across samples and load increments
with open(os.path.join("inputs", "_".join((in_file,"output_struct",in_suffix))+'.json'),'r') as f:
    abaqus_dict = json.loads(f.readline())

dt_simulation, n_nodes, n_frames = extract_const_frame(abaqus_dict, disp_str)
n_eta = dt_simulation.shape[0] # total number of output points per simulation

#-------------------------------------------------------------------------------

# Standardise the data

# Normalise inputs such that training data is on the unit hypercube
tc, t_min, t_max = normalise_inputs(XT_sim)
# Normalise test data in the same way as the training data for consistency
t_pred = normalise_inputs(XT_pred, x_min=t_min, x_max=t_max)

# Standardise the outputs to have zero mean (for each row) and unit standard
# deviation (for each displacement)
eta, mu_dt, sd_dt = standardise_vector_output(dt_simulation, q_y = q_y)

#-------------------------------------------------------------------------------

# Perform dimension reduction on data via SVD
K_eta, p_eta = svd_basis(eta, exp_tol = exp_tol, print_output = print_svd_output, csv_label = "nonlinear_"+in_file)

# Output training data mean and basis vectors to json if required
if output_json:
    # Write basis to .json file
    out_dict = basis_mean_to_json(n_frames, n_nodes, K_eta, mu_dt, disp_str)
#    print(out_dict.keys())
#    ASDSAD
#   THIS WAS CREATING A JSON FROM A STRING BEFORE AS I CONVERTED TO JSON TWICE. NEED TO OPEN THIS AND TEST
# ALSO TEST THE VTK PRODUCED USING THE JSON LOOKS OK
# WILL THEN USE THE SAME CODE TO GET THE DICTIONARY IN THE CORRECT FORMAT, THEN PASS
# DIRECT TO THE VTK WRITING CODE WHICH I'LL TURN INTO A HEADER
    with open(os.path.join("outputs","basis_nonlinear_"+in_file+".json"),'w') as f:
        f.write(json.dumps(out_dict))

# ALSO OPTION TO OUTPUT VTK


# Thoughts Add standardistation of the output data here as well. 
# Then do pairs plots after dimension reductions

## plot Design of Experiments and test points
#pairs(tc,col="blue", pch=4, 
#      main = "Normalised Training Data Input values",
##      cex=1.5,
#      cex.labels = 1.75,
#      cex.axis=1.5,
#      cex.lab=1.5)
#pairs(t_pred,col = "blue", pch=4, 
#      main = "Normalised Test Data Input Values",
#      cex=1.5,
#      cex.labels = 1.75,
#      cex.axis=1.5,
#      cex.lab=1.5)
