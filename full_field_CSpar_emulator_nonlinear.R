
# Fit an emulator to the full-field output from a C-Spar Finite Element model
# Input data are the Longitudinal Modulus E11, ply thickness, and the torsional
# stiffness of a spring representing an uncertain boundary condition
# Output data are the axial displacements of the nodes of the FE model across
# a range of load steps of a nonlinear model
# Follows emulator aspect of D. Higdon et al, "Computer Model Calibration Using 
# High-Dimensional Output",Journal of the American Statistical Association,2008.
# This code handles pre and post processing.

library(data.table)
library(rstan)
library(matrixStats)
library(rjson)

# Set current working directory. This should be modified to match the directory
# of the user
setwd("C:/Users/cs2361/Documents/CSpar_Calibration/")

# include functions which are called in this code
source("source/estimate_mode.R")
source("source/covariance_matrices.R")
source("source/abaqus_json.R")

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation

p_eta = 5 # Number of basis functions to be retained for the emulator from SVD
disp_str = "w" # String which identifies the displacement component of interest (u,v, or w)

#-------------------------------------------------------------------------------

# Set up simulation data

# Load in emulator training data input values from Design of Experiments. 
in_file = "LHSDesign40x4" # File identifier string which is in both input and output files
XT_sim = fread(paste("inputs/",in_file,".csv", sep = ""))

# Determine useful quantities from model inputs and outputs. Variable names 
# match the notation of Higdon et al. 2008
q = ncol(XT_sim)          # number of calibration inputs
tc = as.matrix(XT_sim)    # Convert training data input points to a matrix for passing to stan
m = nrow(XT_sim)          # sample size of computer simulation data

# Load in  training data output displacement values from Abaqus. Here I've used
# a similar naming convention to the inputs to automate changes. The file name 
# can be changed manually if need be.
# Here I've stored outputs structured in a json across samples and load increments
abaqus_displacements <- fromJSON(file = paste("inputs/",in_file,"_output_struct.json", sep=""))
# DataFrame might be a more intuitive format for R...
# Extract displacements from json, and store as matrix where each column is a 
# training sample, and the rows are the concatenation of displacement output
# across all output frames
sorted_data = extract_const_frame(abaqus_displacements,"w")
dt_simulation = sorted_data[[1]]
n_nodes = sorted_data[[2]] # Number of nodes
n_frames = sorted_data[[3]] # Number of frames
n_eta = nrow(dt_simulation) # total number of output points per simulation

#-------------------------------------------------------------------------------