import sys
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import os.path

header_dir = "C:\\Users\\cs2361\\Documents\\Bayesian_Model_Calibration\\source" # directory of sampling headers
sys.path.append(header_dir) # add header directory to path
from LHS_Design import transformed_LHS  # Import Latin Hypercube module

# Sample prior distributions using Latin Hypercube Sampling. Samples are weighted to encourage a space-filling design across intervals of equal probability
# Define a list of inputs upon which priors are defined, with entries composed of a list of name, distribution types and
# parameters (mostly taken from NIAR AS4/8552 statistical report unless otherwise specified)
# Units kN, mm
#inputs = [
#    ['E11', 'Gaussian', 115.6, 6.0, True],    # NIAR data, E11c, RTD 
#    ['E22', 'Gaussian', 9.24, 6.0, True],     # NIAR data, E22t, RTD
#    ['nu12','Gaussian', 0.335, 12.123, True], # NIAR data, nu_12,c RTD. Not adjusting CoV as the observed value is quite high. Also consider using tensile measurement but adjusting up to 6% 
#    ['nu23','Gaussian', 0.487, 12.0, True],   # "Modelling and simulation methodology for unidirectional composite laminates in a Virtual Test Lab framework". Large CoV assumed as unspecified in the paper. Consider using smaller CoV
#    ['G12', 'Gaussian', 4.826, 6.0, True],    # NIAR data, G12,s RTD.
#    ['t_ply','Gaussian', 0.196, 5.0, True]     # Mean taken from past data of material used in Bath, mean chosen from Engineering judgement of a conservative but reasonable spread in values
#    ]

# 3D Example
#inputs = [
#    ['E11', 'Gaussian', 115.6, 6.0, True],    # NIAR data, E11c, RTD 
#    ['t_ply', 'Gaussian', 0.196, 5.0, True],  # Mean taken from past data of material used in Bath, mean chosen from Engineering judgement of a conservative but reasonable spread in values
#    ['K', 'Loguniform', 100.0, 1.0e9, True]   # Bounds taken from rough paramteric study to check for switch from Simply-Supported to clamped behaviour
#    # ['K', 'Lognormal', 11.0, 3.25, True]   # Parameters taken from rough paramteric study to check for switch from Simply-Supported to clamped behaviour. Bit of a fudge
#    ]

# 7D Example
#inputs = [
#    ['E11', 'Gaussian', 115.6, 6.0, True],    # NIAR data, E11c, RTD 
#    ['E22', 'Gaussian', 9.24, 6.0, True],     # NIAR data, E22t, RTD
#    ['nu12','Gaussian', 0.335, 12.123, True], # NIAR data, nu_12,c RTD. Not adjusting CoV as the observed value is quite high. Also consider using tensile measurement but adjusting up to 6% 
#    ['nu23','Gaussian', 0.487, 12.0, True],   # "Modelling and simulation methodology for unidirectional composite laminates in a Virtual Test Lab framework". Large CoV assumed as unspecified in the paper. Consider using smaller CoV
#    ['G12', 'Gaussian', 4.826, 6.0, True],    # NIAR data, G12,s RTD.
#    ['t_ply','Gaussian', 0.196, 5.0, True],     # Mean taken from past data of material used in Bath, mean chosen from Engineering judgement of a conservative but reasonable spread in values
#    ['K', 'Loguniform', 100.0, 1.0e9, True]   # Bounds taken from rough paramteric study to check for switch from Simply-Supported to clamped behaviour
#]

# 8D Example 
#inputs = [
#    ['E11', 'Gaussian', 115.6, 6.0, True],    # NIAR data, E11c, RTD 
#    ['E22', 'Gaussian', 9.24, 6.0, True],     # NIAR data, E22t, RTD
#    ['nu12','Gaussian', 0.335, 12.123, True], # NIAR data, nu_12,c RTD. Not adjusting CoV as the observed value is quite high. Also consider using tensile measurement but adjusting up to 6% 
#    ['nu23','Gaussian', 0.487, 12.0, True],   # "Modelling and simulation methodology for unidirectional composite laminates in a Virtual Test Lab framework". Large CoV assumed as unspecified in the paper. Consider using smaller CoV
#    ['G12', 'Gaussian', 4.826, 6.0, True],    # NIAR data, G12,s RTD.
#    ['t_ply','Gaussian', 0.196, 5.0, True],     # Mean taken from past data of material used in Bath, mean chosen from Engineering judgement of a conservative but reasonable spread in values
#    ['K', 'Loguniform', 100.0, 1.0e9, True],   # Bounds taken from rough paramteric study to check for switch from Simply-Supported to clamped behaviour
#    ['x_spring','Uniform', -10.0, 65.0, True]
#]

# 4D Example
# inputs = [
#    ['E11', 'Gaussian', 115.6, 6.0, True],    # NIAR data, E11c, RTD 
#    ['t_ply','Gaussian', 0.196, 5.0, True],     # Mean taken from past data of material used in Bath, mean chosen from Engineering judgement of a conservative but reasonable spread in values
#    ['K', 'Loguniform', 100.0, 1.0e9, True],   # Bounds taken from rough paramteric study to check for switch from Simply-Supported to clamped behaviour
#    ['x_spring','Uniform', -10.0, 65.0, True]
#]

# Another 4D Example
# inputs = [
#    ['E11', 'Gaussian', 115.6, 6.0, True], # NIAR data, E11c, RTD 
#    ['t_ply','Gaussian', 0.196, 5.0, True],     # Mean taken from past data of material used in Bath, mean chosen from Engineering judgement of a conservative but reasonable spread in values
#    ['LFlange_theta', 'Gaussian', 0.0, 4.0, True],# +/- 3 standard deviations within 12 degrees
#    ['RFlange_theta', 'Gaussian', 0.0, 4.0, True] # +/- 3 standard deviations within 12 degrees
#]

# Another 7D Example
# inputs = [
#    ['E11', 'Gaussian', 115.6, 6.0, True], # NIAR data, E11c, RTD 
#    ['E22', 'Gaussian', 9.24, 6.0, True],  # NIAR data, E22t, RTD
#    ['nu12','Gaussian', 0.335, 12.123, True], # NIAR data, nu_12,c RTD. Not adjusting CoV as the observed value is quite high. Also consider using tensile measurement but adjusting up to 6% 
#    ['G12', 'Gaussian', 4.826, 6.0, True],    # NIAR data, G12,s RTD.
#    ['t_ply','Gaussian', 0.196, 5.0, True],     # Mean taken from past data of material used in Bath, mean chosen from Engineering judgement of a conservative but reasonable spread in values
#    ['LFlange_theta', 'Gaussian', 0.0, 4.0, True], # +/- 3 standard deviations within 12 degrees
#    ['RFlange_theta', 'Gaussian', 0.0, 4.0, True] # +/- 3 standard deviations within 12 degrees
#    ]

# Full example
inputs = [
    ['E11', 'Gaussian', 115.6, 6.0, True], # NIAR data, E11c, RTD
    ['E22', 'Gaussian', 9.24, 6.0, True],  # NIAR data, E22t, RTD
    ['nu12','Gaussian', 0.335, 12.123, True], # NIAR data, nu_12,c RTD. Not adjusting CoV as the observed value is quite high. Also consider using tensile measurement but adjusting up to 6% 
    ['nu23','Gaussian', 0.487, 12.0, True],   # "Modelling and simulation methodology for unidirectional composite laminates in a Virtual Test Lab framework". Large CoV assumed as unspecified in the paper. Consider using smaller CoV
    ['G12', 'Gaussian', 4.826, 6.0, True],    # NIAR data, G12,s RTD.
    ['t_ply','Gaussian', 0.196, 4.0, True],   # Mean taken from past data of material used in Bath, mean chosen from Engineering judgement of a conservative but reasonable spread in values
    ['LFlange_theta', 'Gaussian', 0.0, 5.0/3.0, True], # +/- 3 standard deviations within 5 degrees
    ['RFlange_theta', 'Gaussian', 0.0, 5.0/3.0, True], # +/- 3 standard deviations within 5 degrees
    ['K', 'Lognormal', 16.0, 1.0, True],   # Bounds taken from rough paramteric study to check for switch from Simply-Supported to clamped behaviour
    ]

# Example for new spar (IM7 properties)
# Simple example assuming no misalignment in x direction
inputs = [
    ['E11', 'Gaussian', 140.9, 6.0, True], # NIAR, RTD
    #['G12', 'Gaussian', 4.69, 6.0, True], # NIAR, RTD
    ['t_ply', 'Gaussian', 0.125, 4.0, True], #  Meng Yi's Refs
    ['K', 'Loguniform', 10.0, 25000.0, True] # Initial test case - switch probably to log-gamma once I've had time to think about it
    #['x_spring_error', 'Uniform', -1.0, 1.0, True]  # Halfway between two settings
    #['pivot_offset_error', 'Uniform', -5.0, 5.0, True] # Rig is precision machined, pick this value to incorporate uncertainty in trimmed spar length, and amount embedded in end caps. Basically this changes the effective length
]

# inputs = [
#    ['E11', 'Gaussian', 140.9, 6.0, True], # NIAR, RTD
#    ['G12', 'Gaussian', 4.69, 6.0, True], # NIAR, RTD
#    ['t_ply', 'Gaussian', 0.125, 4.0, True], #  Meng Yi's Refs
#    ['K', 'Loggamma', 1.0, 0.3, True], # Initial test case - switch probably to log-gamma once I've had time to think about it
#    ['x_spring_error', 'Uniform', -1.0, 1.0, True],  # Halfway between two settings
#    ['pivot_offset_error', 'Uniform', -5.0, 5.0, True] # Rig is precision machined, pick this value to incorporate uncertainty in trimmed spar length, and amount embedded in end caps. Basically this changes the effective length
#]

# Alternative with higher modulus
# inputs = [
#    ['E11', 'Gaussian', 150.0, 6.0, True], # NIAR, RTD
#    ['G12', 'Gaussian', 4.69, 6.0, True], # NIAR, RTD
#    ['t_ply', 'Gaussian', 0.125, 4.0, True], #  Meng Yi's Refs
#    ['K', 'Loggamma', 1.0, 0.3, True], # Initial test case - switch probably to log-gamma once I've had time to think about it
#    ['x_spring_error', 'Uniform', -1.0, 1.0, True],  # Halfway between two settings
#    ['pivot_offset_error', 'Uniform', -5.0, 5.0, True] # Rig is precision machined, pick this value to incorporate uncertainty in trimmed spar length, and amount embedded in end caps. Basically this changes the effective length
#]

inputs = [
    ['E11', 'Gaussian', 140.9, 6.0, True], # NIAR, RTD
    #['G12', 'Gaussian', 4.69, 6.0, True], # NIAR, RTD
    ['t_ply', 'Gaussian', 0.125, 4.0, True], #  Meng Yi's Refs
    ['K', 'Loggamma', 1.0, 0.3, True] # Initial test case - switch probably to log-gamma once I've had time to think about it
    #['x_spring_error', 'Uniform', -1.0, 1.0, True]  # Halfway between two settings
    #['pivot_offset_error', 'Uniform', -5.0, 5.0, True] # Rig is precision machined, pick this value to incorporate uncertainty in trimmed spar length, and amount embedded in end caps. Basically this changes the effective length
]


# Write code for outputting prior info
inputs = [
    ['E11', 'Gaussian', 140.9, 6.0, True], # NIAR, RTD
    ['G12', 'Gaussian', 4.69, 6.0, True], # NIAR, RTD
    ['t_ply', 'Gaussian', 0.125, 4.0, True], #  Meng Yi's Refs
    ['K', 'Loggamma', 1.0, 0.3, True], # Initial test case - switch probably to log-gamma once I've had time to think about it
    ['x_spring_error', 'Uniform', -1.0, 1.0, True],  # Halfway between two settings
    ['pivot_offset_error', 'Uniform', -5.0, 5.0, True], # Rig is precision machined, pick this value to incorporate uncertainty in trimmed spar length, and amount embedded in end caps. Basically this changes the effective length
    ['x_misalign_slope', 'Gaussian', 0.0, 1.0/370.0, True] # Such that one standard deviation results in a misalignment of 1mm at each end (of the spar)
]

inputs = [
    ['E11', 'Gaussian', 140.9, 6.0, True], # NIAR, RTD
    ['G12', 'Gaussian', 4.69, 6.0, True], # NIAR, RTD
    ['t_ply', 'Gaussian', 0.125, 4.0, True], #  Meng Yi's Refs
    ['x_spring_error', 'Uniform', -1.0, 1.0, True],  # Halfway between two settings
    ['Flange_theta', 'Gaussian', 0.0, 1.0, True]
]



# What about misalignment in y direction?

N = 50 # Number of samples to be generated
xLHS = transformed_LHS(inputs, N, sampler_package="scikit-optimize", sampler_kwargs={"lhs_type":"classic","criterion":"maximin", "iterations":10000})
d = len(inputs) # Number of inputs

# If a log distribution, output the natural log as this is a more natural scale for the Gaussian process
# This would be less messy if using Pandas, or if Priors object were returned/computed separately
for i, input in enumerate(inputs):
    if input[1] == "Lognormal" or input[1] == "Loguniform":
        input[0] = "log_" + input[0]
        xLHS[:,i] = np.log(xLHS[:,i])

# At least convert this bit to Pandas next...
# Write Design of Experiments to a csv file
head_string = ','.join([input[0] for input in inputs])

# Check if output file exists to avoid over-writing
filename = "LHSDesign" + str(N) + "x" + str(d) + ".csv"
i = 1 # counter to be appended to filename if it already exists
# Check if the filename is already in use, and if so append an integer to filename until a non-existing file is found
while os.path.exists(filename):
    filename = "LHSDesign" + str(N) + "x" + str(d) + "_" + str(i) + ".csv"
    i = i + 1

np.savetxt(filename, xLHS, delimiter=",", header = head_string, comments = "")
# When using as in a top-level file it might be nice for each column to be attributed to a named variable. 
# Don't bother for now as this will probably be easier when I know how to use Pandas
# Consider scikit-opt for optimised Latin Hypercubes
# Write the calibration parameters file
for input in inputs:
    if input[1] == "Loguniform":
        input[2] = np.log(input[2])
        input[3] = np.log(input[3])
    elif input[1] == "Gaussian" and input[2]!= 0.0:
        input[3] = input[2]*input[3]/100.0

print(inputs)

out_frame = pd.DataFrame({"distribution" : [input[1] for input in inputs], 
 "param_1" : [input[2] for input in inputs],
 "param_2" : [input[3] for input in inputs]}, 
 index = [input[0] for input in inputs])

out_frame.to_csv(filename.strip(".csv")+"_tf_param.csv")