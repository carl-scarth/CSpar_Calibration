import mogp_emulator as mogp
from scipy.stats import norm, uniform, loguniform, lognorm
from scipy.stats.qmc import LatinHypercube

import numpy as np
import os.path

# Consider creating a folder called utils and importing the below classes separately as Jean does

class Prior:
# Define class of prior distribution objects containing type of distribution all important statistical attributes
    def __init__(self,name, distribution, param_1, param_2, weight_LHS):
        # Initialise attributes of prior based upon the type of distribution
        self.name = name
        self.distribution = distribution
        self.weight_LHS = weight_LHS
        if (self.distribution == 'Gaussian'):
            # If Gaussian, the prior is described by a mean and Coefficient of Variation
            self.mu = param_1
            self.COV = param_2
        elif (self.distribution == 'Lognormal'):
            self.mu = param_1
            self.s = param_2
        elif (self.distribution == 'Uniform') or (self.distribution == 'Loguniform'):
            # If Uniform or loguniform, the prior is described by a lower and upper bound
            self.lb = param_1
            self.ub = param_2

    @property
    def sigma(self):
        if self.distribution == 'Gaussian':
            if self.mu != 0:
                return self.mu*self.COV/100.0
            else:
                # If the input has zero mean, then the standard deviation is specified directly rather than via the COV
                return self.COV
        else:
            return None
            # Shouldn't be used, but if necessary could actually specify this
        
    @property
    def min(self):
        if self.distribution == 'Gaussian':
            return self.mu - 3.0*self.sigma
        elif self.distribution == 'Lognormal':
            # Defined by +/- 3 standard deviations of log(x)
            return np.exp(self.mu - 3.0*self.s)
        elif (self.distribution == 'Uniform') or (self.distribution == 'Loguniform'):
            return self.lb
        
    @property
    def max(self):
        if self.distribution == 'Gaussian':
            return self.mu + 3.0*self.sigma
        elif self.distribution == 'Lognormal':
            return np.exp(self.mu + 3.0*self.s)
        elif (self.distribution == 'Uniform') or (self.distribution == 'Loguniform'):
            return self.ub

# Sample prior distributions using Latin Hypercube Sampling. Samples are weighted to encourage a space-filling design across intervals of equal probability
# Eventually we'd define the inputs via a separate file and run the DoE as a function.

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
inputs = [
    ['E11', 'Gaussian', 115.6, 6.0, True], # NIAR data, E11c, RTD 
    ['t_ply','Gaussian', 0.196, 5.0, True],     # Mean taken from past data of material used in Bath, mean chosen from Engineering judgement of a conservative but reasonable spread in values
    ['LFlange_theta', 'Gaussian', 0.0, 4.0, True],# +/- 3 standard deviations within 12 degrees
    ['RFlange_theta', 'Gaussian', 0.0, 4.0, True] # +/- 3 standard deviations within 12 degrees
]

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

# Also define a log-normal option as may be more appropriate for K given the sudden nature of the switch

N = 75 # Number of samples to be generated
d = len(inputs) # Number of inputs

# Create a list of prior objects containing details of priors for each input
priors = []
for item in inputs:
    priors.append(Prior(item[0],item[1],item[2],item[3],item[4]))      

# Create a Latin Hypercube Sample (LHS) object defined on the unit hypercube [0, 1]^d
LHSobj = mogp.LatinHypercubeDesign([(0.0, 1.0) for entry in inputs])
# Generate a set of N samples using the LHS object
FLHS = LHSobj.sample(N)

# Loop over each of the columns of the Design of Experiments and convert from interval [0, 1] onto the 
# actual input space using the appropriate PPF (inverse of the CDF) for each variable
# Maybe tidy later to get rid of conditional statement. A good way of doing it would be to encode the 
# 'distribution' attribute consistent with the scipy objects or even renaming these during the import
xLHS = np.empty([N, d])
for i, item in enumerate(priors):
    if not item.weight_LHS:
        xLHS[:,i] = [uniform.ppf(xij, loc = item.min, scale = (item.max - item.min)) for xij in FLHS[:,i]]
    else:
        if item.distribution == 'Gaussian':
            xLHS[:,i] = [norm.ppf(xij, loc = item.mu, scale = item.sigma) for xij in FLHS[:,i]]
        elif item.distribution == 'Lognormal':
            # scipy respresents lognormal distribution in an unusual way
            xLHS[:,i] = [lognorm.ppf(xij, item.s, loc = 0.0, scale = np.exp(item.mu)) for xij in FLHS[:,i]]
        elif item.distribution == 'Uniform':
            # Note, the scipy uniform distribution is bounded by [loc, loc + scale]
            xLHS[:,i] = [uniform.ppf(xij, loc = item.lb, scale = (item.ub - item.lb)) for xij in FLHS[:,i]]
        elif item.distribution == 'Loguniform':
            xLHS[:,i] = [loguniform.ppf(xij, item.lb, item.ub) for xij in FLHS[:,i]]

# Write Design of Experiments to a csv file
head_string = ','.join([item.name for item in priors])

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