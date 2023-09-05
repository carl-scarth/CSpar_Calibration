
# This code is an application of the multivariate calibration formulation 
# proposed in 'Computer Model Calibration using High-dimensional Output', by
# Higdon et al, JASA, to a CSpar finite element model with uncertain inputs 
# using DIC data from one experiment. A simplified version of Higdon et al. 
# with n = 1 is implemented here. Sampling is undertaken in stan. This code 
# handles pre and post processing.

library(data.table)
library(rstan)
library(matrixStats)
# library(MASS)
# library(colormap)

# Set current working directory. This should be modified to match the directory
# of the user
setwd("C:/Users/cs2361/Documents/CSpar_Calibration/")

# include functions which are called in this code
source("source/transform_input_output.R")
# source("source/utils.R")
# source("source/dimension_reduction.R")
# source("source/prior_posterior_plots.R")
# source("source/gp_predictions.R")

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation
# Might be able to delete some of these later
p_eta = 7 # Number of basis functions retained for the emulator from SVD
# exp_tol = 1e-6 # Tolerance variance fraction used to assess SVD convergence
disp_str = "w" # String which identifies the displacement component of interest (u,v, or w)
# Define parameters of the gamma prior on the error associated with truncating
# the series expansion for the model output
# a_eta = 1.0     # Shape parameter for the lambda_eta prior
# b_eta = 0.0001  # Rate parameter for the lambda_eta prior 
iter = 4000 # Number of samples per chain
chains = 3 # Number of chains for simulation
# print_svd_output = TRUE # Print diagnostic output of svd to the terminal?

p_eta = 5 # Number of basis functions to be retained for the emulator from SVD

#-------------------------------------------------------------------------------

# Define prior distribution parameters for passing to stan

# Pre-processing for BC example (Mean and coefficients of variation for Gaussian
# inputs)
E11_mu = 115.6
t_ply_mu = 0.196
E11_cov = 6.0 
t_ply_cov = 5.0
# Bounds for log-uniform inputs
#K_lb = 100.0
#K_ub = 1.0e9

# Define data_frame of prior parameters (this could be done via csv?)
# tf_param <- data.frame(distribution = c("Gaussian","Gaussian","Loguniform"),
#                       param_1 = c(E11_mu, t_ply_mu, log(K_lb)),
#                       param_2 = c(E11_mu*E11_cov/100, t_ply_mu*t_ply_cov/100, log(K_ub)))
# row.names(tf_param) <- c("E11","t_ply","log_K")

# Additional pre-processing for flange-rotation example
flange_theta_mu = 0.0
flange_theta_sigma = 4.0
# Define data_frame of prior parameters
tf_param <- data.frame(distribution = c("Gaussian","Gaussian","Gaussian","Gaussian"),
                      param_1 = c(E11_mu, t_ply_mu, flange_theta_mu, flange_theta_mu ),
                      param_2 = c(E11_mu*E11_cov/100, t_ply_mu*t_ply_cov/100, flange_theta_sigma, flange_theta_sigma))
row.names(tf_param) <- c("E11","t_ply","LFlange_theta","RFlange_theta")

#-------------------------------------------------------------------------------

# Set up simulation data (need to retain this for normalisation of inputs)

# Load in emulator training data input values from Design of Experiments. 
# in_file = "LHSDesign50x3" # File identifier string for input and output csvs
in_file = "LHSDesign40x4" # File identifier string for input and output csvs
XT_sim = fread(paste("inputs/",in_file,".csv", sep = ""))

# In this example I fit the emulator to the log of spring stiffness K, which is 
# a more natural choice of values across which outputs are expected for
# variations in this input
# XT_sim$K = log(XT_sim$K)
# colnames(XT_sim)[3] = "log_K"

# Determine useful quantities from model inputs and outputs. Variable names 
# match the notation of Higdon et al. 2008
q = ncol(XT_sim)          # number of calibration inputs
tc = as.matrix(XT_sim)    # Convert to a matrix for passing to stan
m = nrow(XT_sim)          # sample size of computer simulation data

# Each row of XT_sim corresponds to a block of three columns of displacement 
# data, with a column for each component u,v,w 
# abaqus_displacements = fread(paste("inputs/",in_file,"_displacements.csv", sep=""))
abaqus_displacements = fread(paste("inputs/",in_file,"_fixed_100kN.csv", sep=""))
n_eta = nrow(abaqus_displacements) # number of output points per simulation

# Extract the displacement for the component of interest and store in a matrix
dt_simulation = matrix(NA,nrow = n_eta, ncol = m)
for (i in 1:m){
  dt_simulation[,i] = abaqus_displacements[[as.name(paste(disp_str,sprintf("_%d", i),sep=""))]]
}

#-------------------------------------------------------------------------------

# All training data points are normalised onto the unit hypercube [0,1]^q before
# being passed to stan. The same transformation is applied to the test points

# Determine the maximum and minimum value of each input within the training data
t_min = colMins(tc)
t_max = colMaxs(tc)
# Normalise the inputs using these maximum and minimum values
tc = normalise_inputs(tc, t_min, t_max)

# For all implemented distributions the transformation of the 1st parameter is the same
tf_param$p1_trans = normalise_inputs(tf_param$param_1, t_min, t_max)
# The transformation of the 2nd parameter is different if this is a standard deviation
tf_param$p2_trans = NA
for (i in 1:q){
  # Second parameter of a Gaussian is a standard deviation
  if (tf_param$distribution[i] == "Gaussian"){
    tf_param$p2_trans[i] = normalise_inputs(tf_param$param_2[i], t_min[i], t_max[i], sd = TRUE)
  } else if ((tf_param$distribution[i] == "Uniform") | (tf_param$distribution[i] == "Loguniform")){
    tf_param$p2_trans[i] = normalise_inputs(tf_param$param_2[i], t_min[i], t_max[i])
  } else {
    stop("Error: Non-implemented distribution")
  }
}

# plot Design of Experiments and test points
pairs(tc,col="blue", pch=4, 
      main = "Normalised Training Data Input values",
      cex=1.5,
      cex.labels = 1.75,
      cex.axis=1.5,
      cex.lab=1.5)

#-------------------------------------------------------------------------------

# WORK FROM HERE!!! CONSIDER STORING USEFUL QUANTITIES, E.G. NORMALISATION 
# QUANTITIES, BASIS FUNCTIONS ETC IN A CSV FILE AS IN EMULATOR MODES TO REDUCE 
# DUPLICATION - REPEAT ABOVE FOR NONLINEAR CODE AND IN EXISTING EMULATOR CODE

# This portion of code performs SVD on the model data, which is standardised 
# such that the weight coefficients of the resulting expansion will have zero 
# mean and unit variance

# Centre the simulation output for each element. This guarantees that the weights (w)
# will have zero (sample) mean. Note that this same transformation must also be 
# applied to experimental data y, likely requiring some interpolation. 
# An alternative approach would be to centre the data using the overall mean, 
# calculating the sample mean of the transformed output, then specify this 
# (constant) mean value to each Gaussian process. This is arguably easier and 
# might reduce this possibility of errors being introduced.
mu_dt = rowMeans(dt_all_simulation)
dt_all_cen = sweep(dt_all_simulation,1,mu_dt,"-")

# Divide by the overall standard deviation of the simulation data to standardise
# to unit variance
sd_dt = sd(as.matrix(dt_all_cen))
dt_all_cen = dt_all_cen/sd_dt
# Store standardised simulation output as eta and convert to matrix to pass to stan
eta = as.matrix(dt_all_cen)

# Perform SVD on centred data
dt_svd = svd(dt_all_cen)

# Extract the first p_eta basis functions from the svd. 
# standardised such that the weights (columns of sqrt(m-1)*v) have unit variance
# K_eta = dt_svd$u[,1:p_eta]*matrix(rep(dt_svd$d[1:p_eta],2),nrow = n_eta, ncol = p_eta, byrow = TRUE)/sqrt(m-1)
K_eta = dt_svd$u[,1:p_eta]*matrix(dt_svd$d[1:p_eta],nrow = n_eta, ncol = p_eta, byrow = TRUE)/sqrt(m-1)
K_eta = dt_svd$u


# Sanity check that weights of SVD have zero mean and unit variance
print("mean of reduced dimension output w = ")
print(colMeans(dt_svd$v))
print("standard deviations of reduced dimension output w = ")
print(colSds(dt_svd$v*sqrt(m-1)))

# write basis functions to file for external plotting
write.csv(K_eta, "outputs/basis.csv", row.names = FALSE)
# Also write the mean vector
write.csv(mu_dt, "outputs/model_mean.csv", row.names = FALSE)

# Plot magnitude of singular value d with increasing number of basis functions. 
# This gives an indication of how much each base contributes to the output variance.
# I've used this plot to select which value of p_eta to choose with this method
d_r = dt_svd$d[1:p_eta]
d_r_norm = d_r^2/sum(d_r^2)

par(mfrow = c(1,2))
plot(1:p_eta,d_r_norm,"type"="p","col"="red","pch"=4,"lwd"=3,cex=1.5,'xlab' = "Feature",'ylab'="normalised d_i",cex.axis=1.75,cex.lab=1.75)
plot(2:p_eta,d_r_norm[-1],"type"="p","col"="red","pch"=4,"lwd"=3,cex=1.5,'xlab' = "Feature",'ylab'="normalised d_i",cex.axis=1.75,cex.lab=1.75)
# How many basis functions are required to get within a tolerance fraction of the
# variance
tol = 1e-6
cumsum(d_r_norm) > (1-tol)
basis_tol = 1:p_eta
basis_tol = basis_tol[cumsum(d_r_norm) > (1-tol)]
basis_tol = basis_tol[1]
print("Number of basis functions required to represent output within tolerance = ")
print(basis_tol)

#-------------------------------------------------------------------------------

# Load in the experimental data, and centre using the same method as used to 
# stanrdardise the model output. This requires some interpolation to determine the
# mean model output at each measurement location

# Read in DIC data. Each row has entries for which element the measurement has 
# been mapped to, and the natural coordinates for that point within the element
# experimental_data = read.table("inputs/Ext_LCorner_Image_0115_0.tiff_aligned_adjusted.csv", sep = ',', header = TRUE)
experimental_data = read.table("inputs/Ext_LCorner_Image_0115_0.tiff_nat_coord_rad_trim.csv", sep = ',', header = TRUE)
n_y = nrow(experimental_data)# Number of observations
y_element = experimental_data$Element
hr = as.matrix(experimental_data[,c("h","r")])
exp_displacement = experimental_data$W

# I also need to load mesh connectivity information about the outer surface of
# the spar. This must be consistent across all of the training data for this 
# method to work.
connectivity = read.table("inputs/outer_surface_elements.csv", sep = ',', header = TRUE)

# Define the position of the nodes in natural coordinates (this is 4,1,5,8) as
# this is what I used to determine the coordinates.
# It would probably make more sense to use standard ordering for a quad element,
# but I would need to re-do the mapping
HR = matrix(c(1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0),nrow=2,ncol=4,byrow = TRUE)

# For consistency with simulation output, the DIC data must be centred
# using with the mean model output at their locations. As this is  unknown 
# precisely, this must be approximated through interpolation. Here I use the
# same linear interpolation functions of the isoparametric elements used by
# Abaqus
mu_y = rep(NA,n_y)
for (i in 1:n_y) {
  # Might be clearer just to write out basis equations in full... 
  bases = 1.0 + matrix(hr[i,],2,4)*HR
  # The complete basis functions for the quad element are given by the product
  # of those in h and r, contained in the rows of "bases"
  # Note that the element index is in Python indexing convention, but connectivities are in the Abaqus convention
  # plus 2 as the first two nodes are the reference points, which are numbered
  # according to a different system as they are defined directly on the assembly
  mu_y[i] = sum(bases[1,]*bases[2,]*mu_dt[as.numeric(connectivity[y_element[i]+1,])+2])/4.0
}

# Centre the experimental data using the interpolated mean model output
# output mean at data point, residual and relative error (with mean) across data points (consider other full-field metrics)
residual = abs(exp_displacement - mu_y)
rel_error = (residual/abs(mu_y))*100
write.csv(cbind(experimental_data[c("X","Y","Z")],mu_y,residual,rel_error), "outputs/mean_error.csv", row.names = FALSE)
exp_displacement_cen = (exp_displacement - mu_y)/sd_dt
# Store experimental data as y and convert to vector to pass to stan
y = as.vector(exp_displacement_cen)
 
#-------------------------------------------------------------------------------

# In this Section, W_y, the (prior) precision of the observation error is specified
# Here we give a value to the expected standard deviation of the observation error,
# which will later be weighted by parameter lambda_y, which will be given a 
# strong prior centred roughly around a value of 1. Here I assumed that error at
# each of the DIC observations is iid with standard deviation sigma_error. A 
# value of 0.01 has been chosen considering samples loaded at 0kN.

# In practice there will also be an error correlated across all observations due 
# to the fact that I had to define a zero datum for the DIC data. This can be
# added to every entry of the covariance matrix as sigma_shift, but if not
# careful this will dominate. I'm unsure of a precise value, but I suggest
# investigating differing values of 0, 0.0033, 0.01, and 0.025

# From the data, the DIC error appears correlated, with regions of higher noise
# appearing in clusters rather than as white noise. Ultimately, it would be good
# to represent these as random fields about which we infer the amplitude and
# correlation length parameters, but this requires a different statistical model

# Finally, there is an error which I do not account for, due to misalignment of
# the DIC data and model. This error would obviously be correlated, applying to
# all data points, but would not take the form of a uniform value. In some places
# such a shift would result in higher values than the "truth", in some cases
# lower. I'm unsure of how to account for this in practice.

# There will also be a "time" correlated error in-line with Higdon et al, but 
# there is no need to account for this here for one load value

# sigma_error = 0.01 # standard deviation associated with noise (a choice of 0.005 would also be reasonable if this is too large)
# sigma_shift = 0.0033 # standard deviation due to a shift in the zero value of the DIC data
# sigma_shift = 0
# Sigma_y = diag(rep(sigma_error^2,n_y)) + matrix(sigma_shift^2, nrow = n_y, ncol = n_y)

# Because the experimental data has been standardised by dividing through by the 
# standard deviation of the model output, the covariance of the measurement error
# must also be divided through by the variance of the model output for consistency
# Sigma_y = Sigma_y/(sd_dt^2)
# Convert covariance matrix to a precision by taking the inverse, as this is what
# is specified to stan as in Higdon et al.
# W_y = solve(Sigma_y) # There are almost certainly more efficient ways of implementing this...
W_y = diag(rep(1,n_y))

#-------------------------------------------------------------------------------

# The statistical model for the experimental data is y = Kw(theta) + Dv + e. For this 
# model it is necessary to evaluate the basis fuctions, K (as determined above 
# using SVD) at the physical locations of the experimental measurements in y. 
# To do this is is necessary do interpolate the basis functions, as was the case
# with the mean model output.

# Firstly interpolate emulator basis vectors, K.
K_y = matrix(NA,n_y,p_eta)
for (i in 1:n_y) {
  # Might be clearer just to write out basis equations in full... 
  bases = 1.0 + matrix(hr[i,],2,4)*HR
  # The complete basis functions for the quad element are given by the product
  # of those in h and r, contained in the rows of "bases"
  # Note that the element index is in Python indexing convention, but connectivities are in the Abaqus convention
  K_y[i,] = colSums(matrix(bases[1,]*bases[2,],4,p_eta)*K_eta[as.numeric(connectivity[y_element[i]+1,])+2,])/4.0
  connectivity[y_element[i]+1,]
}

output_frame = cbind(experimental_data[c("X","Y","Z")],K_y)
for (i in 1:p_eta) {
  colnames(output_frame)[3+i] = sprintf("K_y,%d",i)
}
# output interpolated bases for plotting
write.csv(output_frame, "outputs/interpolated_basis.csv", row.names = FALSE)


# Next specify the basis vectors for the discrepancy D.
# In the formulation proposed by Higdon et al. it is possible to used different
# bases to model discrepancy as opposed to model output. This is an opportunity
# to express prior belief about model accuracy, for instance, a linear model 
# will have higher discrepancy at higher applied loads.

# Define set of orthogonal polynomials for model discrepancy bases
z_norm = (experimental_data$Z/420)
D_y = matrix(0,n_y,10)
# Could possible automate Gram-Schmidt to do this in a loop...
D_y[,1] = sqrt(3)*z_norm
D_y[,2] = sqrt(80)*(z_norm^2 - 3*z_norm/4)
D_y[,3] = 15*sqrt(7)*(z_norm^3 - 4*z_norm^2/3 + 2*z_norm/5)
D_y[,4] = 168*z_norm^4 - 315*z_norm^3 + 180*z_norm^2 - 30*z_norm
D_y[,5] = 210*sqrt(11)*(z_norm^5 - 12*z_norm^4/5 + 2*z_norm^3 - 2*z_norm^2/3 + z_norm/14)
D_y[,6] = 792*sqrt(13)*(z_norm^6 - 35*z_norm^5/12 + 35*z_norm^4/11 - 35*z_norm^3/22 + 35*z_norm^2/99  - 7*z_norm/264)
D_y[,7] = 3003*sqrt(15)*(z_norm^7 - 24*z_norm^6/7 + 60*z_norm^5/13 - 40*z_norm^4/13 + 150*z_norm^3/143 - 24*z_norm^2/143 + 4*z_norm/429)
D_y[,8] = 11440*sqrt(17)*(z_norm^8 - 63*z_norm^7/16 + 63*z_norm^6/10 - 21*z_norm^5/4 + 63*z_norm^4/26 -63*z_norm^3/104 + 21*z_norm^2/286 - 9*z_norm/2860)
D_y[,9] = 43758*sqrt(19)*(z_norm^9 - 40*z_norm^8/9 + 140*z_norm^7/17 - 140*z_norm^6/17 + 245*z_norm^5/51 - 28*z_norm^4/17 + 70*z_norm^3/221 - 20*z_norm^2/663 + 5*z_norm/4862)
D_y[,10] = 167960*sqrt(21)*(z_norm^10 - 99*z_norm^9/20 + 198*z_norm^8/19 - 231*z_norm^7/19 + 2772*z_norm^6/323 - 4851*z_norm^5/1292 + 1617*z_norm^4/1615 - 99*z_norm^3/646 + 99*z_norm^2/8398 - 11*z_norm/33592)

# Read in nominal value of node coordinates for determining basis functions for
# predictions
D_eta = matrix(0,n_eta,10)
node_coords = as.matrix(fread("inputs/CSpar_sam_mesh_nodes.csv"))
z_norm = node_coords[,4]/420
D_eta[,1] = sqrt(3)*z_norm
D_eta[,2] = sqrt(80)*(z_norm^2 - 3*z_norm/4)
D_eta[,3] = 15*sqrt(7)*(z_norm^3 - 4*z_norm^2/3 + 2*z_norm/5)
D_eta[,4] = 168*z_norm^4 - 315*z_norm^3 + 180*z_norm^2 - 30*z_norm
D_eta[,5] = 210*sqrt(11)*(z_norm^5 - 12*z_norm^4/5 + 2*z_norm^3 - 2*z_norm^2/3 + z_norm/14)
D_eta[,6] = 792*sqrt(13)*(z_norm^6 - 35*z_norm^5/12 + 35*z_norm^4/11 - 35*z_norm^3/22 + 35*z_norm^2/99  - 7*z_norm/264)
D_eta[,7] = 3003*sqrt(15)*(z_norm^7 - 24*z_norm^6/7 + 60*z_norm^5/13 - 40*z_norm^4/13 + 150*z_norm^3/143 - 24*z_norm^2/143 + 4*z_norm/429)
D_eta[,8] = 11440*sqrt(17)*(z_norm^8 - 63*z_norm^7/16 + 63*z_norm^6/10 - 21*z_norm^5/4 + 63*z_norm^4/26 -63*z_norm^3/104 + 21*z_norm^2/286 - 9*z_norm/2860)
D_eta[,9] = 43758*sqrt(19)*(z_norm^9 - 40*z_norm^8/9 + 140*z_norm^7/17 - 140*z_norm^6/17 + 245*z_norm^5/51 - 28*z_norm^4/17 + 70*z_norm^3/221 - 20*z_norm^2/663 + 5*z_norm/4862)
D_eta[,10] = 167960*sqrt(21)*(z_norm^10 - 99*z_norm^9/20 + 198*z_norm^8/19 - 231*z_norm^7/19 + 2772*z_norm^6/323 - 4851*z_norm^5/1292 + 1617*z_norm^4/1615 - 99*z_norm^3/646 + 99*z_norm^2/8398 - 11*z_norm/33592)

D_y = D_y[,1:p_delta]
D_eta = D_eta[,1:p_delta]

# Add in extra terms to account for change in slope at ends, and a mis-located
# pivot. Increase p_delta accordingly
#####D_y = D_y[,1:(p_delta-3)]
# D_y = D_y[,1:(p_delta-2)]
##### D_y = cbind(D_y, 4*(experimental_data$X/55 - 0.5)*(experimental_data$Z/420 - 0.5))
# D_y = cbind(D_y, 2*experimental_data$X/55*(experimental_data$Z/420 - 0.5))
# D_y = cbind(D_y, 2*(experimental_data$X/55-1)*(experimental_data$Z/420 - 0.5))
# qr(D_y)$rank # calculate the rank of B
####D_eta = D_eta[,1:(p_delta-3)]
# D_eta = D_eta[,1:(p_delta-2)]
####D_eta = cbind(D_eta, 4*(node_coords[,2]/55 - 0.5)*(node_coords[,4]/420 - 0.5))
# D_eta = cbind(D_eta, 2*node_coords[,2]/55*(node_coords[,4]/420 - 0.5))
# D_eta = cbind(D_eta, 2*(node_coords[,2]/55-1)*(node_coords[,4]/420 - 0.5))



# Write functions to csv files for plotting
write.csv(D_eta, "outputs/discrepancy_basis.csv", row.names = FALSE)
output_frame = cbind(experimental_data[c("X","Y","Z")],D_y)
for (i in 1:p_delta) {
  colnames(output_frame)[3+i] = sprintf("D_y,%d",i)
}
write.csv(output_frame, "outputs/discrepancy_basis_interp.csv", row.names = FALSE)

# The formulation proposed by Higdon et al combines matrices K_y and D_y into a
# single matrix B, which is subsequently used to calculate some pseudo-inverses
# to calculate combined reduced-dimensional variable z_hat. If there are
# duplicated columns in B, as will be the case when using identical bases for 
# the emulator and discrepancy, this will result in these pseudo-inverses not 
# being possible due to B not being full-rank. An solution was proposed by
# specifying B = B_tilde * L, where B_tilde is full rank such that B * z can be 
# substituted with B_tilde * z_tilde, where z_tilde = L*z. The pseudo-inverses
# are then carried out for z_tilde.

#-------------------------------------------------------------------------------

# Code for if identical bases are used for emulator and discrepancy

# In this case I use identical bases for both the emulator and discrepancy. It
# seems reasonable that in the absence of prior information the same bases could
# be used. While a simpler formulation could be used in this case, I chose to 
# stick with that of Higdon et al. for the sake of generality.
#D_y = K_y[,1:p_delta]

# Not used here, but for predictions. Define Discrepancy basis matrix at points
# at which predictions are made
#D_eta = K_eta[,1:p_delta]

# In this case it is easy to specify B_tilde as equal to K_y, as this contains
# all of the unique bases, then define L as a matrix mapping the corresponding 
# columns of this matrix to also represent the discrepancy bases
#L = matrix(0,nrow = p_eta, ncol = p_eta + p_delta)
#L[1:p_delta,1:p_delta] = diag(p_delta)
#L[,(p_delta+1):(p_delta+p_eta)] = diag(p_eta)
# B_tilde = K_y

# In general, if similar but not identical bases are used for D_y and K_y it can
# result in a matrix which is numerically singular. In such an instance L and 
# B_tilde should be determined numerically.

#-------------------------------------------------------------------------------

# This portion of the code deals with the matrix algebra from Section 2.2.4 of 
# Higdon et al., calculating the necessary quantities for passing to stan, where
# the sampling is undertaken

# Assemble B matrix, which maps between the experimental data points,y, and the
# reduced-dimension coefficients, v (discrepancy) and u (emulator) corresponding 
# to these points. For n = 1 experiment this is simply the concatenation of D_y 
# and K_y
B = cbind(D_y,K_y)

# As stated above, B may not be full rank, which can make it impossible to
# calculating the inverse of B'*W_y*B. Below is the implementation of two 
# methods for accounting for this possibility, as well as the standard method
# from Higdon et al. The other possibilities are:
  # i) Do sampling using alternative full-rank B_tilde as described above
# ii) Add a small diagonal ridge term to the B'*W_y*B matrix which is to be
#     inverted
# Note that I define the inverse here rather than solving the equations via 
# other more efficient means, as the inverse is used for multiple calculations
ridge = 1e-4 # magnitude of the ridge term if used
rank_B = qr(B)$rank # calculate the rank of B
# Calculate B'*W_y*B
BTWyB = t(B)%*%W_y%*%B
# Is B full rank?
if (rank_B < (p_eta+p_delta)){
  if (exists("B_tilde")){
    # If B is not full rank, but matrix B_tilde is pre-defined, then use this 
    # full-rank matrix instead of B
    print("B is not full column rank, using B_tilde instead")
    # Calculate B_tilde'*W_y*B_tilde
    BTWyB = t(B_tilde)%*%W_y%*%B_tilde
    BTWyBinv = solve(BTWyB) # Determine inverse
    # Update B with B_tilde for all subsequent calculations
    B = B_tilde
  } else {
    # If B is not full rank, but an alternative full-rank matrix has not been 
    # provided, add a ridge to B'*W_y*B
    print("Adding ridge as B is not full column rank and a reduced-rank alternative has not been specified")
    BTWyB = BTWyB + diag(p_eta+p_delta)*ridge
    BTWyBinv = solve(BTWyB) # Determine inverse
  }
} else {
  # If B is full rank then B'*W_y*B may be solved directly
  print("B is full column rank, solve directly")
  BTWyBinv = solve(BTWyB)
}

# Calculate inverse of K'*K. This is used for multiple calculations and so it is
# more efficient to store the inverse than to solve the equations via other means
# Here K is the matrix of emulator basis functions evaluated for full model 
# output, arranged as specified at the end of Section 2.2.2 of Higdon et al.
# Note that this is a very large matrix (m x n_eta) by (m x p_eta), and as such
# calculating K'*K is computationally expensive, and requires a prohibitively 
# large amount of memory. To make this code possible I've calculated closed-form
# expressions for this matrix product, and inputted K'*K directly.

# Populate dense matrix of products of basis vectors, 
# with KTK_dense[i,j] = k_i'*k_j
# These product of individual basis vectors make up the entries of K'*K
KTK_dense = t(K_eta)%*%K_eta
# Note that the below code is the general expression for this matrix product if
# basis vectors are not orthogonal. As the basis vectors obtained by SVD are 
# orthogonal such that k_i'*k_j = 0 if i not equal to j, the code for populating
# off-diagonal terms has been commented out. This may be adapted to non-
# orthogonal bases by un-commenting the nested for loop.
KTK = matrix(0,m*p_eta,m*p_eta)
for (i in 1:p_eta){
  # Populate diagonal terms
  KTK[((i-1)*m+1):(i*m),((i-1)*m+1):(i*m)] = diag(m)*KTK_dense[i,i]
  # Note, Need to use seq_len as colon operator in R does not assume the list of 
  # indices must always increase, and j in 1:(i+1) throws an error when i = p_eta
  # for (j in (i + seq_len(p_eta-i))){
  #  # Populate off-diagonal terms. Note that these will be zero if the ks are orthogonal
  #  KTK[((i-1)*m+1):(i*m),((j-1)*m+1):(j*m)] = diag(m)*KTK_dense[i,j]
  #  KTK[((j-1)*m+1):(j*m),((i-1)*m+1):(i*m)] = diag(m)*KTK_dense[i,j]
  #}
}
# Determine inverse of K'*K. Note that as this matrix is diagonal for
# orthogonal basis functions in such cases it would be more efficient to 
# directly input the inverse, by inputting the reciprocal of the entries in the 
# above loop.
KTKinv = solve(KTK)

# Calculate K'*eta, which is used to calculate OLS solution w_hat. This 
# calculation has also been implemented by directly inputting the closed-form
# expression for this product, due to the large size of K.
# Note that in Higdon et al, eta is reshaped from a n_eta x m matrix into a 
# (m x n_eta) x 1 vector, as eta = [eta_1',eta_2',...,eta_m']'.
# K'*eta therefore results in a stack of
# [k_1'*eta_1;...;k_1'*eta_m;k_2'*eta_1;k_2'*eta_m;...;k_p_eta'*eta_m]
# It is more efficient to do this product by keeping K_eta stored as a matrix
# and likewise retaining eta as a matrix (for the purposes of this comment let's
# call this matrix eta_mat), then calculating K'*eta by obtaining eta_mat'*K_eta
# reshaped column-wise into a vector, as below
KTeta = t(eta)%*%(K_eta)
KTeta = as.vector(KTeta)

# Determine pseudo-inverses of basis matrices to get reduced-dimensional 
# representation of model output and experimental data
w_hat = KTKinv%*%KTeta # reduced-dimensional representation of emulator training data
# Calculate reduced dimensional representation of discrepancy and emulator 
# corresponding to experimental data.
# Note. I have not implemented an efficient version of the below matrix algebra
# as the quantity of experimental data is small in this application. When using
# DIC data, I will likely have to implement more closed-form expressions for the
# matrix products as in the case of K.I've implemented the more efficient version
# for b_y_dash - just not the pseudo-inverse
BTWyy = t(B)%*%W_y%*%y
# vu_hat = BTWyBinv%*%t(B)%*%W_y%*%y 
vu_hat = BTWyBinv%*%BTWyy



# Group together all samples from experimental data and model output into a 
# single vector
z_hat = as.vector(rbind(vu_hat,w_hat))

# Define parameters for (gamma) prior distributions on the emulator error and 
# observation error terms
a_eta = 1.0     # Shape parameter for the lambda_eta prior
b_eta = 0.0001  # Rate parameter for the lambda_eta prior 
a_y = 5.0       # Shape parameter for the lambda_y prior
# b_y = 5.0       # Rate parameter for the lambda_y prior
# b_y = 0.05
b_y = 0.5

a_eta_dash = a_eta+(0.5*(m*(n_eta-p_eta))) # Adjusted shape parameter for the lambda_eta prior (Eq. 11 Higdon et al.)
a_y_dash = a_y+(0.5*(n_y-rank_B)) # Adjusted shape parameter for the lambda_y prior.

# Stack reshaped model output eta into a single vector. This is the format used 
# in Higdon et al., although it wasn't necessary for the above calculations it
# is used to adjust the prior parameters.
eta_vec = as.vector(eta)
# Re-arraged version of b_eta_dash from that in Eq. (11) for the sake of 
# computational efficiency, using the fact that:
# eta'*(I - K*(K'*K)^-1*K')*eta = eta'*eta - (K'*eta)'*w_hat
b_eta_dash = as.numeric(b_eta + (0.5*(t(eta_vec)%*%eta_vec - t(KTeta) %*% w_hat)))
#b_y_dash = as.numeric(b_y + (0.5*(t(y)%*%(W_y - W_y%*%B%*%BTWyBinv%*%t(B)%*%W_y)%*%y)))
b_y_dash = as.numeric(b_y + (0.5*(t(y)%*%W_y%*%y - t(BTWyy) %*% vu_hat)))

#-------------------------------------------------------------------------------


# This segment of code deals with setting up the environment for stan, passing 
# data to stan, and running the correct stan code depending upon which method is
# required for sampling

# Set up the environment for the Stan model to run in parallel. Taken from:
# https://betanalpha.github.io/assets/case_studies/gaussian_processes.html#21_Simulating_From_A_Gaussian_Process
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
util = new.env()
par(family="CMU Serif", las=1, bty="l", cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 5))

# If basis matrix B is not full column rank, and a full-rank equivalent B_tilde 
# is used in its place, a different stan code is required to calculate as the
# dimension of the covariance matrix will be smaller, and will consequently 
# require some adjustment. This calculation must happen within stan as the 
# matrix in question depends upon uncertain parameters.


# Sort out conditionals later. Having just one general purpose code might be
# better, at least for the discrepancy bit


# Code for a single set of discrepancy parameters, using the full rank code
stan_data = list(m=m, q=q, n_eta=n_eta, n_y=n_y, p_eta=p_eta, p_delta=p_delta, a_eta_dash = a_eta_dash, a_y_dash = a_y_dash, b_eta_dash = b_eta_dash, b_y_dash = b_y_dash, z_hat = z_hat, tf_param_1=tf_param_1, tf_param_2=tf_param_2, tc = tc, KTKinv = KTKinv, BTWyBinv = BTWyBinv)
fit = stan(file = "source/Full_Field_Calibration_Higdon_Prior_Spec_one_disc.stan",
           data = stan_data,
           iter = 4000,
           chains = 3)


# Code for a single set of discrepancy parameters, using the reduced rank code
stan_data = list(m=m, q=q, n_eta=n_eta, n_y=n_y, p_eta=p_eta, p_delta=p_delta, rank_B=rank_B,a_eta_dash = a_eta_dash, a_y_dash = a_y_dash, b_eta_dash = b_eta_dash, b_y_dash = b_y_dash, z_hat = z_hat, tf_param_1=tf_param_1, tf_param_2=tf_param_2, tc = tc, KTKinv = KTKinv, BTWyBinv = BTWyBinv, L = L)
fit = stan(file = "source/Full_Field_Calibration_Higdon_Red_Rank_one_disc.stan",
           data = stan_data,
           iter = 4000,
           chains = 3)


# Conditional statement for testing which stan code to run
if ((rank_B < (p_eta+p_delta)) && exists("B_tilde")){
  # B is not full rank, but alternative B_tilde has been provided
  
  # Create list of relevant data for input to stan.
  # Note that it is necessary to pass L to stan in this case for adjusting the 
  # covariance matrix
  stan_data = list(m=m, q=q, n_eta=n_eta, n_y=n_y, p_eta=p_eta, p_delta=p_delta, rank_B=rank_B,a_eta_dash = a_eta_dash, a_y_dash = a_y_dash, b_eta_dash = b_eta_dash, b_y_dash = b_y_dash, z_hat = z_hat, tf_param_1=tf_param_1, tf_param_2=tf_param_2, tc = tc, KTKinv = KTKinv, BTWyBinv = BTWyBinv, L = L)
  # Run stan code to sample from posterior distributions
  fit = stan(file = "source/Full_Field_Calibration_Higdon_Red_Rank.stan",
               data = stan_data,
               iter = 4000,
               chains = 3)
} else {
  # If B is of full column rank, or a ridge is used to force B to be of full rank
  # then the code may be used as described directly in Higdon et al.
  
  # Create list of relevent data for input to stan
  stan_data = list(m=m, q=q, n_eta=n_eta, n_y=n_y, p_eta=p_eta, p_delta=p_delta, a_eta_dash = a_eta_dash, a_y_dash = a_y_dash, b_eta_dash = b_eta_dash, b_y_dash = b_y_dash, z_hat = z_hat, tf_param_1=tf_param_1, tf_param_2=tf_param_2, tc = tc, KTKinv = KTKinv, BTWyBinv = BTWyBinv)
  # Run stan code to sample from posterior distributions
  fit = stan(file = "source/Full_Field_Calibration_Higdon_Prior_Spec.stan",
             data = stan_data,
             iter = 4000,
             chains = 3)
}

#-------------------------------------------------------------------------------

# This section of code deals with post-processing of the data coming out of the 
# simulations

# plot trace plots of simulation samples
# dev.new(noRStudioGD = TRUE)  # generate plots in separate window
stan_trace(fit, pars = c("tf", "rho_w", "lambda_w", "lambda_v","lambda_y","lambda_eta"))

# Summarise results to check convergence
print(fit, pars = c("tf", "rho_w", "lambda_w", "lambda_v","lambda_y","lambda_eta"))

# extract samples from stan output

samples <- extract(fit)
N_samples = dim(samples$rho_w)[1] # get total number of samples

# Extract label of inputs for plots
labels = colnames(XT_sim)

# Produce plots of posterior and prior distributions of correlation parameters
# (rho) for emulator. Here the rows corresponding to the different principal
# components, whereas the columns correspond to the different calibration inputs
# A values of rho close to 1 implies that an input does strongly affect the 
# model output
dev.new(noRStudioGD = TRUE)
par(mfrow = c(p_eta,q))
rho_plot = seq(0,1, length.out = 100)
for (i in 1:p_eta){
  for (j in 1:q){
    # Plot histogram of posterior distribution
    hist(samples$rho_w[,(i-1)*q+j],
         main = labels[j],
         #main = paste("rho_w,",as.character(i),",",as.character(j)),
         xlab = paste("rho_w,",as.character(i),",",as.character(j)),
         col = "firebrick1",
         breaks = 25,
         freq = FALSE,
         xlim = c(0,1), 
         cex.lab=1.5,
         cex.axis=1.5)
    # Plot prior distribution
    prior_plot = dbeta(rho_plot,shape1=1,shape2=0.1)
    lines(rho_plot,prior_plot,lwd=3,col="blue")
  }
}

# Plot precision parameters for emulator. Each plot corresponds to a different
# principal component
dev.new(noRStudioGD = TRUE)
par(mfrow = c(1,p_eta))
lambda_plot = seq(0,2.5, length.out = 100)
for (i in 1:p_eta) {
  # plot posterior
  hist(samples$lambda_w[,i],
       main = paste("lambda_w",as.character(i)),
       xlab = paste("lambda_w",as.character(i)),
       col = "firebrick1",
       breaks = 25,
       freq = FALSE,
       xlim = c(0,2.5),
       cex.axis=1.5,
       cex.lab=1.5)
  # plot prior
  prior_plot = dgamma(lambda_plot,shape=5,rate=5)
  lines(lambda_plot,prior_plot,lwd=3,col="blue")
}

# Plot precision parameters for the discrepancy. Each plot corresponds to a 
# different basis function
# dev.new(noRStudioGD = TRUE)
# par(mfrow = c(1,p_delta))
# lambda_plot = seq(0,1E5, length.out = 100)
# for (i in 1:p_delta){
#   # plot posterior
#  hist(samples$lambda_v[,i],
#       main = paste("lambda_v",as.character(i)),
#       xlab = paste("lambda_v",as.character(i)),
#       col = "firebrick1",
#       breaks = 25,
#       freq = FALSE,
#       #xlim = c(0,75),
#       cex.axis=1.5,
#       cex.lab=1.5)
#  # plot prior
#  prior_plot = dgamma(lambda_plot,shape=1.0,rate=0.0001)
#  lines(lambda_plot,prior_plot,lwd=3,col="blue")
#}

# Version for if just one set of discrepancy parameters
dev.new(noRStudioGD = TRUE)
lambda_plot = seq(0,1E5, length.out = 100)
# plot posterior
hist(samples$lambda_v,
  main = "lambda_v",
  xlab = "lambda_v",
  col = "firebrick1",
  breaks = 25,
  freq = FALSE,
  # xlim = c(0,1),
  #ylim = c(0,0.001),
  cex.axis=1.5,
  cex.lab=1.5)
# plot prior
prior_plot = dgamma(lambda_plot,shape=1.0,rate=0.0001)
lines(lambda_plot,prior_plot,lwd=3,col="blue")



# Plot precision of the observation error and PCA truncation error
lambda_plot = seq(0,10.0, length.out = 100)
dev.new(noRStudioGD = TRUE)
par(mfrow = c(1,2))
#plot posterior
hist(samples$lambda_y,
     main = "lambda_y",
     xlab = "lambda_y",
     col = "firebrick1",
     breaks = 25,
     freq = FALSE,
#     xlim = c(0,10.0),
     cex.axis=1.5,
     cex.lab=1.5)
# plot prior
prior_plot = dgamma(lambda_plot,shape=a_y,rate=b_y)
lines(lambda_plot,prior_plot,lwd=3,col="blue")
# plot of prior adjusted for equivalent, reduced-dimension normal-gamma model
adj_prior_plot = dgamma(lambda_plot,shape=a_y_dash,rate=b_y_dash)
lines(lambda_plot,adj_prior_plot,lwd=3,col="green")

hist(samples$lambda_eta,
     main = "lambda_eta",
     xlab = "lambda_eta",
     col = "firebrick1",
     breaks = 25,
     freq = FALSE,
     xlim = c(0,5e6),
     cex.axis=1.5,
     cex.lab=1,5)
# plot prior
lambda_plot = seq(0,5E6, length.out = 1000)
prior_plot = dgamma(lambda_plot,shape=a_eta,rate=b_eta)
lines(lambda_plot,prior_plot,lwd=3,col="blue")
adj_prior_plot = dgamma(lambda_plot,shape=a_eta_dash,rate=b_eta_dash)
lines(lambda_plot,adj_prior_plot,lwd=3,col="green")

# Before plotting their distributions, the calibration parameters must be 
# transformed back onto their original scale, recalling that they were initially
# transformed onto [0,1]
tf_trans = samples$tf
tf_trans = tf_trans*matrix(rep(t_max-t_min,N_samples),ncol=q,byrow=TRUE) + matrix(rep(t_min,N_samples),ncol=q,byrow=TRUE)

# Plot prior and posterior distribution of the calibration parameters.
# Note that the calibration parameters were initially defined as deviations from
# their nominal values in this application. Here I've added these nominal values
# back onto the prior and posterior plots, such that they are more readable by 
# engineers
dev.new(noRStudioGD = TRUE) # plot in new window
par(mfrow = c(1,q))
# plot histogram of posterior
for (i in 1:q){
  hist(tf_trans[,i], 
       main = labels[i],
       xlab = labels[i],
       col = "brown1",
       breaks = 25,
       freq = FALSE,
       xlim = c(t_min[i],t_max[i]),
       cex.lab = 1.25,
       cex.axis = 1.25)
  # overlay plot of prior distribution
  t_plot = seq(t_min[i],t_max[i], length.out = 100)
  if (i < q){
    print(i)
    prior_plot = dnorm(t_plot,mean = tf_param_1[i]*(t_max[i]-t_min[i]) + t_min[i], sd = tf_param_2[i]*(t_max[i]-t_min[i]))
    } else {
    prior_plot = dunif(t_plot, min = (tf_param_1[i]*(t_max[i]-t_min[i]) + t_min[i]), max = (tf_param_2[i]*(t_max[i]-t_min[i]) + t_min[i]))
  }
  lines(t_plot,prior_plot,lwd=3,"col"="blue")
}


# estimate mode of the posterior distribution
modes = rep(0,q)
for (i in 1:q){
  modes[i] = estimate_mode(tf_trans[,i])
}
print("calibration parameter modes = ")
print(modes)

print("calibration parameter means = ")
print(colMeans(tf_trans))

#-------------------------------------------------------------------------------


# Code for making predictions - here I'm assuming we're using the transformed 
# version for matrices which are not full rank

# N_sam_plot = 1000
N_sam_plot = 1
# Pick N_sam_plot random samples without repetition
rand_ind = sample.int(N_samples,N_sam_plot)
#N_sam_plot = N_samples
#rand_ind = 1:N_samples


v_star_mu = matrix(0, p_delta, N_sam_plot)
w_star_mu = matrix(0, p_eta, N_sam_plot)
v_star_sigma = array(0, c(p_delta, p_delta, N_sam_plot))
w_star_sigma = array(0, c(p_eta, p_eta, N_sam_plot))

v_star = matrix(0, p_delta, N_sam_plot)
w_star = matrix(0, p_eta, N_sam_plot)

delta_sam = matrix(0, n_eta, N_sam_plot)
delta_y = matrix(0, n_y, N_sam_plot)
eta_sam = matrix(0, n_eta, N_sam_plot)
eta_y = matrix(0, n_y, N_sam_plot)


sam_ind = 1
# Sort out predictions for full-rank code from here!!
for (i in 1:N_sam_plot){
  sam_ind = rand_ind[i]
  # Construct covariance matrix for the current sample
  
  # I could probably package the below code into an external function
  
  # Extract the inferred parameters for the current sample
  beta_w_i = samples$beta_w[sam_ind,]
  #lambda_v_i = samples$lambda_v[sam_ind,]
  # Code for if there are only one set of discrepancy parameters to learn about
  # Probably need some kind of conditional statement
  lambda_v_i = samples$lambda_v[sam_ind]
  
  lambda_w_i = samples$lambda_w[sam_ind,]
  tf_i = samples$tf[sam_ind, ]
  lambda_y_i = samples$lambda_y[i]
  lambda_eta_i = samples$lambda_eta[i]
  
  # THE ABOVE SHOULD WORK FOR ALL PREDICTION METHODS CHECK
  
  # Firstly we need to reconstruct the covariance matrix for the training data points,
  # Sigma_w_hat
  
  # Covariance matrix discrepancy weights.Here I assume that each weight is an
  # independent Gaussian process (i.e. F = p_delta, |G_i| = 1 for all i = 1,...p_delta)
  # Code for set of discrepancy parameters per basis functio
  # sigma_v = diag(1/lambda_v_i)
  # Code for if there is only one set of discrepancy parameters
  sigma_v = diag(rep(1/lambda_v_i,p_delta))
  
  # Define the covariance matrix for the emulator weights, evaluated at the
  # experimental data points. As we only have one experiment, and
  # R(theta,theta) = 1, this also reduces to a diagonal matrix.
  sigma_u = diag(1/lambda_w_i)

  # Combine contributions from sigma_u and sigma_v into a single block-diagonal
  # matrix for subsequent transformation onto their full-rank equivalent combining
  # terms from both
  sigma_vu_join = matrix(0, p_eta+p_delta, p_eta+p_delta)
  sigma_vu_join[1:p_delta, 1:p_delta] = sigma_v
  sigma_vu_join[-(1:p_delta), -(1:p_delta)] = sigma_u
  
  # Define the covariance matrix for the emulator weights, evaluated at the
  # training data points.
  sigma_w = matrix(0, m*p_eta, m*p_eta)
  
  for (j in 1:p_eta) {
    # Calculate the covariance matrix for the ith emulator weight
    sigma_w[((j-1)*m+1):(j*m), ((j-1)*m+1):(j*m)] = ARD_SE_cov(tc, lambda_w_i[j], beta_w_i[((j-1)*q+1):(j*q)], 0)
  }
  
  # Evaluate cross-covariance of emulator weights between training data and experimental data points
  sigma_uw = matrix(0, p_eta, m*p_eta)
  for (j in 1:p_eta) {
    sigma_uw[j,((j-1)*m+1):(j*m)] = ARD_SE_cov_non_sym(t(as.matrix(tf_i)), tc, lambda_w_i[j], beta_w_i[((j-1)*q+1):(j*q)])
  }
  
  # The covariance matrices used in predictions will differ depending on whether
  # or not B is full rank
  if ((rank_B < (p_eta+p_delta)) && exists("B_tilde")){
    # Transform contributions from sigma uv to into their full rank equivalent
    Lsigma_uw = L %*% rbind(matrix(0, p_delta, m*p_eta), sigma_uw)
    # Assemble covariance matrix for joint experimental and model data
    sigma_z = matrix(0,rank_B+m*p_eta, rank_B+m*p_eta)
    sigma_z[1:rank_B,1:rank_B] = L %*% sigma_vu_join %*% t(L)
    sigma_z[-(1:rank_B),-(1:rank_B)] = sigma_w
    sigma_z[1:rank_B,-(1:rank_B)] = Lsigma_uw
    sigma_z[-(1:rank_B),1:rank_B] = t(Lsigma_uw)
  } else {
    # Assemble covariance matrix for joint experimental and model data
    sigma_z = matrix(0,p_delta+(m+1)*p_eta, p_delta+(m+1)*p_eta)
    sigma_z[1:p_delta,1:p_delta] = sigma_v
    sigma_z[(p_delta+1):(p_delta+p_eta),(p_delta+1):(p_delta+p_eta)] = sigma_u
    sigma_z[-(1:(p_delta+p_eta)),-(1:(p_delta+p_eta))] = sigma_w
    sigma_z[(p_delta+1):(p_delta+p_eta),-(1:(p_delta+p_eta))] = sigma_uw
    sigma_z[-(1:(p_delta+p_eta)),(p_delta+1):(p_delta+p_eta)] = t(sigma_uw)
  }
  # Adjust the covariance matrix sigma_z to include transformed emulator and 
  # experimental error terms
  sigma_z_hat = matrix(0,rank_B+m*p_eta, rank_B+m*p_eta)
  sigma_z_hat[1:rank_B,1:rank_B] = BTWyBinv/lambda_y_i
  sigma_z_hat[-(1:rank_B),-(1:rank_B)] = KTKinv/lambda_eta_i
  sigma_z_hat = sigma_z_hat + sigma_z

  # Likewise covariance matrix with predictions will depend upon whether full-
  # rank format is used
  if ((rank_B < (p_eta+p_delta)) && exists("B_tilde")){
    # Define cross-correlation terms for z with discrepancy prediction, v_star. Note
    # That this does not need to be adjusted to as sigma_z_hat was, as the error 
    # terms do not apply to predictions
    sigma_z_v_star = rbind(sigma_v,matrix(0,p_eta,p_delta))
    # Need to multiply by L to obtain full-rank equivalent
    Lsigma_z_v_star = L %*% sigma_z_v_star
    # The resulting matrix is then  padded with zeros for correlation with w
    Lsigma_z_v_star = rbind(Lsigma_z_v_star,matrix(0,m*p_eta,p_delta))
  
    # Define cross-correlation of emulator predictions w_star with emulator training
    # data. Note that in this application, because we are making a prediction at the
    # same value of controlled input x, and uncontrolled inputs are set to their
    # calibrated values as was the case for the experimental data, then this is 
    # given by sigma_u_w', Note that this won't be the case if we want to make a
    # prediction at a different x, or want to evaluate just the emulator across a 
    # range of t values
    sigma_w_w_star = t(sigma_uw)

    # Define cross-correlation terms for z with emulator prediction w_star
    sigma_z_w_star = rbind(matrix(0,p_delta,p_eta), sigma_u)
    Lsigma_z_w_star = L %*% sigma_z_w_star
    # Combine with the correlation matrix for w to get the full correlation with z
    Lsigma_z_w_star = rbind(Lsigma_z_w_star,sigma_w_w_star)
  } else {
    # Define cross-correlation terms for z with discrepancy prediction, v_star. 
    # Note that this does not need to be adjusted to as sigma_z_hat was, as the 
    # error terms do not apply to predictions
    sigma_z_v_star = rbind(sigma_v,matrix(0,(m+1)*p_eta,p_delta))
    
    # Define cross-correlation of emulator predictions w_star with emulator training
    # data. Note that in this application, because we are making a prediction at the
    # same value of controlled input x, and uncontrolled inputs are set to their
    # calibrated values as was the case for the experimental data, then this is 
    # given by sigma_u_w', Note that this won't be the case if we want to make a
    # prediction at a different x, or want to evaluate just the emulator across a 
    # range of t values
    sigma_w_w_star = t(sigma_uw)
    # Define cross-correlation terms for z with emulator prediction w_star
    sigma_z_w_star = rbind(matrix(0,p_delta,p_eta), sigma_u, sigma_w_w_star)
  }
  
  # Define correlation of discrepancy prediction v_star with itself
  sigma_v_star = sigma_v
  
  # Define correlation of emulator prediction with itself
  sigma_w_star = sigma_u
  
  # Consider adding a nugget...
  
  # Two different methods for making predictions. The first method is quicker, but 
  # I think the second is more numerically stable, which seems to make a different
  # when the variance is small. Consider using this if the simulation is taking 
  # too long
  
  # Explicitly calculating inverse, then reusing for all predictions. 
  # Ainv = solve(sigma_z_hat)
  # v_mu = t(Lsigma_z_v_star) %*% Ainv %*% z_hat
  # w_mu = t(Lsigma_z_w_star) %*% Ainv %*% z_hat
  # v_sigma = sigma_v_star - (t(Lsigma_z_v_star) %*% Ainv %*% Lsigma_z_v_star)
  # w_sigma = sigma_w_star - (t(Lsigma_z_w_star) %*% Ainv %*% Lsigma_z_w_star)

  # Solving using solve
  Ainv_z_hat = solve(sigma_z_hat,z_hat)
  # Store mean and covariance matrices of discrepancy and adjusted prediction Gaussian processes
  if ((rank_B < (p_eta+p_delta)) && exists("B_tilde")){
    v_star_mu[,i] = t(Lsigma_z_v_star) %*% Ainv_z_hat
    w_star_mu[,i] = t(Lsigma_z_w_star) %*% Ainv_z_hat
    v_star_sigma[,,i] = sigma_v_star - (t(Lsigma_z_v_star) %*% solve(sigma_z_hat,Lsigma_z_v_star))
    w_star_sigma[,,i] = sigma_w_star - (t(Lsigma_z_w_star) %*% solve(sigma_z_hat,Lsigma_z_w_star))
  } else {
    v_star_mu[,i] = t(sigma_z_v_star) %*% Ainv_z_hat
    w_star_mu[,i] = t(sigma_z_w_star) %*% Ainv_z_hat
    v_star_sigma[,,i] = sigma_v_star - (t(sigma_z_v_star) %*% solve(sigma_z_hat,sigma_z_v_star))
    w_star_sigma[,,i] = sigma_w_star - (t(sigma_z_w_star) %*% solve(sigma_z_hat,sigma_z_w_star))
  }

  # Sample from Gaussian process (might not be necessary - could just cheat and 
  # look at the mean
  v_star[,i] = mvrnorm(n = 1, v_star_mu[,i], v_star_sigma[,,i])
  w_star[,i] = mvrnorm(n = 1, w_star_mu[,i], w_star_sigma[,,i])
  
  # Generate individual Samples
  delta_sam[,i] = (D_eta %*% v_star[,i])*sd_dt
  delta_y[,i] = (D_y %*% v_star[,i])*sd_dt
  eta_sam[,i] = (K_eta %*% w_star[,i])*sd_dt + mu_dt
  eta_y[,i] = (K_y %*% w_star[,i])*sd_dt + mu_y
}

# Good to here!!!

# Plot histogram of reduced coefficients
dev.new(noRStudioGD = TRUE) # plot in new window
par(mfrow = c(1, p_eta))
for (i in 1:p_eta){
  hist(w_star[i,], 
       main =  paste("Feature",as.character(i)),
       xlab = paste("w_star", as.character(i)),
       col = "brown1",
       breaks = 25,
       freq = FALSE,
       xlim = c(-4,4),
       cex.lab = 1.25,
       cex.axis = 1.25)
  # overlay plot of prior distribution
  w_plot = seq(-4,4, length.out = 100)
  prior_plot = dnorm(w_plot,mean = 0, sd = 1)
  lines(w_plot,prior_plot,lwd=3,"col"="blue")
}

dev.new(noRStudioGD = TRUE) # plot in new window
par(mfrow = c(1, p_delta))
for (i in 1:p_delta){
  hist(w_star[i,], 
       main =  paste("Feature",as.character(i)),
       xlab = paste("v_star", as.character(i)),
       col = "brown1",
       breaks = 25,
       freq = FALSE,
       xlim = c(-4,4),
       cex.lab = 1.25,
       cex.axis = 1.25)
  # overlay plot of prior distribution
  w_plot = seq(-4,4, length.out = 100)
  prior_plot = dnorm(w_plot,mean = 0, sd = 1)
  lines(w_plot,prior_plot,lwd=3,"col"="blue")
}

# Also look at samples of the mean for output
# Multiply stored mean values by the appropriate basis matrix and transform back
# onto their correct scales
delta_mu = (D_eta %*% v_star_mu)*sd_dt
delta_y_mu = (D_y %*% v_star_mu)*sd_dt
eta_mu = (K_eta %*% w_star_mu)*sd_dt + mu_dt
eta_y_mu = (K_y %*% w_star_mu)*sd_dt + mu_y

# Calculate the diagonal terms of the covariance matrix. Just do this as the
# full covariance is too big. Take the square root (standard deviation) as this
# is more meaningful
delta_sigma = matrix(0, n_eta, N_sam_plot)
delta_y_sigma = matrix(0, n_y, N_sam_plot)
eta_sigma = matrix(0, n_eta, N_sam_plot)
eta_y_sigma = matrix(0, n_y, N_sam_plot)
for (i in 1:N_sam_plot){
  delta_sigma[,i] = sqrt(as.matrix(rowSums((D_eta %*% v_star_sigma[,,i]) * D_eta)))*sd_dt
  delta_y_sigma[,i] = sqrt(as.matrix(rowSums((D_y %*% v_star_sigma[,,i]) * D_y)))*sd_dt
  eta_sigma[,i] = sqrt(as.matrix(rowSums((K_eta %*% w_star_sigma[,,i]) * K_eta)))*sd_dt
  eta_y_sigma[,i] = sqrt(as.matrix(rowSums((K_y %*% w_star_sigma[,,i]) * K_y)))*sd_dt
}

# Integrate uncertainty out of mean and standard deviation by taking average
eta_mu_mu = rowMeans(eta_mu)
delta_mu_mu = rowMeans(delta_mu)
eta_y_mu_mu = rowMeans(eta_y_mu)
delta_y_mu_mu = rowMeans(delta_y_mu)

eta_sigma_mu = rowMeans(eta_sigma)
delta_sigma_mu = rowMeans(delta_sigma)
eta_y_sigma_mu = rowMeans(eta_y_sigma)
delta_y_sigma_mu = rowMeans(delta_y_sigma)



# Also separate emulator

# Here it would be good to take some averages!! We haven't averaged out the calibration parameters!!



# Write all output to text files for plotting
write.csv(eta_mu, "outputs/eta_mu.csv", row.names = FALSE)
write.csv(delta_mu, "outputs/delta_mu.csv", row.names = FALSE)
write.csv(eta_y_mu, "outputs/eta_y_mu.csv", row.names = FALSE)
write.csv(delta_y_mu, "outputs/delta_y_mu.csv", row.names = FALSE)
write.csv(delta_sam, "outputs/delta_sam.csv", row.names = FALSE)
write.csv(eta_sam, "outputs/eta_sam.csv", row.names = FALSE)
write.csv(delta_y, "outputs/delta_y.csv", row.names = FALSE)
write.csv(eta_y, "outputs/eta_y.csv", row.names = FALSE)
write.csv(delta_sigma, "outputs/delta_sigma.csv", row.names = FALSE)
write.csv(delta_y_sigma, "outputs/delta_y_sigma.csv", row.names = FALSE)
write.csv(eta_sigma, "outputs/eta_sigma.csv", row.names = FALSE)
write.csv(eta_y_sigma, "outputs/eta_y_sigma.csv", row.names = FALSE)

write.csv(eta_mu_mu, "outputs/eta_mu_mu.csv", row.names = FALSE)
write.csv(delta_mu_mu, "outputs/delta_mu_mu.csv", row.names = FALSE)
write.csv(eta_y_mu_mu, "outputs/eta_y_mu_mu.csv", row.names = FALSE)
write.csv(delta_y_mu_mu, "outputs/delta_y_mu_mu.csv", row.names = FALSE)

write.csv(eta_sigma_mu, "outputs/eta_sigma_mu.csv", row.names = FALSE)
write.csv(delta_sigma_mu, "outputs/delta_sigma_mu.csv", row.names = FALSE)
write.csv(eta_y_sigma_mu, "outputs/eta_y_sigma_mu.csv", row.names = FALSE)
write.csv(delta_y_sigma_mu, "outputs/delta_y_sigma_mu.csv", row.names = FALSE)
