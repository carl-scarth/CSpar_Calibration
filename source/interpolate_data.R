# Standalone script for interpolation using the ABAQUS natural coordinates and 
# interpolation functions

# Set working directory
setwd("C:/Users/cs2361/Documents/CSpar_Calibration/")

# Load in the experimental data which contains the natural coordinate definition
experimental_data = read.table("inputs/Ext_LCorner_Image_0115_0.tiff_nat_coord_rad_trim.csv", sep = ',', header = TRUE)
n_y = nrow(experimental_data)# Number of observations
y_element = experimental_data$Element # element index of each point
hr = as.matrix(experimental_data[,c("h","r")]) # natural coordinates of each point within the element
exp_displacement = experimental_data$W

# load mesh connectivity information about the outer surface of the spar. 
connectivity = read.table("inputs/outer_surface_elements.csv", sep = ',', header = TRUE)

# load in nodal output for quantity which is to be interpolated
f_node = read.table("inputs/nominal_inputs_displacements.csv", sep = ',', header = TRUE)$w_1


# Define the position of the nodes in natural coordinates (this is 4,1,5,8) as
# this is what I used to determine the coordinates.
# It would probably make more sense to use standard ordering for a quad element,
# but I would need to re-do the mapping
HR = matrix(c(1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0),nrow=2,ncol=4,byrow = TRUE)

# Create vector of interpolated output. Loop over each data point and perform
# the interpolation
f_y = rep(NA,n_y)
for (i in 1:n_y) {
  # Might be clearer just to write out basis equations in full... 
  bases = 1.0 + matrix(hr[i,],2,4)*HR
  # The complete basis functions for the quad element are given by the product
  # of those in h and r, contained in the rows of "bases"
  # Note that the element index is in Python indexing convention, but connectivities are in the Abaqus convention
  # plus 2 as the first two nodes are the reference points, which are numbered
  # according to a different system as they are defined directly on the assembly
  f_y[i] = sum(bases[1,]*bases[2,]*f_node[as.numeric(connectivity[y_element[i]+1,])+2])/4.0
}

# Centre the experimental data using the interpolated mean model output
# output mean at data point, residual and relative error (with mean) across data points (consider other full-field metrics)
residual = exp_displacement - f_y
abs_residual = abs(residual)
rel_error = (residual/abs(f_y))*100
write.csv(cbind(experimental_data[c("X","Y","Z")],f_y,residual,abs_residual,rel_error), "outputs/interp_error.csv", row.names = FALSE)
