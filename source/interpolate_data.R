# Scripts for interpolating model and emulator input/output quantities to the 
# coordinates of a DIC point cloud given surface mesh connectivity info, a 
# n_nodes x n_quantities matrix of quantities to be interpolated, and
# a list of matched elements and natural coordinates for each point in the cloud
library(data.table)

intp_nodes_to_cloud <- function(el_ind, gh, f_node, connectivity = NULL, conn_file = NULL, GH = NULL, skip_nodes = 0){
  # Interpolate nodal quantities f_node from nodes of a quad element mesh to 
  # n_cloud point cloud coordinates. Element connectivity can be specified 
  # either via a matrix or csv file.
  # el_ind = n_cloud vector of integer element indices, matching an element in 
  # the connectivity matrix with a point in the cloud
  # gh = n_cloud x 2 matrix of element natural coordinates, with columns g and h
  # f_node = n_nodes x n_f matrix of nodal values of n_f quantities of interest
  # to be interpolated to point cloud coordinates
  # connectivity = n_nodes x 4 matrix containing element connectivities
  # NOTE: Connectivities are provided in the Abaqus index convention, not Python
  # conn_file = string name of csv file containing element connectivities
  # GH = 2 x 4 matrix with natural coordinates of the 4 nodes connected to each
  # element with rows g and h. If not specified default Abaqus values are used.
  # skip_nodes = number of nodes to skip at the beginning of f_node This is 
  # useful when the first set of points in the output are not referenced by the 
  # connectivity file
  
  print(skip_nodes)
  if (skip_nodes > 0){
    print(nrow(f_node))
    f_node = f_node[-(1:skip_nodes),,drop=FALSE]
    # Need extra argument to prevent R converting to vector when only 1D
    print(nrow(f_node))
  }

  n_cloud = nrow(gh) # Number of experimental data points
  n_f = ncol(f_node) # Number of output quantities to interpolate
  print(n_cloud)
  print(n_f)
  
  # Load connectivity matrix from file if it hasn't been provided directly, or
  # throw and error if conn_file also hasn't been provided
  if (is.null(connectivity)){
    if (is.null(conn_file)){
      stop("Please provide connectvitiy matrix \"connectivity\", or csv file \"conn_file\" containing element connectivities")
    } else {
      connectivity = fread(conn_file)
    }
  }

  # Define default Abaqus values for nodal natural coordinates, if not inputted
  if (is.null(GH)){
    GH = matrix(c(-1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0),nrow=2,ncol=4,byrow = TRUE)
  }

  # Create vector of interpolated output. Loop over each data point and perform
  # the interpolation
  # f_cloud = rep(NA,n_y)
  #for (i in 1:n_y) {
    # Might be clearer just to write out basis equations in full... 
    #bases = 1.0 + matrix(hr[i,],2,4)*HR
    # The complete basis functions for the quad element are given by the product
    # of those in h and r, contained in the rows of "bases"
    # Note that the element index is in Python indexing convention, but connectivities are in the Abaqus convention
    # plus 2 as the first two nodes are the reference points, which are numbered
    # according to a different system as they are defined directly on the assembly
    #f_cloud[i] = sum(bases[1,]*bases[2,]*f_node[as.numeric(connectivity[y_element[i]+1,])+2])/4.0
  #}
  
  # Interpolate using the abaus interpolation functions for isoparametric quad
  # elements
  f_cloud = matrix(NA,n_cloud,n_f)
  
  #CHECK IF WORKS FOR 1D - HAD ISSUES ABOVE WITH R CONVERTING BACK TO VECTOR
  # USE DROP= FALSE KEYWORD AS ABOVE IF ISSUES
  for (i in 1:n_cloud) {
    bases = 1.0 + matrix(gh[i,],2,4)*GH # 1-dimensional basis functions
    f_cloud[i,] = colSums(matrix(bases[1,]*bases[2,],4,n_f)*f_node[as.numeric(connectivity[el_ind[i],]),])/4.0
  }

  return(f_cloud)
}