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
  # the connectivity matrix
  # gh = n_cloud x 2 matrix of element natural coordinates, with columns g and h
  # f_node = n_nodes x n_f matrix of nodal values of n_f quantities of interest
  # to be interpolated to point cloud coordinates
  # connectivity = n_nodes x 4 matrix containing element connectivities
  # conn_file = string name of csv file containing element connectivities
  # GH = 2 x 4 matrix with natural coordinates of the 4 nodes connected to each
  # element with rows g and h. If not specified default Abaqus values are used.
  # skip_nodes = number of nodes to skip at the beginning of f_node This is 
  # useful when the first set of points in the output are not referenced by the 
  # connectivity file
  # NOTE: All element and node indices are passed in the Abaqus (and R) 
  # convention, not Python
  
  if (skip_nodes > 0){
    # Need extra argument to prevent R converting to vector when only 1D
    f_node = f_node[-(1:skip_nodes),,drop=FALSE]
  }

  n_cloud = nrow(gh) # Number of experimental data points
  n_f = ncol(f_node) # Number of output quantities to interpolate
  
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
  
  # Interpolate using the Abaqus interpolation functions for isoparametric quad
  # elements
  f_cloud = matrix(NA,n_cloud,n_f)
  for (i in 1:n_cloud) {
    bases = 1.0 + matrix(gh[i,],2,4)*GH # 1-dimensional basis functions
    f_cloud[i,] = colSums(matrix(bases[1,]*bases[2,],4,n_f)*f_node[as.numeric(connectivity[el_ind[i],]),])/4.0
  }

  return(f_cloud)
}

# Wrapper function for interpolating incremented model data (i.e. multiple 
# applied loads) to DIC point cloud data at matched increments
intp_nodes_to_cloud_inc <- function(el_ind, gh, f_node, inc_ind, n_incs, connectivity = NULL, conn_file = NULL, GH = NULL, skip_nodes = 0) {
  # n_incs = integer number of increments
  # el_ind = vector of integer element indices for each point in the cloud 
  # across all increments, matching an element in the connectivity matrix
  # gh = n_points x 2 matrix of element natural coordinates, with columns g and h
  # Both this and el_ind are the concatenation of points across all increments
  # f_node = (n_nodes*n_incs) x n_f matrix of nodal values of n_f quantities of 
  # interest, concatenated across all increments
  # inc ind = vector of integers indicating which increment each point in the
  # cloud belongs to
  # All other inputs are as defined above for intp_nodes_to_cloud
  # NOTE: It  isn't necessary that there are the same number of points in the 
  # cloud across all increments. It is assumed that there are the same number of
  # nodes across all increments.
  f_cloud = c()
  for (i in 1:n_incs){
    print(i)
    f_i = as.matrix(f_node[((i-1)*n_nodes+1):(i*n_nodes),,drop=FALSE])
    f_cloud = c(f_cloud, intp_nodes_to_cloud(el_ind[inc_ind == i],gh[inc_ind == i,], f_i, connectivity = connectivity, conn_file = conn_file, skip_nodes=skip_nodes))
  }

  return(f_cloud)

}