import numpy as np

def intp_nodes_to_cloud(el_ind, gh, f_node, conn, GH = [], skip_nodes = 0):
    # Interpolate nodal quantities f_node from nodes of a quad element mesh to 
    # n_cloud point cloud coordinates. 
    # el_ind = n_cloud vector of integer element indices, matching an element in
    # the connectivity matrix
    # gh = n_cloud x 2 matrix of element natural coordinates, with columns g and h
    # f_node = n_nodes x n_f matrix of nodal values of n_f quantities of interest
    # to be interpolated to point cloud coordinates
    # conn = n_nodes x 4 matrix containing element connectivities
    # GH = 2 x 4 matrix with natural coordinates of the 4 nodes connected to each
    # element with rows g and h. If not specified default Abaqus values are used.
    # skip_nodes = number of nodes to skip at the beginning of f_node This is 
    # useful when the first set of points in the output are not referenced by the 
    # connectivity file

    if (skip_nodes > 0):
        f_node = f_node[skip_nodes:,:]
    
    n_cloud = gh.shape[0] # Number of experimental data points
    n_f = f_node.shape[1] # Number of output quantities to interpolate

    # Define default Abaqus values for nodal natural coordinates, if not inputted
    if len(GH) < 1:
        #GH = np.array([[-1.0,1.0,1.0,-1.0],[-1.0,-1.0,1.0,1.0]])
        GH = np.array([[-1.0,-1.0],[1.0,-1.0],[1.0,1.0],[-1.0,1.0]])
    
    # Interpolate using the functions for Abaqus isoparametric quad elements
    f_cloud = np.empty((n_cloud,n_f))
    for i, el in enumerate(el_ind):
        # bases = 1.0 + gh[i,:]*GH.T
        bases = 1.0 + gh[i,:]*GH
        prod_bases = np.prod(bases, axis = 1, keepdims = True)*f_node[conn[el,:],:]
        # bases = 1.0 + gh[i,:].T*GH # 1-dimensional basis functions
        # prod_bases = np.prod(bases, axis=0)*f_node[conn[el,:],:]
        # f_cloud[i,:] = np.sum(prod_bases, axis=1)/4.0
        f_cloud[i,:] = np.sum(prod_bases, axis=0)/4.0

    return(f_cloud)
