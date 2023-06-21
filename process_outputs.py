from odbAccess import *
import numpy as np

# output to command prompt
sys.stdout = sys.__stdout__

# Potentially pass the analysis name as a variable at the command line.
analysis = sys.argv[-1]

# Open odb file. Assumes the scrip is run from the parent directory, and not from the Abaqus subfolder
o1 = session.openOdb(name = analysis + '.odb')
# Get last frame of the load step.
lastFrame = o1.steps['Step-1'].frames[-1]
# Extract displacements from this frame
displacements = lastFrame.fieldOutputs['U']

# Consider just outputting displacement on surfaces
# this should be easy enough to get all nodes in the 
# surface plies thanks to Jean's nodesets, ideally
# I'd be able to filter further... Use the below code to find all nodesets in the assembly
# print(o1.rootAssembly.instances['M'].nodeSets)
# There aren't any appropriate node sets defined. Will need to do this in the mesh generation.
# Have a look at this later/ask Jean how this works

# Loop over all of the displacement values and store both the displacement vector and nodal coordiantes
# to numpy arrays
disp_out = np.empty([len(displacements.values),3])
nodes_out = np.empty([len(displacements.values),3])
#print(o1.rootAssembly.getNodeFromLabel(1))
for i, value in enumerate(displacements.values):
    # Reference points are defined directly on the assembly, not on an instance
    if (value.instance is None):
        disp_out[i,:] = value.data
        nodes_out[i,:] = [node.coordinates for node in o1.rootAssembly.nodes if node.label == value.nodeLabel][0]
    else:
        disp_out[i,:] = value.data
        nodes_out[i,:] = value.instance.getNodeFromLabel(value.nodeLabel).coordinates

# Write displacements and nodal coordinates to a csv file
np.savetxt(analysis + "_displacement.csv", disp_out, delimiter=",")
np.savetxt(analysis + "_nodes.csv", nodes_out, delimiter=",")