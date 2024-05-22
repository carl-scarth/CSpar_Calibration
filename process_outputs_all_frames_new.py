from odbAccess import *
import numpy as np
import json

# output to command prompt
sys.stdout = sys.__stdout__

# Potentially pass the analysis name as a variable at the command line.
analysis = sys.argv[-1]
fix_to_ground = bool(sys.argv[-2])

# Open odb file. Assumes the scrip is run from the parent directory, and not from the Abaqus subfolder
o1 = session.openOdb(name = analysis + '.odb')

# Get last frame of the load step.
# lastFrame = o1.steps['Step-1'].frames[-1]

# Loop over all of the displacement values and store both the displacement vector and nodal coordiantes
# to numpy arrays. Only need to output one frame for nodes (if that)

n_nodes = len(o1.steps['Step-1'].frames[0].fieldOutputs['U'].values) # Number of nodes
n_frames = len(o1.steps['Step-1'].frames) # Number of Frames

# The structure would be nicer if stored as json files. Maybe also numpy has an option for outputting 3D arrays to csv. Might be worth exploring...
disp_out = np.empty([n_nodes,3*n_frames]) 
nodes_out = np.empty([n_nodes,3])
nodes_list = []
RFs_out = np.empty([n_frames,3])
increments = []
# out_dict = {"Frame" : [], "Increment" : [], "Reaction_Forces" : [], "Displacements" : [], "Nodes" : []}
out_dict = {"Frame" : [], "Nodes" : []}
# Will need to change below if adding in more steps
for i, frame in enumerate(o1.steps['Step-1'].frames):
    # Extract displacements from this frame
    displacements = frame.fieldOutputs['U']
    disp_i = []
    for j, value in enumerate(displacements.values):
        disp_out[j,3*i:3*(i+1)] = value.data
        disp_i.append(value.data.tolist())
        if i == 0:
            nodes_out[j,:] = value.instance.getNodeFromLabel(value.nodeLabel).coordinates
            nodes_list.append(nodes_out[j,:].tolist())

    # Store reaction forces at fixed end of the spar
    if fix_to_ground:
        RFs = frame.fieldOutputs['RF'].getSubset(region = o1.rootAssembly.instances['M'].nodeSets['FIXED_SUP_RP'])
    else:
        RFs = frame.fieldOutputs['RF'].getSubset(region = o1.rootAssembly.instances['M'].nodeSets['GROUND'])
    RFs_out[i,:] = RFs.values[0].data
    increments.append(frame.frameValue)

    # Write dictionary for output in json format
    out_dict["Frame"].append({"Index" : i,
                              "Increment" : frame.frameValue,
                              "RFs" : RFs.values[0].data.tolist(),
                              "Displacements" : disp_i})
    out_dict["Nodes"] = nodes_list
print(RFs_out)
# Write displacements and nodal coordinates to a csv file
np.savetxt(analysis + "_displacement.csv", disp_out, delimiter=",")
np.savetxt(analysis + "_nodes.csv", nodes_out, delimiter=",")
np.savetxt(analysis + "_RFs.csv", RFs_out, delimiter=",")
np.savetxt(analysis + "_increments.csv", increments, delimiter=",")

# Is there a nicer way of formatting these strings for readability? 
with open(analysis + "_output.json",'w') as f:
    f.write(json.dumps(out_dict))