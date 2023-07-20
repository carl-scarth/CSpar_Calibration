# Plots the model output used to train the emulator
import numpy as np
import meshio

shell_mesh = True
# Open the file in question
output_file = 'LHSDesign40x4_1_fixed_100kN.csv'
# output_file = 'LHSDesign50x3_1_displacements.csv'
# output_file = 'spring_study_displacements.csv'
# Open the output files
if shell_mesh:
    node_file = "CSpar_sam_shell_mesh_nodes.csv" # Nodes of nominal input (ignores geometric uncertainty)
    element_file = "CSpar_sam_shell_mesh_elements.csv" # Element connectivity
else:
    node_file = "CSpar_sam_mesh_nodes.csv" # Nodes of nominal input (ignores geometric uncertainty)
    element_file = "CSpar_sam_mesh_elements.csv" # Element connectivity


# Read in the element and node definitions
elements = np.loadtxt(element_file, dtype = int, delimiter = ',')
nodes = np.loadtxt(node_file, delimiter = ',')
# first column is just an index, which I don't need here as everything is (I think) in ascending order
# with no gaps
nodes = nodes[:,1:4]

# Read in displacement data
#displacement = np.loadtxt(output_file, delimiter = ',', skiprows = 0)
# Now I've added headers to the simulation data - could do something which spots the error perhaps
displacement = np.loadtxt(output_file, delimiter = ',', skiprows = 1)
with open(output_file,'r') as f:
    output_header = f.readline().strip().split(',')

# Define list of indice of the nodes which define each face of the brick (i.e. which column of the connectivity)
# Note that first column of elements is just the element number, hence indexing starts at 1
face_nodes = [[1,2,3,4],[5,6,7,8],[1,2,6,5],[2,3,7,6],[3,4,8,7],[4,1,5,8]]
faces = np.empty([0,4], dtype = int)
for face in face_nodes:
    elements_face = elements[:,face]
    faces = np.concatenate((faces,elements_face))


# Re-writing from abaqus to python indexing for meshio
faces = faces - 1
# Delete first two points which are reference points - these are defined directly on the assembly, rather than
# the part instance, and so have a different numbering system to that stored in the connectivity file
nodes = nodes[2:,:]

# Create dictionary using basis data
output_dict = {}
#for i in range(displacement.shape[1]//3):    
    # output_dict["u_" + str(i+1)] = displacement[:,3*i]
    # output_dict["v_" + str(i+1)] = displacement[:,3*i+1]
    # output_dict["w_" + str(i+1)] = displacement[:,3*i+2]

for i, entry in enumerate(output_header):
    output_dict[entry] = displacement[2:,i]

print(output_dict)
# Get data in the correct format for meshio
meshio.Mesh(points = nodes, cells = [("quad",faces)], point_data = output_dict).write("DoE_Displacements_40x4_1.vtk", file_format="vtk")