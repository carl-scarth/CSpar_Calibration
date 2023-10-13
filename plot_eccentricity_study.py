import json
import matplotlib.pyplot as plt
import pandas as pd

# Load in simulation data from json
with open("eccentricity_study_json_-350.0_max_inc=0.05.json",'r') as f:
    data = f.readline()
displacements = json.loads(data)

# Load in input values
inputs = pd.read_csv("eccentricity_study.csv")
x_spring = inputs["x_spring"]

# Plot force-displacement of reference point at the spar tip. Second node is 
# the reference point. Longitudinal displacement is the 3rd coordinate u_3.
# Negate as negative displacement
ref_w = [[-frame["Displacements"][1][2] for frame in sample["Frame"]]for sample in displacements["Sample"]]
# Get the reaction force
RF_z = [[frame["RFs"][2] for frame in sample["Frame"]]for sample in displacements["Sample"]]

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)
legend_label = ["x_spring = "+str(eccentricity)+"mm" for eccentricity in x_spring]
lines = []
for i, sample in enumerate(ref_w):
    ax.plot(sample, RF_z[i], label = legend_label[i])

ax.legend()
label_font = {'family': 'serif', 'size': 16,}
ax.set_ylabel("Compressive Force (kN)", fontdict = label_font)
ax.set_xlabel("End displacement (mm)", fontdict = label_font)

# ax.legend(lines)
plt.show()