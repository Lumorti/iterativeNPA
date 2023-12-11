
# We have an input file in the format
# Angle: 0.123   L1: 5.5   L2: 6.7   L3: 7.8

# Parse the input into a list
data = []
import sys
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
    for line in lines:
        if "Angle:" in line:
            split_line = line.split()
            data.append([float(split_line[1]), float(split_line[3]), float(split_line[5]), float(split_line[7])])
        elif "Vars:" in line and len(data) > 10:
            break

print(data)

# Plot the data
import matplotlib.pyplot as plt
import numpy as np

# for each angle, plot the level 1 and level 2 values
# this should be done around some central point, such that it creates two circles
for i in range(len(data)):
    angle = data[i][0]
    l1 = data[i][1]
    l2 = data[i][2]
    l3 = data[i][3]
    plt.plot(l1*np.cos(angle), l1*np.sin(angle), 'rx')
    plt.plot(l2*np.cos(angle), l2*np.sin(angle), 'bx')
    plt.plot(l3*np.cos(angle), l3*np.sin(angle), 'gx')
plt.show()






