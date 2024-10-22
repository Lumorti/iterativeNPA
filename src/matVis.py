#!/usr/env/bin python3
import sys
import matplotlib.pyplot as plt
import numpy as np

# Read the file given by the arg
lines = []
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()

# Search for the line "Moment matrix:"
mat = []
for i, line in enumerate(lines):
    if line.startswith("Moment matrix"):
        j = i + 1
        while j < len(lines) and len(lines[j].strip()) > 0 and not lines[j].startswith("Solving"):
            lines[j] = lines[j].replace("+", "").replace("-", "").replace("<", "").replace(">", "").strip()
            mat.append(lines[j].split())
            j += 1
        break

print(mat)
topRow = mat[0]

# Count how many times each element appears
repeated = np.zeros((len(mat), len(mat[0])))
for i in range(len(mat)):
    for j in range(i+1, len(mat[0])):
        if mat[i][j] in ["A1", "A2", "B1", "B2", "A1B1", "A1B2", "A2B1", "A2B2", "A1B3", "A2B3", "A3B1", "A3B2"]:
            for k in range(len(mat)):
                for l in range(k+1, len(mat[0])):
                    if mat[i][j] == mat[k][l]:
                        repeated[i][j] += 1
                        repeated[j][i] += 1

# cost is such that we want as much grouping as possible in the top left
def cost(m):
    c = 0
    # halfSize = int(len(m)/2)
    # for i in range(halfSize):
        # for j in range(i+1, halfSize):
            # if m[i][j] > 0:
                # c -= 1
    # for i in range(halfSize, len(m)):
        # for j in range(i+1, len(m)):
            # if m[i][j] > 0:
                # c -= 1
    matSize = 60
    for i in range(matSize):
        for j in range(i+1, matSize):
            if m[i][j] > 0:
                c -= matSize
    return c

# Scatter plot the above matrix
# fig, ax = plt.subplots()
# ax.scatter(range(len(mat)), range(len(mat[0])), s=repeated.flatten())
# plt.show()

print(mat[0])

# Do many swaps trying to minimize the cost
if len(sys.argv) >= 3:
    prevCost = cost(repeated)
    for itNum in range(10000):
        i = np.random.randint(1, len(repeated))
        j = np.random.randint(1, len(repeated))
        repeated[[i, j]] = repeated[[j, i]]
        repeated[:, [i, j]] = repeated[:, [j, i]]
        topRow[i], topRow[j] = topRow[j], topRow[i]
        currentCost = cost(repeated)
        print(itNum, currentCost)
        if currentCost > prevCost:
            repeated[[i, j]] = repeated[[j, i]]
            repeated[:, [i, j]] = repeated[:, [j, i]]
            topRow[i], topRow[j] = topRow[j], topRow[i]
        else:
            prevCost = currentCost
    print(" ".join(topRow[:60]))

# Plot the new
X,Y = np.meshgrid(np.arange(repeated.shape[1]), np.arange(repeated.shape[0]))
plt.scatter(X.flatten(), Y.flatten(), c=repeated.flatten())
plt.colorbar()
plt.gca().invert_yaxis()
plt.show()
