
# Plot the semidefinite cone defined by
# X = [[1, x, y], [x, 1, z], [y, z, 1]]

import numpy as np
import matplotlib.pyplot as plt
from cvxpy import *

# Pick a random objective function
np.random.seed(1)

# Optimize to find points on the edge of the cone
points = []
for i in range(1000):
    c = np.random.randn(3, 1)
    x = Semidef(3)
    obj = Maximize(c.T*x)
    constraints = [x[0, 0] == 1, x[1, 1] == 1, x[2, 2] == 1]
    prob = Problem(obj, constraints)
    prob.solve()
    x = x.value
    points.append(x)

# Plot the cone
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_xlim3d(-1, 1)
ax.set_ylim3d(-1, 1)
ax.set_zlim3d(-1, 1)

# Plot the points
points = np.array(points)
ax.scatter(points[:, 0, 0], points[:, 1, 1], points[:, 2, 2], s=1)

# Show the plot
plt.show()








