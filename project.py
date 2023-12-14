import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cp

# Plot the semidefinite cone defined by
# X = [[1, x, y], [x, 1, z], [y, z, 1]]
np.random.seed(1)

# Optimize to find points on the edge of the cone
# points = []
# samples = 30
# for c1 in np.linspace(-1, 1, samples):
    # print(c1)
    # for c2 in np.linspace(-1, 1, samples):
        # for c3 in np.linspace(-1, 1, samples):
            # c = np.array([c1, c2, c3])
            # x = cp.Variable((3,3), symmetric=True)
            # obj = cp.Minimize(c[0]*x[0,1] + c[1]*x[0,2] + c[2]*x[1,2])
            # constraints = [
                    # x >> 0,
                    # x[0,0] == 1, 
                    # x[1,1] == 1, 
                    # x[2,2] == 1,
                # ]
            # prob = cp.Problem(obj, constraints)
            # prob.solve()
            # points.append(x.value)

# Sample random points, check if they are on the cone
points = []
maxTries = 100000000000
samples = 1000
for i in range(maxTries):

    # Check if the matrix is positive semidefinite
    # x = np.random.randn(3)
    # X = np.array([[1, x[0], x[1]], [x[0], 1, x[2]], [x[1], x[2], 1]])
    x = np.random.randn(10)
    x[8] = 1
    X = np.array([[1, x[0], x[1], x[2], x[3]], 
                  [x[0], 1, x[4], x[5], x[6]], 
                  [x[1], x[4], 1, x[7], x[8]],
                  [x[2], x[5], x[7], 1, x[9]],
                  [x[3], x[6], x[8], x[9], 1]])
    vals, vecs = np.linalg.eig(X)
    if np.all(vals >= 0):
        points.append(x)
        if len(points) > samples:
            break

# Plot the cone
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_xlim3d(-2, 2)
ax.set_ylim3d(-2, 2)
ax.set_zlim3d(-2, 2)

# Plot the points
points = np.array(points)
# ax.scatter(points[:,0,1], points[:,0,2], points[:,1,2], color='r')
xPoints = []
yPoints = []
zPoints = []
for i in range(len(points)):
    xPoints.append(points[i,5])
    yPoints.append(points[i,6])
    zPoints.append(points[i,7])
# ax.scatter(points[:,5], points[:,6], points[:,7], color='r')
ax.scatter(xPoints, yPoints, zPoints, color='r')

# Show the plot
plt.show()








