import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cp
from matplotlib.widgets import Slider
from multiprocessing import Pool
import os

# Fixed seed
np.random.seed(1)

# Sample random points, check if they are on the cone
def getPoints(a):
    points = []
    maxTries = 100000000000
    samples = 10000
    for i in range(maxTries):
        x = np.random.randn(10)
        X = np.array([[1, x[0], x[1], x[2], x[3]], 
                      [x[0], 1, x[4], x[5], x[6]], 
                      [x[1], x[4], 1, x[7], x[8]],
                      [x[2], x[5], x[7], 1, x[9]],
                      [x[3], x[6], x[8], x[9], 1]])
        vals, vecs = np.linalg.eig(X)
        while np.any(vals < 0):
            vals[vals < 0] = 0
            X = vecs @ np.diag(vals) @ vecs.T
            for i in range(X.shape[0]):
                X[i,i] = 1
            vals, vecs = np.linalg.eig(X)
        x = [X[0,1], X[0,2], X[0,3], X[0,4], X[1,2], X[1,3], X[1,4], X[2,3], X[2,4], X[3,4]]
        points.append(x)
        if len(points) >= samples:
            break
        if a == 0 and i % 1000 == 0:
            print(float(i)/maxTries*100, '%')
    return points

# Get the points in parallel
threads = 8
pool = Pool(threads)
points = pool.map(getPoints, range(threads))
points = np.array(points)
points = points.reshape(-1, 10)
print(points.shape)

# Plot the cone
fig = plt.figure()
ax = plt.gca()
ax.set_xlabel('x')
ax.set_ylabel('y')
axSlider1 = fig.add_axes([0.2, 0.95, 0.65, 0.03])
axSlider2 = fig.add_axes([0.2, 0.90, 0.65, 0.03])
# ax = fig.gca(projection='3d')
# ax.set_zlabel('z')
# ax.set_xlim3d(-2, 2)
# ax.set_ylim3d(-2, 2)
# ax.set_zlim3d(-2, 2)

# Plot the points
val1 = 0
val2 = 0
points = np.array(points)
def update():
    global val1, val2
    xPoints = []
    yPoints = []
    for i in range(len(points)):
        if abs(points[i,7]-val1) <= 1e-1 and abs(points[i,8]-val2) <= 1e-1:
            xPoints.append(points[i,5])
            yPoints.append(points[i,6])
    ax.clear()
    ax.scatter(xPoints, yPoints, color='r', s=2)
def setVal1(val):
    global val1
    val1 = val
    update()
def setVal2(val):
    global val2
    val2 = val
    update()

slider1 = Slider(axSlider1, 'Val 1', -1.5, 1.5, valinit=1, valfmt='%f')
slider2 = Slider(axSlider2, 'Val 2', -1.5, 1.5, valinit=1, valfmt='%f')
slider1.on_changed(setVal1)
slider2.on_changed(setVal2)
update()

# Show the plot
plt.show()








