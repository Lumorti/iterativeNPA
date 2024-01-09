import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cp
from matplotlib.widgets import Slider
from multiprocessing import Pool
import os
import math
import random
import matplotlib
from scipy.spatial import ConvexHull
from itertools import product, combinations
from functools import reduce
matplotlib.use("GTK4Agg")

# Fixed seed
np.random.seed(1)

# 3 or 10
draw = True
thingsToDraw = [
        {"name": "box", "draw": True, "regen": False},
        {"name": "data", "draw": True, "regen": False},
        {"name": "sphube", "draw": False, "regen": False},
        {"name": "cone", "draw": False, "regen": False},
        {"name": "test", "draw": True, "regen": True},
    ]
fourthVal = 0.0
threads = 10
thresh = 0.30
tol = 0.01

# Sample random points, check if they are on the cone
def getPoints(a):
    points = []
    samples = 1000000
    while len(points) < samples:
        X = np.random.randn(5,5)
        for i in range(X.shape[0]):
            X[i,:] /= np.linalg.norm(X[i,:])
        X = X @ X.T
        if abs(X[1,3]) > thresh or abs(X[1,4]) > thresh or abs(X[2,3]) > thresh or abs(X[2,4]) > thresh:
            points.append([X[1,3], X[1,4], X[2,3], X[2,4]])
        if a == 0 and len(points) % (samples//10) == 0:
            print(100.0*float(len(points))/samples, "%")
    return points

# Check if the four variable point is in the full 10D cone
def checkPoint(var1, var2, var3, var4):
    X = cp.Variable((5,5), symmetric=True)
    cons = [
        X[0,0] == 1,
        X[1,1] == 1,
        X[2,2] == 1,
        X[3,3] == 1,
        X[4,4] == 1,
        X[1,3] == var1,
        X[1,4] == var2,
        X[2,3] == var3,
        X[2,4] == var4,
        X >> 0
    ]
    prob = cp.Problem(cp.Minimize(0), cons)
    try:
        prob.solve(solver=cp.MOSEK, mosek_params={"MSK_IPAR_NUM_THREADS": 1})
        return prob.status == "optimal"
    except:
        return False

# Same as above, but sample uniformly over the target vars
def getPointsUniform(a):
    points = []
    var4 = fourthVal
    pointsPer = 50
    count = 0
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in np.linspace(-1, 1, pointsPer):
            for var3 in np.linspace(-1, 1, pointsPer):
                count += 1
                if checkPoint(var1, var2, var3, var4):
                    points.append([var1, var2, var3])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    points = np.array(points)
    hull = ConvexHull(points)
    points = points[hull.vertices,:]
    return points

# Check if a point is in the squircle
def checkPointSquircle(x, y, z, s, r):
    s2r2 = s**2 / float(r**2)
    s4r4 = s**4 / float(r**4)
    x2 = x**2
    y2 = y**2
    z2 = z**2
    if x2 + y2 + z2 - s2r2*(x2*y2 + x2*z2 + y2*z2) + s4r4*x2*y2*z2 <= r**2:
        return True
    else:
        return False

# Get the point cloud representation of the squircle
def getPointsUniformSquircle(a):
    points = []
    pointsPer = 60
    count = 0
    s = 0.970
    r = 1.09
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                if checkPointSquircle(var1, var2, var3, s, r):
                    points.append([var1, var2, var3])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    points = np.array(points)
    hull = ConvexHull(points)
    points = points[hull.vertices,:]
    return points

# Check if a point is in the 3D SDP cone
def checkPointCone(x, y, z):
    X = [[1, x, y], 
         [0, 1, z], 
         [0, 0, 1]]
    X = np.array(X)
    for i in range(3):
        for j in range(0, i):
            X[i,j] = X[j,i]
    vals, vecs = np.linalg.eig(X)
    return np.all(vals > 0)

# Get the point cloud representation of the 3D SDP cone
def getPointsUniformCone(a):
    points = []
    pointsPer = 60
    count = 0
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                if checkPointCone(var1, var2, var3):
                    points.append([var1, var2, var3])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    points = np.array(points)
    hull = ConvexHull(points)
    points = points[hull.vertices,:]
    return points

# Check if a point is in the test domain TODO
def checkPointTest(x, y, z):
    # term = x**2 + y**2 + z**2 + 2*x*y*z
    # return term <= 1
    # X = [[1, term, term], 
         # [0, 1, term], 
         # [0, 0, 1]]
    x2 = x**2
    y2 = y**2
    z2 = z**2
    prod = x*y*z
    prod2 = prod**2
    sumPrd = x2*y2 + x2*z2 + y2*z2
    term = x2 + y2 + z2
    # X0 = [[1, term], 
         # [0, 1+sumPrd]]
    X0 = [[2, term, term], 
         [0, 1+sumPrd, term], 
         [0, 0, 1+sumPrd]]
    # X0 = [[1, x2, y2, z2], 
         # [0, 1+x2, 0, 0], 
         # [0, 0, 1+y2, 0],
         # [0, 0, 0, 1+z2]]
    X1 = [[1, x, y], 
         [0, 1, fourthVal*z], 
         [0, 0, 1]]
    X0 = np.array(X0)
    X1 = np.array(X1)
    X = fourthVal**2*X1 + (1-fourthVal**2)*X0
    # X = [[1, x, y, z], 
         # [0, 1, z, y], 
         # [0, 0, 1, x],
         # [0, 0, 0, 1]]
    for i in range(X.shape[0]):
        for j in range(0, i):
            X[i,j] = X[j,i]
    vals, vecs = np.linalg.eig(X)
    return np.all(vals > 0)

# Get the point cloud representation of the test region
def getPointsUniformTest(a):
    points = []
    pointsPer = 60
    count = 0
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                if checkPointTest(var1, var2, var3):
                    points.append([var1, var2, var3])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    points = np.array(points)
    hull = ConvexHull(points)
    points = points[hull.vertices,:]
    return points

# Sample random points in the cone
def writePoints(points, filename):
    print("Writing ", len(points), " points to file " + filename + "...")
    with open(filename, "w") as f:
        for p in points:
            for i in range(len(p)):
                f.write(str(p[i]))
                if i < len(p)-1:
                    f.write(",")
            f.write("\n")

# Load points from a file
def readPoints(filename):
    print("Loading points from file " + filename + "...")
    with open(filename, "r") as f:
        points = []
        for line in f:
            x = [float(i) for i in line.split(",")]
            points.append(x)
        print("Loaded ", len(points), " points")
        return np.array(points)

# For each thing to each draw or recalulate
pointArray = []
nameArray = []
for thingToDraw in thingsToDraw:

    # If told to draw the box from -1 to 1
    if thingToDraw["name"] == "box" and thingToDraw["draw"]:
        points = [
            [-1, -1, -1],
            [-1, -1,  1],
            [-1,  1, -1],
            [-1,  1,  1],
            [ 1, -1, -1],
            [ 1, -1,  1],
            [ 1,  1, -1],
            [ 1,  1,  1]
        ]
        pointArray.append(np.array(points))
        nameArray.append("Box")

    # If told to draw a standard-ish sphube
    elif thingToDraw["name"] == "sphube" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            pool2 = Pool(threads)
            points2 = pool2.map(getPointsUniformSquircle, range(threads))
            points2 = reduce(lambda a, b: np.concatenate((a,b)), points2)
            points2 = points2.reshape(-1, 3)
            writePoints(points2, "pointsSphube.txt")
        else:
            points2 = readPoints("pointsSphube.txt")
        if thingToDraw["draw"]:
            pointArray.append(points2)
            nameArray.append("Sphube")

    # If told to draw the 3D SDP cone
    elif thingToDraw["name"] == "cone" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            pool3 = Pool(threads)
            points3 = pool3.map(getPointsUniformCone, range(threads))
            points3 = reduce(lambda a, b: np.concatenate((a,b)), points3)
            points3 = points3.reshape(-1, 3)
            writePoints(points3, "pointsCone.txt")
        else:
            points3 = readPoints("pointsCone.txt")
        if thingToDraw["draw"]:
            pointArray.append(points3)
            nameArray.append("Cone")

    # If told to draw the test region
    elif thingToDraw["name"] == "test" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            pool3 = Pool(threads)
            points3 = pool3.map(getPointsUniformTest, range(threads))
            points3 = reduce(lambda a, b: np.concatenate((a,b)), points3)
            points3 = points3.reshape(-1, 3)
            writePoints(points3, "pointsTest.txt")
        else:
            points3 = readPoints("pointsTest.txt")
        if thingToDraw["draw"]:
            pointArray.append(points3)
            nameArray.append("Test Region")

    # If told to draw the data for a specific fourth value
    elif thingToDraw["name"] == "data" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            pool = Pool(threads)
            points = pool.map(getPointsUniform, range(threads))
            points = reduce(lambda a, b: np.concatenate((a,b)), points)
            points = points.reshape(-1, 3)
            writePoints(points, "points" + str(fourthVal) + ".txt")
        else:
            points = readPoints("points" + str(fourthVal) + ".txt")
        if thingToDraw["draw"]:
            pointArray.append(points)
            nameArray.append("Data")

# If testing the points
if not draw:
    pool = Pool(threads)
    chunkSize = len(points)//threads
    chunks = [points[i:i+chunkSize] for i in range(0, len(points), chunkSize)]
    numsNotPSD = pool.map(testPoints, chunks)
    totalChecks = sum([x[1] for x in numsNotPSD])
    totalNotPSD = sum([x[0] for x in numsNotPSD])
    if totalChecks == 0:
        print("No checks done")
    else:
        percentPSD = 100.0*(1.0-float(totalNotPSD)/totalChecks)
        print("Total checks: ", totalChecks)
        print("Percent correct: ", percentPSD)

# If drawing the points
else:
    print("Setting up plot...")
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("y")
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(-1.2, 1.2)
    ax.set_proj_type("ortho")
    styles = ["r-", "b-", "g-", "y-", "c-", "m-", "k-"]
    for i, points in enumerate(pointArray):
        style = styles[i]
        hull = ConvexHull(points[:,:])
        for s in hull.simplices:
            s = np.append(s, s[0])
            ax.plot(points[s, 0], points[s, 1], points[s, 2], style)
        ax.plot([], [], style, label=nameArray[i])
    fig.tight_layout()
    fig.legend()
    plt.show()

