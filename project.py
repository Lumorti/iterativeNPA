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
thingsToDraw = [
        {"name": "box",    "draw": False,  "regen": False, "check": False},
        {"name": "data",   "draw": True,  "regen": False, "check": False},
        {"name": "sphube", "draw": False, "regen": False, "check": False},
        {"name": "cone",   "draw": False, "regen": False, "check": False},
        {"name": "test",   "draw": True,  "regen": True,  "check": False},
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
def checkPointSphube(x, y, z):
    s = 0.970
    r = 1.09
    s2r2 = s**2 / float(r**2)
    s4r4 = s**4 / float(r**4)
    x2 = x**2
    y2 = y**2
    z2 = z**2
    if x2 + y2 + z2 - s2r2*(x2*y2 + x2*z2 + y2*z2) + s4r4*x2*y2*z2 <= r**2:
        return True
    else:
        return False

# Check if a point is in the box
def checkPointBox(x, y, z):
    return abs(x) <= 1 and abs(y) <= 1 and abs(z) <= 1

# Get the point cloud representation of the squircle
def getPointsUniformSphube(a):
    points = []
    pointsPer = 60
    count = 0
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                if checkPointSphube(var1, var2, var3):
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
    a = fourthVal
    x2 = x**2
    y2 = y**2
    z2 = z**2
    a2 = a**2
    s = x2*y2 + z2*(x2 + y2) + a2*(x2 + y2 + z2)
    t = x2 + y2 + z2 + a2
    # X0 = [[1, term], 
         # [0, 1+sumPrd]]
    X0 = [[2, t, t], 
         [0, 1+s, t], 
         [0, 0, 1+s]]
    # X0 = [[1, x2, y2, z2], 
         # [0, 1+x2, 0, 0], 
         # [0, 0, 1+y2, 0],
         # [0, 0, 0, 1+z2]]
    X1 = [[1, x, y], 
         [0, 1, a*z], 
         [0, 0, 1]]
    X0 = np.array(X0)
    X1 = np.array(X1)
    X = a*a*X1 + (1-a*a)*X0
    # X = X0
    # b = x*y*z*a
    # c = x*x + y*y + z*z + a*a
    X = [[1, 0, x, y], 
         [0, 1, z, a], 
         [0, 0, 1, 0],
         [0, 0, 0, 1]]
    X = np.array(X)
    for i in range(X.shape[0]):
        for j in range(0, i):
            if X[i,j] == 0:
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
            points2 = pool2.map(getPointsUniformSphube, range(threads))
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

# Sample many points and see which regions they are in
def checkAll(a):
    counts = {}
    for thing in thingsToDraw:
        for thing2 in thingsToDraw:
            if thing["check"] and thing2["check"]:
                counts[(thing["name"], thing2["name"])] = 0
    asDict = {}
    for thing in thingsToDraw:
        asDict[thing["name"]] = thing["check"]
    pointsPer = 20
    count = 0
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                results = {}
                if asDict["box"]:
                    results["box"] = checkPointBox(var1, var2, var3)
                if asDict["test"]:
                    results["test"] = checkPointTest(var1, var2, var3)
                if asDict["cone"]:
                    results["cone"] = checkPointCone(var1, var2, var3)
                if asDict["sphube"]:
                    results["sphube"] = checkPointSphube(var1, var2, var3)
                if asDict["data"]:
                    results["data"] = checkPoint(var1, var2, var3, fourthVal)
                for thing in thingsToDraw:
                    for thing2 in thingsToDraw:
                        if thing["check"] and thing2["check"]:
                            counts[(thing["name"], thing2["name"])] += (results[thing["name"]] == results[thing2["name"]])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    for thing in thingsToDraw:
        for thing2 in thingsToDraw:
            if thing["check"] and thing2["check"]:
                counts[(thing["name"], thing2["name"])] *= 100.0 / float(pointsPer**3)
    return counts

# If we need to do any analysis
checkAny = False
for thing in thingsToDraw:
    if thing["check"]:
        checkAny = True
        break
if checkAny:
    pool = Pool(threads)
    results = pool.map(checkAll, range(threads))
    results = reduce(lambda a, b: {k: a.get(k, 0) + b.get(k, 0) for k in set(a) | set(b)}, results)
    print(results)

# Draw everything in 3D if there's something to draw
if len(pointArray) > 0:
    print("Setting up plot...")
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
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

