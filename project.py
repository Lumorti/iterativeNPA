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
        ("name": "box", "draw": True, "regen": False)
        ("name": "sphube", "draw": True, "regen": False)
        ("name": "cone", "draw": True, "regen": False)
        ("name": "data", "draw": True, "regen": False)
    ]
mode = 20
fourthVal = 1.0
threads = 10
thresh = 0.30
tol = 0.01

# Sample random points, check if they are on the cone
def getPoints(a):
    points = []
    if mode == 3:
        samples = 10000
    elif mode == 10:
        samples = 1000000
    while len(points) < samples:
        if mode == 3:
            X = np.random.randn(3,3)
        elif mode == 10:
            X = np.random.randn(5,5)
        for i in range(X.shape[0]):
            X[i,:] /= np.linalg.norm(X[i,:])
        X = X @ X.T
        if mode == 3:
            if abs(X[0,1]) > thresh or abs(X[0,2]) > thresh or abs(X[1,2]) > thresh:
                points.append([X[0,1], X[0,2], X[1,2]])
        elif mode == 10:
            if abs(X[1,3]) > thresh or abs(X[1,4]) > thresh or abs(X[2,3]) > thresh or abs(X[2,4]) > thresh:
                points.append([X[1,3], X[1,4], X[2,3], X[2,4]])
        if a == 0 and len(points) % (samples//10) == 0:
            print(100.0*float(len(points))/samples, "%")
    return points

# Form an SDP with cvxpy and see if its feasible
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

def checkPointCone(var1, var2, var3):
    X = cp.Variable((3,3), symmetric=True)
    cons = [
        X[0,0] == 1,
        X[1,1] == 1,
        X[2,2] == 1,
        X[0,1] == var1,
        X[0,2] == var2,
        X[1,2] == var3,
        X >> 0
    ]
    prob = cp.Problem(cp.Minimize(0), cons)
    try:
        prob.solve(solver=cp.MOSEK, mosek_params={"MSK_IPAR_NUM_THREADS": 1})
        return prob.status == "optimal"
    except:
        return False

# Check if a point is in the squircle
def checkPointSquircle(var1, var2, var3, p, r):
    if var1**p + var2**p + var3**p <= r**p:
        return True
    else:
        return False

# Same as above, but sample uniformly over the target vars TODO
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

# Get the point cloud representation of the squircle
def getPointsUniformSquircle(a):
    points = []
    pointsPer = 60
    count = 0
    p = 4
    r = 1.15
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                if checkPointSquircle(var1, var2, var3, p, r):
                    points.append([var1, var2, var3])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    points = np.array(points)
    hull = ConvexHull(points)
    points = points[hull.vertices,:]
    return points

# Get the point cloud representation of the squircle
def getPointsUniformCone(a):
    points = []
    pointsPer = 60
    count = 0
    p = 4
    r = 1.15
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

# Test whether my simplification is a subset
def testPoints(points):
    np.random.seed()
    numInequal = 0 
    totalDone = 0
    for j, point in enumerate(points):
        A1B1 = point[0]
        A1B2 = point[1]
        A2B1 = point[2]
        A2B2 = point[3]
        A1B1sq = A1B1**2
        A1B2sq = A1B2**2
        A2B1sq = A2B1**2
        A2B2sq = A2B2**2
        Ys = []
        if abs(A2B1) >= 0.95: 
            Ys.append(np.array([[1,      A1B1,  A1B2],
                                [0,      1,     A2B2*A2B1],
                                [0,      0,     1]]))
        if len(Ys) > 0:
            YIsPSD = True
            for l in range(len(Ys)):
                for i in range(len(Ys[l])):
                    for k in range(i+1, len(Ys[l])):
                        if Ys[l][k][i] == 0:
                            Ys[l][k][i] = Ys[l][i][k]
                vals, vecs = np.linalg.eig(Ys[l])
                YIsPSD = np.all(vals >= 0.0)
                if not YIsPSD:
                    break
            totalDone += 1
            if not YIsPSD:
                numInequal += 1
                print("X: ")
                print(X)
                print(XIsPSD)
                print(min(val))
                print("Ys: ")
                for Y in Ys:
                    print(Y)
                print(YIsPSD)
                print(min(vals))
                break
        if j % (len(points)//10) == 0:
            print(100.0*float(j)/len(points), "%")
    return numInequal, totalDone

# TODO simplify everything
pointArray = []
for thingToDraw in thingsToDraw:

    # If told to draw the box from -1 to 1
    if thingToDraw.name == "box" and thingToDraw.draw:
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

    # If told to draw a standard-ish sphube
    elif thingToDraw.name == "sphube" and thingToDraw.draw:
        if thingToDraw.regen:
            pool2 = Pool(threads)
            points2 = pool2.map(getPointsUniformSquircle, range(threads))
            points2 = reduce(lambda a, b: np.concatenate((a,b)), points2)
            points2 = points2.reshape(-1, 3)
            pointArray.append(points2)
            writePoints(points2, "pointsSphube.txt")
        else:
            points2 = readPoints("pointsSphube.txt")
            pointArray.append(points2)

    # If told to draw the 3D SDP cone
    elif thingToDraw.name == "cone" and thingToDraw.draw:
        if thingToDraw.regen:
            pool3 = Pool(threads)
            points3 = pool3.map(getPointsUniformCone, range(threads))
            points3 = reduce(lambda a, b: np.concatenate((a,b)), points3)
            points3 = points3.reshape(-1, 3)
            pointArray.append(points3)
            writePoints(points3, "pointsCone.txt")
        else:
            points3 = readPoints("pointsCone.txt")
            pointArray.append(points3)

    # If told to draw the data for a specific fourth value
    elif thingToDraw.name == "data" and thingToDraw.draw:
        if thingToDraw.regen:
            pool = Pool(threads)
            points = pool.map(getPointsUniform, range(threads))
            points = reduce(lambda a, b: np.concatenate((a,b)), points)
            points = points.reshape(-1, 3)
            writePoints(points, "points" + str(fourthVal) + ".txt")
        else:
            points = readPoints("points" + str(fourthVal) + ".txt")
            pointArray.append(points)

# If told to regenerate the point list
if regeneratePointList:

    # If asking for either
    if mode == 3 or mode == 10:

        # Get the points in parallel
        print("Generating points...")
        pool = Pool(threads)
        points = pool.map(getPoints, range(threads))
        points = np.array(points)
        if mode == 3:
            points = points.reshape(-1, 3)
        elif mode == 10:
            points = points.reshape(-1, 4)

        # Write them to file
        writePoints(points, "points"+str(mode)+".txt")

    # If asking for both
    elif mode == 30:

        # Same as above, but switch the global
        print("Generating points for mode 3...")
        mode = 3
        pool = Pool(threads)
        points3 = pool.map(getPoints, range(threads))
        points3 = np.array(points3)
        points3 = points3.reshape(-1, 3)
        writePoints(points3, "points3.txt")
        print("Generating points for mode 10...")
        mode = 10
        pool2 = Pool(threads)
        points10 = pool2.map(getPoints, range(threads))
        points10 = np.array(points10)
        points10 = points10.reshape(-1, 4)
        writePoints(points10, "points10.txt")
        mode = 30

    elif mode == 20:

        # Uniformly sample over the projected space, then see if they are in the OG cone
        pool = Pool(threads)
        points = pool.map(getPointsUniform, range(threads))
        points = reduce(lambda a, b: np.concatenate((a,b)), points)
        points = points.reshape(-1, 3)
        writePoints(points, "points" + str(fourthVal) + ".txt")

        # Same as above but for squircle
        if fourthVal >= 0.5:
            pool2 = Pool(threads)
            points2 = pool2.map(getPointsUniformCone, range(threads))
            points2 = reduce(lambda a, b: np.concatenate((a,b)), points2)
            points2 = points2.reshape(-1, 3)
            writePoints(points2, "pointsCone.txt")
        else:
            pool2 = Pool(threads)
            points2 = pool2.map(getPointsUniformSquircle, range(threads))
            points2 = reduce(lambda a, b: np.concatenate((a,b)), points2)
            points2 = points2.reshape(-1, 3)

# Otherwise load them from file
else:

    # If just asking for one type
    if mode == 3 or mode == 10:
        points = readPoints("points"+str(mode)+".txt")

    elif mode == 20:
        points = readPoints("points" + str(fourthVal) + ".txt")
        if fourthVal >= 0.5:
            points2 = readPoints("pointsCone.txt")
        else:
            pool2 = Pool(threads)
            points2 = pool2.map(getPointsUniformSquircle, range(threads))
            points2 = reduce(lambda a, b: np.concatenate((a,b)), points2)
            points2 = points2.reshape(-1, 3)

    # Or if asking for both
    elif mode == 30:
        points3 = readPoints("points3.txt")
        points10 = readPoints("points10.txt")

# If testing the points
if not draw:

    # Test the points in parallel
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

    # 3D
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
    if mode == 3 or mode == 10:
        print("Processing points...")
        points = np.array(points)
        if mode == 10:
            newPoints = []
            for i in range(len(points)):
                if abs(points[i,3]-fourthVal) <= tol:
                    newPoints.append([points[i,0], points[i,1], points[i,2]])
            points = np.array(newPoints)
        print("Calculating convex hull...")
        hull = ConvexHull(points[:,:])
        print("Drawing...")
        for s in hull.simplices:
            s = np.append(s, s[0])
            ax.plot(points[s, 0], points[s, 1], points[s, 2], "r-")
    elif mode == 20:
        print("Processing points...")
        points = np.array(points)
        if mode == 10:
            newPoints = []
            for i in range(len(points)):
                if abs(points[i,3]-fourthVal) <= tol:
                    newPoints.append([points[i,0], points[i,1], points[i,2]])
            points = np.array(newPoints)
        print("Calculating convex hull...")
        hull = ConvexHull(points[:,:])
        hull2 = ConvexHull(points2[:,:])
        print("Drawing...")
        for s in hull.simplices:
            s = np.append(s, s[0])
            ax.plot(points[s, 0], points[s, 1], points[s, 2], "r-")
        for s in hull2.simplices:
            s = np.append(s, s[0])
            ax.plot(points2[s, 0], points2[s, 1], points2[s, 2], "b-")
    elif mode == 30:
        print("Processing points...")
        points3 = np.array(points3)
        points10 = np.array(points10)
        newPoints = []
        for i in range(len(points10)):
            if abs(points10[i,3]-fourthVal) <= tol:
                newPoints.append([points10[i,0], points10[i,1], points10[i,2]])
        points10 = np.array(newPoints)
        print("Calculating convex hull...")
        hull3 = ConvexHull(points3[:,:])
        hull10 = ConvexHull(points10[:,:])
        print("Drawing...")
        for s in hull3.simplices:
            s = np.append(s, s[0])
            ax.plot(points3[s, 0], points3[s, 1], points3[s, 2], "r-")
        for s in hull10.simplices:
            s = np.append(s, s[0])
            ax.plot(points10[s, 0], points10[s, 1], points10[s, 2], "b-")
    r = [-1, 1]
    for s, e in combinations(np.array(list(product(r,r,r))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s,e), color="g")
    fig.tight_layout()
    plt.show()









