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
from scipy.optimize import minimize
from itertools import product, combinations
from functools import reduce
matplotlib.use("GTK4Agg")

# Fixed seed
np.random.seed(1)
random.seed(0)

# 3 or 10
thingsToDraw = [
        {"name": "box",    "draw": False,  "regen": False, "check": False},
        {"name": "data",   "draw": True,  "regen": False, "check": False},
        {"name": "sphube", "draw": False, "regen": False, "check": False},
        {"name": "cone",   "draw": False, "regen": False, "check": False},
        {"name": "test",   "draw": False,  "regen": False,  "check": False},
    ]
optimizeSDP = True
fourthVal = 1.0
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
    print("Generating points for fourth variable = ", var4)
    pointsPer = 40
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

# Check if a point is in the test domain
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

# Write points to a file
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

pointArray = []
nameArray = []

# Check if a point is in the a custom cone 
def checkPointCustom(vec, point, matSize):
    X = np.eye(matSize)
    vecInd = 0
    for i in range(matSize):
        for j in range(0, i):
            X[j,i] = vec[vecInd]
            vecInd += 1
            for el in point:
                X[j,i] += vec[vecInd]*el
                vecInd += 1
            X[i,j] = X[j,i]
    vals, vecs = np.linalg.eig(X)
    return np.all(vals > 0)

# Cost function for fitting the cone
def cost(vec, points, matSize):
    dist = 0
    for point in points:
        X = np.eye(matSize)
        vecInd = 0
        for i in range(matSize):
            for j in range(0, i):
                X[j,i] = vec[vecInd]
                vecInd += 1
                for el in point:
                    X[j,i] += vec[vecInd]*el
                    vecInd += 1
                X[i,j] = X[j,i]
        vals, vecs = np.linalg.eig(X)
        dist += abs(np.min(vals))
    dist /= len(points)
    print(dist)
    return dist

# Get the point cloud representation of the test region
def getPointsUniformCustom(args):
    coeffs, matSize, a = args
    points = []
    pointsPer = 60
    count = 0
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                if checkPointCustom(coeffs, [var1, var2, var3], matSize):
                    points.append([var1, var2, var3])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    points = np.array(points)
    hull = ConvexHull(points)
    points = points[hull.vertices,:]
    return points

# Get the point cloud representation of the test region
def getPointsUniformCustom4(args):
    coeffs, matSize, a = args
    points = []
    pointsPer = 60
    count = 0
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                if checkPointCustom(coeffs, [var1, var2, var3, fourthVal], matSize):
                    points.append([var1, var2, var3])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    points = np.array(points)
    hull = ConvexHull(points)
    points = points[hull.vertices,:]
    return points

# Brute force a section of the search space
def bruteForce(args):
    matSize, numVars, points, thread = args
    res = 10000000
    perThread = (2**numVars)//threads
    start = thread*perThread
    end = (thread+1)*perThread
    for i in range(start, end):
        vec = [0 for i in range(numVars)]
        for j in range(numVars):
            vec[j] = (i >> j) % 2
        res = min(res, cost(vec, points, matSize))
        if thread == 0 and i % 100 == 0:
            print(100.0*float(i)/perThread, "%, min so far = ", res)
    return res

# If told to optimize to minimize the SDP hull TODO
if optimizeSDP:
    matSize = 5
    allPoints = readPoints("data/pointsAll.csv")
    points = []
    pointsToSample = 1000000000
    numVars = (matSize*(matSize+1)//2-matSize) * (allPoints.shape[1]+1)
    vec = [1 for i in range(numVars)]
    # vec = np.random.rand(numVars)
    # vec = np.array([-0.21806701, -0.30024651, -0.2024949 , -0.15697767, -0.12269202,
       # -0.482077  ,  0.26369564,  0.28087138, -0.21673742, -0.02869441,
       # -0.22670792,  0.08536931, -0.09920705, -0.0567977 , -0.23431801,
       # -0.07499742, -0.50900859,  0.13968817,  0.15751262,  0.40836086,
        # 0.50056849,  0.09140297, -0.13158019,  0.47424429,  0.12773484,
       # -0.22482994, -0.38876055,  0.31887217, -0.06631954, -0.01572061,
       # -0.27017813, -0.09920275,  0.24522633, -0.15417556,  0.32584317,
       # -0.09420153, -0.02195686,  0.27172932, -0.00135903,  0.05110155])
    pointsToSample = min(pointsToSample, len(allPoints))
    for i in range(pointsToSample):
        rand = random.randint(0, len(allPoints)-1)
        points.append(allPoints[rand])
        np.delete(allPoints, rand)

    # x = cp.Variable(numVars)
    # mats = [cp.Variable((matSize,matSize), symmetric=True) for i in range(len(points))]
    # cons = []
    # lam = cp.Variable()
    # for i, mat in enumerate(mats):
        # cons.append(mat >> 0)
        # varInd = 0
        # for j in range(matSize):
            # cons.append(mat[j,j] == 1)
            # for k in range(0, j):
                # term = 0
                # for el in points[i]:
                    # term += el*x[varInd]
                    # varInd += 1
                # cons.append(mat[k,j] == term)
    # for i in range(numVars):
        # cons.append(x[i] >= 0)
        # cons.append(x[i] <= 1)
    # obj = cp.sum(x)
    # prob = cp.Problem(cp.Maximize(obj), cons)
    # print(prob)
    # prob.solve(solver=cp.MOSEK, mosek_params={"MSK_IPAR_NUM_THREADS": 1})
    # print("lam = ", lam.value)
    # print("x = ", x.value)
    # cost(x.value, points, matSize)

    # res = minimize(cost, vec, args = (points, matSize), method = "nelder-mead", tol = 1e-8, options = {"maxiter": 10000})
    # res = minimize(cost, vec, args = (points, matSize), method = "BFGS", tol = 1e-8, options = {"maxiter": 10000})
    # res = minimize(cost, vec, args = (points, matSize), method = "L-BFGS-B", tol = 1e-2, options = {"maxiter": 10000})
    res = minimize(cost, vec, args = (points, matSize), method = "COBYLA", tol = 1e-6, options = {"maxiter": 10000})
    finalVec = res.x
    # finalVec = vec

    # pool = Pool(threads)
    # toSplit = [(matSize, numVars, points, i) for i in range(threads)]
    # results = pool.map(bruteForce, toSplit)
    # res = 10000000
    # for result in results:
        # res = min(res, result)

    # print("result = ", res)
    # print("num vars = ", numVars)
    # print("final vector = ", finalVec)

    if points[0].shape[0] == 3:
        pool2 = Pool(threads)
        toSplit2 = [(finalVec, matSize, i) for i in range(threads)]
        points2 = pool2.map(getPointsUniformCustom, toSplit2)
        points2 = reduce(lambda a, b: np.concatenate((a,b)), points2)
        points2 = points2.reshape(-1, 3)
        pointArray.append(points2)
        nameArray.append("Custom")
    elif points[0].shape[0] == 4:
        pool2 = Pool(threads)
        toSplit2 = [(finalVec, matSize, i) for i in range(threads)]
        points2 = pool2.map(getPointsUniformCustom4, toSplit2)
        points2 = reduce(lambda a, b: np.concatenate((a,b)), points2)
        points2 = points2.reshape(-1, 3)
        pointArray.append(points2)
        nameArray.append("Custom")
    else:
        exit()

# For each thing to each draw or recalulate
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
            writePoints(points2, "data/pointsSphube.csv")
        else:
            points2 = readPoints("data/pointsSphube.csv")
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
            writePoints(points3, "data/pointsCone.csv")
        else:
            points3 = readPoints("data/pointsCone.csv")
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
            writePoints(points3, "data/pointsTest.csv")
        else:
            points3 = readPoints("data/pointsTest.csv")
        if thingToDraw["draw"]:
            pointArray.append(points3)
            nameArray.append("Test Region")

    # If told to draw the data for a specific fourth value
    elif thingToDraw["name"] == "data" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            for val in np.linspace(-1, 1, 30):
                fourthVal = val
                pool = Pool(threads)
                points = pool.map(getPointsUniform, range(threads))
                points = reduce(lambda a, b: np.concatenate((a,b)), points)
                points = points.reshape(-1, 3)
                writePoints(points, "data/points" + str(fourthVal) + ".csv")
        else:
            points = readPoints("data/points" + str(fourthVal) + ".csv")
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

