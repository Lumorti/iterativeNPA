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
thingsToDraw = {
        "box"      : {"draw": False,  "regen": False, "check": False},
        "data"     : {"draw": True,  "regen": False, "check": False},
        "sphube"   : {"draw": False, "regen": False, "check": False},
        "cone"     : {"draw": False, "regen": False, "check": False},
        "test"     : {"draw": True,  "regen": True,  "check": False},
        "optimized": {"draw": False,  "regen": False,  "check": False},
    }
# fourthVal = 0.0
# fourthVal = 0.24137931034482762
# fourthVal = 0.5172413793103448
fourthVal = 0.7931034482758621
# fourthVal = 1.0
threads = 10
thresh = 0.30
tol = 0.01

# Join two numpy arrays together if they aren't empty
def concat(a, b):
    if a.shape[0] == 0:
        return b
    elif b.shape[0] == 0:
        return a
    else:
        return np.concatenate((a,b))

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
    return points

# Check if a point is in the test domain
# TODO found a pretty good approx
def checkPointTest(x, y, z):
    a = fourthVal
    x2 = x**2
    y2 = y**2
    z2 = z**2
    a2 = a**2
    m = x + y + z + a
    s = x2*y2 + z2*(x2 + y2) + a2*(x2 + y2 + z2)
    t = x2 + y2 + z2 + a2
    # s = x2*y2 + z2*(x2 + y2)
    # t = x2 + y2 + z2
    X0 = [[2, t, t], 
         [0, 1+s, t], 
         [0, 0, 1+s]]
    X1 = [[1, x, y], 
         [0, 1, z], 
         [0, 0, 1]]
    X0 = np.array(X0)
    X1 = np.array(X1)
    for i in range(X0.shape[0]):
        for j in range(0, i):
            X0[i,j] = X0[j,i]
            X1[i,j] = X1[j,i]
    b = 1 - a
    I = np.eye(3)
    X = np.zeros((6,6))
    X[:3,:3] = X0 + a*I
    X[3:,3:] = X1 + b*I
    # X = (a*X1 + b*I)@(b*X0 + a*I)
    # X = a*(X0@X1) + b*(X1@X0)
    # print(X1)
    # print(np.power(X1,a))
    # print(X)
    # X = X0
    # b = x*y*z*a
    # c = x*x + y*y + z*z + a*a
    # X = [[1, x, 0, x], 
         # [0, 1, y, 0], 
         # [0, 0, 1, z],
         # [0, 0, 0, 1]]
    # X = [[2, t], 
         # [0, 1+s]]
    # X = [[2, t, t, t], 
         # [0, 1+s, t, t], 
         # [0, 0, 1+s, t],
         # [0, 0, 0, 1+s]]
    # X = [[1, x, y, 0, 0], 
         # [0, 2, z, 0, 0], 
         # [0, 0, 1, 0, 0],
         # [0, 0, 0, 1, 0],
         # [0, 0, 0, 0, 1]]
    # e = 3 - s
    # X = [[e, x, 0, 0], 
         # [0, e, y, 0], 
         # [0, 0, e, z],
         # [0, 0, 0, e]]
    # X = np.array(X)
    # for i in range(X.shape[0]):
        # for j in range(0, i):
            # if X[i,j] == 0:
                # X[i,j] = X[j,i]
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
def checkPointCustom(vec, point, numMats, matSize):
    X = np.zeros((matSize, matSize))
    vecInd = 0
    allGood = True
    for k in range(numMats):
        for i in range(matSize):
            for j in range(0, i+1):
                X[j,i] = vec[vecInd]
                vecInd += 1
                for el in point:
                    X[j,i] += vec[vecInd]*el
                    vecInd += 1
                for el in point:
                    X[j,i] += vec[vecInd]*el**2
                    vecInd += 1
                X[i,j] = X[j,i]
        vals, vecs = np.linalg.eig(X)
        if not np.all(vals > 0):
            allGood = False
            break
    return allGood

# Cost function for fitting the cone
def cost(vec, points, pointsNo, numMats, matSize):
    dist = 0
    numYes = 0
    numNo = 0
    for point in points:
        vecInd = 0
        overallMin = 1000000
        for k in range(numMats):
            X = np.zeros((matSize, matSize))
            for i in range(matSize):
                for j in range(0, i+1):
                    X[j,i] = vec[vecInd]
                    vecInd += 1
                    for el in point:
                        X[j,i] += vec[vecInd]*el
                        vecInd += 1
                    for el in point:
                        X[j,i] += vec[vecInd]*el**2
                        vecInd += 1
                    X[i,j] = X[j,i]
            vals, vecs = np.linalg.eig(X)
            minVal = np.min(vals)
            overallMin = min(overallMin, minVal)
        if overallMin < 0:
            dist += 0.1*abs(overallMin)
        else:
            dist += abs(overallMin)
            # dist = max(dist, abs(np.min(vals)))
            # if np.all(vals >= 0):
                # numYes += 1
    # for point in pointsNo:
        # X = np.zeros((matSize, matSize))
        # vecInd = 0
        # for i in range(matSize):
            # for j in range(0, i+1):
                # X[j,i] = vec[vecInd]
                # vecInd += 1
                # for el in point:
                    # X[j,i] += vec[vecInd]*el
                    # vecInd += 1
                # X[i,j] = X[j,i]
        # vals, vecs = np.linalg.eig(X)
        # if np.any(vals < 0):
            # numNo += 1
        # dist -= np.min(vals)
    dist /= len(points)
    cost = dist**2 - numNo - numYes
    print("no = ", numNo, ", yes = ", numYes, ", dist = ", dist, ", cost = ", cost)
    # cost = dist**2
    # print("avg dist =", dist, ", cost =", cost)
    return cost

# Get the point cloud representation of the test region
def getPointsUniformCustom(args):
    coeffs, numMats, matSize, a = args
    points = []
    pointsPer = 60
    count = 0
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                if checkPointCustom(coeffs, [var1, var2, var3], numMats, matSize):
                    points.append([var1, var2, var3])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    points = np.array(points)
    return points

# Get the point cloud representation of the test region
def getPointsUniformCustom4(args):
    coeffs, numMats, matSize, a = args
    points = []
    pointsPer = 60
    count = 0
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                if checkPointCustom(coeffs, [var1, var2, var3, fourthVal], numMats, matSize):
                    points.append([var1, var2, var3])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    points = np.array(points)
    return points

# If told to optimize to minimize the SDP hull
if thingsToDraw["optimized"]["regen"]:
    matSize = 3
    numMats = 1
    # allPoints = readPoints("data/pointsAll.csv")
    allPoints = readPoints("data/points" + str(fourthVal) + ".csv")
    points = []
    pointsToSample = 1000000000
    numVars = numMats * ((matSize*(matSize+1)//2) * (2*allPoints.shape[1]+1))
    # vec = np.array([1 for i in range(numVars)])
    # vec = np.zeros(numVars)
    vec = np.random.rand(numVars)
    # vec = np.array([-0.21806701, -0.30024651, -0.2024949 , -0.15697767, -0.12269202,
       # -0.482077  ,  0.26369564,  0.28087138, -0.21673742, -0.02869441,
       # -0.22670792,  0.08536931, -0.09920705, -0.0567977 , -0.23431801,
       # -0.07499742, -0.50900859,  0.13968817,  0.15751262,  0.40836086,
        # 0.50056849,  0.09140297, -0.13158019,  0.47424429,  0.12773484,
       # -0.22482994, -0.38876055,  0.31887217, -0.06631954, -0.01572061,
       # -0.27017813, -0.09920275,  0.24522633, -0.15417556,  0.32584317,
       # -0.09420153, -0.02195686,  0.27172932, -0.00135903,  0.05110155])
    # 4x4 with diags, 1+x+y+z+a, 0.08
    # vec = np.array([2.50204461, 0.46072533, 0.54954365, 0.60548681, 0.89528628,
       # 1.65300713, 0.63412646, 0.98233975, 0.89811249, 0.91462576,
       # 2.30651061, 1.0603816 , 0.83007196, 0.7812964 , 0.63145827,
       # 1.21029969, 0.62853171, 0.75748408, 0.66273371, 0.6183474 ,
       # 1.3480396 , 0.71374652, 0.71004322, 0.75188841, 0.80662968,
       # 1.49870684, 0.66534408, 0.66269906, 0.70416922, 0.72317737,
       # 1.30149258, 0.55528545, 0.83009985, 0.69863133, 0.72194382,
       # 1.60931153, 0.79127652, 0.75554167, 0.71476768, 0.78762482,
       # 1.43447154, 0.66580072, 0.67826922, 0.69606782, 0.76085554,
       # 1.4696805 , 0.68817808, 0.70538556, 0.67340316, 0.78624887])
    # vec = np.array([3.69819251, 0.31402345, 0.49310877, 0.44285137, 0.7493408 ,
       # 1.66339902, 0.51581952, 0.96991807, 0.87361816, 0.92034605,
       # 2.3898406 , 1.03594344, 0.86404434, 0.73793836, 0.69029193,
       # 1.47243881, 0.6213243 , 0.72831116, 0.66841164, 0.6077565 ,
       # 1.38059935, 0.68324497, 0.66682239, 0.73401832, 0.76269395,
       # 1.53973629, 0.64654947, 0.65076476, 0.68499182, 0.69765895,
       # 1.31958472, 0.55029784, 0.79959957, 0.67877301, 0.71408517,
       # 1.6276709 , 0.76874916, 0.73333071, 0.69923296, 0.75784035,
       # 1.47260647, 0.64364378, 0.64613486, 0.66667042, 0.74234814,
       # 1.50784845, 0.67151023, 0.66831342, 0.64819513, 0.76928514])
    pointsToSample = min(pointsToSample, len(allPoints))
    for i in range(pointsToSample):
        rand = random.randint(0, len(allPoints)-1)
        points.append(allPoints[rand])
        np.delete(allPoints, rand)
    points = np.array(points)

    center = np.zeros(points.shape[1])
    for point in points:
        center += point
    center /= len(points)
    pointsNo = []
    extraDistance = 0.1
    for point in points:
        pointsNo.append(point + extraDistance*(point-center))
    pointsNo = np.array(pointsNo)

    # print("Generating SDP...")
    # x = cp.Variable(numVars)
    # mats = [cp.Variable((matSize,matSize), symmetric=True) for i in range(len(points))]
    # cons = []
    # for i, mat in enumerate(mats):
        # cons.append(mat >> 0)
        # varInd = 0
        # for j in range(matSize):
            # for k in range(0, j+1):
                # term = x[varInd]
                # varInd += 1
                # for el in points[i]:
                    # term += el*x[varInd]
                    # varInd += 1
                # cons.append(mat[k,j] == term)
    # for i in range(numVars):
        # cons.append(x[i] >= -1)
        # cons.append(x[i] <= 1)
    # obj = 0
    # for i in range(numVars):
        # obj += x[i]
    # prob = cp.Problem(cp.Maximize(obj), cons)
    # print("Solving SDP...")
    # prob.solve(solver=cp.MOSEK, verbose=True, mosek_params={"MSK_IPAR_NUM_THREADS": 1})
    # finalVec = np.array(x.value)

    converge = 1e-6
    maxIters = 1000
    # res = minimize(cost, vec, args = (points, matSize), method = "nelder-mead", tol = 1e-8, options = {"maxiter": 10000})
    # res = minimize(cost, vec, args = (points, matSize), method = "BFGS", tol = 1e-8, options = {"maxiter": 10000})
    # res = minimize(cost, vec, args = (points, matSize), method = "L-BFGS-B", tol = 1e-2, options = {"maxiter": 10000})
    res = minimize(cost, vec, args = (points, pointsNo, numMats, matSize), method = "COBYLA", tol = converge, options = {"maxiter": maxIters})
    print("result = ", res)
    finalVec = np.array(res.x)

    print("final vec = ", finalVec)
    cost(finalVec, points, pointsNo, numMats, matSize)

    if points[0].shape[0] == 3:
        pool2 = Pool(threads)
        toSplit2 = [(finalVec, numMats, matSize, i) for i in range(threads)]
        points2 = pool2.map(getPointsUniformCustom, toSplit2)
        points2 = reduce(concat, points2)
        points2 = points2.reshape(-1, 3)
        hull = ConvexHull(points2)
        points2 = points2[hull.vertices,:]
    elif points[0].shape[0] == 4:
        pool2 = Pool(threads)
        toSplit2 = [(finalVec, numMats, matSize, i) for i in range(threads)]
        points2 = pool2.map(getPointsUniformCustom4, toSplit2)
        points2 = reduce(concat, points2)
        points2 = points2.reshape(-1, 3)
        hull = ConvexHull(points2)
        points2 = points2[hull.vertices,:]

    # Write to file
    writePoints(points2, "data/pointsOptimized.csv")

    # If drawing after regenning
    if thingsToDraw["optimized"]["draw"]:
        pointArray.append(points2)
        nameArray.append("optimized")

# If drawing without regenning
elif thingsToDraw["optimized"]["draw"]:
    points2 = readPoints("data/pointsOptimized.csv")
    pointArray.append(points2)
    nameArray.append("optimized")

# For each thing to each draw or recalulate
for name, thingToDraw in thingsToDraw.items():

    # If told to draw the box from -1 to 1
    if name == "box" and thingToDraw["draw"]:
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
        nameArray.append(name)

    # If told to draw a standard-ish sphube
    elif name == "sphube" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            pool2 = Pool(threads)
            points2 = pool2.map(getPointsUniformSphube, range(threads))
            points2 = reduce(concat, points2)
            points2 = points2.reshape(-1, 3)
            hull = ConvexHull(points2)
            points2 = points2[hull.vertices,:]
            writePoints(points2, "data/pointsSphube.csv")
        else:
            points2 = readPoints("data/pointsSphube.csv")
        if thingToDraw["draw"]:
            pointArray.append(points2)
            nameArray.append(name)

    # If told to draw the 3D SDP cone
    elif name == "cone" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            pool3 = Pool(threads)
            points3 = pool3.map(getPointsUniformCone, range(threads))
            points3 = reduce(concat, points3)
            points3 = points3.reshape(-1, 3)
            hull = ConvexHull(points3)
            points3 = points3[hull.vertices,:]
            writePoints(points3, "data/pointsCone.csv")
        else:
            points3 = readPoints("data/pointsCone.csv")
        if thingToDraw["draw"]:
            pointArray.append(points3)
            nameArray.append(name)

    # If told to draw the test region
    elif name == "test" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            pool3 = Pool(threads)
            points3 = pool3.map(getPointsUniformTest, range(threads))
            points3 = reduce(concat, points3)
            points3 = points3.reshape(-1, 3)
            hull = ConvexHull(points3)
            points3 = points3[hull.vertices,:]
            writePoints(points3, "data/pointsTest.csv")
        else:
            points3 = readPoints("data/pointsTest.csv")
        if thingToDraw["draw"]:
            pointArray.append(points3)
            nameArray.append(name)

    # If told to draw the data for a specific fourth value
    elif name == "data" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            for val in np.linspace(-1, 1, 30):
                fourthVal = val
                pool = Pool(threads)
                points = pool.map(getPointsUniform, range(threads))
                points = reduce(concat, points)
                points = points.reshape(-1, 3)
                hull = ConvexHull(points)
                points = points[hull.vertices,:]
                writePoints(points, "data/points" + str(fourthVal) + ".csv")
        else:
            points = readPoints("data/points" + str(fourthVal) + ".csv")
        if thingToDraw["draw"]:
            pointArray.append(points)
            nameArray.append(name)

# Sample many points and see which regions they are in
def checkAll(a):
    counts = {}
    for name, thing in thingsToDraw.items():
        for name2, thing2 in thingsToDraw.items():
            if thing["check"] and thing2["check"]:
                counts[(name, name2)] = 0
    pointsPer = 20
    count = 0
    fullRegion = np.linspace(-1, 1, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                results = {}
                if thingsToDraw["box"]["check"]:
                    results["box"] = checkPointBox(var1, var2, var3)
                if thingsToDraw["test"]["check"]:
                    results["test"] = checkPointTest(var1, var2, var3)
                if thingsToDraw["cone"]["check"]:
                    results["cone"] = checkPointCone(var1, var2, var3)
                if thingsToDraw["sphube"]["check"]:
                    results["sphube"] = checkPointSphube(var1, var2, var3)
                if thingsToDraw["data"]["check"]:
                    results["data"] = checkPoint(var1, var2, var3, fourthVal)
                for name, thing in thingsToDraw.items():
                    for name2, thing2 in thingsToDraw.items():
                        if thing["check"] and thing2["check"]:
                            counts[(name, name2)] += (results[name] == results[name2])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    for name, thing in thingsToDraw.items():
        for name2, thing2 in thingsToDraw.items():
            if thing["check"] and thing2["check"]:
                counts[(name, name2)] *= 100.0 / float(pointsPer**3)
    return counts

# If we need to do any analysis
checkAny = False
for thing in thingsToDraw.values():
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

