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
np.random.seed(0)
random.seed(0)

# Which things to draw, regenerate, and check
onlyCheckEdge = True
thingsToDraw = {
        "box"         : {"draw": False,  "regen": False, "check": False},
        "data"        : {"draw": True,  "regen": False, "check": True},
        "dataPartial" : {"draw": False, "regen": False, "check": False},
        "sphube"      : {"draw": False, "regen": False, "check": False},
        "cone"        : {"draw": False, "regen": False, "check": False},
        "optimized"   : {"draw": False,  "regen": False,  "check": False},
        "test"        : {"draw": True,  "regen": True,  "check": True},
        "point"       : {"draw": True,  "regen": True,  "check": False},
    }
fourthVal = 0.0
# fourthVal = 0.25
# fourthVal = 0.5
# fourthVal = 0.75
# fourthVal = 1.0
pointsPer = 120
pointsPerCheck = 30
limMin = -1.05
limMax = 1.05
threads = 8
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

# Check if the four variable point is in the partial <10D cone
def checkPointPartial(var1, var2, var3, var4):
    X = cp.Variable((4,4), symmetric=True)
    cons = [
        X[0,0] == 1,
        X[1,1] == 1,
        X[2,2] == 1,
        X[3,3] == 1,
        X[0,2] == var1,
        X[0,3] == var2,
        X[1,2] == var3,
        X[1,3] == var4,
        X >> 0
    ]
    prob = cp.Problem(cp.Minimize(0), cons)
    try:
        prob.solve(solver=cp.MOSEK, mosek_params={"MSK_IPAR_NUM_THREADS": 1})
        return prob.status == "optimal"
    except:
        return False

# Check if the four variable point is in the full 10D cone
def checkPoint(var1, var2, var3, var4, minEig=False):
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
        if minEig:
            vals, vecs = np.linalg.eig(X.value)
            return min(vals)
        else:
            return prob.status == "optimal"
    except:
        return False

# Same as above, but sample uniformly over the target vars
def getPointsUniform(a):
    points = []
    var4 = fourthVal
    print("Generating points for fourth variable = ", var4)
    count = 0
    fullRegion = np.linspace(limMin, limMax, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in np.linspace(limMin, limMax, pointsPer):
            for var3 in np.linspace(limMin, limMax, pointsPer):
                count += 1
                if checkPoint(var1, var2, var3, var4):
                    points.append([var1, var2, var3])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPer**3), "%")
    points = np.array(points)
    return points

# Same as above, but sample uniformly over the target vars, but over a reduced region
def getPointsUniformPartial(a):
    points = []
    var4 = fourthVal
    print("Generating points for fourth variable = ", var4)
    count = 0
    fullRegion = np.linspace(limMin, limMax, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in np.linspace(limMin, limMax, pointsPer):
            for var3 in np.linspace(limMin, limMax, pointsPer):
                count += 1
                if checkPointPartial(var1, var2, var3, var4):
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
    count = 0
    fullRegion = np.linspace(limMin, limMax, pointsPer)
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
    return np.all(vals >= -1e-5)

# Get the point cloud representation of the 3D SDP cone
def getPointsUniformCone(a):
    points = []
    count = 0
    fullRegion = np.linspace(limMin, limMax, pointsPer)
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

def Sqrt(x):
    return math.sqrt(x)

# Check if a point is in the test domain
def checkPointTest(x, y, z):

    a = fourthVal
    b = 3*(1-a)

    x2 = x**2
    y2 = y**2
    z2 = z**2
    a2 = a**2
    x4 = x**4
    y4 = y**4
    z4 = z**4
    a4 = a**4
    x6 = x**6
    y6 = y**6
    z6 = z**6
    a6 = a**6
    m = x + y + z + a
    s1 = x2 + y2 + z2 + a2
    s2 = x2*y2 + z2*(x2 + y2) + a2*(x2 + y2 + z2)
    s3 = x2*y2*z2 + a2*(x2*y2 + x2*z2 + y2*z2)
    s4 = x2*y2*z2*a2
    p1 = x4 + y4 + z4 + a4
    # return s1 - s2 + s3 <= 4 and abs(x) <= 1 and abs(y) <= 1 and abs(z) <= 1
    # return s1 - 1.0*s2 + 1.5*s3 <= 1
    rt3o2 = math.sqrt(3.0)/2.0
    # return x**4+y**4+z**4 <= 3.0*(rt3o2**4)
    # v = 1 - 3*Sqrt(3)*x - (28*x**2)/3.0 + 3*Sqrt(3)*x**3 + (25*x**4)/3.0 - y**2 + x**2*y**2 - 3*Sqrt(3)*z - 2*x*y**2*z - z**2 + x**2*z**2 + y**2*z**2 + 3*Sqrt(3)*z**3 
    # print(x,y,z,v)
    # return (abs(v) <= 1e-1)
    # return (-1+y**2)*(-1+z**2)-2*x*y*z*a-a**2+x**2*(-1+x**2) >= 0

    # TODO 
    if abs(x) > 1 or abs(y) > 1 or abs(z) > 1:
        return False
    # coeff = (9.0/4.0 - 2.0) / ((1-rt3o2**2)**3)
    # return 2 - s1 + 2*x*y*z*a + coeff*(1-x2)*(1-y2)*(1-z2)*(1-a2) >= 0
    # return 1.0-x2-y2+0.5*x2*y2 >= 0.0
    coeff = 16
    return 1.0000000000000002 + 0.03927473781378879*y - 0.03927473780993075*z - 0.32289196022724304*x**2 - 0.17648042579037382*x*y + 0.17648042579314377*x*z - 0.3228919602289153*y**2 - 0.2628353304819463*y*z - 0.3228919602250197*z**2 + 0.052771513067287154*x**2*y - 0.052771513071611424*x**2*z - 0.17248880040359144*y**3 - 0.05277151306658405*y**2*z + 0.05277151307088638*y*z**2 + 0.17248880040228845*z**3 + 0.15638462705278344*x**4 + 0.1297898530953499*x**3*y - 0.12978985309921556*x**3*z - 0.022706693472113598*x**2*y**2 + 0.1812009895350408*x**2*y*z - 0.02270669347489869*x**2*z**2 + 0.12978985309617294*x*y**3 - 0.1216671584000173*x*y**2*z + 0.1216671584038794*x*y*z**2 - 0.12978985310059352*x*z**3 + 0.1563846270554527*y**4 + 0.193298258312934*y**3*z - 0.022706693472616973*y**2*z**2 + 0.19329825831456743*y*z**3 + 0.15638462705350464*z**4 + 0.00902684799430618*x**4*y - 0.009026847993394482*x**4*z - 0.031775558874438184*x**2*y**3 + 0.29215751482644636*x**2*y**2*z - 0.2921575148227299*x**2*y*z**2 + 0.031775558868669604*x**2*z**3 + 0.13919787114263796*y**5 - 0.00902684799886655*y**4*z - 0.03177555887302186*y**3*z**2 + 0.03177555887270084*y**2*z**3 + 0.009026847991736288*y*z**4 - 0.13919787114681342*z**5 - 0.9953957475353586*x**6 - 0.10615321402496394*x**5*y + 0.10615321402548367*x**5*z + 0.19720806994310314*x**4*y**2 - 0.17093530619639055*x**4*y*z + 0.19720806994416706*x**4*z**2 + 0.18510021136546054*x**3*y**3 - 0.18638660331346713*x**3*y**2*z + 0.18638660331376594*x**3*y*z**2 - 0.18510021136455174*x**3*z**3 + 0.19720806994207155*x**2*y**4 + 0.2775887708738374*x**2*y**3*z + 0.7800068327826088*x**2*y**2*z**2 + 0.2775887708721724*x**2*y*z**3 + 0.19720806994421503*x**2*z**4 - 0.10615321402800677*x*y**5 + 0.11477427925945088*x*y**4*z + 0.18638660331256401*x*y**3*z**2 - 0.18638660331209198*x*y**2*z**3 - 0.1147742792626624*x*y*z**4 + 0.1061532140279548*x*z**5 - 0.9953957475362483*y**6 - 0.15809580560984002*y**5*z + 0.19720806994404338*y**4*z**2 + 0.2756729252440131*y**3*z**3 + 0.19720806994416745*y**2*z**4 - 0.15809580560606187*y*z**5 - 0.9953957475388879*z**6 + 0.008171988987177772*x**6*y - 0.008171988987099449*x**6*z - 0.01056025677271797*x**4*y**3 - 0.056469429832836385*x**4*y**2*z + 0.05646942982812289*x**4*y*z**2 + 0.010560256777932625*x**4*z**3 - 0.027634534401614494*x**2*y**5 - 0.05646942982876589*x**2*y**4*z + 0.24361900884602045*x**2*y**3*z**2 - 0.24361900884770213*x**2*y**2*z**3 + 0.05646942983144403*x**2*y*z**4 + 0.02763453440550253*x**2*z**5 - 0.00598380855283516*y**7 - 0.008171988986532014*y**6*z - 0.027634534401813335*y**5*z**2 + 0.010560256776598359*y**4*z**3 - 0.010560256774937882*y**3*z**4 + 0.02763453440268286*y**2*z**5 + 0.008171988987150769*y*z**6 + 0.005983808554456169*z**7 + 0.1619030807098169*x**8 - 0.011669743700786106*x**7*y + 0.011669743701014866*x**7*z + 0.32744890874042315*x**6*y**2 - 0.009362041245048053*x**6*y*z + 0.3274489087398257*x**6*z**2 - 0.016276788156376246*x**5*y**3 + 0.012513531918895843*x**5*y**2*z - 0.01251353191966579*x**5*y*z**2 + 0.016276788157527946*x**5*z**3 - 0.026607263894940406*x**4*y**4 - 0.014835000159162115*x**4*y**3*z - 0.8613312609537082*x**4*y**2*z**2 - 0.01483500015972758*x**4*y*z**3 - 0.026607263894937335*x**4*z**4 - 0.01627678815685482*x**3*y**5 + 0.009960940715245961*x**3*y**4*z - 0.46004370685195106*x**3*y**3*z**2 + 0.46004370685104773*x**3*y**2*z**3 - 0.009960940714893049*x**3*y*z**4 + 0.016276788158926043*x**3*z**5 + 0.32744890874145743*x**2*y**6 - 0.01863661809657871*x**2*y**5*z - 0.8613312609547316*x**2*y**4*z**2 - 0.6851509972439722*x**2*y**3*z**3 - 0.8613312609563593*x**2*y**2*z**4 - 0.01863661809798224*x**2*y*z**5 + 0.32744890874163024*x**2*z**6 - 0.011669743699621024*x*y**7 + 0.006286129883556896*x*y**6*z - 0.01251353191675433*x*y**5*z**2 + 0.009960940716903677*x*y**4*z**3 - 0.009960940715541142*x*y**3*z**4 + 0.012513531915474385*x*y**2*z**5 - 0.006286129883742414*x*y*z**6 + 0.011669743700309525*x*z**7 + 0.16190308070971005*y**8 - 0.017379949806884454*y**7*z + 0.32744890873863075*y**6*z**2 - 0.024241300278130734*y**5*z**3 - 0.026607263895209066*y**4*z**4 - 0.024241300279991246*y**3*z**5 + 0.32744890874098587*y**2*z**6 - 0.0173799498086612*y*z**7 + 0.16190308071040235*z**8 >= 0
    # return 2 - s1 + 2*x*y*z*a + coeff*(1-a2)*(1-x2)*(1-y2)*(1-z2) >= 0
    # coeff0 = 2.77186
    # coeff1 = 2.85679
    # coeff2 = -2.27138
    # coeff3 = -1.01515
    # coeff4 = -1.00968
    # return 2 - s1 + 2*x*y*z*a + (1-a2)*(
                        # 1
                        # +coeff0*(x4+y4+z4)
                        # +coeff1*(x2*y2*z2)
                        # +coeff2*(x2+y2+z2)
                        # +coeff3*(x4*y4*z4)
                        # +coeff4*(x6*y6*z6)
                        # -(x6+y6+z6)
                    # ) >= 0 
    # coeff = 6.15477
    # return 2 - s1 + 2*x*y*z*a + coeff*(1-a2)*(1-0.5*x2-0.5*x4)*(1-0.5*y2-0.5*y4)*(1-0.5*z2-0.5*z4) >= 0

    # X0 = [[2, s1], 
          # [0, 1+s2]]
    X0 = [[2, s1, s1], 
          [0, 1+s2, s1], 
          [0, 0, 1+s2]]
    # X0 = [[1, x, y, z], 
          # [x, 2, z, y], 
          # [y, z, 2, x],
          # [z, y, x, 1]]
    # X0 = [[1, x, y, z, 0], 
          # [0, 1, 0, 0, z], 
          # [0, 0, 1, 0, y],
          # [0, 0, 0, 1, x],
          # [0, 0, 0, 0, 1]]
    X1 = [[1, x, y], 
          [0, 1, z], 
          [0, 0, 1]]
    X0 = np.array(X0)
    X1 = np.array(X1)
    for i in range(X0.shape[0]):
        for j in range(0, i):
            X0[i,j] = X0[j,i]
    for i in range(X1.shape[0]):
        for j in range(0, i):
            X1[i,j] = X1[j,i]
    # b = 1 - a
    I = np.eye(3)
    X = np.zeros((6,6))
    X[:3,:3] = X0 + a*I
    X[3:,3:] = X1 + b*I
    # X = X0
    # b = 0.9
    # X[:3,:3] = I
    # X[3:,3:] = X1
    # X[3,3] -= b
    # X[4,4] += b
    # X[5,5] += b
    # X[3:,3:] = I
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
    return np.all(vals >= -1e-5)

# Get the point cloud representation of the test region
def getPointsUniformTest(a):
    points = []
    count = 0
    fullRegion = np.linspace(limMin, limMax, pointsPer)
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
            if x[0] >= limMin and x[0] <= limMax and x[1] >= limMin and x[1] <= limMax and x[2] >= limMin and x[2] <= limMax:
                points.append(x)
        print("Loaded ", len(points), " points")
        return np.array(points)

pointArray = []
nameArray = []

# Check if a point is in the a custom cone 
def checkPointCustom(vec, point, numMats, matSize):
    # X = np.zeros((matSize, matSize))
    # vecInd = 0
    # allGood = True
    # for k in range(numMats):
        # for i in range(matSize):
            # for j in range(0, i+1):
                # X[j,i] = vec[vecInd]
                # vecInd += 1
                # for el in point:
                    # X[j,i] += vec[vecInd]*el
                    # vecInd += 1
                # for el in point:
                    # X[j,i] += vec[vecInd]*el**2
                    # vecInd += 1
                # X[i,j] = X[j,i]
        # vals, vecs = np.linalg.eig(X)
        # if not np.all(vals >= -1e-5):
            # allGood = False
            # break
    # return allGood
    x = point[0]
    y = point[1]
    z = point[2]
    if len(point) >= 4:
        a = point[3]
    else:
        a = fourthVal
    x2 = x*x
    y2 = y*y
    z2 = z*z
    a2 = a*a
    x4 = x2*x2
    y4 = y2*y2
    z4 = z2*z2
    a4 = a2*a2
    s0 = x + y + z + a
    s1 = x2 + y2 + z2 + a2
    s2 = x2*y2 + z2*(x2 + y2) + a2*(x2 + y2 + z2)
    s3 = x2*y2*z2 + a2*(x2*y2 + z2*(x2 + y2))
    s4 = x2*y2*z2*a2
    t1 = x4 + y4 + z4 + a4
    t2 = x4*y4 + z4*(x4 + y4) + a4*(x4 + y4 + z4)
    t3 = x4*y4*z4 + a4*(x4*y4 + z4*(x4 + y4))
    t4 = x4*y4*z4*a4
    f = 2 - x2 - y2 - z2 - a2 + 2*x*y*z*a + (1-a2)*(vec[0] + vec[1]*s1 + vec[2]*s2 + vec[3]*s3 + vec[4]*s4 + vec[5]*t1 + vec[6]*t2 + vec[7]*t3 + vec[8]*t4)
    dist = abs(f)
    if abs(x) > 1 or abs(y) > 1 or abs(z) > 1:
        return -1
    return dist
    # X = [[vec[0], vec[1]*s0, vec[2]*s1, vec[3]*s2],
         # [0, 1, 0, 0],
         # [0, 0, 1, 0],
         # [0, 0, 0, 1]]
    # X = np.array(X)
    # for i in range(X.shape[0]):
        # for j in range(0, i):
            # if X[i,j] == 0:
                # X[i,j] = X[j,i]
    # vals, vecs = np.linalg.eig(X)
    # dist = min(vals)
    # return dist

# Cost function for fitting the cone
def cost(vec, points, numMats, matSize):
    dist = 0
    for point in points:
        dist += checkPointCustom(vec, point, numMats, matSize)
    cost = dist**2
    print("dist = ", dist, ", cost = ", cost)
    return cost

# Get the point cloud representation of the test region
def getPointsUniformCustom(args):
    coeffs, numMats, matSize, a = args
    points = []
    count = 0
    fullRegion = np.linspace(limMin, limMax, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                if abs(checkPointCustom(coeffs, [var1, var2, var3], numMats, matSize)) <= 1e-3:
                    points.append([var1, var2, var3])
        if a == 0:
            print(100.0*threads*float(count)/(pointsPer**3), "%")
    points = np.array(points)
    return points

# Get the point cloud representation of the test region
def getPointsUniformCustom4(args):
    coeffs, numMats, matSize, a = args
    points = []
    count = 0
    fullRegion = np.linspace(limMin, limMax, pointsPer)
    localRegion = fullRegion[a*pointsPer//threads:(a+1)*pointsPer//threads]
    for var1 in localRegion:
        for var2 in fullRegion:
            for var3 in fullRegion:
                count += 1
                if checkPointCustom(coeffs, [var1, var2, var3, fourthVal], numMats, matSize) >= 1e-3:
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
    numVars = 9
    # numVars = numMats * ((matSize*(matSize+1)//2) * (2*allPoints.shape[1]+1))
    pointsToSample = min(pointsToSample, len(allPoints))
    for i in range(pointsToSample):
        rand = random.randint(0, len(allPoints)-1)
        points.append(allPoints[rand])
        np.delete(allPoints, rand)
    points.append([math.sqrt(3)/2, math.sqrt(3)/2, math.sqrt(3)/2])
    points.append([-math.sqrt(3)/2, math.sqrt(3)/2, math.sqrt(3)/2])
    points.append([math.sqrt(3)/2, -math.sqrt(3)/2, math.sqrt(3)/2])
    points.append([math.sqrt(3)/2, math.sqrt(3)/2, -math.sqrt(3)/2])
    points = np.array(points)

    # Optimization
    converge = 1e-6
    maxIters = 10000
    vec = np.random.rand(numVars)
    res = minimize(cost, vec, args = (points, numMats, matSize), method = "L-BFGS-B", tol = converge, options = {"maxiter": maxIters})
    # res = minimize(cost, vec, args = (points, numMats, matSize), method = "COBYLA", tol = converge, options = {"maxiter": maxIters})
    print("result = ", res)
    finalVec = np.array(res.x)

    print("eval at [1,1,0] = ", checkPointCustom(finalVec, [1,1,0], numMats, matSize))
    print("eval at [sqrt(3)/2,sqrt(3)/2,sqrt(3)/2] = ", checkPointCustom(finalVec, [math.sqrt(3)/2,math.sqrt(3)/2,math.sqrt(3)/2], numMats, matSize))

    print("final vec = ", finalVec)
    cost(finalVec, points, numMats, matSize)

    if points[0].shape[0] == 3:
        with Pool(threads) as pool:
            toSplit = [(finalVec, numMats, matSize, i) for i in range(threads)]
            points = pool.map(getPointsUniformCustom, toSplit)
            points = reduce(concat, points)
            points = points.reshape(-1, 3)
            hull = ConvexHull(points)
            points = points[hull.vertices,:]
    elif points[0].shape[0] == 4:
        with Pool(threads) as pool:
            toSplit = [(finalVec, numMats, matSize, i) for i in range(threads)]
            points = pool.map(getPointsUniformCustom4, toSplit)
            points = reduce(concat, points)
            points = points.reshape(-1, 3)
            hull = ConvexHull(points)
            points = points[hull.vertices,:]

    # Write to file
    writePoints(points, "data/pointsOptimized.csv")

    # If drawing after regenning
    if thingsToDraw["optimized"]["draw"]:
        pointArray.append(points)
        nameArray.append("optimized")

# If drawing without regenning
elif thingsToDraw["optimized"]["draw"]:
    points = readPoints("data/pointsOptimized.csv")
    pointArray.append(points)
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
            with Pool(threads) as pool:
                points = pool.map(getPointsUniformSphube, range(threads))
                points = reduce(concat, points)
                points = points.reshape(-1, 3)
                hull = ConvexHull(points)
                points = points[hull.vertices,:]
                writePoints(points, "data/pointsSphube.csv")
        else:
            points = readPoints("data/pointsSphube.csv")
        if thingToDraw["draw"]:
            pointArray.append(points)
            nameArray.append(name)

    # If told to draw the 3D SDP cone
    elif name == "cone" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            with Pool(threads) as pool:
                points = pool.map(getPointsUniformCone, range(threads))
                points = reduce(concat, points)
                points = points.reshape(-1, 3)
                hull = ConvexHull(points)
                points = points[hull.vertices,:]
                writePoints(points, "data/pointsCone.csv")
        else:
            points = readPoints("data/pointsCone.csv")
        if thingToDraw["draw"]:
            pointArray.append(points)
            nameArray.append(name)

    # If told to draw the test region
    elif name == "test" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            with Pool(threads) as pool:
                points = pool.map(getPointsUniformTest, range(threads))
                points = reduce(concat, points)
                points = points.reshape(-1, 3)
                hull = ConvexHull(points)
                points = points[hull.vertices,:]
                writePoints(points, "data/pointsTest.csv")
        else:
            points = readPoints("data/pointsTest.csv")
        if thingToDraw["draw"]:
            pointArray.append(points)
            nameArray.append(name)

    # If told to draw a test point
    elif name == "point" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["draw"]:
            pointArray.append(np.array([[0.866,0.866,0.866]]))
            nameArray.append(name)
            pointArray.append(np.array([[math.sqrt(2.0+math.sqrt(2.0)) / 2.0,math.sqrt(2.0+math.sqrt(2.0)) / 2.0,1.0 / math.sqrt(2.0)]]))
            nameArray.append(name)
            pointArray.append(np.array([[
                         math.sqrt((1.0+math.sqrt(391.0)/20.0)/2.0),
                         math.sqrt((1.0+math.sqrt(391.0)/20.0)/2.0),
                         3.0/20.0
                     ]]))
            nameArray.append(name)
            pointArray.append(np.array([[
                         1.0/(2*math.sqrt(2-math.sqrt(3))),
                         1.0/(2*math.sqrt(2-math.sqrt(3))),
                         1.0/2.0
                     ]]))
            nameArray.append(name)

    # If told to draw the data for a specific fourth value
    elif name == "data" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            # for val in np.linspace(limMin, limMax, 30):
            for val in [fourthVal]:
                fourthVal = val
                with Pool(threads) as pool:
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

    # If told to draw the data for a specific fourth value
    elif name == "dataPartial" and (thingToDraw["draw"] or thingToDraw["regen"]):
        if thingToDraw["regen"]:
            with Pool(threads) as pool:
                points = pool.map(getPointsUniformPartial, range(threads))
                points = reduce(concat, points)
                points = points.reshape(-1, 3)
                hull = ConvexHull(points)
                points = points[hull.vertices,:]
                writePoints(points, "data/pointsPartial" + str(fourthVal) + ".csv")
        else:
            points = readPoints("data/pointsPartial" + str(fourthVal) + ".csv")
            print("loaded ", len(points), " points")
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
    count = 0
    fullRegion = np.linspace(limMin, limMax, pointsPerCheck)
    localRegion = fullRegion[a*pointsPerCheck//threads:(a+1)*pointsPerCheck//threads]
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
                if thingsToDraw["dataPartial"]["check"]:
                    results["data"] = checkPointPartial(var1, var2, var3, fourthVal)
                for name, thing in thingsToDraw.items():
                    for name2, thing2 in thingsToDraw.items():
                        if thing["check"] and thing2["check"]:
                            counts[(name, name2)] += (results[name] == results[name2])
                            if results[name] != results[name2]:
                                print(var1, var2, var3)
                                if thingsToDraw["data"]["check"]:
                                    minEig = checkPoint(var1, var2, var3, fourthVal, True)
                                    print("minEig = ", minEig)
                                print(name, results[name])
                                print(name2, results[name2])
            if a == 0:
                print(100.0*threads*float(count)/(pointsPerCheck**3), "%")
    for name, thing in thingsToDraw.items():
        for name2, thing2 in thingsToDraw.items():
            if thing["check"] and thing2["check"]:
                counts[(name, name2)] *= 100.0 / float(pointsPerCheck**3)
    return counts

# Same as above, but just the points exactly at the edge of the convex hull of the data
# TODO
def checkAllEdge(a):
    counts = {}
    for name, thing in thingsToDraw.items():
        for name2, thing2 in thingsToDraw.items():
            if thing["check"] and thing2["check"]:
                counts[(name, name2)] = 0
    count = 0
    pointsToCheck = []
    for i in range(len(pointArray)):
        if nameArray[i] == "data" and thingsToDraw["data"]["check"]:
            pointsToCheck = pointArray[i]
            break
        if nameArray[i] == "dataPartial" and thingsToDraw["dataPartial"]["check"]:
            pointsToCheck = pointArray[i]
            break
    localRegion = pointsToCheck[a*pointsPerCheck//threads:(a+1)*pointsPerCheck//threads]
    for (var1, var2, var3) in localRegion:
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
            results["data"] = True
        if thingsToDraw["dataPartial"]["check"]:
            results["dataPartial"] = True
        for name, thing in thingsToDraw.items():
            for name2, thing2 in thingsToDraw.items():
                if thing["check"] and thing2["check"]:
                    counts[(name, name2)] += (results[name] == results[name2])
                    if results[name] != results[name2]:
                        minEig = checkPoint(var1, var2, var3, fourthVal, True)
                        print(var1, var2, var3)
                        print("minEig = ", minEig)
                        print(name, results[name])
                        print(name2, results[name2])
        if a == 0:
            print(100.0*threads*float(count)/(pointsPerCheck**3), "%")
    for name, thing in thingsToDraw.items():
        for name2, thing2 in thingsToDraw.items():
            if thing["check"] and thing2["check"]:
                counts[(name, name2)] *= 100.0 / float(pointsPerCheck**3)
    return counts

# If we need to do any analysis
checkAny = False
for thing in thingsToDraw.values():
    if thing["check"]:
        checkAny = True
        break
if checkAny:
    with Pool(threads) as pool:
        if onlyCheckEdge:
            results = pool.map(checkAllEdge, range(threads))
        else:
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
    ax.set_xlim(limMin, limMax)
    ax.set_ylim(limMin, limMax)
    ax.set_zlim(limMin, limMax)
    ax.margins(x=0, y=0, z=0)
    ax.set_aspect('auto')
    ax.set_proj_type("ortho")
    styles = ["r-", "b-", "g-", "y-", "c-", "m-", "k-"]
    for i, points in enumerate(pointArray):
        style = styles[i]
        if points.shape[0] > 1:
            hull = ConvexHull(points[:,:])
            for s in hull.simplices:
                s = np.append(s, s[0])
                ax.plot(points[s, 0], points[s, 1], points[s, 2], style)
        else:
            ax.scatter(points[0,0], points[0,1], points[0,2], style)
        ax.plot([], [], style, label=nameArray[i])
    fig.tight_layout()
    fig.legend()
    plt.show()

