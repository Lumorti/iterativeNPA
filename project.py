import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cp
from matplotlib.widgets import Slider
from multiprocessing import Pool
import os
import random

# Fixed seed
np.random.seed(1)

# 3 or 10
draw = False
mode = 10
fourthVal = 1.0
threads = 1

# Sample random points, check if they are on the cone
def getPoints(a):
    points = []
    if mode == 3:
        samples = 1000
    elif mode == 10:
        samples = 10000
    for j in range(samples):
        if mode == 3:
            X = np.random.randn(3,3)
        elif mode == 10:
            X = np.random.randn(5,5)
        for i in range(X.shape[0]):
            X[i,:] /= np.linalg.norm(X[i,:])
        X = X @ X.T
        if mode == 10:
            X[2,3] = 1.0
            X[3,2] = 1.0
        val, vec = np.linalg.eig(X)
        if np.all(val >= 0):
            if mode == 3:
                x = [X[0,1], X[0,2], X[1,2]]
            elif mode == 10:
                x = [X[0,1], X[0,2], X[0,3], X[0,4], X[1,2], X[1,3], X[1,4], X[2,3], X[2,4], X[3,4]]
            points.append(x)
        if a == 0 and j % (samples//10) == 0:
            print(100.0*float(j)/samples, '%')
    return points

# Test whether my simplification is a subset
def testPoints(a):
    np.random.seed(a)
    numInequal = 0 
    totalDone = 0
    samples = 10000
    for j in range(samples):
        X = np.random.randn(5,5)
        for i in range(X.shape[0]):
            X[i,:] /= np.linalg.norm(X[i,:])
        X = X @ X.T
        val, vec = np.linalg.eig(X)
        XIsPSD = np.all(val >= 0)
        A1B1 = X[1,3]
        A1B2 = X[1,4]
        A2B1 = X[2,3]
        A2B2 = X[2,4]
        A1B1sq = A1B1**2
        A1B2sq = A1B2**2
        A2B1sq = A2B1**2
        A2B2sq = A2B2**2
        Ys = []
        # Ys.append([[1,      A1B1,       A1B2,       A2B1,       A2B2     ],
                   # [0,      1,          A2B1*A2B2,  A1B2*A2B2,  A1B2*A2B1],
                   # [0,      0,          1,          A1B1*A2B2,  A1B1*A2B1],
                   # [0,      0,          0,          1,          A1B1*A1B2],
                   # [0,      0,          0,          0,          1        ]])
        # f1 = A1B2**2 * A2B1**2 * A2B2**2
        # f2 = A1B1**2 * A1B2**2 * A2B2**2
        # f3 = A1B1**2 * A1B2**2 * A2B1**2
        # f4 = A1B1**2 * A2B1**2 * A2B2**2
        # f1 = A2B2**2
        # f2 = A2B1**2
        # f3 = A1B2**2
        # f4 = A1B1**2
        # Ys.append([[1,      A1B1,    A1B2,       A2B1,       A2B2 ],
                   # [0,      1,       A2B1*A2B2,  A1B2*A2B2,  A1B2*A2B1],
                   # [0,      0,       1,          A1B1*A2B2,  A1B1*A2B1],
                   # [0,      0,       0,          1,          A1B1*A1B2],
                   # [0,      0,       0,          0,          1        ]])
        # Ys.append(np.array([[1,      A1B1,    A1B2,       A2B1,       A2B2],
                           # [0,      1,       0,          0,          0],
                           # [0,      0,       1,          0,          0],
                           # [0,      0,       0,          1,          0],
                           # [0,      0,       0,          0,          1]]))
        # Ys.append(np.array([[1,      A1B1,   A1B2],
                           # [A1B1,    1,      0],
                           # [A1B2,       0,      1]]))
        # f = 0.935
        # fac1 = abs(A2B1)*abs(A2B2)
        # fac2 = abs(A1B1)*abs(A1B2)
        # Ys.append(np.array([[f,      fac1*A1B1,  fac1*A1B2],
                            # [0,      f,     A2B2*A2B1],
                            # [0,      0,     f]]))
        # Ys.append(np.array([[f,      fac2*A2B2,  fac2*A2B1],
                            # [0,      f,     A1B2*A1B1],
                            # [0,      0,     f]]))
        # Ys.append(np.array([[1,      A2B1sq*A2B2sq*A1B1+A1B2sq*A1B1sq*A2B2,  A2B1sq*A2B2sq*A1B2+A1B2sq*A1B1sq*A2B1 ],
                          # [0,      1,     A2B2*A2B1+A1B2*A1B1 ],
                          # [0,      0,     1,   ]]))
        # if abs(A2B1) >= 0.9 or abs(A2B2) >= 0.9: 
        if abs(A2B1) >= 0.95: 
            # f = 1.00
            f = 0.90
            Ys.append(np.array([[f,      A1B1,  A1B2],
                                [0,      f,     A2B2*A2B1],
                                [0,      0,     f]]))
        # if abs(A1B2) >= 0.1 or abs(A1B1) >= 0.1: 
            # Ys.append(np.array([[1,      A1B2sq*A1B1sq*A2B2,  A1B2sq*A1B1sq*A2B1 ],
                              # [0,      1,     A1B2*A1B1 ],
                              # [0,      0,     1,   ]]))
        # if True:
            # Ys.append([[1,      (A2B1**2)*A1B1,  (A2B1**2)*A1B2 ],
                       # [0,      1,               A2B1*A2B2      ],
                       # [0,      0,               1              ]])
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
            if XIsPSD and not YIsPSD:
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
        if a == 0 and j % (samples//10) == 0:
            print(100.0*float(j)/samples, '%')
    return numInequal, totalDone

if not draw:

    # Get the points in parallel
    pool = Pool(threads)
    numsNotPSD = pool.map(testPoints, range(threads))
    totalChecks = sum([x[1] for x in numsNotPSD])
    totalNotPSD = sum([x[0] for x in numsNotPSD])
    if totalChecks == 0:
        print("No checks done")
    else:
        percentPSD = 100.0*(1.0-float(totalNotPSD)/totalChecks)
        print("Total checks: ", totalChecks)
        print("Percent correct: ", percentPSD)
    exit()

# Get the points in parallel
pool = Pool(threads)
points = pool.map(getPoints, range(threads))
points = np.array(points)
print("Points shape: ", points.shape)
if mode == 3:
    points = points.reshape(-1, 3)
elif mode == 10:
    points = points.reshape(-1, 10)

# 2D with sliders
# fig = plt.figure()
# ax = plt.gca()
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# axSlider1 = fig.add_axes([0.2, 0.95, 0.65, 0.03])
# axSlider2 = fig.add_axes([0.2, 0.90, 0.65, 0.03])
# val1 = 1
# val2 = 1
# points = np.array(points)
# def update():
    # global val1, val2
    # xPoints = []
    # yPoints = []
    # for i in range(len(points)):
        # if abs(points[i,7]-val1) <= 1e-1 and abs(points[i,8]-val2) <= 1e-1:
            # xPoints.append(points[i,5])
            # yPoints.append(points[i,6])
    # ax.clear()
    # ax.set_xlim(-1, 1)
    # ax.set_ylim(-1, 1)
    # ax.scatter(xPoints, yPoints, color='r', s=2)
# def setVal1(val):
    # global val1
    # val1 = val
    # update()
# def setVal2(val):
    # global val2
    # val2 = val
    # update()
# slider1 = Slider(axSlider1, 'Val 1', -1.5, 1.5, valinit=1, valfmt='%f')
# slider2 = Slider(axSlider2, 'Val 2', -1.5, 1.5, valinit=1, valfmt='%f')
# slider1.on_changed(setVal1)
# slider2.on_changed(setVal2)
# update()
# plt.show()

# 3D
fig = plt.figure()
ax = plt.gca(projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('y')
points = np.array(points)
if mode == 3:
    xInd = 0
    yInd = 1
    zInd = 2
elif mode == 10:
    xInd = 5
    yInd = 6
    zInd = 8
    newPoints = []
    for i in range(len(points)):
        if abs(points[i,7]-fourthVal) <= 1e-1:
            newPoints.append(points[i,:])
    points = np.array(newPoints)
def update():
    ax.clear()
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(-1.2, 1.2)
    ax.scatter(points[:,xInd], points[:,yInd], points[:,zInd], color='r', s=2)
update()
plt.show()









