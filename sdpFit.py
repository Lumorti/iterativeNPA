import sympy as sp
import numpy as np
import cvxpy as cp
import random
import math
import itertools

pointsOnEdge = [
         [1.0/(2*math.sqrt(2-math.sqrt(3))), 1.0/(2*math.sqrt(2-math.sqrt(3))), 1.0/2.0],
         [math.sqrt((1.0+math.sqrt(391.0)/20.0)/2.0), math.sqrt((1.0+math.sqrt(391.0)/20.0)/2.0), 3.0/20.0],
         [math.sqrt(2.0+math.sqrt(2.0)) / 2.0, math.sqrt(2.0+math.sqrt(2.0)) / 2.0, 1.0 / math.sqrt(2.0)],
         [0.866, 0.866, 0.866],
         [1.0, 1.0, 0.0],
         [1.0, 0.0, 0.0],
    ]

print("Generating edge points...")
allEdgePoints = set([])
for i in range(len(pointsOnEdge)):

    # for each point, consider each term being negative or not
    pointWithNegatives = []
    for j in range(2**3):
        point = []
        for k in range(3):
            if j & (1 << k) and abs(pointsOnEdge[i][k]) >= 1e-6:
                point.append(-pointsOnEdge[i][k])
            else:
                point.append(pointsOnEdge[i][k])
        
        # Generate permutations of this point
        perms = set(itertools.permutations(point))
        allEdgePoints = allEdgePoints.union(perms)
        
print("Generated " + str(len(allEdgePoints)) + " edge points")
print("Constructing problem...")
matSize = 6
numVars = int(4*matSize*(matSize+1)/2)
coeffs = [cp.Variable() for i in range(numVars)]
lam = 1e-5
identity = np.identity(matSize)
mats = []
cons = []
for coeff in coeffs:
    cons.append(coeff >= -1)
    cons.append(coeff <= 1)
for point in allEdgePoints:
    newMat = cp.Variable((matSize,matSize), symmetric=True)
    cons.append(newMat >> 0)
    cons.append(lam*identity - newMat >> 0)
    index = 0
    for i in range(matSize):
        for j in range(i, matSize):
            cons.append(newMat[i,j] == coeffs[index] + point[0]*coeffs[index+1] + point[1]*coeffs[index+2] + point[2]*coeffs[index+3])
            index += 4
    mats.append(newMat)

print("Solving problem...")
prob = cp.Problem(cp.Maximize(sum(coeffs)), cons)
prob.solve(solver=cp.MOSEK, verbose=True, mosek_params={"MSK_IPAR_NUM_THREADS": 8})
print(prob.status)
print(prob.value)
print([coeff.value for coeff in coeffs])
exit()


