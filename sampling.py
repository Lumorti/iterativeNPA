
import sympy as sp
import numpy as np
import cvxpy as cp
import random


def checkPoint(x, y, z, a):
    X = cp.Variable((4,4), symmetric=True)
    cons = [
        X[0,0] == 1,
        X[1,1] == 1,
        X[2,2] == 1,
        X[3,3] == 1,
        X[0,2] == x,
        X[0,3] == y,
        X[1,2] == z,
        X[1,3] == a,
        X >> 0
    ]
    prob = cp.Problem(cp.Minimize(0), cons)
    try:
        prob.solve(solver=cp.MOSEK, mosek_params={"MSK_IPAR_NUM_THREADS": 1})
        return prob.status == "optimal"
    except:
        return False

def sample(x, y, z, a, numToSample):
    numValid = 0
    for i in range(numToSample):
        m = 2.0
        b = random.uniform(-x, x)
        c = random.uniform(-a, a)
        mat = np.array([[1, b, x, y], [b, 1, z, a], [x, z, 1, c], [y, a, c, 1]])
        eig = np.linalg.eigvals(mat)
        if np.all(eig > 0):
            numValid += 1
    return numValid

numGood = 0
numGoodWhenValid = 0
numValid = 0
pointsToTest = 200
numToSample = 100
for i in range(pointsToTest):
    x = random.uniform(-1, 1)
    y = random.uniform(-1, 1)
    z = random.uniform(-1, 1)
    a = random.uniform(-1, 1)
    pointIsValid = checkPoint(x, y, z, a)
    numSamplesValid = sample(x, y, z, a, numToSample)
    isGood = (pointIsValid and numSamplesValid > 0) or (not pointIsValid and numSamplesValid == 0)
    if pointIsValid:
        print("Point is valid, number of valid samples: ", numSamplesValid, " ", isGood)
        numValid += 1
        if isGood:
            numGoodWhenValid += 1
    else:
        print("Point is inval, number of valid samples: ", numSamplesValid, " ", isGood)
    if isGood:
        numGood += 1
numGoodWhenInvalid = numGood - numGoodWhenValid
numInvalid = pointsToTest - numValid
print(numGood, "/", pointsToTest, "overall")
print(numGoodWhenValid, "/", numValid, "when valid")
print(numGoodWhenInvalid, "/", numInvalid, "when invalid")

