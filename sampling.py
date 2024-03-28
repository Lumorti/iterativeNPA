
import sympy as sp
import numpy as np
import cvxpy as cp
import random
import math
import itertools

# X = cp.Variable((5,5), symmetric=True)
# A1 = cp.Variable()
# A2 = cp.Variable()
# B1 = cp.Variable()
# B2 = cp.Variable()
# A1B1 = cp.Variable()
# A1B2 = cp.Variable()
# A2B1 = cp.Variable()
# A2B2 = cp.Variable()
# A1A2 = cp.Variable()
# B1B2 = cp.Variable()
# cons = [
    # X[0,0] == 1,
    # X[1,1] == 1,
    # X[2,2] == 1,
    # X[3,3] == 1,
    # X[4,4] == 1,
    # X[0,1] == A1,
    # X[0,2] == A2,
    # X[0,3] == B1,
    # X[0,4] == B2,
    # X[1,2] == A1A2,
    # X[1,3] == A1B1,
    # X[1,4] == A1B2,
    # X[2,3] == A2B1,
    # X[2,4] == A2B2,
    # X[3,4] == B1B2,
    # # c == math.sqrt(3)/2,
    # # x == 1.0/(2*math.sqrt(2-math.sqrt(3))),
    # # y == 1.0/(2*math.sqrt(2-math.sqrt(3))),
    # # y == 1.0,
    # X >> 0
# ]
# prob = cp.Problem(cp.Maximize(A1B1+A1B2+A2B1-A2B2), cons)
# prob.solve(solver=cp.MOSEK, mosek_params={"MSK_IPAR_NUM_THREADS": 1})
# np.set_printoptions(edgeitems=30, linewidth=100000, 
    # formatter=dict(float=lambda x: "%.3f" % x))
# print(prob.status)
# print(prob.value)
# print(X.value)
# asNum = np.array(X.value)
# LL = np.linalg.cholesky(asNum + 1e-6*np.eye(5))
# print(LL)
# print(LL@LL.T)
# exit()

X = cp.Variable((4,4), symmetric=True)
x = cp.Variable()
y = cp.Variable()
z = cp.Variable()
a = cp.Variable()
b = cp.Variable()
c = cp.Variable()
cons = [
    X[0,0] == 1,
    X[1,1] == 1,
    X[2,2] == 1,
    X[3,3] == 1,
    X[0,2] == x,
    X[0,3] == y,
    X[1,2] == z,
    X[1,3] == a,
    X[2,3] == c,
    X[0,1] == b,
    # a == 0,
    # z == (2+math.sqrt(2))/4.0,
    # c == math.sqrt(3)/2,
    # x == 1.0/(2*math.sqrt(2-math.sqrt(3))),
    # y == 1.0/(2*math.sqrt(2-math.sqrt(3))),
    a == 0.5,
    x == 1.0,
    y == 1.0,
    z == 1.0,
    X >> 0
]
prob = cp.Problem(cp.Maximize(x+y+z+a), cons)
prob.solve(solver=cp.MOSEK, mosek_params={"MSK_IPAR_NUM_THREADS": 1})
print(prob.status)
print(prob.value)
print(x.value, y.value, z.value)
print(X.value)

A = np.array(X.value)
print(np.linalg.det(A))
print(np.linalg.det(A[0:3,0:3]))
print(np.linalg.det(A[0:2,0:2]))
print(np.linalg.det(A[0:1,0:1]))
print("a = ", a.value)
print("b = ", b.value)
print("c = ", c.value)
print("x = ", x.value)
print("y = ", y.value)
print("z = ", z.value)
b = b.value
c = c.value
a = a.value
x = x.value
y = y.value
z = z.value

print(1 - a**2 - 2*b**2 + a**2*b**2 + b**4 - x**2 + a**2*x**2 + b**2*x**2 + 2*a*b*y - 2*a*b**3*y - 2*a*b*x**2*y - y**2 + b**2*y**2 + x**2*y**2 + 2*b*x*z - 2*a**2*b*x*z - 2*b**3*x*z + 4*a*b**2*x*y*z - 2*b*x*y**2*z - z**2 + a**2*z**2 + b**2*z**2 - 2*a*b*y*z**2 + y**2*z**2) 
exit()

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

