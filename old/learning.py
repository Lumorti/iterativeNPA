
import numpy as np
import cvxpy as cp
import random
import math
import itertools
import random

numPoints = 1000

inputs = []
validity = []

for i in range(numPoints):

    X = cp.Variable((4,4), symmetric=True)
    # x = cp.Variable()
    # y = cp.Variable()
    # z = cp.Variable()
    # a = cp.Variable()
    b = cp.Variable()
    c = cp.Variable()

    x = random.uniform(-1,1)
    y = random.uniform(-1,1)
    z = random.uniform(-1,1)
    a = random.uniform(-1,1)

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
        X >> 0
    ]

    prob = cp.Problem(cp.Maximize(0), cons)   
    prob.solve(solver=cp.MOSEK, mosek_params={"MSK_IPAR_NUM_THREADS": 1})
    print((100.0*i) / numPoints, x, y, z, prob.status)

    feasability = 0.0
    if str(prob.status) == "optimal":
        feasability = 1.0

    inputs.append([x,y,z,a])
    validity.append(feasability)

