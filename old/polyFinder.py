
import sympy as sp
import numpy as np
from sympy import sqrt
import itertools
import math

rt3o2 = sp.sqrt(3)/2

pointsOnEdge = [
         [1.0/(2*sqrt(2-sqrt(3))), 1.0/(2*sqrt(2-sqrt(3))), 1.0/2.0],
         [sqrt((1.0+sqrt(391.0)/20.0)/2.0), sqrt((1.0+sqrt(391.0)/20.0)/2.0), 3.0/20.0],
         [sqrt(2.0+sqrt(2.0)) / 2.0, sqrt(2.0+sqrt(2.0)) / 2.0, 1.0 / sqrt(2.0)],
         [rt3o2, rt3o2, rt3o2],
         [1.0, 1.0, 0.0],
         [1.0, 0.0, 0.0],
    ]
for i in range(len(pointsOnEdge)):
    pointsOnEdge[i] = [float(x) for x in pointsOnEdge[i]]
print("Generating edge points...")
allEdgePoints = set([])
for i in range(len(pointsOnEdge)):
    for j in range(2**3):
        point = []
        for k in range(3):
            if j & (1 << k) and abs(pointsOnEdge[i][k]) >= 1e-6:
                point.append(-pointsOnEdge[i][k])
            else:
                point.append(pointsOnEdge[i][k])
        perms = set(itertools.permutations(point))
        allEdgePoints = allEdgePoints.union(perms)
print("Generated", len(allEdgePoints), "edge points")

degree = 8
numVars = 3
monoms = list(itertools.combinations_with_replacement([i for i in range(numVars+1)], degree))
numCoeffs = len(monoms)
print("Monoms: ", monoms)
print("Number of coefficients: ", numCoeffs)

allEdgePoints = list(allEdgePoints)
negPoints = [
        [1, 1, 1],
        [-1, 1, 1],
        [-1, -1, 1],
        [-1, -1, -1],
        [1, -1, -1],
    ]
posPoints = [
        [0, 0, 0],
    ]
allEdgePoints.extend(negPoints)
allEdgePoints.extend(posPoints)
A = np.zeros((len(allEdgePoints), numCoeffs))
b = np.zeros((len(allEdgePoints), 1))
b[-len(negPoints)-1:-1] = -1
b[-1] = 1
print(b)

for row, point in enumerate(allEdgePoints):
    for col, monom in enumerate(monoms):
        el = 1
        for i in range(degree):
            if monom[i] != 0:
                el *= point[monom[i]-1]
        A[row, col] = el

sol = np.linalg.lstsq(A, b, rcond=None)
x = sol[0]
coeffs = x.transpose()[0]

print("largest residual: ", np.max(np.abs(A @ x - b)))

polyAsString = ""
for i, coeff in enumerate(coeffs):
    if abs(coeff) > 1e-6:
        if i > 0:
            if coeff > 0:
                polyAsString += " + " + str(coeff)
            else:
                polyAsString += " - " + str(-coeff)
        else:
            polyAsString += str(coeff)

        numX = 0
        numY = 0
        numZ = 0
        for j in range(degree):
            if monoms[i][j] == 1:
                numX += 1
            elif monoms[i][j] == 2:
                numY += 1
            elif monoms[i][j] == 3:
                numZ += 1
        if numX == 1:
            polyAsString += "*x"
        elif numX > 1:
            polyAsString += "*x**" + str(numX)
        if numY == 1:
            polyAsString += "*y"
        elif numY > 1:
            polyAsString += "*y**" + str(numY)
        if numZ == 1:
            polyAsString += "*z"
        elif numZ > 1:
            polyAsString += "*z**" + str(numZ) 

print()
print(polyAsString + " >= 0")
print()
print(polyAsString.replace("**", "^").replace("e", "*10^"))

exit()

# x = sp.symbols('x')
# y = sp.symbols('y')
# z = sp.symbols('z')
# a = sp.symbols('a0:36')
# a1 = sp.symbols('a1')
# a2 = sp.symbols('a2')
# a3 = sp.symbols('a3')
# a4 = sp.symbols('a4')
# a5 = sp.symbols('a5')
# a6 = sp.symbols('a6')
# a7 = sp.symbols('a7')
# a8 = sp.symbols('a8')
# a9 = sp.symbols('a9')
# a10 = sp.symbols('a10')
# a11 = sp.symbols('a11')
# a12 = sp.symbols('a12')
# a13 = sp.symbols('a13')
# a14 = sp.symbols('a14')
# a15 = sp.symbols('a15')
# a16 = sp.symbols('a16')
# a17 = sp.symbols('a17')
# a18 = sp.symbols('a18')
# a19 = sp.symbols('a19')
# a20 = sp.symbols('a20')
# a21 = sp.symbols('a21')
# a22 = sp.symbols('a22')
# a23 = sp.symbols('a23')
# a24 = sp.symbols('a24')
# a25 = sp.symbols('a25')
# a26 = sp.symbols('a26')
# a27 = sp.symbols('a27')
# a28 = sp.symbols('a28')
# a29 = sp.symbols('a29')
# a30 = sp.symbols('a30')
# a31 = sp.symbols('a31')
# a32 = sp.symbols('a32')
# a33 = sp.symbols('a33')
# a34 = sp.symbols('a34')
# a35 = sp.symbols('a35')
# poly = 1+x*a[2]+y*a[3]+z*a[4]+x**2*a[5]+x*y*a[6]+x*z*a[7]+y**2*a[8]+y*z*a[9]+z**2*a[10]+x**3*a[11]+x**2*y*a[12]+x**2*z*a[13]+x*y**2*a[14]+x*y*z*a[15]+x*z**2*a[16]+y**3*a[17]+y**2*z*a[18]+y*z**2*a[19]+z**3*a[20]+x**4*a[21]+x**3*y*a[22]+x**3*z*a[23]+x**2*y**2*a[24]+x**2*y*z*a[25]+x**2*z**2*a[26]+x*y**3*a[27]+x*y**2*z*a[28]+x*y*z**2*a[29]+x*z**3*a[30]+y**4*a[31]+y**3*z*a[32]+y**2*z**2*a[33]+y*z**3*a[34]+z**4*a[35]
# poly = a10*z**2 + a22*x**3*y + a23*x**3*z + a27*x*y**3 + a28*x*y**2*z + a3*x*y*z/3 - a3*y**3 + a3*y + a31*y**4 + a32*y**3*z + a33*x**2*y**2 + a33*x**2*z**2 + a33*y**2*z**2 + a33 + a5*x**2 + a9*y*z - x**4*(a33 + a5) + sqrt(3)*x**3*(-3*a10 + 3*a31 - 13*a33 - 3*a5)/6 + x**2*y*z*(2*sqrt(3)*a3 - 3*a9)/9 + x*y*z**2*(3*a22 + 3*a27 + 2*sqrt(3)*a3)/9 - x*y*(a22 + a27) - x*z**3*(3*a10 + 3*a23 - 9*a28 - 3*a31 + 13*a33 + 3*a5)/3 - x*z*(-3*a10 + 9*a28 + 3*a31 - 13*a33 - 3*a5)/3 + sqrt(3)*x*(3*a10 - 3*a31 + 13*a33 + 3*a5)/6 - y**2*(a31 + a33) - y*z**3*(a32 + a9) - z**4*(a10 + a33) + sqrt(3)*z**3*(-3*a10 + 3*a31 - 13*a33 - 3*a5)/6 + sqrt(3)*z*(3*a10 - 3*a31 + 13*a33 + 3*a5)/6


print("Generated", len(allEdgePoints), "edge points")
print("Generating equations...")
for point in allEdgePoints:
    eqns.append(poly.subs({x:point[0], y:point[1], z:point[2]}))
print("Generated", len(eqns), "equations")
        
eqns = []
equalitiesUsed = []
for j in range(len(eqns)-1, 0, -1):
    eqn = sp.simplify(eqns[j])
    if len(eqn.free_symbols) > 0:
        toRemove = next(iter(eqn.free_symbols))
        print("For equation:", eqn)
        solved = sp.solve(eqn, toRemove)
        if len(solved) > 0:
            print("Setting:", toRemove, "->", solved[0])
            equalitiesUsed.append([toRemove, solved[0]])
            for i in range(len(eqns)):
                eqns[i] = sp.simplify(eqns[i].subs({toRemove:solved[0]}))
            poly = sp.simplify(poly.subs({toRemove:solved[0]}))

print()
for eqn in equalitiesUsed:
    print(eqn)
print()
print(poly)
print(poly.free_symbols)
print(len(poly.free_symbols))
b = sp.symbols('b0:'+str(len(poly.free_symbols)-3))
subDict = {}
symbols = list(poly.free_symbols)
bInd = 0
for i in range(len(symbols)):
    if symbols[i] != x and symbols[i] != y and symbols[i] != z:
        subDict[symbols[i]] = b[bInd]
        bInd += 1
print(subDict)
poly = poly.subs(subDict)
print(poly)
# print(str(poly).replace('**', '^').replace('sqrt', '\[Sqrt]'))
exit()

# Loop through all possible combinations of 1, -1, 0 for the coefficients b
choices = [1, -1, 0]
shouldBePositive = [
        {x: 0, y: 0.5, z: 0},
        {x: 0.2, y: 0.0, z: 0},
        {x: 0.2, y: 0.0, z: 0.2},
        {x: 0.7, y: -0.3, z: 0.2},
        {x: -0.2, y: -0.3, z: 0.5},
        {x: -0.2, y: 0.3, z: -0.5}
    ]
shouldBeNegative = [
        {x: 0, y: 2, z: 0},
        {x: 2, y: 0, z: 0},
        {x: 0, y: 0, z: 2},
        {x: -1, y: -2, z: 0},
        {x: -2, y: 0, z: 1},
        {x: 1, y: 1, z: 1},
        {x: -1, y: 1, z: 1},
        {x: -1, y: -1, z: 1},
        {x: -1, y: -1, z: -1},
    ]
for i in range(len(choices)**len(b)):
    if i % 1000 == 0:
        print(i, "/", len(choices)**len(b))
    bVals = []
    for j in range(len(b)):
        bVals.append(choices[(i//(len(choices)**j))%len(choices)])
    polyWithoutB = poly.subs({b[i]:bVals[i] for i in range(len(b))})
    allGood = True
    for point in shouldBePositive:
        if polyWithoutB.subs(point) < 0:
            allGood = False
            break
    for point in shouldBeNegative:
        if polyWithoutB.subs(point) > 0:
            allGood = False
            break
    if allGood:
        print(bVals)
        print(polyWithoutB)
        print()
        break
