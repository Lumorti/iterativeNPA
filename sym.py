from sympy.solvers import solve
from sympy import Symbol, Matrix, symbols, pprint, Poly
from math import prod
import itertools
from more_itertools import chunked, distribute
from multiprocessing import Pool

numVars = 10
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
a = symbols('a0:%d' % numVars)
x2 = x**2
y2 = y**2
z2 = z**2
s1 = x2 + y2 + z2
s2 = x2*y2 + z2*(x2 + y2)
N = Matrix([[2, s1],
            [s1, 1+s2]])
# N = Matrix([[2, s1, s1],
            # [s1, 1+s2, s1],
            # [s1, s1, 1+s2]])
desired = N.det()
# desired = 2 - x**4 - y**4 - z**4
# desired = 1 - x**2 - y**2
# M = Matrix([[a[6], a[0], a[1], a[2]], 
            # [a[0], a[7], a[3], a[4]],
            # [a[1], a[3], a[8], a[5]],
            # [a[2], a[4], a[5], a[9]]])
# M = Matrix([[2, a[0], a[1], a[2]], 
            # [a[0], 1, a[3], a[4]],
            # [a[1], a[3], 1, a[5]],
            # [a[2], a[4], a[5], 1]])
M = Matrix([[1, a[0], a[1], a[2], a[3]], 
            [a[0], 1, a[4], a[5], a[6]],
            [a[1], a[4], 1, a[7], a[8]],
            [a[2], a[5], a[7], 1, a[9]],
            [a[3], a[6], a[8], a[9], 1]])

def getMonoms(thing):
    asPoly = Poly(thing, [x, y, z])
    monoms = [prod(x**k for x, k in zip(asPoly.gens, mon)) for mon in asPoly.monoms()]
    # monoms.sort(key=lambda x: str(x))
    return set(monoms)

desiredMonoms = getMonoms(desired)
det = M.det()

print(desired)
print(M.det())
print(desiredMonoms)

posibles = [0, x, y, z, 1]
threads = 8

def perCore(things):
    firstCore = False
    count = 0
    total = len(posibles)**numVars / threads
    for p in things:
        if list(p) == [posibles[0] for i in range(numVars)]:
            firstCore = True
        break
    for p in things:
        if firstCore:
            count += 1
            if count % 100 == 0:
                print(100*count/total, "%             \r", end="")
        newDet = det.subs(zip(a, p))
        newDetMonoms = getMonoms(newDet)
        # if newDet == desired:
        if newDetMonoms == desiredMonoms:
        # if desiredMonoms.issubset(newDetMonoms):
            print()
            print("found = ", p)
            print("det = ", newDet)
            print("M = ")
            pprint(M.subs(zip(a, p)))
            print()

with Pool(threads) as pool:
    things = distribute(threads, itertools.product(posibles, repeat=numVars))
    pool.map(perCore, things)



