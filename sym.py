

from sympy.solvers import solve
from sympy import Symbol, Matrix, symbols, pprint

numVars = 10
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
a = symbols('a0:%d' % numVars)
desired = 1 - x**2 - y**2 + z*x**2*y**2
# desired = 1 - x**2 - y**2
# M = Matrix([[1, a[0], a[1], a[2]], 
            # [a[0], 1, a[3], a[4]],
            # [a[1], a[3], 1, a[5]],
            # [a[2], a[4], a[5], 1]])
M = Matrix([[1, a[0], a[1], a[2], a[3]], 
            [a[0], 1, a[4], a[5], a[6]],
            [a[1], a[4], 1, a[7], a[8]],
            [a[2], a[5], a[7], 1, a[9]],
            [a[3], a[6], a[8], a[9], 1]])

# get the determinant of the matrix
print(desired)
print(M.det())

posibles = [1, x, y, 0, z]
threads = 8
import itertools
from more_itertools import chunked, distribute
from multiprocessing import Pool

det = M.det()

def perCore(things):
    firstCore = False
    count = 0
    total = len(posibles)**numVars / threads
    for p in things:
        if list(p) == [1 for i in range(numVars)]:
            firstCore = True
        break
    for p in things:
        newDet = det.subs(zip(a, p))
        if firstCore:
            count += 1
            if count % 100 == 0:
                print(100*count/total, "%             \r", end="")
        if newDet == desired:
            print()
            print("found = ", p)
            print("det = ", newDet)
            print("M = ")
            pprint(M.subs(zip(a, p)))
            print()

with Pool(threads) as pool:
    things = distribute(threads, itertools.product(posibles, repeat=numVars))
    pool.map(perCore, things)




