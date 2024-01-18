

from sympy.solvers import solve
from sympy import Symbol, Matrix, symbols

x = Symbol('x')
y = Symbol('y')
numVars = 10
a = symbols('a0:%d' % numVars)
desired = 1 - x**2 - y**2 + 0.5*x**2*y**2
M = Matrix([[a[0], a[1], a[2], a[3]], 
            [a[1], a[4], a[5], a[6]],
            [a[2], a[5], a[7], a[8]],
            [a[3], a[6], a[8], a[9]]])

# get the determinant of the matrix
print(desired)
print(M.det())

posibles = [1, x, y, 0]
threads = 10
import itertools
from more_itertools import chunked, distribute
from multiprocessing import Pool

det = M.det()

def perCore(things):
    firstCore = False
    count = 0
    total = len(posibles)**numVars
    for p in things:
        if p == (1, 1, 1, 1, 1, 1, 1, 1, 1, 1):
            firstCore = True
        newDet = det.subs(zip(a, p))
        if firstCore:
            count += 1
            if count % 100 == 0:
                print(100*count/total, "%             \r", end="")
        if newDet == desired:
            print("found = ", p)
            break

with Pool(threads) as pool:
    things = distribute(threads, itertools.product(posibles, repeat=numVars))
    pool.map(perCore, things)




