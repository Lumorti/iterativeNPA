
import sympy as sp
import numpy as np
from sympy import sqrt

x = sp.symbols('x')
y = sp.symbols('y')
z = sp.symbols('z')
a = sp.symbols('a0:36')
a1 = sp.symbols('a1')
a2 = sp.symbols('a2')
a3 = sp.symbols('a3')
a4 = sp.symbols('a4')
a5 = sp.symbols('a5')
a6 = sp.symbols('a6')
a7 = sp.symbols('a7')
a8 = sp.symbols('a8')
a9 = sp.symbols('a9')
a10 = sp.symbols('a10')
a11 = sp.symbols('a11')
a12 = sp.symbols('a12')
a13 = sp.symbols('a13')
a14 = sp.symbols('a14')
a15 = sp.symbols('a15')
a16 = sp.symbols('a16')
a17 = sp.symbols('a17')
a18 = sp.symbols('a18')
a19 = sp.symbols('a19')
a20 = sp.symbols('a20')
a21 = sp.symbols('a21')
a22 = sp.symbols('a22')
a23 = sp.symbols('a23')
a24 = sp.symbols('a24')
a25 = sp.symbols('a25')
a26 = sp.symbols('a26')
a27 = sp.symbols('a27')
a28 = sp.symbols('a28')
a29 = sp.symbols('a29')
a30 = sp.symbols('a30')
a31 = sp.symbols('a31')
a32 = sp.symbols('a32')
a33 = sp.symbols('a33')
a34 = sp.symbols('a34')
a35 = sp.symbols('a35')
poly = 1+x*a[2]+y*a[3]+z*a[4]+x**2*a[5]+x*y*a[6]+x*z*a[7]+y**2*a[8]+y*z*a[9]+z**2*a[10]+x**3*a[11]+x**2*y*a[12]+x**2*z*a[13]+x*y**2*a[14]+x*y*z*a[15]+x*z**2*a[16]+y**3*a[17]+y**2*z*a[18]+y*z**2*a[19]+z**3*a[20]+x**4*a[21]+x**3*y*a[22]+x**3*z*a[23]+x**2*y**2*a[24]+x**2*y*z*a[25]+x**2*z**2*a[26]+x*y**3*a[27]+x*y**2*z*a[28]+x*y*z**2*a[29]+x*z**3*a[30]+y**4*a[31]+y**3*z*a[32]+y**2*z**2*a[33]+y*z**3*a[34]+z**4*a[35]
# poly = a10*z**2 + a22*x**3*y + a23*x**3*z + a27*x*y**3 + a28*x*y**2*z + a3*x*y*z/3 - a3*y**3 + a3*y + a31*y**4 + a32*y**3*z + a33*x**2*y**2 + a33*x**2*z**2 + a33*y**2*z**2 + a33 + a5*x**2 + a9*y*z - x**4*(a33 + a5) + sqrt(3)*x**3*(-3*a10 + 3*a31 - 13*a33 - 3*a5)/6 + x**2*y*z*(2*sqrt(3)*a3 - 3*a9)/9 + x*y*z**2*(3*a22 + 3*a27 + 2*sqrt(3)*a3)/9 - x*y*(a22 + a27) - x*z**3*(3*a10 + 3*a23 - 9*a28 - 3*a31 + 13*a33 + 3*a5)/3 - x*z*(-3*a10 + 9*a28 + 3*a31 - 13*a33 - 3*a5)/3 + sqrt(3)*x*(3*a10 - 3*a31 + 13*a33 + 3*a5)/6 - y**2*(a31 + a33) - y*z**3*(a32 + a9) - z**4*(a10 + a33) + sqrt(3)*z**3*(-3*a10 + 3*a31 - 13*a33 - 3*a5)/6 + sqrt(3)*z*(3*a10 - 3*a31 + 13*a33 + 3*a5)/6

# rt3o2 = np.sqrt(3)/2
rt3o2 = sp.sqrt(3)/2
eqns = []
eqns.append(poly.subs({x:rt3o2, y:rt3o2, z:rt3o2}))
eqns.append(poly.subs({x:-rt3o2, y:rt3o2, z:rt3o2}))
eqns.append(poly.subs({x:-rt3o2, y:-rt3o2, z:rt3o2}))
eqns.append(poly.subs({x:rt3o2, y:rt3o2, z:-rt3o2}))
eqns.append(poly.subs({x:-rt3o2, y:rt3o2, z:-rt3o2}))
eqns.append(poly.subs({x:rt3o2, y:-rt3o2, z:-rt3o2}))
eqns.append(poly.subs({x:-rt3o2, y:-rt3o2, z:-rt3o2}))
eqns.append(poly.subs({x:1, y:1, z:0}))
eqns.append(poly.subs({x:0, y:1, z:1}))
eqns.append(poly.subs({x:1, y:0, z:1}))
eqns.append(poly.subs({x:-1, y:1, z:0}))
eqns.append(poly.subs({x:0, y:-1, z:1}))
eqns.append(poly.subs({x:-1, y:0, z:1}))
eqns.append(poly.subs({x:1, y:-1, z:0}))
eqns.append(poly.subs({x:0, y:1, z:-1}))
eqns.append(poly.subs({x:1, y:0, z:-1}))
eqns.append(poly.subs({x:-1, y:-1, z:0}))
eqns.append(poly.subs({x:0, y:-1, z:-1}))
eqns.append(poly.subs({x:-1, y:0, z:-1}))
eqns.append(poly.subs({x:1, y:0, z:0}))
eqns.append(poly.subs({x:0, y:1, z:0}))
eqns.append(poly.subs({x:0, y:0, z:1}))
eqns.append(poly.subs({x:-1, y:0, z:0}))
eqns.append(poly.subs({x:0, y:-1, z:0}))
eqns.append(poly.subs({x:0, y:0, z:-1}))
eqns.append(poly.subs({x:0, y:0, z:0})-1)
eqns.append(poly.subs({x:1, y:1, z:1})+1)

equalitiesUsed = []
# for j in range(len(eqns)):
for j in range(len(eqns)-1, 0, -1):
    eqn = sp.simplify(eqns[j])
    if len(eqn.free_symbols) > 0:
        toRemove = next(iter(eqn.free_symbols))
        print("for equation:", eqn)
        solved = sp.solve(eqn, toRemove)
        if len(solved) > 0:
            print("setting:", toRemove, "->", solved[0])
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
print(str(poly).replace('**', '^').replace('sqrt', '\[Sqrt]'))

