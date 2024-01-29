from sympy import *

n = Symbol('n')
x = Symbol('x')
y = Symbol('y')
a = symbols('a0:4')
alpha = Symbol('alpha')
M = Matrix([[1,    x,    a[0], a[1]], 
            [x,    1,    a[2], a[3]],
            [a[0], a[2], 1,    y],
            [a[1], a[3], y,    1]])

det = M.det()
print(det)
print()

maxLogSumExp = exp(alpha*det)
# maxLogSumExp = series(exp(det), x)
print(maxLogSumExp)
print()

integratedX = integrate(maxLogSumExp, (x, -1, 1))
# integratedX = Integral(maxLogSumExp, (x, -1, 1))
# integratedX = integratedX.as_sum(method='midpoint', n=100).expand()
integratedX = simplify(integratedX)
print(integratedX)
print()

# integratedY = integrate(integratedX, (y, -1, 1))
# integratedY = simplify(integratedY)
# print(integratedY)
# print()

# logged = log(integratedY)
# print(logged)
# print()
