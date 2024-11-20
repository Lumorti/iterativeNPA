from sympy.solvers import solve
from sympy import Symbol, Matrix, symbols, pprint, Poly, sqrt

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
b = Symbol('b')
c = Symbol('c')
a = Symbol('a')
M = Matrix([[1, b, x, y], 
            [b, 1, z, 0],
            [x, z, 1, c],
            [y, 0, c, 1]])

posCons = []
posCons.append(M.det())
for i in range(4):
    newCon = M.minor_submatrix(i, i).det()
    if newCon not in posCons:
        posCons.append(newCon)
    for j in range(3):
        newCon = M.minor_submatrix(i, i).minor_submatrix(j, j).det()
        if newCon not in posCons:
            posCons.append(newCon)

for con in posCons:
    print(con)

varsToElim = [c, b]

# For each variable to eliminate
posConSets = [posCons]
for var in varsToElim:

    print("")
    print("Trying to eliminate " + str(var))
    print("")

    # For each branching path
    newBranches = []
    for posCons2 in posConSets:
    
        # Consider when each constraint is zero
        for con in posCons2:
            if con.has(var):

                # Solve for the variable
                solveFor = solve(con, var)

                # Replace the variable in the other constraints
                newCons = []
                for i in range(len(posCons)):
                    replaced = posCons[i].subs(var, solveFor[0])
                    simplified = replaced.simplify()
                    newCons.append(simplified)

                print("")
                print("New Constraints")
                for con2 in newCons:
                    print(con2)
                print("")

                newBranches.append(newCons)

    posConSets = newBranches

print("")
print("Final Constraints")
print("")
for posCons2 in posConSets:
    print("")
    for con in posCons2:
        print(con)


