
from ncpol2sdpa import *

level = 1
# I = [[ 0,   -1,    0 ],
     # [-1,    1,    1 ],
     # [ 0,    1,   -1 ]]
I = [
        [  0,   0,   -1,   0, -1,  -2,   0],
        [  0,   0,    0,   0,  1,   1,   -1],
        [ -1,   0,    0,   0,  0,   1,   1],
        [  0,   0,    0,   0,  -1,   1,   0],
        [ -1,   1,    0,  -1,  0,   0,   0],
        [ -2,   1,    1,   1,  0,   0,   0],
        [  0,  -1,    1,   0,  0,   0,   0],
    ]

P = Probability([3, 3], [3, 3])
objective = -2*define_objective_with_I(I, P)

sdp = SdpRelaxation(P.get_all_operators(), verbose=2)
sdp.get_relaxation(level, objective=objective, 
                   substitutions=P.substitutions)
sdp.solve()
print(sdp.primal)

sdp.write_to_file('sdp.csv')


