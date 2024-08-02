#pragma once
#include <vector>
#include "poly.h"

double solveOptim(Poly& objective, std::vector<Poly>& cons, std::vector<std::vector<std::vector<Poly>>>& momentMatrix, std::map<Mon, std::complex<double>>& startVals, int verbosity=1, int maxIters=1000000, int numExtra=0, double distance=1, double tolerance=1e-8);
