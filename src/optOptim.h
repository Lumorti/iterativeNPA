#pragma once
#include <vector>
#include "poly.h"

double solveOptim(Poly objective, std::vector<Poly> cons, std::vector<std::vector<Poly>> momentMatrix, int verbosity);
