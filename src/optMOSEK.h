#pragma once
#include <vector>
#include "poly.h"

// Convert to MOSEK form and solve
double solveMOSEK(Poly obj, std::vector<std::vector<std::vector<Poly>>>& psd, std::vector<Poly> constraintsZero, std::vector<Poly> constraintsPositive, int verbosity, std::pair<int,int> varBounds = {-1, 1}, std::map<Mon, std::complex<double>>* solution = nullptr);


