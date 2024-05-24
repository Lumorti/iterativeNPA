#pragma once
#include <vector>
#include "poly.h"

// Convert to MOSEK form and solve
std::pair<double,double> solveMOSEK(Poly obj, std::vector<std::vector<std::vector<Poly>>>& psd, std::vector<Poly> constraintsZero, int verbosity);

// Convert to MOSEK form and solve
double maximizeMOSEK(Poly obj, std::vector<std::vector<std::vector<Poly>>>& psd, std::vector<Poly> constraintsZero, std::vector<Poly> constraintsPositive, int verbosity);


