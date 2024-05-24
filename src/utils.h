#pragma once
#include <vector>
#include <complex>
#include <set>
#include <Eigen/Dense>
#include "poly.h"

// Given a matrix and a variable list, return the matrix with the variables replaced
Eigen::MatrixXcd replaceVariables(std::vector<std::vector<Poly>>& momentMatrix, const std::vector<Mon>& variables, const std::vector<std::complex<double>>& varVals);

// Get the eigenvalues and vectors of a matrix after replacement
void getEigens(std::vector<std::vector<Poly>>& momentMatrix, const std::vector<Mon>& variables, const std::vector<std::complex<double>>& varVals, std::vector<std::vector<std::complex<double>>>& eigenvectors, std::vector<std::complex<double>>& eigenvalues);

// Add all single order monomials from a functional to a list of variables
void addSingleMonomials(std::vector<Mon>& variables, Poly functional);

// Generate a moment matrix given the top row as polynomials
std::vector<std::vector<Poly>> generateFromTopRow(std::vector<Poly> monomsInTopRow, int verbosity);

// Generate a list of monomials of a certain level
std::vector<Poly> generateMonomials(std::vector<Mon> variables, int level, int verbosity);

// Generate all moment matrices for a given level given a polynomial
std::vector<std::vector<std::vector<Poly>>> generateAllMomentMatrices(const Poly& functional, std::vector<Poly> zeroCons, int level, int verbosity);

// Add variables from a moment matrix
void addVariables(std::set<Mon>& variables, std::vector<std::vector<Poly>> toAdd);

// Add variables from a polynomial
void addVariables(std::set<Mon>& variables, Poly toAdd);

// Convert from a matrix location to an svec location
int matLocToVecLoc(int i, int j, int n);

// Take the trace of a matrix, assuming it's real
double tr(Eigen::MatrixXcd A);

// Convert a set to a vector
std::vector<Mon> toVector(std::set<Mon> s);

// Generate random number between min and max
double rand(double min, double max);

