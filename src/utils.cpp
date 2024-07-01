#include "utils.h"
#include <iostream>
#include "printing.h"
#include <Eigen/Sparse>
#include <omp.h>

// Given a matrix and a variable list, return the matrix with the variables replaced
Eigen::MatrixXcd replaceVariables(std::vector<std::vector<Poly>>& momentMatrix, const std::map<Mon, std::complex<double>>& varVals) {

    // Replace each variable with its value
    Eigen::MatrixXcd momentMatrixEigen = Eigen::MatrixXcd::Zero(momentMatrix.size(), momentMatrix.size());
    for (size_t i=0; i<momentMatrix.size(); i++) {
        for (size_t j=i; j<momentMatrix[i].size(); j++) {
            for (auto& term : momentMatrix[i][j]) {
                if (term.first.isConstant()) {
                    momentMatrixEigen(i, j) += term.second;
                } else {
                    momentMatrixEigen(i, j) += term.second * varVals.at(term.first);
                }
            }
            momentMatrixEigen(j, i) = std::conj(momentMatrixEigen(i, j));
        }
    }
    return momentMatrixEigen;

}

// Given a matrix and a variable list, return the matrix with the variables replaced
Eigen::MatrixXcd replaceVariables(std::vector<std::vector<Poly>>& momentMatrix, const std::vector<Mon>& variables, const std::vector<std::complex<double>>& varVals) {

    // Replace each variable with its value
    Eigen::MatrixXcd momentMatrixEigen = Eigen::MatrixXcd::Zero(momentMatrix.size(), momentMatrix.size());
    for (size_t i=0; i<momentMatrix.size(); i++) {
        for (size_t j=0; j<momentMatrix[i].size(); j++) {
            for (auto& term : momentMatrix[i][j]) {
                for (size_t l=0; l<variables.size(); l++) {
                    if (term.first == variables[l]) {
                        if (i > j) {
                            momentMatrixEigen(i, j) += term.second * varVals[l];
                        } else {
                            momentMatrixEigen(i, j) += term.second * std::conj(varVals[l]);
                        }
                        break;
                    }
                }
            }
        }
    }
    return momentMatrixEigen;

}

// Get the eigenvalues and vectors of a matrix after replacement
void getEigens(std::vector<std::vector<Poly>>& momentMatrix, const std::vector<Mon>& variables, const std::vector<std::complex<double>>& varVals, std::vector<std::vector<std::complex<double>>>& eigenvectors, std::vector<std::complex<double>>& eigenvalues) {

    // Replace each variable with its value
    Eigen::MatrixXcd momentMatrixEigen = replaceVariables(momentMatrix, variables, varVals);

    // Get the eigenvalues and vectors of this
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(momentMatrixEigen);
    Eigen::MatrixXcd eigenVectorsEigen = es.eigenvectors();
    Eigen::VectorXcd eigenValuesEigen = es.eigenvalues();

    // Copy into the output vectors
    eigenvalues.clear();
    eigenvectors.clear();
    for (int i=0; i<eigenVectorsEigen.cols(); i++) {
        std::vector<std::complex<double>> eigenVector;
        for (int j=0; j<eigenVectorsEigen.rows(); j++) {
            eigenVector.push_back(eigenVectorsEigen(j, i));
        }
        eigenvectors.push_back(eigenVector);
    }
    for (long int i=0; i<eigenValuesEigen.size(); i++) {
        eigenvalues.push_back(eigenValuesEigen(i));
    }

}

// Add all single order monomials from a functional to a list of variables
void addSingleMonomials(std::vector<Mon>& variables, Poly functional) {

    // Iterate through the functional
    for (auto& term : functional) {

        // Iterate through the monomial
        for (size_t j=0; j<term.first.size(); j++) {
            Mon currentMonomial(term.first[j]);

            // Check if this moment is already in the list
            bool found = false;
            for (size_t k=0; k<variables.size(); k++) {
                if (variables[k] == currentMonomial) {
                    found = true;
                    break;
                }
            }

            // If not, add it in the correct place
            if (!found) {
                bool hasBeenAdded = false;
                for (size_t k=0; k<variables.size(); k++) {
                    if (variables[k][0] > currentMonomial[0]) {
                        variables.insert(variables.begin()+k, currentMonomial);
                        hasBeenAdded = true;
                        break;
                    }
                }
                if (!hasBeenAdded) {
                    variables.push_back(currentMonomial);
                }
            }

        }

    }

}

// Generate a moment matrix given the top row as polynomials
std::vector<std::vector<Poly>> generateFromTopRow(std::vector<Poly> monomsInTopRow, int verbosity) {

    // Generate all combinations of the top row
    std::vector<std::vector<Poly>> matrixToReturn = std::vector<std::vector<Poly>>(monomsInTopRow.size(), std::vector<Poly>(monomsInTopRow.size()));
    for (size_t i=0; i<monomsInTopRow.size(); i++) {
        for (size_t j=i; j<monomsInTopRow.size(); j++) {

            // Form the new polynomial
            Poly newPolynomial = monomsInTopRow[j].dagger() * monomsInTopRow[i];
            newPolynomial.reduce();

            // Set the matrix elements
            matrixToReturn[i][j] = newPolynomial;
            matrixToReturn[j][i] = newPolynomial.conj();

        }
    }

    // Return the matrix
    return matrixToReturn;

}

// Generate a list of monomials of a certain level
std::vector<Poly> generateMonomials(std::vector<Mon> variables, int level, int verbosity) {

    // Generate all monomials up to the given level
    std::set<Mon> monomsInTopRow = {};
    if (level >= 1) {
        for (size_t i=0; i<variables.size(); i++) {
            Mon currentMonomial = variables[i];
            //std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
            std::pair<std::complex<double>, Mon> monomCoeff = {0.0, currentMonomial};
            if (!monomsInTopRow.count(monomCoeff.second)) {
                monomsInTopRow.insert(monomCoeff.second);
            }
        }
    }
    if (level >= 2) {
        for (size_t i=0; i<variables.size(); i++) {
            for (size_t j=0; j<variables.size(); j++) {
                Mon currentMonomial = variables[i] * variables[j];
                //std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                std::pair<std::complex<double>, Mon> monomCoeff = {0.0, currentMonomial};
                if (!monomsInTopRow.count(monomCoeff.second)) {
                    monomsInTopRow.insert(monomCoeff.second);
                }
            }
        }
    }
    if (level >= 3) {
        for (size_t i=0; i<variables.size(); i++) {
            for (size_t j=0; j<variables.size(); j++) {
                for (size_t k=0; k<variables.size(); k++) {
                    Mon currentMonomial = variables[i] * variables[j] * variables[k];
                    //std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                    std::pair<std::complex<double>, Mon> monomCoeff = {0.0, currentMonomial};
                    if (!monomsInTopRow.count(monomCoeff.second)) {
                        monomsInTopRow.insert(monomCoeff.second);
                    }
                }
            }
        }
    }
    if (level >= 4) {
        for (size_t i=0; i<variables.size(); i++) {
            for (size_t j=0; j<variables.size(); j++) {
                for (size_t k=0; k<variables.size(); k++) {
                    for (size_t l=0; l<variables.size(); l++) {
                        Mon currentMonomial = variables[i] * variables[j] * variables[k] * variables[l];
                        //std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                        std::pair<std::complex<double>, Mon> monomCoeff = {0.0, currentMonomial};
                        if (!monomsInTopRow.count(monomCoeff.second)) {
                            monomsInTopRow.insert(monomCoeff.second);
                        }
                    }
                }
            }
        }
    }
    if (level >= 5) {
        for (size_t i=0; i<variables.size(); i++) {
            for (size_t j=0; j<variables.size(); j++) {
                for (size_t k=0; k<variables.size(); k++) {
                    for (size_t l=0; l<variables.size(); l++) {
                        for (size_t m=0; m<variables.size(); m++) {
                            Mon currentMonomial = variables[i] * variables[j] * variables[k] * variables[l] * variables[m];
                            //std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                            std::pair<std::complex<double>, Mon> monomCoeff = {0.0, currentMonomial};
                            if (!monomsInTopRow.count(monomCoeff.second)) {
                                monomsInTopRow.insert(monomCoeff.second);
                            }
                        }
                    }
                }
            }
        }
    }
    if (level >= 6) {
        for (size_t i=0; i<variables.size(); i++) {
            for (size_t j=0; j<variables.size(); j++) {
                for (size_t k=0; k<variables.size(); k++) {
                    for (size_t l=0; l<variables.size(); l++) {
                        for (size_t m=0; m<variables.size(); m++) {
                            for (size_t n=0; n<variables.size(); n++) {
                                Mon currentMonomial = variables[i] * variables[j] * variables[k] * variables[l] * variables[m] * variables[n];
                                //std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                                std::pair<std::complex<double>, Mon> monomCoeff = {0.0, currentMonomial};
                                if (!monomsInTopRow.count(monomCoeff.second)) {
                                    monomsInTopRow.insert(monomCoeff.second);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if  (level >= 7) {
        for (size_t i=0; i<variables.size(); i++) {
            for (size_t j=0; j<variables.size(); j++) {
                for (size_t k=0; k<variables.size(); k++) {
                    for (size_t l=0; l<variables.size(); l++) {
                        for (size_t m=0; m<variables.size(); m++) {
                            for (size_t n=0; n<variables.size(); n++) {
                                for (size_t o=0; o<variables.size(); o++) {
                                    Mon currentMonomial = variables[i] * variables[j] * variables[k] * variables[l] * variables[m] * variables[n] * variables[o];
                                    //std::pair<std::complex<double>, Mon> monomCoeff = currentMonomial.reduce();
                                    std::pair<std::complex<double>, Mon> monomCoeff = {0.0, currentMonomial};
                                    if (!monomsInTopRow.count(monomCoeff.second)) {
                                        monomsInTopRow.insert(monomCoeff.second);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Remove the 1 if there is
    std::vector<Poly> asPolys;
    for (auto& monom : monomsInTopRow) {
        std::pair<std::complex<double>, Mon> monomCoeff = monom.reduce();
        if (monomCoeff.second.size() != 0) {
            asPolys.push_back(Poly(monom));
        }
    }

    // Return the monomials (in poly form)
    return asPolys;

}

// Generate all moment matrices for a given level given a polynomial
std::vector<std::vector<std::vector<Poly>>> generateAllMomentMatrices(const Poly& functional, std::vector<Poly> zeroCons, int level, int verbosity) {

    // First get the list of all monomials used by iterating through the polynomial
    std::vector<Mon> variables;
    addSingleMonomials(variables, functional);
    for (size_t i=0; i<zeroCons.size(); i++) {
        addSingleMonomials(variables, zeroCons[i]);
    }

    // Generate all monomials up to the given level
    std::vector<Poly> monomsInTopRow = generateMonomials(variables, level, verbosity);

    // Create the index array to be used for combinations
    std::vector<int> indices;
    for (size_t i=0; i<monomsInTopRow.size(); i++) {
        indices.push_back(i+1);
    }
    std::vector<std::vector<int>> combinations;

    // In this code we just want one big moment matrix
    combinations = {indices};

    // Each moment mat should start with 1
    monomsInTopRow.insert(monomsInTopRow.begin(), Poly(1));
    for (size_t i=0; i<combinations.size(); i++) {
        combinations[i].insert(combinations[i].begin(), 0);
    }

    // Form the moment matrices
    std::vector<std::vector<std::vector<Poly>>> matricesToReturn;
    for (size_t k=0; k<combinations.size(); k++) {

        // Get this top row combination
        std::vector<Poly> monomsInTopRowComb;
        for (size_t i=0; i<combinations[k].size(); i++) {
            monomsInTopRowComb.push_back(monomsInTopRow[combinations[k][i]]);
        }

        // Form that matrix
        std::vector<std::vector<Poly>> newMatrix = generateFromTopRow(monomsInTopRowComb, verbosity);
        matricesToReturn.push_back(newMatrix);

    }

    // Return the moment matrix
    return matricesToReturn;

}

// Add variables from a moment matrix
void addVariables(std::set<Mon>& variables, std::vector<std::vector<Poly>> toAdd) {

    // Iterative through the matrix
    for (size_t i=0; i<toAdd.size(); i++) {
        for (size_t j=0; j<toAdd[i].size(); j++) {

            // Iterate through the polynomial
            for (auto& term : toAdd[i][j]) {
                Mon currentMonomial(term.first);

                // If it's not in the list, add it
                if (!variables.count(currentMonomial)) {
                    variables.insert(currentMonomial);
                }

            }

        }
    }

}

// Add variables from a polynomial
void addVariables(std::set<Mon>& variables, Poly toAdd) {

    // Iterate through the polynomial
    for (auto& term : toAdd) {
        Mon currentMonomial(term.first);

        // If it's not in the list, add it
        if (!variables.count(currentMonomial)) {
            variables.insert(currentMonomial);
        }

    }

}

// Convert from a matrix location to an svec location
int matLocToVecLoc(int i, int j, int n) {
    return i*n + j - i*(i+1)/2;
}

// Take the trace of a matrix, assuming it's real
double tr(Eigen::MatrixXcd A) {
    std::complex<double> trace = A.trace();
    if (std::abs(std::imag(trace)) > 1e-5) {
        std::cout << "WARNING - trace of matrix has non-zero imaginary part" << std::endl;
    }
    return std::real(trace);
}

// Convert a set to a vector
std::vector<Mon> toVector(std::set<Mon> s) {
    std::vector<Mon> v;
    for (auto& elem : s) {
        v.push_back(elem);
    }
    return v;
}

// Generate a random double between min and max
double rand(double min, double max) {
    return min + (max - min) * (double)rand() / RAND_MAX;
}

// Convert a primal SDP problem to a dual problem
void primalToDual(Poly& objective, std::vector<std::vector<std::vector<Poly>>>& momentMatrices, std::vector<Poly>& constraintsZero, std::vector<Poly>& constraintsPositive, bool variableObjective) {

    // First: put into into standard SDP form
    // min C.X s.t. A.X = b, X >= 0
        
    // Define C
    int matSize = momentMatrices[0].size();
    std::vector<std::vector<double>> C = std::vector<std::vector<double>>(matSize, std::vector<double>(matSize, 0));
    std::set<Mon> monsUsed;
    for (int j=0; j<momentMatrices[0].size(); j++) {
        for (int k=0; k<j; k++) {
            Mon key = momentMatrices[0][j][k].getKey();
            if (!monsUsed.count(key)) {
                C[j][k] = -std::real(objective[key]) / 2;
                C[k][j] = C[j][k];
                monsUsed.insert(key);
            }
        }
    }

    // Define the A matrices and b vector
    std::vector<std::vector<std::vector<double>>> As;
    std::vector<double> b;
    for (int i=0; i<momentMatrices[0].size(); i++) {
        for (int j=i; j<momentMatrices[0][i].size(); j++) {

            // Trying to find the matrix relating this to other elements
            std::vector<std::vector<double>> A(matSize, std::vector<double>(matSize, 0));
            if (i == j) { 
                A[i][j] = 1;
            } else {
                A[i][j] = 0.5;
            }
            double bVal = 0;
            bool somethingFound = false;

            // For each term in this poly, try to find an early term 
            // that would reduce this to a constant
            for (auto& term : momentMatrices[0][i][j]) {

                // If it's a constant
                if (term.first.isConstant()) {
                    bVal += std::real(term.second);
                    somethingFound = true;

                // Otherwise find an earlier variable
                } else {
                    Mon monomToFind = term.first;
                    bool found = false;
                    for (int l=0; l<matSize; l++) {
                        for (int m=l; m<matSize; m++) {
                            if (l*matSize + m >= i*matSize + j) {
                                break;
                            }
                            if (momentMatrices[0][l][m].size() == 1 && momentMatrices[0][l][m].getKey() == monomToFind) {
                                if (l == m) { 
                                    A[l][m] = -std::real(term.second);
                                } else {
                                    A[l][m] = -std::real(term.second) / 2;
                                }
                                somethingFound = true;
                                found = true;
                                break;
                            }
                        }
                        if (found) {
                            break;
                        }
                    }

                }

            }

            // Add this matrix
            if (somethingFound) {
                As.push_back(A);
                b.push_back(bVal);
            }

        }
    }

    // Convert the constraints to the dual objective b.y
    Poly newObjective;
    for (int i=0; i<b.size(); i++) {
        if (std::abs(b[i]) > 1e-10) {
            newObjective[Mon("<D" + std::to_string(i) + ">")] = b[i];
        }
    }

    // Convert the objective to the dual constraints C-\sum_i y_i A_i >= 0
    std::vector<std::vector<Poly>> newMomentMat(matSize, std::vector<Poly>(matSize, Poly()));
    int newVarInd = b.size();
    for (int i=0; i<matSize; i++) {
        for (int j=i; j<matSize; j++) {
            if (std::abs(C[i][j]) > 1e-10) {
                // TODO new vars
                if (variableObjective) {
                    Mon newMon("<E" + std::to_string(newVarInd) + ">");
                    newMomentMat[i][j] = Poly(newMon);
                    newVarInd++;
                } else {
                    newMomentMat[i][j][Mon()] = C[i][j];
                }
            }
            for (int k=0; k<As.size(); k++) {
                if (std::abs(As[k][i][j]) > 1e-10) {
                    newMomentMat[i][j][Mon("<D" + std::to_string(k) + ">")] -= As[k][i][j];
                }
            }
        }
    }

    // We don't have any linear constraints
    constraintsZero = {};
    constraintsPositive = {};

    // Symmetrize the matrix
    for (int i=0; i<matSize; i++) {
        for (int j=i+1; j<matSize; j++) {
            newMomentMat[j][i] = newMomentMat[i][j];
        }
    }

    // Change the vars
    momentMatrices = {newMomentMat};
    objective = newObjective;

}


// Convert from sparse to compressed column format
// https://people.sc.fsu.edu/~jburkardt/data/cc/cc.html
void toCC(std::vector<int>& ARows, std::vector<int>& ACols, std::vector<double>& AVals, int numCols) {

    // If empty, do nothing
    if (ARows.size() == 0) {
        return;
    }

    // The arrays to fill
    std::vector<int> colPtrs(numCols+1, AVals.size());
    std::vector<int> rowIndices;
    std::vector<double> values;

    // Sort the vals and rows by columns
    std::vector<std::tuple<int, int, double>> vals;
    for (size_t i=0; i<ARows.size(); i++) {
        vals.push_back(std::make_tuple(ACols[i], ARows[i], AVals[i]));
    }
    std::sort(vals.begin(), vals.end());

    // Fill the arrays
    for (int i=0; i<int(vals.size()); i++) {
        rowIndices.push_back(std::get<1>(vals[i]));
        values.push_back(std::get<2>(vals[i]));
        if (i < colPtrs[std::get<0>(vals[i])]) {
            colPtrs[std::get<0>(vals[i])] = i;
        }
    }

    // Anything not set should be the max
    for (int i=0; i<numCols+1; i++) {
        if (colPtrs[i] == -1) {
            colPtrs[i] = rowIndices.size();
        }
    }

    // Set the output
    ARows = rowIndices;
    ACols = colPtrs;
    AVals = values;

}

// Solve a linear system with Eigen
double solveEigen(Poly& objective, std::vector<Poly>& constraintsZero, int verbosity, int numCores) {

    // Get the list of variables
    std::set<Mon> variableSet;
    variableSet.insert(Mon());
    for (size_t i=0; i<constraintsZero.size(); i++) {
        addVariables(variableSet, constraintsZero[i]);
    }
    addVariables(variableSet, objective);
    variableSet.erase(Mon());
    std::vector<Mon> variables = toVector(variableSet);

    // Cache the variable locations
    std::map<Mon, int> variableLocs;
    for (size_t i=0; i<variables.size(); i++) {
        variableLocs[variables[i]] = i;
    }

    // Set up a sparse Eigen matrix
    Eigen::SparseMatrix<double,Eigen::RowMajor> A(constraintsZero.size(), variables.size());
    Eigen::VectorXd b = Eigen::VectorXd::Zero(constraintsZero.size());
    std::vector<Eigen::Triplet<double>> coefficients;
    for (int i=0; i<constraintsZero.size(); i++) {
        for (auto& term : constraintsZero[i]) {
            if (term.first.isConstant() && std::abs(term.second) > 1e-14) {
                b[i] = -std::real(term.second);
            } else {
                coefficients.push_back(Eigen::Triplet<double>(i, variableLocs[term.first], std::real(term.second)));
            }
        }
    }
    A.setFromTriplets(coefficients.begin(), coefficients.end());
    A.makeCompressed();

    // Output the matrix if verbose
    if (verbosity >= 3) {
        Eigen::MatrixXd denseA = A;
        std::cout << "A:" << std::endl;
        std::cout << denseA << std::endl;
        std::cout << "b:" << std::endl;
        std::cout << b << std::endl;
    }

    // Solve the system
    //Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    //std::cout << "Analyzing pattern..." << std::endl;
    //solver.analyzePattern(A); 
    //std::cout << "Factoring matrix..." << std::endl;
    //solver.factorize(A); 
    //if (solver.info() != Eigen::Success) {
        //std::cout << "ERROR - Failed to factorize matrix" << std::endl;
    //}
    //std::cout << "Solving..." << std::endl;
    //Eigen::VectorXd x = solver.solve(b);

    // Solve the system
    //Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    //std::cout << "Analyzing pattern..." << std::endl;
    //solver.analyzePattern(A); 
    //std::cout << "Factoring matrix..." << std::endl;
    //solver.factorize(A); 
    //if (solver.info() != Eigen::Success) {
        //std::cout << "ERROR - Failed to factorize matrix" << std::endl;
    //}
    //std::cout << "Solving..." << std::endl;
    //Eigen::VectorXd x = solver.solve(b);
    
    omp_set_num_threads(numCores);
    Eigen::setNbThreads(numCores);

    // Solve the system
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
    //solver.setMaxIterations(1000);
    //solver.compute(A);
    //Eigen::VectorXd x = solver.solve(b);
    //std::cout << "Num iterations:  " << solver.iterations() << std::endl;
    //std::cout << "Estimated error: " << solver.error()      << std::endl;

    // Solve the system
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd x = solver.solve(b);
    std::cout << "Num iterations:  " << solver.iterations() << std::endl;
    std::cout << "Estimated error: " << solver.error()      << std::endl;

    // Compute the objective
    std::cout << "Computing objective..." << std::endl;
    double objVal = 0;
    for (auto& term : objective) {
        objVal += std::real(term.second) * x[variableLocs[term.first]];
    }

    return objVal;

}



