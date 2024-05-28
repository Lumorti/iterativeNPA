#include "mosek.h"
#include "utils.h"
#include "printing.h"
#include <iostream>

// Import MOSEK
#include "fusion.h"

// Convert to MOSEK form and solve
double maximizeMOSEK(Poly obj, std::vector<std::vector<std::vector<Poly>>>& psd, std::vector<Poly> constraintsZero, std::vector<Poly> constraintsPositive, int verbosity, std::map<Mon, std::complex<double>>* xMap) {

    // Get the list of variables
    int oneIndex = 0;
    std::set<Mon> variableSet;
    variableSet.insert(Mon());
    for (size_t i=0; i<psd.size(); i++) {
        addVariables(variableSet, psd[i]);
    }
    for (size_t i=0; i<constraintsZero.size(); i++) {
        addVariables(variableSet, constraintsZero[i]);
    }
    for (size_t i=0; i<constraintsPositive.size(); i++) {
        addVariables(variableSet, constraintsPositive[i]);
    }
    addVariables(variableSet, obj);
    std::vector<Mon> variables = toVector(variableSet);

    // Output the variable list
    if (verbosity >= 2) {
        std::cout << "Num variables: " << variables.size() << std::endl;
    }
    if (verbosity >= 3) {
        std::cout << "Variables:" << std::endl;
        for (size_t i=0; i<variables.size(); i++) {
            std::cout << i << " " << variables[i] << std::endl;
        }
        std::cout << std::endl;
    }

    // Cache the variable locations
    std::map<Mon, int> variableLocs;
    for (size_t i=0; i<variables.size(); i++) {
        variableLocs[variables[i]] = i;
    }

    // The c vector defining the objective
    std::vector<double> c(variables.size());
    for (auto& term : obj.polynomial) {

        // Find the location of this variable
        int varLoc = variableLocs[term.first];

        // Add the real and imaginary terms
        c[varLoc] += std::real(term.second);

    }

    // The A matrix defining the equality constraints
    // (c_r + i c_i)*x_r 
    std::vector<int> ARows;
    std::vector<int> ACols;
    std::vector<double> AVals;
    int numCons = 0;
    for (size_t i=0; i<constraintsZero.size(); i++) {

        // For each term in the constraint
        bool realConNonZero = false;
        for (auto& term : constraintsZero[i]) {

            // Find the location of this variable
            int realLoc = variableLocs[term.first];
            double realVal = std::real(term.second);

            // c_r*x_r
            if (std::abs(realVal) > 1e-14) {
                ARows.push_back(numCons);
                ACols.push_back(realLoc);
                AVals.push_back(realVal);
                realConNonZero = true;
            }

        }
        if (realConNonZero) {
            numCons++;
        }
    }

    // Verbose output of the A matrix
    if (verbosity >= 3) {
        std::cout << "ARows: " << ARows << std::endl;
        std::cout << "ACols: " << ACols << std::endl;
        std::cout << "AVals: " << AVals << std::endl;
    }

    // The B matrix defining the positivity constraints
    std::vector<int> BRows;
    std::vector<int> BCols;
    std::vector<double> BVals;
    int numPosCons = 0;
    for (size_t i=0; i<constraintsPositive.size(); i++) {

        // For each term in the constraint
        bool realConNonZero = false;
        for (auto& term : constraintsPositive[i]) {

            // Find the location of this variable
            int realLoc = variableLocs[term.first];
            double realVal = std::real(term.second);

            // c_r*x_r
            if (std::abs(realVal) > 1e-14) {
                BRows.push_back(numPosCons);
                BCols.push_back(realLoc);
                BVals.push_back(realVal);
                realConNonZero = true;
            }

        }
        if (realConNonZero) {
            numPosCons++;
        }

    }

    // The vectors defining the PSD constraints
    std::vector<std::shared_ptr<monty::ndarray<int,1>>> indicesPSDPerMat;
    std::vector<std::shared_ptr<monty::ndarray<double,1>>> coeffsPSDPerMat;
    std::vector<std::pair<int,int>> matDims;
    for (size_t k=0; k<psd.size(); k++) {

        // If the matrix is just one
        if (psd[k].size() == 1) {
            continue;
        }

        // The indices and coefficients for the svec
        int fullMatSize = psd[k].size();
        int sVecSize = fullMatSize * (fullMatSize + 1) / 2;
        int imagOffset = psd[k].size();
        std::vector<int> indicesPSD(sVecSize);
        std::vector<double> coeffsPSD(sVecSize);
        for (size_t i=0; i<psd[k].size(); i++) {
            for (size_t j=i; j<psd[k][i].size(); j++) {

                // Find this in the variable list
                int realLoc = variableLocs[psd[k][i][j].getKey()];

                // Locations in the svec for the real and imaginary parts
                int realInd = matLocToVecLoc(i, j, fullMatSize);

                // If it's a real coeff
                std::complex<double> val = psd[k][i][j].getValue();
                double coeff = std::real(val);
                indicesPSD[realInd] = realLoc;
                coeffsPSD[realInd] = coeff;

                // The off-diagonals are multiplied by sqrt(2)
                if (i != j) {
                    coeffsPSD[realInd] *= std::sqrt(2.0);
                }

            }

        }

        // Construct the dense matrix for debugging
        if (verbosity >= 3) {
            Eigen::MatrixXi indMat = Eigen::MatrixXi::Zero(fullMatSize, fullMatSize);
            Eigen::MatrixXd coeffMat = Eigen::MatrixXd::Zero(fullMatSize, fullMatSize);
            int nextX = 0;
            int nextY = 0;
            for (int i=0; i<sVecSize; i++) {
                indMat(nextX, nextY) = indicesPSD[i];
                coeffMat(nextX, nextY) = coeffsPSD[i];
                nextY++;
                if (nextY == fullMatSize) {
                    nextX++;
                    nextY = nextX;
                }
            }
            std::cout << "Indices matrix: " << std::endl;
            std::cout << indMat << std::endl;
            std::cout << "Coeffs matrix: " << std::endl;
            std::cout << coeffMat << std::endl;
        }

        // Convert to MOSEK form
        auto indicesPSDM = monty::new_array_ptr<int>(indicesPSD);
        auto coeffsPSDM = monty::new_array_ptr<double>(coeffsPSD);

        // Add to the list
        indicesPSDPerMat.push_back(indicesPSDM);
        coeffsPSDPerMat.push_back(coeffsPSDM);
        matDims.push_back({1, sVecSize});

    }

    // Convert to MOSEK form
    auto cM = monty::new_array_ptr<double>(c);
    auto AM = mosek::fusion::Matrix::sparse(numCons, variables.size(), monty::new_array_ptr<int>(ARows), monty::new_array_ptr<int>(ACols), monty::new_array_ptr<double>(AVals));
    auto BM = mosek::fusion::Matrix::sparse(numPosCons, variables.size(), monty::new_array_ptr<int>(BRows), monty::new_array_ptr<int>(BCols), monty::new_array_ptr<double>(BVals));

    // Create a model
    mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});
    if (verbosity >= 3) {
        M->setLogHandler([=](const std::string & msg) {std::cout << msg << std::flush;});
    }

    // Create the main variable vector
    mosek::fusion::Variable::t xM = M->variable(variables.size(), mosek::fusion::Domain::inRange(-1.0, 1.0));

    // The objective function
    M->objective(mosek::fusion::ObjectiveSense::Maximize, mosek::fusion::Expr::dot(cM, xM));

    // The one variable should be fixed
    M->constraint(xM->index(oneIndex), mosek::fusion::Domain::equalsTo(1.0));

    // The matrix of this should be PSD
    for (size_t k=0; k<indicesPSDPerMat.size(); k++) {
        M->constraint(
            mosek::fusion::Expr::sum(
                mosek::fusion::Expr::reshape(
                    mosek::fusion::Expr::mulElm(
                        coeffsPSDPerMat[k], 
                        xM->pick(indicesPSDPerMat[k])
                    ),
                    matDims[k].first,
                    matDims[k].second
                ),
                0
            ), 
            mosek::fusion::Domain::inSVecPSDCone()
        );
    }

    // Linear equality constraints
    if (numCons > 0) {
        M->constraint(mosek::fusion::Expr::mul(AM, xM), mosek::fusion::Domain::equalsTo(0.0));
    }

    // Linear positivity constraints
    if (numPosCons > 0) {
        M->constraint(mosek::fusion::Expr::mul(BM, xM), mosek::fusion::Domain::greaterThan(0.0));
    }

    // Solve the problem
    M->solve();

    // Output the primal objective value
    double objPrimal = M->primalObjValue();
    double upperBound = M->dualObjValue();

    // Get all of the variable values
    auto xMLevel = *(xM->level());
    std::vector<std::complex<double>> variableValues(variables.size(), 0);
    for (size_t i=0; i<variables.size(); i++) {
        variableValues[i] = std::complex<double>(xMLevel[i], 0);
    }

    // Check the eigenvalues of each moment matrix
    if (verbosity >= 2) {
        std::cout << std::endl;
        std::cout << "Objective value: " << objPrimal << std::endl;
        for (size_t i=0; i<psd.size(); i++) {
            std::vector<std::vector<std::complex<double>>> eigenvectors;
            std::vector<std::complex<double>> eigenvalues;
            getEigens(psd[i], variables, variableValues, eigenvectors, eigenvalues);
            double minEig = 1e10;
            for (size_t j=0; j<eigenvalues.size(); j++) {
                if (std::imag(eigenvalues[j]) > 1e-5) {
                    std::cout << "ERROR - Eigenvalue " << j << " has imaginary part: " << std::imag(eigenvalues[j]) << std::endl;
                }
                if (std::real(eigenvalues[j]) < minEig) {
                    minEig = std::real(eigenvalues[j]);
                }
            }
            std::cout << "Min eigenvalue: " << minEig << std::endl;
        }
    }

    // If xMap is not a null pointer
    if (xMap) {
        for (size_t i=0; i<variables.size(); i++) {
            (*xMap)[variables[i]] = variableValues[i];
        }
    }

    // If superverbose, output all monomials
    if (verbosity >= 3) {
        std::cout << "Solution: " << std::endl;
        for (size_t i=0; i<variables.size(); i++) {
            std::cout << variables[i] << ": " << variableValues[i] << std::endl;
        }
        std::cout << "There are " << variables.size() << " variables." << std::endl;

        // Some matrix stuff
        if (psd.size() > 0) {
            Eigen::MatrixXcd A = replaceVariables(psd[0], variables, variableValues);
            std::cout << "Moment matrix:" << std::endl;
            std::cout << psd[0] << std::endl;
            std::cout << "Moment matrix with vars replaced:" << std::endl;
            std::cout << A << std::endl;
            std::vector<std::complex<double>> eigVals;
            std::vector<std::vector<std::complex<double>>> eigVecs;
            getEigens(psd[0], variables, variableValues, eigVecs, eigVals);
            std::cout << "Eigenvalues: " << eigVals << std::endl;
            int numZeroEig = 0;
            for (size_t i=0; i<eigVals.size(); i++) {
                if (std::abs(eigVals[i]) < 1e-8) {
                    numZeroEig++;
                }
            }
            int numPosEig = eigVals.size() - numZeroEig;
            std::cout << "Num zero eig: " << numZeroEig << std::endl;
            std::cout << "Num pos eig: " << numPosEig << std::endl;
        }

        // Eval all of the constraints
        std::map<Mon, std::complex<double>> vals;
        for (size_t i=0; i<variables.size(); i++) {
            vals[variables[i]] = variableValues[i];
        }
        for (size_t i=0; i<constraintsZero.size(); i++) {
            std::cout << "Constraint " << i << ": " << constraintsZero[i].eval(vals) << std::endl;
        }
        for (size_t i=0; i<constraintsPositive.size(); i++) {
            std::cout << "Pos Constraint " << i << ": " << constraintsPositive[i].eval(vals) << std::endl;
        }

    // If verbose, just output monomials that are in the objective
    } else if (verbosity >= 2) {
        std::cout << "Solution: " << std::endl;
        std::set<Mon> variablesInObjSet = {};
        addVariables(variablesInObjSet, obj);
        std::vector<Mon> variablesInObj = toVector(variablesInObjSet);
        for (size_t i=0; i<variablesInObj.size(); i++) {
            for (size_t j=0; j<variables.size(); j++) {
                if (variablesInObj[i] == variables[j]) {
                    std::cout << variablesInObj[i] << ": " << variableValues[j] << std::endl;
                    break;
                }
            }
        }
    } 

    // Return the objective
    return upperBound;

}

// Convert to MOSEK form and solve
std::pair<double,double> solveMOSEK(Poly obj, std::vector<std::vector<std::vector<Poly>>>& psd, std::vector<Poly> constraintsZero, std::vector<Poly> constraintsPositive, int verbosity) {

    // Get the list of variables
    int oneIndex = 0;
    std::set<Mon> variableSet;
    variableSet.insert(Mon());
    for (size_t i=0; i<psd.size(); i++) {
        addVariables(variableSet, psd[i]);
    }
    for (size_t i=0; i<constraintsZero.size(); i++) {
        addVariables(variableSet, constraintsZero[i]);
    }
    for (size_t i=0; i<constraintsPositive.size(); i++) {
        addVariables(variableSet, constraintsPositive[i]);
    }
    addVariables(variableSet, obj);
    std::vector<Mon> variables = toVector(variableSet);

    // Output the variable list
    if (verbosity >= 2) {
        std::cout << "Num variables: " << variables.size() << std::endl;
    }
    if (verbosity >= 3) {
        std::cout << "Variables:" << std::endl;
        for (size_t i=0; i<variables.size(); i++) {
            std::cout << i << " " << variables[i] << std::endl;
        }
        std::cout << std::endl;
    }

    // Cache the variable locations
    std::map<Mon, int> variableLocs;
    for (size_t i=0; i<variables.size(); i++) {
        variableLocs[variables[i]] = i;
    }

    // The c vector defining the objective
    std::vector<double> c(variables.size());
    for (auto& term : obj.polynomial) {

        // Find the location of this variable
        int varLoc = variableLocs[term.first];

        // Add the real and imaginary terms
        c[varLoc] += std::real(term.second);

    }

    // The A matrix defining the equality constraints
    // (c_r + i c_i)*x_r 
    std::vector<int> ARows;
    std::vector<int> ACols;
    std::vector<double> AVals;
    int numCons = 0;
    for (size_t i=0; i<constraintsZero.size(); i++) {

        // For each term in the constraint
        bool realConNonZero = false;
        bool imagConNonZero = false;
        for (auto& term : constraintsZero[i]) {

            // Find the location of this variable
            int realLoc = variableLocs[term.first];
            double realVal = std::real(term.second);
            double imagVal = std::imag(term.second);

            // c_r*x_r
            // c_i*x_r
            if (std::abs(realVal) > 1e-14) {
                ARows.push_back(numCons);
                ACols.push_back(realLoc);
                AVals.push_back(realVal);
                realConNonZero = true;
            }
            if (std::abs(imagVal) > 1e-14) {
                ARows.push_back(numCons+1);
                ACols.push_back(realLoc);
                AVals.push_back(imagVal);
                imagConNonZero = true;
            }

        }
        if (realConNonZero) {
            numCons++;
        }
        if (imagConNonZero) {
            numCons++;
        }
    }

    // Verbose output of the A matrix
    if (verbosity >= 3) {
        std::cout << "ARows: " << ARows << std::endl;
        std::cout << "ACols: " << ACols << std::endl;
        std::cout << "AVals: " << AVals << std::endl;
    }

    // The B matrix defining the positivity constraints
    // x_r^2 + x_i^2
    std::vector<int> BRows;
    std::vector<int> BCols;
    std::vector<double> BVals;
    int numPosCons = 0;
    for (size_t i=0; i<constraintsPositive.size(); i++) {

        // For each term in the constraint
        bool realConNonZero = false;
        bool imagConNonZero = false;
        for (auto& term : constraintsPositive[i]) {

            // Find the location of this variable
            int realLoc = variableLocs[term.first];
            double realVal = std::real(term.second);
            double imagVal = std::imag(term.second);

            // c_r*x_r
            // c_i*x_r
            if (std::abs(realVal) > 1e-14) {
                BRows.push_back(numPosCons);
                BCols.push_back(realLoc);
                BVals.push_back(realVal);
                realConNonZero = true;
            }
            if (std::abs(imagVal) > 1e-14) {
                BRows.push_back(numPosCons+1);
                BCols.push_back(realLoc);
                BVals.push_back(imagVal);
                imagConNonZero = true;
            }

        }
        if (realConNonZero) {
            numPosCons++;
        }
        if (imagConNonZero) {
            numPosCons++;
        }
    }

    // The vectors defining the PSD constraints
    std::vector<std::shared_ptr<monty::ndarray<int,1>>> indicesPSDPerMat;
    std::vector<std::shared_ptr<monty::ndarray<double,1>>> coeffsPSDPerMat;
    std::vector<std::pair<int,int>> matDims;
    for (size_t k=0; k<psd.size(); k++) {

        // If the matrix is just one
        if (psd[k].size() == 1) {
            continue;
        }

        // The indices and coefficients for the svec
        int fullMatSize = 2*psd[k].size();
        int sVecSize = fullMatSize * (fullMatSize + 1) / 2;
        int imagOffset = psd[k].size();
        std::vector<int> indicesPSD(sVecSize);
        std::vector<double> coeffsPSD(sVecSize);
        for (size_t i=0; i<psd[k].size(); i++) {
            for (size_t j=i; j<psd[k][i].size(); j++) {

                // Find this in the variable list
                int realLoc = variableLocs[psd[k][i][j].getKey()];

                // Locations in the svec for the real and imaginary parts
                int realInd = matLocToVecLoc(i, j, fullMatSize);
                int realInd2 = matLocToVecLoc(i+imagOffset, j+imagOffset, fullMatSize);
                int imagInd = matLocToVecLoc(i, j+imagOffset, fullMatSize);
                int imagInd2 = matLocToVecLoc(j, i+imagOffset, fullMatSize);

                // If it's a real coeff
                std::complex<double> val = psd[k][i][j].getValue();
                if (std::imag(val) == 0) {
                    double coeff = std::real(val);
                    indicesPSD[realInd] = realLoc;
                    coeffsPSD[realInd] = coeff;
                    indicesPSD[realInd2] = realLoc;
                    coeffsPSD[realInd2] = coeff;

                // If it's an imaginary coeff
                } else {
                    double coeff = std::imag(val);
                    indicesPSD[imagInd] = realLoc;
                    coeffsPSD[imagInd] = -coeff;
                    indicesPSD[imagInd2] = realLoc;
                    coeffsPSD[imagInd2] = coeff;

                }

                // The off-diagonals are multiplied by sqrt(2)
                if (i != j) {
                    coeffsPSD[realInd] *= std::sqrt(2.0);
                    coeffsPSD[realInd2] *= std::sqrt(2.0);
                    coeffsPSD[imagInd] *= std::sqrt(2.0);
                    coeffsPSD[imagInd2] *= std::sqrt(2.0);
                }


            }

        }

        // Construct the dense matrix for debugging
        if (verbosity >= 3) {
            Eigen::MatrixXi indMat = Eigen::MatrixXi::Zero(fullMatSize, fullMatSize);
            Eigen::MatrixXd coeffMat = Eigen::MatrixXd::Zero(fullMatSize, fullMatSize);
            int nextX = 0;
            int nextY = 0;
            for (int i=0; i<sVecSize; i++) {
                indMat(nextX, nextY) = indicesPSD[i];
                coeffMat(nextX, nextY) = coeffsPSD[i];
                nextY++;
                if (nextY == fullMatSize) {
                    nextX++;
                    nextY = nextX;
                }
            }
            std::cout << "Indices matrix: " << std::endl;
            std::cout << indMat << std::endl;
            std::cout << "Coeffs matrix: " << std::endl;
            std::cout << coeffMat << std::endl;
        }

        // Convert to MOSEK form
        auto indicesPSDM = monty::new_array_ptr<int>(indicesPSD);
        auto coeffsPSDM = monty::new_array_ptr<double>(coeffsPSD);

        // Add to the list
        indicesPSDPerMat.push_back(indicesPSDM);
        coeffsPSDPerMat.push_back(coeffsPSDM);
        matDims.push_back({1, sVecSize});

    }

    // Convert to MOSEK form
    auto cM = monty::new_array_ptr<double>(c);
    auto AM = mosek::fusion::Matrix::sparse(numCons, variables.size(), monty::new_array_ptr<int>(ARows), monty::new_array_ptr<int>(ACols), monty::new_array_ptr<double>(AVals));
    auto BM = mosek::fusion::Matrix::sparse(numPosCons, variables.size(), monty::new_array_ptr<int>(BRows), monty::new_array_ptr<int>(BCols), monty::new_array_ptr<double>(BVals));

    // Create a model
    mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});
    if (verbosity >= 3) {
        M->setLogHandler([=](const std::string & msg) {std::cout << msg << std::flush;});
    }

    // Create the main variable vector
    mosek::fusion::Variable::t xM = M->variable(variables.size(), mosek::fusion::Domain::inRange(-1.0, 1.0));

    // The objective function
    M->objective(mosek::fusion::ObjectiveSense::Maximize, mosek::fusion::Expr::dot(cM, xM));

    // The one variable should be fixed
    M->constraint(xM->index(oneIndex), mosek::fusion::Domain::equalsTo(1.0));

    // The matrix of this should be PSD
    for (size_t k=0; k<indicesPSDPerMat.size(); k++) {
        M->constraint(
            mosek::fusion::Expr::sum(
                mosek::fusion::Expr::reshape(
                    mosek::fusion::Expr::mulElm(
                        coeffsPSDPerMat[k], 
                        xM->pick(indicesPSDPerMat[k])
                    ),
                    matDims[k].first,
                    matDims[k].second
                ),
                0
            ), 
            mosek::fusion::Domain::inSVecPSDCone()
        );
    }

    // Linear equality constraints
    M->constraint(mosek::fusion::Expr::mul(AM, xM), mosek::fusion::Domain::equalsTo(0.0));

    // Linear positivity constraints
    M->constraint(mosek::fusion::Expr::mul(BM, xM), mosek::fusion::Domain::greaterThan(0.0));

    // Solve the problem
    M->solve();

    // Output the primal objective value
    double objPrimal = M->primalObjValue();
    double upperBound = M->dualObjValue();

    // Get all of the variable values
    auto xMLevel = *(xM->level());
    std::vector<std::complex<double>> variableValues(variables.size(), 0);
    for (size_t i=0; i<variables.size(); i++) {
        variableValues[i] = std::complex<double>(xMLevel[i], 0);
    }

    // Check the eigenvalues of each moment matrix
    if (verbosity >= 2) {
        std::cout << std::endl;
        std::cout << "Objective value: " << objPrimal << std::endl;
        for (size_t i=0; i<psd.size(); i++) {
            std::vector<std::vector<std::complex<double>>> eigenvectors;
            std::vector<std::complex<double>> eigenvalues;
            getEigens(psd[i], variables, variableValues, eigenvectors, eigenvalues);
            double minEig = 1e10;
            for (size_t j=0; j<eigenvalues.size(); j++) {
                if (std::imag(eigenvalues[j]) > 1e-5) {
                    std::cout << "ERROR - Eigenvalue " << j << " has imaginary part: " << std::imag(eigenvalues[j]) << std::endl;
                }
                if (std::real(eigenvalues[j]) < minEig) {
                    minEig = std::real(eigenvalues[j]);
                }
            }
            std::cout << "Min eigenvalue: " << minEig << std::endl;
        }
    }

    // If superverbose, output all monomials
    if (verbosity >= 3) {
        std::cout << "Solution: " << std::endl;
        for (size_t i=0; i<variables.size(); i++) {
            std::cout << variables[i] << ": " << variableValues[i] << std::endl;
        }
        std::cout << "There are " << variables.size() << " variables." << std::endl;
        Eigen::MatrixXcd A = replaceVariables(psd[0], variables, variableValues);
        std::cout << "Moment matrix:" << std::endl;
        std::cout << psd[0] << std::endl;
        std::cout << "Moment matrix with vars replaced:" << std::endl;
        std::cout << A << std::endl;

        // Eval all of the constraints
        std::map<Mon, std::complex<double>> vals;
        for (size_t i=0; i<variables.size(); i++) {
            vals[variables[i]] = variableValues[i];
        }
        for (size_t i=0; i<constraintsZero.size(); i++) {
            std::cout << "Constraint " << i << ": " << constraintsZero[i].eval(vals) << std::endl;
        }

    // If verbose, just output monomials that are in the objective
    } else if (verbosity >= 2) {
        std::cout << "Solution: " << std::endl;
        std::set<Mon> variablesInObjSet = {};
        addVariables(variablesInObjSet, obj);
        std::vector<Mon> variablesInObj = toVector(variablesInObjSet);
        for (size_t i=0; i<variablesInObj.size(); i++) {
            for (size_t j=0; j<variables.size(); j++) {
                if (variablesInObj[i] == variables[j]) {
                    std::cout << variablesInObj[i] << ": " << variableValues[j] << std::endl;
                    break;
                }
            }
        }
    } 

    M->objective(mosek::fusion::ObjectiveSense::Minimize, mosek::fusion::Expr::dot(cM, xM));
    M->solve();
    double objPrimal2 = M->primalObjValue();
    double lowerBound = M->dualObjValue();

    // Return the objective
    return {lowerBound, upperBound};

}

