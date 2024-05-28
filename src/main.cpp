// Standard includes
#include <iostream>
#include <vector>

// Import MOSEK
#include "fusion.h"

// Import Eigen
#include <Eigen/Dense>

// Optim
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"

// Local libs
#include "poly.h"
#include "utils.h"
#include "mon.h"
#include "mosek.h"
#include "optim.h"

// Generic entry function
int main(int argc, char* argv[]) {

    // Define the Bell scenario
    int numOutcomesA = 2;
    int numOutcomesB = 2;
    int level = 1;
    int subLevel = -1;
    bool useDual = false;
    Poly bellFunc("<A1B1>+<A1B2>+<A2B1>-<A2B2>");
    int testing = 0;
    int numToSample = 1000;
    int verbosity = 1;
    int maxIters = 1000;
    double startingPenalty = 1e5;
    int numPenalties = 5;
    bool use01 = false;
    std::string seed = "";
    std::string problemName = "CHSH";
    std::vector<std::string> extraMoments;

    // Process command-line args
    for (int i=1; i<argc; i++) {
        std::string argAsString = std::string(argv[i]);

        // Setting number of outcomes
        if (argAsString == "-a") {
            numOutcomesA = std::stoi(argv[i+1]);
            i++;
        } else if (argAsString == "-b") {
            numOutcomesB = std::stoi(argv[i+1]);
            i++;

        // CHSH (c = 2, q = 2sqrt(2) w/ level 1)
        } else if (argAsString == "--chsh" || argAsString == "--CHSH") {
            bellFunc = Poly("<A1B1>+<A1B2>+<A2B1>-<A2B2>");
            problemName = "CHSH";

        // If told to use a 0/1 output
        } else if (argAsString == "-01") {
            use01 = true;

        // I3322
        // l1 (7x7) = 5.5
        // l2 (28x28) = 5.00612
        // l3 (88x88) = 5.0035
        // l4 (244x244) = 5.0035
        } else if (argAsString == "--I3322" || argAsString == "--i3322") {
            if (use01) {
                bellFunc = Poly("-<A2>-<B1>-2<B2>+<A1B1>+<A1B2>+<A2B1>+<A2B2>-<A1B3>+<A2B3>-<A3B1>+<A3B2>");
            } else {
                bellFunc = Poly("<A1>-<A2>+<B1>-<B2>-<A1B1>+<A1B2>+<A2B1>-<A2B2>+<A1B3>+<A2B3>+<A3B1>+<A3B2>");
            }
            problemName = "I3322";

        // Set the level
        } else if (argAsString == "-l") {
            level = std::stoi(argv[i+1]);
            i++;

        // If adding an extra moment to the top row
        } else if (argAsString == "-e") {
            extraMoments.push_back(std::string(argv[i+1]));
            i++;

        // Set the iteration limit
        } else if (argAsString == "-i") {
            maxIters = std::stoi(argv[i+1]);
            i++;

        // Set the seed
        } else if (argAsString == "-S") {
            seed = std::string(argv[i+1]);
            i++;

        // Set the sub level
        } else if (argAsString == "-s") {
            subLevel = std::stoi(argv[i+1]);
            i++;

        // Set the number to sample
        } else if (argAsString == "-n") {
            numToSample = std::stoi(argv[i+1]);
            i++;

        // If we're testing
        } else if (argAsString == "-t") {
            testing = std::stoi(argv[i+1]);
			i++;

        // If setting verbosity
        } else if (argAsString == "-v") {
            verbosity = std::stoi(argv[i+1]);
            i++;

        // Output the help
        } else if (argAsString == "-h" || argAsString == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -a <num>        Number of outcomes for Alice" << std::endl;
            std::cout << "  -b <num>        Number of outcomes for Bob" << std::endl;
            std::cout << "  --chsh          Use the CHSH scenario" << std::endl;
            std::cout << "  --I3322         Use the I3322 scenario" << std::endl;
            std::cout << "  -l <num>        Level of the moment matrix" << std::endl;
            std::cout << "  -s <num>        Sub-level of the moment matrix" << std::endl;
            std::cout << "  -i <num>        Iteration limit" << std::endl;
            std::cout << "  -n <num>        Number of moment matrices to sample" << std::endl;
            std::cout << "  -S <str>        Seed for the random number generator" << std::endl;
            std::cout << "  -v <num>        Set the verbosity level" << std::endl;
            std::cout << "  -p <num>        Set the starting penalty" << std::endl;
            std::cout << "  -P <num>        Set the number of penalties to use" << std::endl;
            std::cout << "  -t              Test the moment matrix" << std::endl;
            std::cout << "  -e <str>        Add an extra moment to the top row" << std::endl;
            std::cout << "  -h              Display this help message" << std::endl;
            std::cout << "  -D              Use the dual of the problem" << std::endl;
            std::cout << "  -01             Use 0/1 output instead of -1/1" << std::endl;
            return 0;

        // If using the dual
        } else if (argAsString == "-D") {
            useDual = true;

        // Set the starting penalty
        } else if (argAsString == "-p") {
            startingPenalty = std::stod(argv[i+1]);
            i++;

        // Set the number of penalty adjustments
        } else if (argAsString == "-P") {
            numPenalties = std::stoi(argv[i+1]);
            i++;

        // Otherwise we don't know what this is
        } else {
            std::cout << "Unknown argument: " << argAsString << std::endl;
            return 1;

        }
    }

    // If the seed isn't set
    if (seed == "") {
        srand(time(NULL));
    } else {
        srand(std::stoi(seed));
    }

    // Define the moment matrix
    std::vector<std::vector<std::vector<Poly>>> momentMatrices = generateAllMomentMatrices(bellFunc, {}, level, verbosity);
    if (verbosity >= 1) {
        int largestMomentMatrix = 0;
        for (int i=0; i<momentMatrices.size(); i++) {
            if (momentMatrices[i].size() > largestMomentMatrix) {
                largestMomentMatrix = momentMatrices[i].size();
            }
        }
        std::cout << "Largest moment matrix has size " << largestMomentMatrix << std::endl;
    }

    // If told to add extra to the top row
    if (extraMoments.size() > 0) {
        std::vector<Poly> topRow = momentMatrices[0][0];
        for (int i=0; i<extraMoments.size(); i++) {
            topRow.push_back(Poly(extraMoments[i]));
        }
        momentMatrices[0] = generateFromTopRow(topRow, use01);
    }

    // If using zero-one output, our diagonals aren't 1
    if (use01) {
        for (int i=0; i<momentMatrices.size(); i++) {
            for (int j=0; j<momentMatrices[i].size(); j++) {
                momentMatrices[i][j][j] = momentMatrices[i][j][0];
            }
        }
    }

    // Define other constraints
    Poly objective = bellFunc;
    std::vector<Poly> constraintsZero;
    std::vector<Poly> constraintsPositive;

    // The idea of using the dual to find cons that aren't in level 1 TODO
    //if (testing == 1) {

        //std::cout << "Primal moment matrix: " << std::endl;
        //std::cout << momentMatrices[0] << std::endl << std::endl;
        //std::cout << "Primal objective: " << std::endl;
        //std::cout << objective << std::endl << std::endl;

        //// Keep iterating, for now just a fixed number of times
        //for (int iter=0; iter<maxIters; iter++) {

            //// Get the L1 solution
            //std::vector<Mon> varNamesL1;
            //std::vector<double> varValsL1;
            //std::vector<std::vector<std::vector<Poly>>> momentMatricesL1 = generateAllMomentMatrices(bellFunc, 1, subLevel, numToSample, use01);
            //std::vector<std::vector<std::vector<Poly>>> momentMatricesL2 = generateAllMomentMatrices(bellFunc, 2, subLevel, numToSample, use01);
            //double resL1 = maximizeMOSEK(objective, momentMatricesL1, constraintsZero, constraintsPositive, verbosity, varNamesL1, varValsL1, use01);
            //Eigen::MatrixXd X = replaceVariables(momentMatricesL1[0], varNamesL1, varValsL1);
            //std::cout << "L1 result: " << resL1 << std::endl;
            //std::cout << "L1 solution: " << std::endl;
            //std::cout << X << std::endl << std::endl;

            //// Then find the dual for the L2
            //momentMatrices = momentMatricesL2;

            //// Define C
            //int matSize = momentMatrices[0].size();
            //std::vector<std::vector<double>> C = std::vector<std::vector<double>>(matSize, std::vector<double>(matSize, 0));
            //for (int i=0; i<objective.size(); i++) {
                //monomial monToFind = objective[i].second;
                //int locX = -1;
                //int locY = -1;
                //for (int j=0; j<momentMatrices[0].size(); j++) {
                    //for (int k=j+1; k<momentMatrices[0][j].size(); k++) {
                        //if (momentMatrices[0][j][k].size() == 1 && momentMatrices[0][j][k][0].second == monToFind) {
                            //locX = j;
                            //locY = k;
                            //break;
                        //}
                    //}
                //}
                //C[locX][locY] = -objective[i].first/2.0;
                //C[locY][locX] = -objective[i].first/2.0;
            //}

            //// Define the A matrices and b vector
            //std::vector<doubleMatrix> As;
            //std::vector<double> b;
            //for (int i=0; i<momentMatrices[0].size(); i++) {
                //for (int j=i; j<momentMatrices[0][i].size(); j++) {

                    //// Trying to find the matrix relating this to other elements
                    //doubleMatrix A(matSize, doubleVector(matSize, 0));
                    //if (i == j) { 
                        //A[i][j] = 1;
                    //} else {
                        //A[i][j] = 0.5;
                    //}
                    //double bVal = 0;
                    //bool somethingFound = false;

                    //// For each term in this poly, try to find an early term 
                    //// that would reduce this to a constant
                    //for (int k=0; k<momentMatrices[0][i][j].size(); k++) {

                        //// If it's a constant
                        //if (momentMatrices[0][i][j][k].second.size() == 0) {
                            //bVal += momentMatrices[0][i][j][k].first;
                            //somethingFound = true;
                            //continue;
                        //}

                        //// Otherwise find an earlier variable
                        //monomial monomToFind = momentMatrices[0][i][j][k].second;
                        //bool found = false;
                        //for (int l=0; l<momentMatrices[0].size(); l++) {
                            //for (int m=l; m<momentMatrices[0][l].size(); m++) {
                                //if (l*momentMatrices[0].size() + m >= i*momentMatrices[0].size() + j) {
                                    //break;
                                //}
                                //if (momentMatrices[0][l][m].size() == 1 && momentMatrices[0][l][m][0].second == monomToFind) {
                                    //if (l == m) { 
                                        //A[l][m] = -momentMatrices[0][i][j][k].first;
                                    //} else {
                                        //A[l][m] = -momentMatrices[0][i][j][k].first / 2.0;
                                    //}
                                    //somethingFound = true;
                                    //found = true;
                                    //break;
                                //}
                            //}
                            //if (found) {
                                //break;
                            //}
                        //}

                    //}

                    //// Add this matrix
                    //if (somethingFound) {
                        //As.push_back(A);
                        //b.push_back(bVal);
                    //}

                //}
            //}
                
            //// Convert the constraints to the dual objective b.y
            //Poly newObjective;
            //for (int i=0; i<b.size(); i++) {
                //if (std::abs(b[i]) > 1e-10) {
                    //newObjective.push_back(std::make_pair(b[i], stringToMonomial("<D" + std::to_string(i) + ">")));
                //}
            //}
            //std::cout << "Dual objective: " << newObjective << std::endl;

            //// Convert the objective to the dual constraints C-\sum_i y_i A_i >= 0
            //std::vector<std::vector<Poly>> newMomentMat = std::vector<std::vector<Poly>>(matSize, std::vector<Poly>(matSize, Poly()));
            //for (int i=0; i<matSize; i++) {
                //for (int j=i; j<matSize; j++) {
                    //newMomentMat[i][j].push_back(std::make_pair(C[i][j], monomial()));
                    //for (int k=0; k<As.size(); k++) {
                        //if (std::abs(As[k][i][j]) > 1e-10) {
                            //newMomentMat[i][j].push_back(std::make_pair(-As[k][i][j], stringToMonomial("<D" + std::to_string(k) + ">")));
                        //}
                    //}
                //}
            //}

            //// Symmetrize the matrix
            //for (int i=0; i<matSize; i++) {
                //for (int j=i+1; j<matSize; j++) {
                    //newMomentMat[j][i] = newMomentMat[i][j];
                //}
            //}
            //std::cout << "Dual moment matrix: " << std::endl;
            //std::cout << newMomentMat << std::endl << std::endl;

            //// Objective being the inner product with the primal solution
            //newObjective = Poly();
            //for (int i=0; i<X.rows(); i++) {
                //for (int j=i; j<X.cols(); j++) {
                    //for (int k=0; k<newMomentMat[i][j].size(); k++) {
                        //if (i == j) {
                            //newObjective.push_back(std::make_pair(X(i,j), newMomentMat[i][j][k].second));
                        //} else {
                            //newObjective.push_back(std::make_pair(2.0*X(i,j), newMomentMat[i][j][k].second));
                        //}
                    //}
                //}
            //}

            //// Solve the dual
            //Poly objectiveDual = newObjective;
            //std::vector<std::vector<std::vector<Poly>>> momentMatricesDual = {newMomentMat};
            //std::vector<Poly> constraintsZeroDual = {};
            //std::vector<Poly> constraintsPositiveDual = {};
            //std::vector<double> varVals;
            //std::vector<Mon> varNames;
            //solveMOSEK(objectiveDual, momentMatricesDual, constraintsZeroDual, constraintsPositiveDual, verbosity, varNames, varVals);

            //return 0;

            //// Using this dual info, generate the sub-problem
            //std::vector<std::vector<std::vector<Poly>>> momentMatricesSub = generateAllMomentMatrices(bellFunc, level+1, subLevel, numToSample, use01);
            //Eigen::MatrixXd Y = replaceVariables(momentMatricesDual[0], varNames, varVals);
            //std::cout << "Y: " << std::endl;
            //std::cout << Y << std::endl << std::endl;
            //Poly objectiveSub;
            //for (int i=0; i<Y.rows()-iter; i++) {
                //for (int j=i; j<Y.cols()-iter; j++) {
                    //for (int k=0; k<momentMatricesSub[0][i][j].size(); k++) {
                        //if (i == j) {
                            //objectiveSub.push_back(std::make_pair(-Y(i,j), momentMatricesSub[0][i][j][k].second));
                        //} else {
                            //objectiveSub.push_back(std::make_pair(-2.0*Y(i,j), momentMatricesSub[0][i][j][k].second));
                        //}
                    //}
                //}
            //}
            //std::cout << "Sub problem objective: " << std::endl;
            //objectiveSub = simplify(objectiveSub);
            //std::cout << objectiveSub << std::endl << std::endl;
            //std::vector<Poly> constraintsZeroSub = {};
            //std::vector<Poly> constraintsPositiveSub = {};

            //// Solve the sub problem
            //std::vector<double> varValsSub;
            //std::vector<monomial> varNamesSub;
            //double valSub = solveMOSEK(objectiveSub, momentMatricesSub, constraintsZeroSub, constraintsPositiveSub, verbosity, varNamesSub, varValsSub);

            //// Check for convergence
            //if (std::abs(valSub) < 1e-5) {
                //std::cout << "Converged due to sub-problem giving zero" << std::endl;
                //break;
            //}

            //// Output the vars
            ////std::cout << "Vars: " << std::endl;
            ////for (int i=0; i<varNamesSub.size(); i++) {
                ////std::cout << varNamesSub[i] << " = " << varValsSub[i] << std::endl;
            ////}

            //// Generate the new constraint 
            //Poly newPositivityCon = objectiveSub;
            //newPositivityCon.push_back(std::make_pair(-valSub, monomial()));
            ////for (int i=0; i<objectiveSub.size(); i++) {
                ////for (int j=0; j<varNamesSub.size(); j++) {
                    ////if (varNamesSub[j] == objectiveSub[i].second) {
                        ////newPositivityCon.push_back(std::make_pair(objectiveSub[i].first, objectiveSub[i].second));
                        ////break;
                    ////}
                ////}
            ////}
            //newPositivityCon = simplify(newPositivityCon);
            //for (int i=0; i<newPositivityCon.size(); i++) {
                //newPositivityCon[i].first *= -1;
            //}
            //std::cout << "New positivity constraint: " << std::endl;
            //std::cout << newPositivityCon << std::endl << std::endl;

            //// Add this to the original main moment matrix, expanding by one
            //momentMatrices[0].push_back(std::vector<Poly>(matSize+1, Poly()));
            //for (int i=0; i<matSize; i++) {
                //momentMatrices[0][i].push_back(Poly());
            //}
            //momentMatrices[0][matSize][matSize] = newPositivityCon;

            //// Output the new moment matrix
            //std::cout << "New moment matrix: " << std::endl;
            //std::cout << momentMatrices[0] << std::endl << std::endl;

        //}

        //return 0;

    //}

    // Problem specific testing
    if (testing == 2) {
        std::vector<int> rowsToRemove = {};
        if (problemName == "CHSH") {
            rowsToRemove = {0};
        }
        for (int i=rowsToRemove.size()-1; i>=0; i--) {
            momentMatrices[0].erase(momentMatrices[0].begin() + rowsToRemove[i]);
            for (int j=0; j<momentMatrices[0].size(); j++) {
                momentMatrices[0][j].erase(momentMatrices[0][j].begin() + rowsToRemove[i]);
            }
        }
        if (problemName == "I3322") {
            constraintsZero.push_back(Poly("<A1>+<A2>"));
            constraintsZero.push_back(Poly("<B1>+<B2>"));

            constraintsZero.push_back(Poly("<A1B1>+<A1B2>"));
            constraintsZero.push_back(Poly("<A2B1>+<A2B2>"));
            constraintsZero.push_back(Poly("<A1B1>+<A2B1>"));

            constraintsZero.push_back(Poly("<A1B3>-<A2B3>"));
            constraintsZero.push_back(Poly("<A3B1>-<A3B2>"));
            constraintsZero.push_back(Poly("<A1B3>-<A3B2>"));
        }
    }

    if (testing == 3) {

        // Remove a random number of rows
        double chanceToRemove = rand() / (double)RAND_MAX;
        std::cout << "Chance to remove: " << chanceToRemove << std::endl;
        for (int i=0; i<momentMatrices[0].size(); i++) {
            if (rand() / (double)RAND_MAX < chanceToRemove) {
                momentMatrices[0].erase(momentMatrices[0].begin() + i);
                for (int j=0; j<momentMatrices[0].size(); j++) {
                    momentMatrices[0][j].erase(momentMatrices[0][j].begin() + i);
                }
                i--;
            }
        }
        std::cout << "Removed " << chanceToRemove << " of the rows" << std::endl;
        std::cout << "Size of moment matrix: " << momentMatrices[0].size() << "x" << momentMatrices[0][0].size() << std::endl;
        if (momentMatrices[0].size() == 0) {
            momentMatrices[0] = {{Poly(1)}};
        }

    }

    // TODO enforce that last row is mix of all the others
    if (testing == 4) {

        // Get the variable list from the moment matrix
        std::set<Mon> variableSet;
        variableSet.insert(Mon());
        for (size_t i=0; i<momentMatrices.size(); i++) {
            addVariables(variableSet, momentMatrices[i]);
        }
        for (size_t i=0; i<constraintsZero.size(); i++) {
            addVariables(variableSet, constraintsZero[i]);
        }
        addVariables(variableSet, objective);
        std::vector<Mon> variables = toVector(variableSet);

        // Turn this into a map
        std::map<Mon, Mon> varMap;
        for (size_t i=0; i<variables.size(); i++) {
            if (variables[i] != 1) {
                varMap[variables[i]] = Mon("<X" + std::to_string(i) + ">");
            }
        }

        // Apply the map
        for (size_t i=0; i<momentMatrices[0].size(); i++) {
            for (size_t j=0; j<momentMatrices[0][i].size(); j++) {
                momentMatrices[0][i][j] = momentMatrices[0][i][j].applyMap(varMap);
            }   
        }

        // Output the new matrix
        std::cout << momentMatrices[0] << std::endl;

        Poly newObj = objective.applyMap(varMap);
        std::cout << "Objective: " << newObj << std::endl;

        // The new constraints
        int numXs = variables.size();
        std::vector<Poly> newCons;
        for (int i=0; i<momentMatrices[0].size(); i++) {
            Poly newCon(momentMatrices[0][momentMatrices[0].size()-1][i]);
            for (int j=momentMatrices[0].size()-2; j>=0; j--) {
                newCon += Mon("<X" + std::to_string(numXs + j) + ">") * momentMatrices[0][j][i];
            }
            std::cout << newCon << " = 0 " << std::endl;
            newCons.push_back(newCon);
        }

        solveOptim(newObj, newCons);

        return 0;

    }

    if (testing >= 5) {

        // Solve just as a linear system
        std::vector<std::vector<Poly>> momentMatricesL1 = generateAllMomentMatrices(bellFunc, {}, 1, verbosity)[0];
        std::vector<std::vector<Poly>> momentMatricesL2 = generateAllMomentMatrices(bellFunc, {}, 2, verbosity)[0];
        std::vector<std::vector<Poly>> momentMatricesL3 = generateAllMomentMatrices(bellFunc, {}, 3, verbosity)[0];
        std::vector<std::vector<Poly>> momentMatCopy = momentMatrices[0];
        momentMatrices = {};

        std::map<Mon, std::complex<double>> xMap;
        xMap[Mon()] = 1;
        for (int i=0; i<momentMatCopy.size(); i++) {
            for (int j=0; j<momentMatCopy[i].size(); j++) {
                xMap[momentMatCopy[i][j].getKey()] = 0;
            }
        }

        double prevBound = 0;
        for (int i=0; i<testing; i++) {

            std::cout << "Linear solving " << i << std::endl;
            double res = maximizeMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, &xMap);

            //if (i == 200) {
                //std::cout << "Switching to L2" << std::endl;
                //momentMatCopy = momentMatricesL2;
            //}
            prevBound = res;

            // Get the resulting matrix
            std::cout << "Getting matrix" << std::endl;
            Eigen::MatrixXd X = Eigen::MatrixXd::Zero(momentMatCopy.size(), momentMatCopy.size());
            for (size_t i=0; i<momentMatCopy.size(); i++) {
                for (size_t j=i; j<momentMatCopy[i].size(); j++) {
                    X(i, j) = xMap[momentMatCopy[i][j].getKey()].real();
                    X(j, i) = X(i, j);
                }
            }

            // Get the eigenvalues and vectors
            std::cout << "Eigensolving" << std::endl;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(X);
            Eigen::VectorXd eigVals = es.eigenvalues().real();
            Eigen::MatrixXd eigVecs = es.eigenvectors().real();

            std::cout << constraintsPositive.size() << " " << res << " " << eigVals.minCoeff() << std::endl;

            // For each negative eigenvalue, store the eigenvector
            std::vector<Eigen::VectorXd> negEigVecs;
            for (int i=0; i<eigVals.size(); i++) {
                if (eigVals(i) < 0) {
                    negEigVecs.push_back(eigVecs.col(i));
                }
            }
            if (negEigVecs.size() == 0) {
                break;
            }

            // Add these as new cons
            std::cout << "Adding constraints" << std::endl;
            for (int i=0; i<negEigVecs.size(); i++) {
                Poly newCon;
                for (int j=0; j<momentMatCopy.size(); j++) {
                    for (int k=0; k<momentMatCopy[j].size(); k++) {
                        newCon[momentMatCopy[j][k].getKey()] += negEigVecs[i](j) * negEigVecs[i](k);
                    }
                }
                constraintsPositive.push_back(newCon);
            }

        }

        return 0;

    }

    //if (testing >= 6) {

        //// We are gonna test converting the moment matrix into a set of linear constraints
        //// A >= 0   becomes x^T A x >= 0 for some random x vectors

        //// Make a bunch of random vectors
        //int numVecs = testing;
        //std::vector<std::vector<double>> randomVectors;
        //for (int i=0; i<numVecs; i++) {
            //std::vector<double> vec(momentMatrices[0].size(), 0);
            //for (int j=0; j<vec.size(); j++) {
                //vec[j] = rand(-0.1, 0.1);
                ////if (std::abs(vec[j]) < 0.2) {
                    ////vec[j] = 0;
                ////} else if (vec[j] > 0) {
                    ////vec[j] = 1;
                ////} else {
                    ////vec[j] = -1;
                ////}
            //}
            //randomVectors.push_back(vec);
        //}

        //// newCon = x^T A x
        //std::vector<Poly> newCons;
        //for (int i=0; i<numVecs; i++) {
            //Poly newCon;
            //for (int j=0; j<momentMatrices[0].size(); j++) {
                //for (int k=0; k<momentMatrices[0][j].size(); k++) {
                    //newCon += momentMatrices[0][j][k] * randomVectors[i][j] * randomVectors[i][k];
                //}
            //}
            //newCons.push_back(newCon);
        //}
        //momentMatrices = {};

    //}

    // If using the dual TODO
    if (useDual) {
        primalToDual(objective, momentMatrices, constraintsZero, constraintsPositive);
    }

    // Output the problem
    if (verbosity >= 2) {
        if (objective.size() > 0) {
            std::cout << "Objective: " << std::endl;
            std::cout << objective << std::endl << std::endl;
        }
        if (momentMatrices.size() > 0) {
            for (int i=0; i<momentMatrices.size(); i++) {
                std::cout << "Moment matrix " << i << ": " << std::endl;
                std::cout << momentMatrices[i] << std::endl << std::endl;
            }
        }
        if (constraintsZero.size() > 0) {
            std::cout << "Zero constraints: " << std::endl;
            std::cout << constraintsZero << std::endl << std::endl;
        }
        if (constraintsPositive.size() > 0) {
            std::cout << "Positive constraints: " << std::endl;
            std::cout << constraintsPositive << std::endl << std::endl;
        }
    }

    // Convert to MOSEK form and solve
    double res = maximizeMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity);
    std::cout << "Result: " << res << std::endl;

    // If I3322, convert to the 0/1 version too
    if (problemName == "I3322" && !use01) {
        std::cout << "Result in 0/1: " << (res/4.0)-1.0 << std::endl;
    }

    // Exit without errors
    return 0;

}
