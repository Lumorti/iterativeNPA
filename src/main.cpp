// Standard includes
#include <iostream>
#include <vector>

// Import Eigen
#include <Eigen/Dense>

// Local libs
#include "poly.h"
#include "utils.h"
#include "mon.h"
#include "mosek.h"
#include "optim.h"

// Generic entry function
int main(int argc, char* argv[]) {

    // Assume random seed unless later set
    srand(time(NULL));

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
    int maxIters = 1;
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
        // l2 (37x37) = 5.00379
        // l3 (253x253) = 5.0035
        } else if (argAsString == "--I3322" || argAsString == "--i3322") {
            if (use01) {
                bellFunc = Poly("-<A2>-<B1>-2<B2>+<A1B1>+<A1B2>+<A2B1>+<A2B2>-<A1B3>+<A2B3>-<A3B1>+<A3B2>");
            } else {
                bellFunc = Poly("<A1>-<A2>+<B1>-<B2>-<A1B1>+<A1B2>+<A2B1>-<A2B2>+<A1B3>+<A2B3>+<A3B1>+<A3B2>");
            }
            problemName = "I3322";

        // Randomized version of the above TODO
        // for -S 1
        // l1: 3.52733
        // l2: 3.37212 (0.06s)
        } else if (argAsString == "--R3322" || argAsString == "--r3322") {
            bellFunc = Poly("<A1>-<A2>+<B1>-<B2>-<A1B1>+<A1B2>+<A2B1>-<A2B2>+<A1B3>+<A2B3>+<A3B1>+<A3B2>");
            for (auto term : bellFunc) {
                bellFunc[term.first] = rand(-1, 1);
            }
            problemName = "R3322";

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
            srand(std::stoi(seed));
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

    // Define the moment matrix
    std::vector<std::vector<std::vector<Poly>>> momentMatrices = generateAllMomentMatrices(bellFunc, {}, level, verbosity);

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

    // If using the dual
    if (useDual) {
        if (verbosity >= 2) {
            std::cout << "Primal objective: " << std::endl;
            std::cout << objective << std::endl << std::endl;
            std::cout << "Primal moment matrix: " << std::endl;
            std::cout << momentMatrices[0] << std::endl << std::endl;
        }
        primalToDual(objective, momentMatrices, constraintsZero, constraintsPositive);
    }

    // The idea of using the dual to find cons that aren't in level 1 TODO
    if (testing == 1) {

        // Generate the moment matrices
        std::vector<std::vector<std::vector<Poly>>> momentMatricesL1 = generateAllMomentMatrices(bellFunc, {}, 1, verbosity);
        std::vector<std::vector<std::vector<Poly>>> momentMatricesL2 = generateAllMomentMatrices(bellFunc, {}, std::max(level, 2), verbosity);

        // Test the level 2 case
        std::map<Mon, std::complex<double>> varValsL2OG;
        double resL2OG = solveMOSEK(objective, momentMatricesL2, constraintsZero, constraintsPositive, verbosity-1, {-1, 1}, &varValsL2OG);

        // Keep iterating, for now just a fixed number of times
        double prevRes = 0;
        for (int iter=0; iter<maxIters; iter++) {

            // Get the L1 primal solution
            std::map<Mon, std::complex<double>> varValsL1;
            double resL1 = solveMOSEK(objective, momentMatricesL1, constraintsZero, constraintsPositive, verbosity-1, {-1, 1}, &varValsL1);
            Eigen::MatrixXcd Y = replaceVariables(momentMatricesL1[0], varValsL1);
            if (verbosity >= 1) {
                std::cout << iter << " " << resL1 << " " << resL2OG << " " << (resL1 < resL2OG - 1e-6) << std::endl;
            }

            // Stop if nothing changes
            if (std::abs(resL1-prevRes) < 1e-10) {
                break;
            }
            prevRes = resL1;

            // Debugging output
            if (verbosity >= 3) {
                std::cout << "Primal moment matrix: " << std::endl;
                std::cout << momentMatricesL1[0] << std::endl << std::endl;
                std::cout << "Primal objective: " << std::endl;
                std::cout << objective << std::endl << std::endl;
                std::cout << "L1 result: " << resL1 << std::endl;
                std::cout << "L1 solution: " << std::endl;
                std::cout << Y.real() << std::endl << std::endl;
                std::cout << "L2 moment matrix: " << std::endl;
                std::cout << momentMatricesL2[0] << std::endl << std::endl;
            }

            // Create a new moment matrix with plain variables
            std::vector<std::vector<std::vector<Poly>>> newMomentMatrices = momentMatricesL2;
            int nextVarInd = 0;
            for (int i=0; i<momentMatricesL2[0].size(); i++) {
                for (int j=i; j<momentMatricesL2[0][i].size(); j++) {
                    newMomentMatrices[0][i][j] = Poly("<M" + std::to_string(nextVarInd) + ">");
                    newMomentMatrices[0][j][i] = newMomentMatrices[0][i][j];
                    nextVarInd++;
                }
            }

            // Get the list of all variables
            std::map<Mon, std::complex<double>> varValsAll = varValsL1;
            for (int i=0; i<momentMatricesL2[0].size(); i++) {
                for (int j=0; j<momentMatricesL2[0][i].size(); j++) {
                    Mon key = momentMatricesL2[0][i][j].getKey();
                    if (varValsAll.find(key) == varValsAll.end()) {
                        varValsAll[key] = 0;
                    }
                }
            }

            // The mapping from expectation values to P variables
            std::map<Mon, int> monToInd;
            std::map<int, Mon> indToMon;
            int nextInd = 0;
            for (auto mon : varValsAll) {
                monToInd[mon.first] = nextInd;
                indToMon[nextInd] = mon.first;
                nextInd++;
            }

            // The new objective, min Y.P
            Poly newObjective;
            std::vector<Poly> newZeroCons(nextInd);
            std::vector<Poly> newPosCons;
            for (auto mon : varValsL1) {
                int pInd = monToInd[mon.first];
                newZeroCons[pInd] = Poly("<P" + std::to_string(pInd) + ">");
                newObjective -= Poly(mon.second, "<P" + std::to_string(pInd) + ">");
            }
            for (int i=0; i<momentMatricesL2[0].size(); i++) {
                for (int j=0; j<momentMatricesL2[0][i].size(); j++) {
                    newZeroCons[monToInd[momentMatricesL2[0][i][j].getKey()]] -= newMomentMatrices[0][i][j];
                }
            }

            // Debugging output
            if (verbosity >= 3) {
                std::cout << "P map: " << std::endl;
                for (auto mon : varValsAll) {
                    std::cout << mon.first << " -> " << monToInd[mon.first] << std::endl;
                }
                std::cout << "Objective sub:" << std::endl;
                std::cout << newObjective << std::endl << std::endl;
                std::cout << "Zero cons sub:" << std::endl;
                std::cout << newZeroCons << std::endl << std::endl;
                std::cout << "Moment matrix sub: " << std::endl;
                std::cout << newMomentMatrices[0] << std::endl << std::endl;
            }

            // Solve the sub problem
            std::map<Mon, std::complex<double>> varValsL2;
            double resL2 = solveMOSEK(newObjective, newMomentMatrices, newZeroCons, newPosCons, verbosity-1, {-100, 100}, &varValsL2);

            // Stop if the new constraint isn't even negative
            if (resL2 < 1e-10) {
                break;
            }

            // Generate the new constraint
            Poly newCon;
            for (auto mon : varValsL1) {
                newCon += Poly(varValsL2[Mon("<P" + std::to_string(monToInd[mon.first]) + ">")], mon.first);
            }
            newCon.clean();
            constraintsPositive.push_back(newCon);
            if (verbosity >= 2) {
                std::cout << "New constraint: " << newCon << std::endl;
                std::cout << "New con evalled: " << newCon.eval(varValsL1) << std::endl;
            }

        }

        return 0;

    }

    if (testing == 3) {

        // Put random values into the moment matrix
        std::map<Mon, std::complex<double>> varVals;
        for (int i=0; i<momentMatrices[0].size(); i++) {
            for (int j=0; j<momentMatrices[0][i].size(); j++) {
                varVals[momentMatrices[0][i][j].getKey()] = std::complex<double>(rand(-1, 1), 0);
            }
        }
        Eigen::MatrixXcd X = replaceVariables(momentMatrices[0], varVals);

        // Check the rank
        Eigen::FullPivLU<Eigen::MatrixXcd> lu(X);
        std::cout << "Rank: " << lu.rank() << " / " << X.rows() << std::endl;

        std::cout << "Moment matrix: " << std::endl;
        std::cout << momentMatrices[0] << std::endl;


        std::cout << "Moment matrix with random values: " << std::endl;
        std::cout << X << std::endl;

        // Do Guassian elimination on matrix X
        std::vector<std::vector<Poly>> newMat = momentMatrices[0];
        for (int i=0; i<X.rows(); i++) {
            for (int j=i+1; j<X.rows(); j++) {
                std::complex<double> toMult = X(j,i) / X(i,i);
                X.row(j) -= X.row(i) * toMult;
                for (int k=0; k<X.cols(); k++) {
                    newMat[j][k] -= newMat[i][k] * toMult;
                }
            }
        }

        std::cout << "After GE: " << std::endl;
        std::cout << X << std::endl;

        std::cout << "After GE: " << std::endl;
        std::cout << newMat << std::endl;

        // Add the diagonals as positivity constraints
        for (int i=0; i<X.rows(); i++) {
            constraintsPositive.push_back(newMat[i][i]);
        }

    }

    if (testing == 4) {

        double sol = solveOptim(objective, constraintsPositive, momentMatrices[0], verbosity);
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
        for (int iter=0; iter<testing; iter++) {

            std::cout << "Linear solving " << iter << std::endl;
            double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, {-1, 1}, &xMap);

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

    if (testing == 2) {

        // Start from all zeros
        std::map<Mon, std::complex<double>> startVals;
        std::map<Mon, std::complex<double>> endVals;
        for (int i=0; i<momentMatrices[0].size(); i++) {
            for (int j=0; j<momentMatrices[0][i].size(); j++) {
                std::cout << momentMatrices[0][i][j] << std::endl;
                if (!momentMatrices[0][i][j].isConstant()) {
                    if (i == j) {
                        //startVals[momentMatrices[0][i][j].getKey()] = std::complex<double>(-1, 0);
                        //endVals[momentMatrices[0][i][j].getKey()] = std::complex<double>(rand(0, 1), 0);
                        startVals[momentMatrices[0][i][j].getKey()] = std::complex<double>(0, 0);
                        endVals[momentMatrices[0][i][j].getKey()] = std::complex<double>(0, 0);
                    } else {
                        startVals[momentMatrices[0][i][j].getKey()] = std::complex<double>(0, 0);
                        endVals[momentMatrices[0][i][j].getKey()] = std::complex<double>(rand(-1, 1), 0);
                    }
                }
            }
        }
        //for (auto mon : objective) {
            //if (!mon.first.isConstant()) {
                //endVals[mon.first] = std::complex<double>(1, 0) * mon.second;
            //}
        //}

        double stepSize = 0.1;
        for (int i=0; i<maxIters; i++) {

            std::map<Mon, std::complex<double>> varVals;
            for (auto& mon : startVals) {
                varVals[mon.first] = (startVals[mon.first] + endVals[mon.first]) / 2.0;
            }

            // Check if X is positive
            Eigen::MatrixXcd X = replaceVariables(momentMatrices[0], varVals);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(X);
            double minEig = es.eigenvalues().minCoeff();
            std::cout << X << std::endl;
            std::cout << "Min eig: " << minEig << std::endl;
            double newObjVal = objective.eval(varVals).real();
            std::cout << "Objective: " << newObjVal << std::endl;

            if (minEig < 0) {
                endVals = varVals;
            } else {
                startVals = varVals;
            }

            //if (minEig < 0) {
                //for (auto& mon : objective) {
                    //varVals[mon.first] -= stepSize * mon.second;
                //}
                //stepSize *= 0.9;
            //}

            //// Travel in the direction of the objective
            //for (auto& mon : objective) {
                //varVals[mon.first] += stepSize * mon.second;
            //}

        }

        return 0;

    }

    // Non-verbose output
    if (verbosity >= 1) {
        int largestMomentMatrix = 0;
        for (int i=0; i<momentMatrices.size(); i++) {
            if (momentMatrices[i].size() > largestMomentMatrix) {
                largestMomentMatrix = momentMatrices[i].size();
            }
        }
        std::cout << "Largest moment matrix has size " << largestMomentMatrix << std::endl;
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
    std::pair<int,int> varBounds = {-1, 1};
    if (useDual) {
        varBounds = {-10000, 10000};
    }
    double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varBounds);
    std::cout << "Result: " << res << std::endl;

    // If I3322, convert to the 0/1 version too
    if (problemName == "I3322" && !use01) {
        std::cout << "Result in 0/1: " << (res/4.0)-1.0 << std::endl;
    }

    // Exit without errors
    return 0;

}
