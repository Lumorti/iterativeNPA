// Standard includes
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>

// Import OpenMP
#include <omp.h>

// Import Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Import Spectra
//#include <Spectra/SymEigsSolver.h>

// Local libs
#include "poly.h"
#include "utils.h"
#include "mon.h"
#include "optMOSEK.h"
#include "optSCS.h"
#include "optOptim.h"

// Define addition and subtraction of maps
std::map<Mon, std::complex<double>> operator+(const std::map<Mon, std::complex<double>>& a, const std::map<Mon, std::complex<double>>& b) {
    std::map<Mon, std::complex<double>> res = a;
    for (auto& pair : b) {
        res[pair.first] += pair.second;
    }
    return res;
}
std::map<Mon, std::complex<double>> operator-(const std::map<Mon, std::complex<double>>& a, const std::map<Mon, std::complex<double>>& b) {
    std::map<Mon, std::complex<double>> res = a;
    for (auto& pair : b) {
        res[pair.first] -= pair.second;
    }
    return res;
}

// Generic entry function
int main(int argc, char* argv[]) {

    // Assume random seed unless later set
    srand(time(NULL));

    // Define the Bell scenario
    int level = 1;
    bool useDual = false;
    Poly bellFunc("<A1B1>+<A1B2>+<A2B1>-<A2B2>");
    int testing = 0;
    int verbosity = 1;
    int maxIters = 10000000;
    int numCores = 1;
    bool use01 = false;
    double stepSize = 100;
    double tolerance = 1e-6;
    std::string solver = "MOSEK";
    std::string seed = "";
    std::string problemName = "CHSH";
    std::vector<std::string> extraMoments;

    // Process command-line args
    for (int i=1; i<argc; i++) {
        std::string argAsString = std::string(argv[i]);

        // CHSH (c = 2, q = 2sqrt(2) w/ level 1)
        if (argAsString == "--chsh" || argAsString == "--CHSH") {
            bellFunc = Poly("<A1B1>+<A1B2>+<A2B1>-<A2B2>");
            problemName = "CHSH";

        // If told to use a 0/1 output
        } else if (argAsString == "-01") {
            use01 = true;

        // If told to use SCS
        } else if (argAsString == "-C") {
            solver = "SCS";

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

        // Randomized version of the above
        // for -S 1
        // l1: 3.52733
        // l2: 3.37212 (0.06s)
        } else if (argAsString == "--R3322" || argAsString == "--r3322") {
            bellFunc = Poly("<A1>-<A2>+<B1>-<B2>-<A1B1>+<A1B2>+<A2B1>-<A2B2>+<A1B3>+<A2B3>+<A3B1>+<A3B2>");
            for (auto term : bellFunc) {
                bellFunc[term.first] = rand(-1, 1);
            }
            problemName = "R3322";

        // Randomized version for arbitrary number of inputs TODO
        //
        // for -S 1 -RXX22 100
        // SCS 34s 1044
        // PRJ 28s 1068 2839i
        // CEN  0s 1137
        //
        // for -S 1 -RXX22 200
        // SCS 61m 3009
        // PRJ  7m 3069 4507i
        // CEN  0m 3274
        //
        } else if (argAsString == "--RXX22") {
            int numInputs = std::stoi(argv[i+1]);
            bellFunc = Poly();
            for (int i=1; i<=numInputs; i++) {
                bellFunc += Poly("<A" + std::to_string(i) + ">");
                bellFunc -= Poly("<B" + std::to_string(i) + ">");
            }
            for (int i=1; i<=numInputs; i++) {
                for (int j=1; j<=numInputs; j++) {
                    bellFunc += Poly("<A" + std::to_string(i) + "B" + std::to_string(j) + ">");
                }
            }
            bellFunc.randomize();
            problemName = "RXX22";
            i++;

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

        // Set the step size
        } else if (argAsString == "-A") {
            stepSize = std::stod(argv[i+1]);
            i++;

        // If we're testing
        } else if (argAsString == "-t") {
            testing = std::stoi(argv[i+1]);
			i++;

        // If setting verbosity
        } else if (argAsString == "-v") {
            verbosity = std::stoi(argv[i+1]);
            i++;

        // If setting core count
        } else if (argAsString == "-c") {
            numCores = std::stoi(argv[i+1]);
            omp_set_num_threads(numCores);
            Eigen::setNbThreads(numCores);
            i++;

        // Output the help
        } else if (argAsString == "-h" || argAsString == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  --chsh          Use the CHSH scenario" << std::endl;
            std::cout << "  --I3322         Use the I3322 scenario" << std::endl;
            std::cout << "  --R3322         Use the randomized I3322 scenario" << std::endl;
            std::cout << "  --RXX22         Use the randomized scenario with arbitrary number of inputs" << std::endl;
            std::cout << "  -c <int>        Number of CPU threads to use" << std::endl;
            std::cout << "  -l <int>        Level of the moment matrix" << std::endl;
            std::cout << "  -i <int>        Iteration limit" << std::endl;
            std::cout << "  -S <str>        Seed for the random number generator" << std::endl;
            std::cout << "  -v <int>        Set the verbosity level" << std::endl;
            std::cout << "  -p <dbl>        Set the starting penalty" << std::endl;
            std::cout << "  -P <int>        Set the number of penalties to use" << std::endl;
            std::cout << "  -T <dbl>        Set the tolerance" << std::endl;
            std::cout << "  -t <int>        Testing some new idea" << std::endl;
            std::cout << "  -e <str>        Add an extra moment to the top row" << std::endl;
            std::cout << "  -h              Display this help message" << std::endl;
            std::cout << "  -D              Use the dual of the problem" << std::endl;
            std::cout << "  -C              Use SCS to solve rather than MOSEK" << std::endl;
            std::cout << "  -A <dbl>        Set the step size" << std::endl;
            std::cout << "  -01             Use 0/1 output instead of -1/1" << std::endl;
            return 0;

        // If setting the tolerance
        } else if (argAsString == "-T") {
            tolerance = std::stod(argv[i+1]);
            i++;

        // If using the dual
        } else if (argAsString == "-D") {
            useDual = true;

        // Otherwise we don't know what this is
        } else {
            std::cout << "Unknown argument: " << argAsString << std::endl;
            return 1;

        }
    }

    // Define the moment matrix
    if (verbosity >= 1) {
        std::cout << "Generating moment matrices..." << std::endl;
    }
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
        if (verbosity >= 3) {
            std::cout << "Primal objective: " << std::endl;
            std::cout << objective << std::endl << std::endl;
            std::cout << "Primal moment matrix: " << std::endl;
            std::cout << momentMatrices[0] << std::endl << std::endl;
        }
        if (verbosity >= 1) {
            std::cout << "Converting to dual..." << std::endl;
        }
        primalToDual(objective, momentMatrices, constraintsZero, constraintsPositive);
        if (verbosity >= 3) {
            std::cout << "Dual objective: " << std::endl;
            std::cout << objective << std::endl << std::endl;
            std::cout << "Dual moment matrix: " << std::endl;
            std::cout << momentMatrices[0] << std::endl << std::endl;
        }
    }

    // The idea of using the dual to find cons that aren't in level 1
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

        //std::cout << "Moment matrix: " << std::endl;
        //std::cout << momentMatrices[0] << std::endl;


        //std::cout << "Moment matrix with random values: " << std::endl;
        //std::cout << X << std::endl;

        // Do Guassian elimination on matrix X
        //std::vector<std::vector<Poly>> newMat = momentMatrices[0];
        //for (int i=0; i<X.rows(); i++) {
            //for (int j=i+1; j<X.rows(); j++) {
                //std::complex<double> toMult = X(j,i) / X(i,i);
                //X.row(j) -= X.row(i) * toMult;
                //for (int k=0; k<X.cols(); k++) {
                    //newMat[j][k] -= newMat[i][k] * toMult;
                //}
            //}
        //}

        //std::cout << "After GE: " << std::endl;
        //std::cout << X << std::endl;

        //std::cout << "After GE: " << std::endl;
        //std::cout << newMat << std::endl;

        //// Add the diagonals as positivity constraints
        //for (int i=0; i<X.rows(); i++) {
            //constraintsPositive.push_back(newMat[i][i]);
        //}

    }

    if (testing == 4) {

        double sol = solveOptim(objective, constraintsPositive, momentMatrices[0], verbosity);
        return 0;

    }

    if (testing == 5) {

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
        for (int iter=0; iter<maxIters; iter++) {

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

            if (eigVals.minCoeff() > tolerance) {
                break;
            }

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

            // Limit the number of constraints, removing old ones
            //if (constraintsPositive.size() > 150) {
                //constraintsPositive.erase(constraintsPositive.begin(), constraintsPositive.begin() + constraintsPositive.size() - 100);
            //}

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

    // Go very far in the objective, then project back onto the set
    if (testing == 2) {

        // Starting output
        if (verbosity >= 1) {
            std::cout << "Moment matrix has size " << momentMatrices[0].size() << std::endl;
            int estimatedIters = int(22.4834 * momentMatrices[0].size() / 2.0 + 283.873);
            std::cout << "Should require on the order of " << estimatedIters << " iterations" << std::endl;
            std::cout << "Starting projection..." << std::endl;
        }

        // Get an easy bound by setting vars to the identity * minimum eigen
        std::map<Mon, std::complex<double>> varVals0;
        if (verbosity >= 1) {
            std::set<Mon> vars0;
            std::set<Mon> varsDiag;
            addVariables(vars0, momentMatrices[0]);
            for (int i=0; i<momentMatrices[0].size(); i++) {
                varsDiag.insert(momentMatrices[0][i][i].getKey());
            }
            for (auto& mon : vars0) {
                varVals0[mon] = 0;
            }
            Eigen::MatrixXd X0 = replaceVariables(momentMatrices[0], varVals0).real();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(X0);
            double minEig = es.eigenvalues().minCoeff();
            for (auto& mon : varsDiag) {
                varVals0[mon] = minEig;
            }
            double easyBound = -objective.eval(varVals0).real();
            std::cout << "Center bound: " << easyBound << std::endl;
            //X0 = replaceVariables(momentMatrices[0], varVals0).real();
            //es.compute(X0);
            //std::cout << "Min eig after: " << es.eigenvalues().minCoeff() << std::endl;

        }

        // Create new moment matrix and then add equalities
        std::vector<std::vector<Poly>> newMomentMatrix = momentMatrices[0];
        int newInd = 0;
        for (int i=0; i<newMomentMatrix.size(); i++) {
            for (int j=i; j<newMomentMatrix[i].size(); j++) {
                newMomentMatrix[i][j] = Poly("<M" + std::to_string(newInd) + ">");
                newMomentMatrix[j][i] = newMomentMatrix[i][j];
                constraintsZero.push_back(momentMatrices[0][i][j] - newMomentMatrix[i][j]);
                newInd++;
            }
        }
        momentMatrices = {newMomentMatrix};

        // Output things if super verbose
        if (verbosity >= 3) {
            std::cout << "New moment matrix: " << std::endl;
            std::cout << momentMatrices[0] << std::endl;
            std::cout << "Zero constraints: " << std::endl;
            for (auto con : constraintsZero) {
                std::cout << con << std::endl;
            }
        }

        // Get the list of variables
        std::set<Mon> vars;
        addVariables(vars, momentMatrices[0]);
        for (auto con : constraintsZero) {
            addVariables(vars, con);
        }
        if (vars.count(Mon()) > 0) {
            vars.erase(Mon());
        }
        std::vector<Mon> varList(vars.begin(), vars.end());
        std::map<Mon, int> varLocs;
        for (int i=0; i<varList.size(); i++) {
            varLocs[varList[i]] = i;
        }

        // Starting point
        double distance = stepSize;
        std::map<Mon, std::complex<double>> varVals;
        for (auto& mon : vars) {
            if (objective.contains(mon)) {
                if (objective[mon].real() > 0) {
                    varVals[mon] = distance;
                } else {
                    varVals[mon] = -distance;
                }
            } else {
                varVals[mon] = 0;
            }
            //if (varVals0.count(mon) > 0) {
                //varVals[mon] = varVals0[mon];
            //}
        }

        // The matrix defining the linear constraints
        Eigen::SparseMatrix<double> A(constraintsZero.size(), varList.size());
        std::vector<Eigen::Triplet<double>> tripletList;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(constraintsZero.size());
        for (int i=0; i<constraintsZero.size(); i++) {
            for (auto& term : constraintsZero[i]) {
                if (term.first.size() != 0) {
                    tripletList.push_back(Eigen::Triplet<double>(i, varLocs[term.first], term.second.real()));
                }
            }
            b(i) = -constraintsZero[i][Mon()].real();
        }
        A.setFromTriplets(tripletList.begin(), tripletList.end());

        // If verbose, output everything
        if (verbosity >= 3) {
            std::cout << "vars: " << std::endl;
            for (int i=0; i<varList.size(); i++) {
                std::cout << varList[i] << " " << varVals[varList[i]] << std::endl;
            }
            std::cout << "A: " << std::endl;
            std::cout << A << std::endl;
            std::cout << "b: " << std::endl;
            std::cout << b << std::endl;
        }

        // Set up the linear solver
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> linSolver;
        linSolver.setTolerance(tolerance/10.0);

        // Precaculate the pseudo-inverse
        std::cout << "Calculating inverse..." << std::endl;
        Eigen::SparseMatrix<double> ASquare = A * A.transpose();
        linSolver.compute(ASquare);
        Eigen::SparseMatrix<double> I(ASquare.cols(), ASquare.cols());
        I.setIdentity();
        Eigen::SparseMatrix<double> APseudoSparse = A.transpose() * linSolver.solve(I);
        //std::cout << "Making sparse..." << std::endl;
        //std::vector<Eigen::Triplet<double>> tripletListPseudo;
        //for (int i=0; i<APseudo.rows(); i++) {
            //for (int j=0; j<APseudo.cols(); j++) {
                //if (std::abs(APseudo(i, j)) > 1e-10) {
                    //tripletListPseudo.push_back(Eigen::Triplet<double>(i, j, APseudo(i, j)));
                //}
            //}
        //}
        //Eigen::SparseMatrix<double> APseudoSparse(APseudo.rows(), APseudo.cols());
        //APseudoSparse.setFromTriplets(tripletListPseudo.begin(), tripletListPseudo.end());
        std::cout << "APseudo has size: " << APseudoSparse.rows() << " " << APseudoSparse.cols() << std::endl;
        std::cout << "APseudo has num nonzeros: " << APseudoSparse.nonZeros() << std::endl;

        // Keep iterating until reaching limit or convergence TODO
        linSolver.compute(A);
        double prevMinEig = -1;
        double prevLinError = 1e-3;
        Eigen::VectorXd x = Eigen::VectorXd::Zero(varList.size());
        std::vector<Eigen::VectorXd> xList;
        std::vector<double> yVals;
        std::vector<int> prevMaxIters;
        for (int iter=0; iter<maxIters; iter++) {

            // Check the objective
            std::chrono::steady_clock::time_point timeStart = std::chrono::steady_clock::now();
            double newObjVal = -objective.eval(varVals).real();

            // SDP projection
            Eigen::MatrixXd X = replaceVariables(momentMatrices[0], varVals).real();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(X);
            Eigen::VectorXd eigVals = es.eigenvalues().real();
            Eigen::MatrixXd eigVecs = es.eigenvectors().real();
            double minEig = eigVals.minCoeff();
            Eigen::MatrixXd diagEigVals = Eigen::MatrixXd::Zero(X.rows(), X.cols());
            for (int i=0; i<eigVals.size(); i++) {
                diagEigVals(i, i) = std::max(eigVals(i), 0.0);
            }
            Eigen::MatrixXd proj = eigVecs * diagEigVals * eigVecs.transpose();
            for (int i=0; i<proj.rows(); i++) {
                for (int j=i; j<proj.cols(); j++) {
                    varVals[momentMatrices[0][i][j].getKey()] = proj(i, j);
                }
            }

            // Convert to a vector
            for (int i=0; i<varList.size(); i++) {
                x(i) = varVals[varList[i]].real();
            }
            Eigen::VectorXd errVec = A*x - b;
            double errorLin = errVec.norm();

            // Linear projection
            //linSolver.setTolerance(std::min(prevLinError / 100.0, 1e-6));
            //linSolver.setTolerance(tolerance/10.0);
            //Eigen::VectorXd projX = linSolver.solveWithGuess(b, x);
            //Eigen::VectorXd projX = x - APseudoDense * errVec;
            Eigen::VectorXd projX = x - APseudoSparse * errVec;
            for (int i=0; i<varList.size(); i++) {
                varVals[varList[i]] = projX(i);
            }

            // Average the two vectors
            Eigen::VectorXd avgX = (x + projX) / 2;

            // Calculate how long each iteration takes
            std::chrono::steady_clock::time_point timeFinished = std::chrono::steady_clock::now();
            int perIter = std::chrono::duration_cast<std::chrono::milliseconds>(timeFinished - timeStart).count();
            double changeInLinError = std::abs(errorLin / prevLinError);
            double changeInMinEig = std::abs(minEig / prevMinEig);
            if (changeInLinError > 1) {
                changeInLinError = 0;
            }
            if (changeInMinEig > 1) {
                changeInMinEig = 0;
            }
            int estimatedMaxIter;
            if (changeInLinError > changeInMinEig) {
                estimatedMaxIter = std::log(tolerance) / std::log(changeInLinError);
            } else {
                estimatedMaxIter = std::log(tolerance) / std::log(changeInMinEig);
            }
            if (estimatedMaxIter > 0) {
                prevMaxIters.push_back(estimatedMaxIter);
            }
            if (prevMaxIters.size() > 5) {
                prevMaxIters.erase(prevMaxIters.begin());
            }
            if (prevMaxIters.size() > 0) {
                estimatedMaxIter = 0;
                for (int i=0; i<prevMaxIters.size(); i++) {
                    estimatedMaxIter += prevMaxIters[i];
                }
                estimatedMaxIter /= prevMaxIters.size();
            }

            // Per iteration output
            if (verbosity == -1) {
                for (size_t i=0; i<varList.size(); i++) {
                    std::cout << x(i);
                    if (i < varList.size()-1) {
                        std::cout << "\t ";
                    } else {
                        std::cout << std::endl;
                    }
                }
                for (size_t i=0; i<varList.size(); i++) {
                    std::cout << projX(i);
                    if (i < varList.size()-1) {
                        std::cout << "\t ";
                    } else {
                        std::cout << std::endl;
                    }
                }
                for (size_t i=0; i<varList.size(); i++) {
                    std::cout << avgX(i);
                    if (i < varList.size()-1) {
                        std::cout << "\t ";
                    } else {
                        std::cout << std::endl;
                    }
                }
            } else if (verbosity >= 1) {
                std::cout << std::setw(8) << iter << "i " << std::setw(12) << newObjVal << " " << std::setw(12) << minEig << " " << std::setw(12) << errorLin << " " << std::setw(6) << perIter << "ms/i " << std::setw(8) << estimatedMaxIter << "i           \r" << std::flush;
            }

            // Stopping condition
            prevLinError = errorLin;
            prevMinEig = minEig;
            if (std::abs(errorLin) < tolerance && minEig > -tolerance) {
                break;
            }

        }

        if (verbosity >= 1) {
            std::cout << std::endl;
            varVals[Mon()] = 1;
            std::cout << "Final objective: " << -objective.eval(varVals).real() << std::endl;
            Eigen::MatrixXcd X = replaceVariables(momentMatrices[0], varVals);
            if (verbosity >= 3) {
                std::cout << "Final moment matrix: " << std::endl;
                std::cout << X << std::endl;
                std::cout << "Final var vals: " << std::endl;
                for (size_t i=0; i<varList.size(); i++) {
                    std::cout << varList[i] << " -> " << varVals[varList[i]] << std::endl;
                }
                for (int i=0; i<X.rows(); i++) {
                    double rowSum = 0;
                    for (int j=0; j<X.cols(); j++) {
                        if (i != j) {
                            rowSum += std::abs(X(i, j));
                        }
                    }
                    std::cout << "Diag = " << std::abs(X(i, i)) << " Row sum = " << rowSum << std::endl;
                }
            }
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(X);
            double minEig = es.eigenvalues().minCoeff();
            std::cout << "Final min eig: " << minEig << std::endl;
            Eigen::VectorXd x = Eigen::VectorXd::Zero(varList.size());
            for (size_t i=0; i<varList.size(); i++) {
                x(i) = varVals[varList[i]].real();
            }
            std::cout << "Final linear error: " << (A*x - b).norm() << std::endl;
        }

        return 0;

    }

    // Non-verbose output
    if (verbosity >= 1) {
        int largestMomentMatrix = 0;
        for (size_t i=0; i<momentMatrices.size(); i++) {
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
        varBounds = {-100000, 100000};
    }
    double res = 0;
    if (solver == "MOSEK") {
        std::cout << "Solving with MOSEK..." << std::endl;
        res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varBounds);
    } else if (solver == "SCS") {
        std::cout << "Solving with SCS..." << std::endl;
        res = solveSCS(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varBounds);
    }
    if (useDual) {
        res = -res;
    }
    std::cout << "Result: " << res << std::endl;

    // If I3322, convert to the 0/1 version too
    if (problemName == "I3322" && !use01) {
        std::cout << "Result in 0/1: " << (res/4.0)-1.0 << std::endl;
    }

    // Exit without errors
    return 0;

}
