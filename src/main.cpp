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

// Import Optim
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"

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

// The structure used to hold the data for the optimisation
struct optimData {
    Poly objective;
    Eigen::SparseMatrix<double> A;
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>* linSolver;
    Eigen::VectorXd b;
    std::vector<std::vector<Poly>> momentMatrix;
    std::vector<Mon> varList;
    int maxIters;
};

// Cost/gradient function for optim
static double gradFunction(const Eigen::VectorXd& x, Eigen::VectorXd* gradOut, void* optData) {
    
    // Recast this generic pointer into the correct format
    optimData* optDataRecast = reinterpret_cast<optimData*>(optData);

    // Convert the x vector into a map
    std::map<Mon, std::complex<double>> xMap;
    for (int i=0; i<x.size(); i++) {
        xMap[optDataRecast->varList[i]] = x(i);
    }

    // Objective value
    double obj = -optDataRecast->objective.eval(xMap).real();

    // Calculate linear error
    double errorLin = (optDataRecast->A*x - optDataRecast->b).norm();

    // Calculate the eigenspectrum
    Eigen::MatrixXd X = replaceVariables(optDataRecast->momentMatrix, xMap).real();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(X);
    Eigen::VectorXd eigVals = es.eigenvalues().real();
    Eigen::MatrixXd eigVecs = es.eigenvectors().real();
    double minEig = eigVals.minCoeff();

    // The cost function to minimize
    double cost = errorLin*errorLin + minEig*minEig;

    // Calculate the gradient
    if (gradOut) {

        // SDP projection
        Eigen::MatrixXd diagEigVals = Eigen::MatrixXd::Zero(X.rows(), X.cols());
        for (int i=0; i<eigVals.size(); i++) {
            diagEigVals(i, i) = std::max(eigVals(i), 0.0);
        }
        Eigen::MatrixXd proj = eigVecs * diagEigVals * eigVecs.transpose();
        for (int i=0; i<proj.rows(); i++) {
            for (int j=i; j<proj.cols(); j++) {
                xMap[optDataRecast->momentMatrix[i][j].getKey()] = proj(i, j);
            }
        }

        // Convert back to a vector
        Eigen::VectorXd xPos = Eigen::VectorXd::Zero(x.size());
        for (int i=0; i<xMap.size(); i++) {
            xPos(i) = xMap[optDataRecast->varList[i]].real();
        }
        Eigen::VectorXd errVec = optDataRecast->A*xPos - optDataRecast->b;
        double errorLin = errVec.norm();

        // Linear projection
        Eigen::VectorXd projX = optDataRecast->linSolver->solveWithGuess(optDataRecast->b, xPos);

        // Reset the gradient
        *gradOut = Eigen::VectorXd::Zero(x.size());

        // Set the gradient
        for (int i=0; i<gradOut->size(); i++) {
            (*gradOut)(i) = x(i) - projX(i);
        }

        // Early convergence
        if (errorLin < 1e-7 && minEig > -1e-7) {
            *gradOut = Eigen::VectorXd::Zero(x.size());
        }

        // Per iteration output
        std::cout << "obj=" << obj << "  lin=" << errorLin << "  eig=" << minEig << "       \r" << std::flush;

    }

    // Return the cost
    return cost;

}

// Cost/gradient function for optim TODO
static double gradFunctionFlipped(const Eigen::VectorXd& x, Eigen::VectorXd* gradOut, void* optData) {
    
    // Recast this generic pointer into the correct format
    optimData* optDataRecast = reinterpret_cast<optimData*>(optData);

    // Convert the x vector into a map
    std::map<Mon, std::complex<double>> xMap;
    for (int i=0; i<x.size(); i++) {
        xMap[optDataRecast->varList[i]] = x(i);
    }

    // Objective value
    double obj = -optDataRecast->objective.eval(xMap).real();

    // Calculate linear error
    double errorLin = (optDataRecast->A*x - optDataRecast->b).norm();

    // Calculate the eigenspectrum
    Eigen::MatrixXd X = replaceVariables(optDataRecast->momentMatrix, xMap).real();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(X);
    Eigen::VectorXd eigVals = es.eigenvalues().real();
    double minEig = eigVals.minCoeff();

    // The cost function to minimize
    double cost = errorLin*errorLin + minEig*minEig;

    // Calculate the gradient
    if (gradOut) {

        // Linear projection
        Eigen::VectorXd projX = optDataRecast->linSolver->solveWithGuess(optDataRecast->b, x);

        // Convert back to a map
        for (int i=0; i<x.size(); i++) {
            xMap[optDataRecast->varList[i]] = projX(i);
        }

        // SDP projection
        X = replaceVariables(optDataRecast->momentMatrix, xMap).real();
        es.compute(X);
        eigVals = es.eigenvalues().real();
        Eigen::MatrixXd eigVecs = es.eigenvectors().real();
        Eigen::MatrixXd diagEigVals = Eigen::MatrixXd::Zero(X.rows(), X.cols());
        for (int i=0; i<eigVals.size(); i++) {
            diagEigVals(i, i) = std::max(eigVals(i), 0.0);
        }
        Eigen::MatrixXd proj = eigVecs * diagEigVals * eigVecs.transpose();
        for (int i=0; i<proj.rows(); i++) {
            for (int j=i; j<proj.cols(); j++) {
                xMap[optDataRecast->momentMatrix[i][j].getKey()] = proj(i, j);
            }
        }

        // Convert back to a vector
        Eigen::VectorXd xPos = Eigen::VectorXd::Zero(x.size());
        for (int i=0; i<xMap.size(); i++) {
            xPos(i) = xMap[optDataRecast->varList[i]].real();
        }

        // Reset the gradient
        *gradOut = Eigen::VectorXd::Zero(x.size());

        // Set the gradient
        for (int i=0; i<gradOut->size(); i++) {
            (*gradOut)(i) = x(i) - xPos(i);
        }

        // Early convergence
        if (errorLin < 1e-7 && minEig > -1e-7) { 
            *gradOut = Eigen::VectorXd::Zero(x.size());
        }

        // Per iteration output
        std::cout << "obj=" << obj << "  lin=" << errorLin << "  eig=" << minEig << "       \r" << std::flush;

    }

    // Return the cost
    return cost;

}

// Generic entry function
int main(int argc, char* argv[]) {

    // Start a timer
    std::chrono::steady_clock::time_point timeStart = std::chrono::steady_clock::now();

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

        // Randomized version for arbitrary number of inputs
        //
        // for -S 1 -RXX22 100
        // SCS 34s 1044
        // PRJ 28s 1068 2839i
        // LBFGS 1s 1064
        // CEN  0s 1137
        //
        // for -S 1 -RXX22 200
        // SCS 61m 3009
        // PRJ  7m 3069 4507i
        // LBFGS 3s 3064
        // CEN  0m 3274
        //
        // time ./run -S 1 --RXX22 100 -l 1 -D
        // MOSEK 1m45s 1044 8GB
        // SCS 57s 1044 38MB
        // LBFGS 0.8s 1064 78MB
        // CEN 0s 1137
        //
        } else if (argAsString == "--RXX22" || argAsString == "--rxx22") {
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

    // Start a timer
    std::chrono::steady_clock::time_point timeFinishedGenerating = std::chrono::steady_clock::now();

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

    // Linearized SDP
    if (testing == 6) {

        // Starting output
        if (verbosity >= 1) {
            std::cout << "Moment matrix has size " << momentMatrices[0].size() << std::endl;
            std::cout << "Starting linearization..." << std::endl;
        }

        // Generate the random vectors
        if (verbosity >= 1) {
            std::cout << "Generating random vectors..." << std::endl;
        }
        int numVecs = int(stepSize);
        Eigen::MatrixXd randomVecs = Eigen::MatrixXd::Random(numVecs, momentMatrices[0].size());

        // Get an easy bound by setting vars to the identity * minimum eigen
        std::map<Mon, std::complex<double>> varVals0;
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
        X0 = replaceVariables(momentMatrices[0], varVals0).real();
        es.compute(X0);
        Eigen::MatrixXd X0EigVecs = es.eigenvectors().real();
        int nonRandomCount = 0;
        for (int i=0; i<std::min(numVecs, int(X0EigVecs.cols())); i++) {
            randomVecs.row(i) = X0EigVecs.col(i).transpose();
            nonRandomCount++;
        }
        if (nonRandomCount < numVecs) {
            for (int i=0; i<momentMatrices[0].size(); i++) {
                randomVecs.row(nonRandomCount) = Eigen::MatrixXd::Zero(1, momentMatrices[0].size());
                randomVecs(nonRandomCount, i) = 1;
                nonRandomCount++;
                if (nonRandomCount >= numVecs) {
                    break;
                }
            }
            
        }

        // For now build the full matrix as the sum of the projectors
        if (verbosity >= 1) {
            std::cout << "Building the full matrix..." << std::endl;
        }
        std::vector<std::vector<Poly>> momentMatEigens(momentMatrices[0].size(), std::vector<Poly>(momentMatrices[0].size(), Poly()));
        for (int i=0; i<numVecs; i++) {
            Eigen::MatrixXd proj = randomVecs.row(i).transpose() * randomVecs.row(i);
            Mon newVar("<L" + std::to_string(i) + ">");
            for (int j=0; j<momentMatrices[0].size(); j++) {
                for (int k=0; k<momentMatrices[0][j].size(); k++) {
                    momentMatEigens[j][k] += newVar * proj(j, k);
                }
            }
            constraintsPositive.push_back(Poly(newVar));
        }

        // Add the equality constraints
        if (verbosity >= 1) {
            std::cout << "Adding the equality constraints..." << std::endl;
        }
        for (int j=0; j<momentMatrices[0].size(); j++) {
            for (int k=0; k<momentMatrices[0][j].size(); k++) {
                constraintsZero.push_back(momentMatEigens[j][k] - momentMatrices[0][j][k]);
            }
        }

        // Set the objective
        //if (verbosity >= 1) {
            //std::cout << "Setting the objective..." << std::endl;
        //}
        //Poly newObj;
        //std::set<Mon> monsUsed;
        //for (int j=0; j<momentMatrices[0].size(); j++) {
            //for (int k=0; k<j; k++) {
                //Mon key = momentMatrices[0][j][k].getKey();
                //if (!monsUsed.count(key)) {
                    //newObj += momentMatEigens[j][k] * momentMatrices[0][j][k].getValue();
                    //monsUsed.insert(key);
                //}
            //}
        //}
        //objective = newObj;

        // Output things if super verbose
        if (verbosity >= 3) {
            std::cout << "New moment matrix: " << std::endl;
            std::cout << momentMatEigens << std::endl;
            std::cout << "New equality constraints: " << std::endl;
            std::cout << constraintsZero << std::endl;
            std::cout << "New positivity constraints: " << std::endl;
            std::cout << constraintsPositive << std::endl;
        }

        // Solve the linear program with MOSEK
        std::map<Mon, std::complex<double>> varVals;
        momentMatrices = {};
        double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, {-100000, 100000}, &varVals);
        std::cout << "Result: " << res << std::endl;

        return 0;

    }

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
        Eigen::MatrixXd X0EigVecs = es.eigenvectors().real();
        if (verbosity >= 1) {
            std::cout << "Center bound: " << easyBound << std::endl;

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
        if (verbosity >= 1) {
            std::cout << "Starting distance: " << distance << std::endl;
        }
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

        // Start a timer
        std::chrono::steady_clock::time_point timeFinishedConverting = std::chrono::steady_clock::now();

        // Set up the linear solver
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> linSolver;
        linSolver.setTolerance(tolerance/100.0);
        linSolver.compute(A);

        // Convert to a vector
        Eigen::VectorXd x = Eigen::VectorXd::Zero(varList.size());
        for (int i=0; i<varList.size(); i++) {
            x(i) = varVals[varList[i]].real();
        }

        // Test With Optim
        optimData optData;
        optData.objective = objective;
        optData.A = A;
        optData.linSolver = &linSolver;
        optData.b = b;
        optData.momentMatrix = momentMatrices[0];
        optData.varList = varList;
        optim::algo_settings_t settings;
        settings.print_level = 0;
        settings.iter_max = maxIters;
        settings.grad_err_tol = 1e-12;
        settings.rel_sol_change_tol = 1e-12;
        settings.rel_objfn_change_tol = 1e-12;
        settings.lbfgs_settings.par_M = 6;
        settings.lbfgs_settings.wolfe_cons_1 = 1e-4;
        settings.lbfgs_settings.wolfe_cons_2 = 0.9;
        settings.gd_settings.method = 0;

        double objPrev = 10000000;
        for (int k=0; k<varList.size(); k++) {

            std::map<Mon, std::complex<double>> varValsToTest;

            // Set this variable to plus distance
            x(k) = distance;
            Eigen::VectorXd xLin1 = linSolver.solveWithGuess(b, x);
            std::cout << "testing " << varList[k] << " = " << distance << std::endl;
            bool successLin1 = optim::lbfgs(xLin1, gradFunction, &optData, settings);
            std::cout << std::endl;
            for (int i=0; i<varList.size(); i++) {
                varValsToTest[varList[i]] = xLin1(i);
            }
            double obj1 = -objective.eval(varValsToTest).real();

            // Set this variable to minus distance
            x(k) = -distance;
            Eigen::VectorXd xLin2 = linSolver.solveWithGuess(b, x);
            std::cout << "testing " << varList[k] << " = " << -distance << std::endl;
            bool successLin2 = optim::lbfgs(xLin2, gradFunction, &optData, settings);
            std::cout << std::endl;
            for (int i=0; i<varList.size(); i++) {
                varValsToTest[varList[i]] = xLin2(i);
            }
            double obj2 = -objective.eval(varValsToTest).real();

            // Set this variable back to zero
            x(k) = 0;
            Eigen::VectorXd xLin3 = linSolver.solveWithGuess(b, x);
            std::cout << "testing " << varList[k] << " = 0" << std::endl;
            bool successLin3 = optim::lbfgs(xLin3, gradFunction, &optData, settings);
            std::cout << std::endl;
            for (int i=0; i<varList.size(); i++) {
                varValsToTest[varList[i]] = xLin3(i);
            }
            double obj3 = -objective.eval(varValsToTest).real();

            double objMin = std::min(obj1, std::min(obj2, obj3));
            if (objMin == obj1) {
                x(k) = distance;
            } else if (objMin == obj2) {
                x(k) = -distance;
            } else {
                x(k) = 0;
            }

        }

        // Travel a bit in the objective direction TODO
        //int numExtra = 1000;
        //double prevObj = 10000000;
        //for (int extraIter=0; extraIter<numExtra; extraIter++) {

            //double distanceToTry = distance;
            //double bestObj = prevObj;
            //double bestDist = distanceToTry;
            //std::map<Mon, std::complex<double>> bestVarVals;
            //for (int j=0; j<10; j++) {

                //// Move a bit in the direction
                //std::map<Mon, std::complex<double>> varValsNew = varVals;
                //for (auto& term : objective) {
                    //if (std::abs(term.second) > 0) {
                        //varValsNew[term.first] += distanceToTry;
                    //} else {
                        //varValsNew[term.first] -= distanceToTry;
                    //}
                //}
                //Eigen::VectorXd xNew = Eigen::VectorXd::Zero(varList.size());
                //for (int i=0; i<varList.size(); i++) {
                    //xNew(i) = varValsNew[varList[i]].real();
                //}

                //// The last iteration should be longer
                //if (extraIter == numExtra-1) {
                    //settings.iter_max = 1000000;
                //} else {
                    //settings.iter_max = maxIters;
                //}

                //// Solve
                //Eigen::VectorXd projX = linSolver.solveWithGuess(b, xNew);
                //bool success = optim::lbfgs(projX, gradFunction, &optData, settings);
                //for (int i=0; i<varList.size(); i++) {
                    //varValsNew[varList[i]] = projX(i);
                //}

                //// Calculate the objective
                //double newObj = -objective.eval(varValsNew).real();
                //std::cout << "distance=" << distanceToTry << " obj=" << newObj << "  best=" << bestObj << std::endl;

                //// See if this is the best so far
                //if (newObj < bestObj) {
                    //bestObj = newObj;
                    //bestVarVals = varValsNew;
                    //bestDist = 2*distanceToTry;
                    //if (j != 0) {
                        //break;
                    //}
                //}

                //// Reduce the distance and try again
                //distanceToTry *= 0.5;
                //if (distanceToTry < 1e-6) {
                    //break;
                //}

            //}
            //std::cout << std::endl;

            ////if (std::abs(bestObj - prevObj) < 1e-4) {
                ////break;
            ////}

            //if (distanceToTry > 1e-6) {
                //varVals = bestVarVals;
                //prevObj = bestObj;
                //distance = bestDist;
            //} else {
                //break;
            //}

            //// Adjust the distance
            ////if (newObj >= prevObj || std::abs(newObj - prevObj) < 1e-4) {
                ////distance *= 0.5;
            ////}
            ////prevObj = newObj;

        //}

        // Verbose output of the final solution
        if (verbosity >= 1) {
            std::cout << std::endl;
            int numIters = settings.opt_iter;
            std::cout << "Iterations needed: " << numIters << std::endl;
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

        // Start a timer
        std::chrono::steady_clock::time_point timeFinishedSolving = std::chrono::steady_clock::now();

        // Output timings
        if (verbosity >= 1) {
            int timeToGen = std::chrono::duration_cast<std::chrono::seconds>(timeFinishedGenerating - timeStart).count();
            int timeToConvert = std::chrono::duration_cast<std::chrono::seconds>(timeFinishedConverting - timeFinishedGenerating).count();
            int timeToSolve = std::chrono::duration_cast<std::chrono::seconds>(timeFinishedSolving - timeFinishedConverting).count();
            if (timeToGen == 0) {
                timeToGen = std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedGenerating - timeStart).count();
                std::cout << "Time to generate: " << timeToGen << "ms" << std::endl;
            } else {
                std::cout << "Time to generate: " << timeToGen << "s" << std::endl;
            }
            if (timeToConvert == 0) {
                timeToConvert = std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedConverting - timeFinishedGenerating).count();
                std::cout << "Time to convert: " << timeToConvert << "ms" << std::endl;
            } else {
                std::cout << "Time to convert: " << timeToConvert << "s" << std::endl;
            }
            if (timeToSolve == 0) {
                timeToSolve = std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedSolving - timeFinishedConverting).count();
                std::cout << "Time to solve: " << timeToSolve << "ms" << std::endl;
            } else {
                std::cout << "Time to solve: " << timeToSolve << "s" << std::endl;
            }
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

    // Start a timer
    std::chrono::steady_clock::time_point timeFinishedSolving = std::chrono::steady_clock::now();

    // Output timings
    if (verbosity >= 1) {
        int timeToGen = std::chrono::duration_cast<std::chrono::seconds>(timeFinishedGenerating - timeStart).count();
        int timeToSolve = std::chrono::duration_cast<std::chrono::seconds>(timeFinishedSolving - timeFinishedGenerating).count();
        if (timeToGen == 0) {
            timeToGen = std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedGenerating - timeStart).count();
            std::cout << "Time to generate: " << timeToGen << "ms" << std::endl;
        } else {
            std::cout << "Time to generate: " << timeToGen << "s" << std::endl;
        }
        if (timeToSolve == 0) {
            timeToSolve = std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedSolving - timeFinishedGenerating).count();
            std::cout << "Time to solve: " << timeToSolve << "ms" << std::endl;
        } else {
            std::cout << "Time to solve: " << timeToSolve << "s" << std::endl;
        }
    }

    // Exit without errors
    return 0;

}
