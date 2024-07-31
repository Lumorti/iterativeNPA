// Standard includesSCS
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>

// Import OpenMP
#include <omp.h>

// Import Eigen
#include <Eigen/Dense>

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
    double tolerance = 1e-8;
    int numExtra = 0;
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

        // If told to use a specific solver
        } else if (argAsString == "-s") {
            if (std::string(argv[i+1]) == "s") {
                solver = "SCS";
            } else if (std::string(argv[i+1]) == "m") {
                solver = "MOSEK";
            } else if (std::string(argv[i+1]) == "o") {
                solver = "Optim";
            }
            i++;

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

        // If setting num extra
        } else if (argAsString == "-E") {
            numExtra = std::stoi(argv[i+1]);
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
            std::cout << "  -E <int>        Number of extra iterations to do" << std::endl;
            std::cout << "  -h              Display this help message" << std::endl;
            std::cout << "  -D              Use the dual of the problem" << std::endl;
            std::cout << "  -s S            Use SCS as the solver" << std::endl;
            std::cout << "  -s M            Use MOSEK as the solver" << std::endl;
            std::cout << "  -s o            Use Optim as the solver" << std::endl;
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

    // Do level 1, then 2, then 3 TODO
    if (testing == 2) {

        std::map<Mon, std::complex<double>> startVals;
        double res = 0;

        // Level 1
        momentMatrices = generateAllMomentMatrices(bellFunc, {}, 1, verbosity);
        objective = bellFunc;
        constraintsZero = {};
        constraintsPositive = {};
        primalToDual(objective, momentMatrices, constraintsZero, constraintsPositive);
        res = solveOptim(objective, constraintsZero, momentMatrices, startVals, verbosity, maxIters, numExtra, stepSize, tolerance);

        // Level 2
        momentMatrices = generateAllMomentMatrices(bellFunc, {}, 2, verbosity);
        objective = bellFunc;
        constraintsZero = {};
        constraintsPositive = {};
        primalToDual(objective, momentMatrices, constraintsZero, constraintsPositive);
        res = solveOptim(objective, constraintsZero, momentMatrices, startVals, verbosity, maxIters, numExtra, stepSize, tolerance);

        // Level 3
        momentMatrices = generateAllMomentMatrices(bellFunc, {}, 3, verbosity);
        objective = bellFunc;
        constraintsZero = {};
        constraintsPositive = {};
        primalToDual(objective, momentMatrices, constraintsZero, constraintsPositive);
        res = solveOptim(objective, constraintsZero, momentMatrices, startVals, verbosity, maxIters, numExtra, stepSize, tolerance);

        return 0;

    }

    // Solve the problem
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
    } else if (solver == "Optim") {
        std::cout << "Solving with Optim..." << std::endl;
        std::map<Mon, std::complex<double>> startVals;
        res = solveOptim(objective, constraintsZero, momentMatrices, startVals, verbosity, maxIters, numExtra, stepSize, tolerance);
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
        if (timeToGen <= 1) {
            timeToGen = std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedGenerating - timeStart).count();
            std::cout << "Time to generate: " << timeToGen << "ms" << std::endl;
        } else {
            std::cout << "Time to generate: " << timeToGen << "s" << std::endl;
        }
        if (timeToSolve <= 1) {
            timeToSolve = std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedSolving - timeFinishedGenerating).count();
            std::cout << "Time to solve: " << timeToSolve << "ms" << std::endl;
        } else {
            std::cout << "Time to solve: " << timeToSolve << "s" << std::endl;
        }
    }

    // Exit without errors
    return 0;

}
