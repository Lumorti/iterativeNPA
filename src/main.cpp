// Standard includesSCS
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>

// Import OpenMP
#include <omp.h>

// Import Eigen
#include <Eigen/Dense>

// PolyNC
#include "../../PolyNC/src/polync.h"

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

// Convert an integer to a character
std::string intToChar(int i) {
    if (i < 26) {
        return std::string(1, 'A' + i);
    } else {
        return std::string(1, 'a' + i - 26);
    }
}

// Get a submatrix
std::vector<std::vector<Poly>> fromIndices(std::vector<std::vector<Poly>>& mat, std::vector<int> indices) {
    std::vector<std::vector<Poly>> res;
    for (int i=0; i<indices.size(); i++) {
        std::vector<Poly> newRow;
        for (int j=0; j<indices.size(); j++) {
            newRow.push_back(mat[indices[i]][indices[j]]);
        }
        res.push_back(newRow);
    }
    return res;
}

// Round to +1 or -1
int roundTo1(double x) {
    return x >= 0 ? 1 : -1;
}

// From https://stackoverflow.com/questions/28711797/generating-n-choose-k-permutations-in-c
template <typename T>
void combinations(const std::vector<T>& v, std::size_t count) {
    assert(count <= v.size());
    std::vector<bool> bitset(v.size() - count, 0);
    bitset.resize(v.size(), 1);
    do {
        for (std::size_t i = 0; i != v.size(); ++i) {
            if (bitset[i]) {
                std::cout << v[i] << " ";
            }
        }
        std::cout << std::endl;
    } while (std::next_permutation(bitset.begin(), bitset.end()));
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
    bool envelopes = false;
    Poly bellFunc("<A1B1>+<A1B2>+<A2B1>-<A2B2>");
    int testing = 0;
    int verbosity = 1;
    int maxIters = 10000000;
    int maxMoments = 10000000;
    int numCores = 1;
    bool use01 = false;
    double stepSize = 100;
    double tolerance = 1e-8;
    int numExtra = 0;
    bool generateData = false;
    std::string solver = "MOSEK";
    std::string seed = "";
    std::string problemName = "CHSH";
    std::vector<std::string> extraMoments;
    std::string forcedTopRow = "";

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

        // If told to generate data
        } else if (argAsString == "-d") {
            generateData = true;
            verbosity = 0;

        // If adding an extra moment to the top row
        } else if (argAsString == "-e") {
            extraMoments.push_back(std::string(argv[i+1]));
            i++;

        // Set the iteration limit
        } else if (argAsString == "-i") {
            maxIters = std::stoi(argv[i+1]);
            i++;

        // Set the max moment limit
        } else if (argAsString == "-m") {
            maxMoments = std::stoi(argv[i+1]);
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
            std::cout << "  --RXX22 <int>   Use the randomized scenario with arbitrary number of inputs" << std::endl;
            std::cout << "  --ising <int>   A randomized fully-connected Ising model" << std::endl;
            std::cout << "  -O <str>        Use a custom objective" << std::endl;
            std::cout << "  -V              Use McCormick envelopes" << std::endl;
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
            std::cout << "  -m <int>        Maximum size of moment matrix" << std::endl;
            std::cout << "  -M <str>        Use a specific top row for the moment matrix" << std::endl;
            std::cout << "  -d              Generate data" << std::endl;
            std::cout << "  -h              Display this help message" << std::endl;
            std::cout << "  -D              Use the dual of the problem" << std::endl;
            std::cout << "  -s S            Use SCS as the solver" << std::endl;
            std::cout << "  -s M            Use MOSEK as the solver" << std::endl;
            std::cout << "  -s o            Use Optim as the solver" << std::endl;
            std::cout << "  -A <dbl>        Set the step size" << std::endl;
            std::cout << "  -01             Use 0/1 output instead of -1/1" << std::endl;
            return 0;

        // If using a custom objective
        } else if (argAsString == "-O") {
            bellFunc = Poly(argv[i+1]);
            problemName = "CUSTOM";
            i++;

        // If using a specific top row
        } else if (argAsString == "-M") {
            forcedTopRow = std::string(argv[i+1]);
            i++;

        // If setting the tolerance
        } else if (argAsString == "-T") {
            tolerance = std::stod(argv[i+1]);
            i++;

        // If using the dual
        } else if (argAsString == "-D") {
            useDual = true;

        // If using envelopes
        } else if (argAsString == "-V") {
            envelopes = true;

        // If told to use the Ising model
        } else if (argAsString == "--ising") {
            int numInputs = std::stoi(argv[i+1]);
            problemName = "ISING";
            if (numInputs > 50) {
                std::cout << "Error - too big a model, can only have up to 50" << std::endl;
                return 1;
            }
            bellFunc = Poly();
            for (int i=0; i<numInputs; i++) {
                bellFunc += Poly("<" + intToChar(i) + "1>");
            }
            for (int i=0; i<numInputs; i++) {
                for (int j=i+1; j<numInputs; j++) {
                    bellFunc += Poly("<" + intToChar(i) + "1" + intToChar(j) + "1>");
                }
            }
            bellFunc.randomize();
            i++;

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

    // If we should use a specific top row
    if (forcedTopRow != "") {

        // Split the string
        std::vector<std::string> monStrings;
        std::string currentMon = "";
        for (int i=0; i<forcedTopRow.size(); i++) {
            if (forcedTopRow[i] == ' ') {
                monStrings.push_back(currentMon);
                currentMon = "";
            } else if (forcedTopRow[i] != '+' && forcedTopRow[i] != '-') {
                currentMon += forcedTopRow[i];
            }
        }
        if (currentMon.size() > 0) {
            monStrings.push_back(currentMon);
        }

        // Add the inverse of each
        int ogSize = monStrings.size();
        for (int i=0; i<ogSize; i++) {
            std::string rev;
            for (int j=monStrings[i].size()-2; j>=0; j-=2) {
                rev += monStrings[i][j];
                rev += monStrings[i][j+1];
            }
            if (rev.size() > 0) {
                if (std::find(monStrings.begin(), monStrings.end(), rev) == monStrings.end()) {
                    monStrings.push_back(rev);
                }
            }
        }

        // Add braket to each
        for (int i=0; i<monStrings.size(); i++) {
            if (monStrings[i][0] != '<' && monStrings[i] != "1") {
                monStrings[i] = "<" + monStrings[i] + ">";
            }
        }

        // Print the top row
        for (int i=0; i<monStrings.size(); i++) {
            std::cout << monStrings[i] << " ";
        }
        std::cout << std::endl;

        // Convert to polynomials
        std::vector<Poly> topRow = {};
        for (int i=0; i<monStrings.size(); i++) {
            topRow.push_back(Poly(monStrings[i]));
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

    // Manual paritioning of the moment matrix 
    if (testing == -1) {
        
        std::vector<int> indexList;

        indexList = {};
        for (int i=0; i<28; i++) {
            indexList.push_back(i);
        }
        std::vector<std::vector<Poly>> mat1 = fromIndices(momentMatrices[0], indexList);

        indexList = {};
        for (int i=0; i<12; i++) {
            indexList.push_back(i);
        }
        for (int i=27; i<43; i++) {
            indexList.push_back(i);
        }
        std::vector<std::vector<Poly>> mat2 = fromIndices(momentMatrices[0], indexList);

        indexList = {};
        for (int i=0; i<12; i++) {
            indexList.push_back(i);
        }
        for (int i=43; i<59; i++) {
            indexList.push_back(i);
        }
        std::vector<std::vector<Poly>> mat3 = fromIndices(momentMatrices[0], indexList);

        indexList = {};
        for (int i=0; i<12; i++) {
            indexList.push_back(i);
        }
        for (int i=59; i<75; i++) {
            indexList.push_back(i);
        }
        std::vector<std::vector<Poly>> mat4 = fromIndices(momentMatrices[0], indexList);

        momentMatrices = {mat1, mat2, mat3, mat4};

    }

    // Define other constraints
    Poly objective = bellFunc;
    std::vector<Poly> constraintsZero;
    std::vector<Poly> constraintsPositive;

    // If using McCormick envelopes TODO
    if (envelopes) {

        // Count how many unique As and Bs there are
        int numA = 0;
        int numB = 0;
        for (auto& term : bellFunc) {
            for (int i=0; i<term.first.size(); i++) {
                if (term.first[i].first == 'A') {
                    numA = std::max(numA, term.first[i].second);
                } else if (term.first[i].first == 'B') {
                    numB = std::max(numB, term.first[i].second);
                }
            }
        }

        // For each pair of operators
        std::vector<Mon> opListA;
        std::vector<Mon> opListB;
        for (int i=1; i<=numA; i++) {
            opListA.push_back(Mon("<A" + std::to_string(i) + ">"));
        }
        for (int i=1; i<=numB; i++) {
            opListB.push_back(Mon("<B" + std::to_string(i) + ">"));
        }
        for (int i=0; i<opListA.size(); i++) {
            for (int j=0; j<opListB.size(); j++) {

                // Get the operators
                Poly A = opListA[i];
                Poly B = opListB[j];
                Poly AB = A * B;

                // (A + I) (B + I) >= 0
                // <AB> + <A> + <B> + 1 >= 0
                constraintsPositive.push_back(AB + A + B + 1);

                // (A + I) (I - B) >= 0
                // <A> - <AB> + 1 - <B> >= 0
                constraintsPositive.push_back(A - AB + 1 - B);

                // (I - A) (B + I) >= 0
                // <B> + 1 - <AB> - <A> >= 0
                constraintsPositive.push_back(B + 1 - AB - A);

                // (I - A) (I - B) >= 0
                // 1 - <A> - <B> + <AB> >= 0
                constraintsPositive.push_back(1 - A - B + AB);

            }
        }

    }

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

    // If told to generate data, keep trying sets of moments TODO
    if (generateData) {

        // Possible list of monoms to choose from
        std::vector<Mon> monomialsInProblemVec;
        addSingleMonomials(monomialsInProblemVec, bellFunc);
        std::vector<Poly> possibleMonomials = generateMonomials(monomialsInProblemVec, level, 0);
        std::vector<Poly> topRow = {Poly(1)};
        for (int i=0; i<possibleMonomials.size(); i++) {
            if (possibleMonomials[i].size() > 0 && possibleMonomials[i].getKey().size() == 1) {
                topRow.push_back(possibleMonomials[i]);
            }
        }
        std::map<Mon, std::complex<double>> varVals;
        std::pair<int,int> varBounds = {-1, 1};
        int amountExtra = std::min(maxMoments, int(possibleMonomials.size()));
        std::vector<Poly> topRowL1 = topRow;

        // Iterate over all sets of moments
        std::vector<Poly> v;
        for (int i=0; i<possibleMonomials.size(); i++) {
            if (std::find(topRow.begin(), topRow.end(), possibleMonomials[i]) == topRow.end()) {
                v.push_back(possibleMonomials[i]);
            }
        }
        int count = amountExtra;
        std::vector<bool> bitset(v.size() - count, false);
        for (int i=0; i<count; i++) {
            bitset.push_back(true);
        }
        do {

            // Generate the moment matrix
            topRow = topRowL1;
            for (std::size_t i = 0; i != v.size(); ++i) {
                if (bitset[i]) {
                    topRow.push_back(v[i]);
                }
            }
            momentMatrices[0] = generateFromTopRow(topRow, use01);

            // Solve the problem
            double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varBounds, 0, &varVals);
            std::cout << objective;
            std::cout << " | ";
            for (int j=0; j<topRow.size(); j++) {
                std::cout << topRow[j] << " ";
            }
            std::cout << "| ";
            std::cout << topRow.size();
            std::cout << " | ";
            std::cout << std::setprecision(10) << res << std::endl;

        } while (std::next_permutation(bitset.begin(), bitset.end()));

        // Iterate over all sets of moments
        //std::vector<bool> mask(possibleMonomials.size(), false);
        //for (int i=0; i<amount; i++) {
            //mask[mask.size()-1-i] = true;
        //}
        //std::vector<int> indices;
        //for (int j=0; j<possibleMonomials.size(); j++) {
            //indices.push_back(j);
        //}
        //do {

            //// Add this permutation to the moment matrix
            //topRow = {Poly(1)};
            //for (int j=0; j<amount; j++) {
                ////topRow.push_back(possibleMonomials[indices[j]]);
                //if (mask[j]) {
                    //topRow.push_back(possibleMonomials[indices[j]]);
                //}
            //}
            //momentMatrices[0] = generateFromTopRow(topRow, use01);

            //// Solve the problem
            //double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varBounds, &varVals);
            //std::cout << objective;
            //std::cout << " | ";
            //for (int j=0; j<topRow.size(); j++) {
                //std::cout << topRow[j] << " ";
            //}
            //std::cout << "| ";
            //std::cout << topRow.size();
            //std::cout << " | ";
            //std::cout << std::setprecision(10) << res << std::endl;

        //// Permute the index list
        //} while (std::next_permutation(indices.begin(), indices.end()));
        //} while (std::prev_permutation(mask.begin(), mask.end()));

        return 0;

    }

    // Represent it as a linear combination of random positive matrices TODO
    if (testing >= 6) {

        // Generate the list of positive matrices
        int numMats = testing;
        int matSize = momentMatrices[0].size();
        std::vector<Eigen::MatrixXd> positiveMats;
        for (int i=0; i<numMats; i++) {
            Eigen::VectorXd v = Eigen::VectorXd::Random(matSize);
            Eigen::MatrixXd mat = v * v.transpose();
            positiveMats.push_back(mat);
        }
        //for (int i=0; i<numMats; i++) {
            //Eigen::VectorXd v = Eigen::VectorXd::Random(matSize);
            //for (int j=0; j<matSize; j++) {
                //if (v(j) < 0.5) {
                    //v(j) = 0;
                //} else {
                    //v(j) = 1;
                //}
            //}
            //Eigen::MatrixXd mat = v * v.transpose();
            //positiveMats.push_back(mat);
        //}
        //for (int i=0; i<matSize; i++) {
            //Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(matSize, matSize);
            //mat(i, i) = 1;
            //positiveMats.push_back(mat);
        //}
        if (verbosity >= 3) {
            for (int i=0; i<positiveMats.size(); i++) {
                std::cout << positiveMats[i] << std::endl;
                std::cout << std::endl;
            }
        }

        // Generate the moment matrix as a linear combination of these
        std::vector<std::vector<Poly>> newMomentMatrix(matSize, std::vector<Poly>(matSize, Poly(0)));
        for (int i=0; i<positiveMats.size(); i++) {
            for (int j=0; j<matSize; j++) {
                for (int k=0; k<matSize; k++) {
                    newMomentMatrix[j][k] += Poly("<L" + std::to_string(i) + ">") * Poly(positiveMats[i](j, k));
                }
            }
        }

        // Equate these two, element-wise, as new linear zero constraints
        for (int i=0; i<matSize; i++) {
            for (int j=i; j<matSize; j++) {
                constraintsZero.push_back(newMomentMatrix[i][j] - momentMatrices[0][i][j]);
            }
        }
        momentMatrices = {};
        if (verbosity >= 3) {
            std::cout << "Linear zero constraints: " << std::endl;
            for (int i=0; i<constraintsZero.size(); i++) {
                std::cout << constraintsZero[i] << std::endl;
            }
        }

        // All the L's are positive
        for (int i=0; i<positiveMats.size(); i++) {
            constraintsPositive.push_back(Poly("<L" + std::to_string(i) + ">"));
        }

        // Solve the now linear program
        std::map<Mon, std::complex<double>> varVals;
        std::pair<int,int> varBounds = {-1, 1};
        double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varBounds, 0, &varVals);
        std::cout << res << std::endl;

        return 0;

    }

    // Iterative scheme TODO
    if (testing == 7) {

        // Solve this problem
        std::map<Mon, std::complex<double>> varVals;
        std::pair<int,int> varBounds = {-1, 1};
        double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varBounds, 0, &varVals);
        std::cout << res << std::endl;

        // Add the constraint that the objective function is bounded
        constraintsPositive.push_back(Poly(res) - objective);
        res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varBounds, 0, &varVals);
        std::cout << res << std::endl;

        // A1, A2, B1, B2
        //constraintsPositive.push_back(Poly(2*sqrt(2)) - Poly("<A1B1>+<A1B2>+<A2B1>-<A2B2>"));
        //constraintsPositive.push_back(Poly("<A1B1>+<A1B2>+<A2B1>-<A2B2>") - Poly(-2*sqrt(2)));

        // A1, A2, B2, B3
        //constraintsPositive.push_back(Poly(2*sqrt(2)) - Poly("<A1B3>+<A1B2>+<A2B3>-<A2B2>"));
        //constraintsPositive.push_back(Poly("<A1B3>+<A1B2>+<A2B3>-<A2B2>") - Poly(-2*sqrt(2)));

        // A1, A2, B1, B3
        //constraintsPositive.push_back(Poly(2*sqrt(2)) - Poly("<A1B1>+<A1B3>+<A2B1>-<A2B3>"));
        //constraintsPositive.push_back(Poly("<A1B1>+<A1B3>+<A2B1>-<A2B3>") - Poly(-2*sqrt(2)));

        // A1, A3, B2, B3
        //constraintsPositive.push_back(Poly(2*sqrt(2)) - Poly("<A1B3>+<A1B2>+<A3B3>-<A3B2>"));
        //constraintsPositive.push_back(Poly("<A1B3>+<A1B2>+<A3B3>-<A3B2>") - Poly(-2*sqrt(2)));

        // A2, A3, B1, B2
        //constraintsPositive.push_back(Poly(2*sqrt(2)) - Poly("<A2B1>+<A2B2>+<A3B1>-<A3B2>"));
        //constraintsPositive.push_back(Poly("<A2B1>+<A2B2>+<A3B1>-<A3B2>") - Poly(-2*sqrt(2)));

        // Individually bound each variable
        for (auto& pair : varVals) {
            if (pair.first.size() == 0) {
                continue;
            }
            std::cout << "For variable " << pair.first << std::endl;
            Poly newObj = Poly(pair.first);
            res = solveMOSEK(newObj, momentMatrices, constraintsZero, constraintsPositive, 0, varBounds, 0, &varVals);
            std::cout << "max: " << res << std::endl;
            Poly newCon = Poly(res) - newObj;
            std::cout << "adding constraint: " << newCon << " >= 0" << std::endl;
            constraintsPositive.push_back(newCon);
            newObj = Poly(-1, pair.first);
            res = solveMOSEK(newObj, momentMatrices, constraintsZero, constraintsPositive, 0, varBounds, 0, &varVals);
            std::cout << "min: " << -res << std::endl;
            newCon = -newObj - Poly(-res);
            constraintsPositive.push_back(newCon);
            std::cout << "adding constraint: " << newCon << " >= 0" << std::endl;
        }

        // Bound each positivity constraint
        for (int i=0; i<constraintsPositive.size(); i++) {
            std::cout << "For constraint " << constraintsPositive[i] << std::endl;
            res = solveMOSEK(constraintsPositive[i], momentMatrices, constraintsZero, constraintsPositive, 0, varBounds, 0, &varVals);
            std::cout << "max: " << res << std::endl;
            Poly newCon = Poly(res) - constraintsPositive[i];
            constraintsPositive.push_back(newCon);
            std::cout << "adding constraint: " << newCon << " >= 0" << std::endl;
            res = solveMOSEK(-constraintsPositive[i], momentMatrices, constraintsZero, constraintsPositive, 0, varBounds, 0, &varVals);
            std::cout << "min: " << -res << std::endl;
            newCon = -constraintsPositive[i] - Poly(-res);
            constraintsPositive.push_back(newCon);
            std::cout << "adding constraint: " << newCon << " >= 0" << std::endl;
        }

        // Solve this problem
        res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varBounds, 0, &varVals);
        std::cout << "new i3322: " << res << std::endl;

        return 0;

    }

    // Try doing smaller submatrices at a time
    if (testing == 4) {

        // Get a list of all the variables
        std::set<Mon> monomialList;
        addVariables(monomialList, momentMatrices[0]);
        std::set<Mon> monomsInObjective;
        addVariables(monomsInObjective, objective);

        // Init the var list
        std::map<Mon, std::complex<double>> currentVals;
        for (auto& mon : monomialList) {
            currentVals[mon] = 0;
        }

        // For some iteration count
        for (int i=0; i<maxIters; i++) {

            // Pick a variables to optimize over
            //std::set<Mon> varsToOptimize = monomialList;
            std::set<Mon> varsToOptimize = monomsInObjective;
            int numVars = 7;
            while (varsToOptimize.size() < numVars + monomsInObjective.size()) {
                Mon var = *std::next(monomialList.begin(), rand(0, monomialList.size()-1));
                if (varsToOptimize.find(var) == varsToOptimize.end()) {
                    varsToOptimize.insert(var);
                }
            }

            // Randomly remove some
            //int toRemove = 10;
            //for (int j=0; j<toRemove; j++) {
                //Mon var = *std::next(varsToOptimize.begin(), rand(0, varsToOptimize.size()-1));
                //varsToOptimize.erase(var);
            //}

            // Fix everything else in the moment matrix
            std::vector<std::vector<Poly>> fixedMomentMatrix = momentMatrices[0];
            for (int i=0; i<fixedMomentMatrix.size(); i++) {
                for (int j=0; j<fixedMomentMatrix[i].size(); j++) {
                    Mon mon = fixedMomentMatrix[i][j].getKey();
                    if (varsToOptimize.find(mon) == varsToOptimize.end()) {
                        fixedMomentMatrix[i][j] = Poly(currentVals[mon]);
                    }
                }
            }

            // Solve this problem
            std::map<Mon, std::complex<double>> varVals;
            std::pair<int,int> varBounds = {-1, 1};
            std::vector<std::vector<std::vector<Poly>>> fixedMats = {fixedMomentMatrix};
            double res = solveMOSEK(objective, fixedMats, constraintsZero, constraintsPositive, verbosity, varBounds, 0, &varVals);
            std::cout << res << std::endl;

            // Set the variables
            for (auto& pair : varVals) {
                currentVals[pair.first] = pair.second;
            }

        }

        return 0;

    }

    // Try to turn the SDP into a linear program
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
            double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, {-1, 1}, 0, &xMap);

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

            if (eigVals.minCoeff() > -tolerance) {
                break;
            }

            // For each negative eigenvalue, store the eigenvector
            std::vector<Eigen::VectorXd> negEigVecs;
            for (int i=0; i<eigVals.size(); i++) {
                if (eigVals(i) < -tolerance) {
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

    // Do level 1, then 2, then 3
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

    // Same as below, but just checking eigenvalues
    if (testing == 8) {

        std::vector<Mon> monomialsInProblemVec;
        addSingleMonomials(monomialsInProblemVec, bellFunc);
        std::vector<Poly> possibleMonomials = generateMonomials(monomialsInProblemVec, level, 0);
        std::vector<Poly> topRow = {Poly(1)};
        for (int i=0; i<possibleMonomials.size(); i++) {
            topRow.push_back(possibleMonomials[i]);
        }
        momentMatrices[0] = generateFromTopRow(topRow, use01);
        std::map<Mon, std::complex<double>> varVals;
        double resOG = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, {-1, 1}, 0, &varVals);
        std::cout << "Initial: " << resOG << std::endl;

        // Move a small distance in the direction of the objective
        std::map<Mon, std::complex<double>> varValsNew = varVals;
        double dist = tolerance;
        for (auto& pair : objective) {
            varValsNew[pair.first] += dist * pair.second;
        }

        std::vector<std::set<int>> queue;
        for (int i=1; i<topRow.size(); i++) {
            queue.push_back({i});
        }

        // For each thing in the queue, try removing those
        std::set<std::set<int>> tried;
        std::vector<Poly> bestTopRow = topRow;
        while (queue.size() > 0) {

            // Get the next thing
            std::cout << "Queue size: " << queue.size() << std::endl;
            std::set<int> toRemove = queue.back();
            queue.pop_back();

            // Generate the new top row
            std::vector<Poly> newTopRow;
            for (int i=0; i<topRow.size(); i++) {
                if (toRemove.find(i) == toRemove.end()) {
                    newTopRow.push_back(topRow[i]);
                }
            }

            // Generate the new moment matrix
            momentMatrices[0] = generateFromTopRow(newTopRow, use01);

            // Get the eigenvalues
            Eigen::MatrixXcd A = replaceVariables(momentMatrices[0], varValsNew);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A);
            Eigen::VectorXd eigVals = es.eigenvalues().real();

            // Solve the problem
            //std::cout << "Removed ";
            //for (auto& i : toRemove) {
                //std::cout << i << " ";
            //}
            double res = eigVals.minCoeff();
            std::cout << "Size: " << momentMatrices[0].size() << ", best: " << bestTopRow.size() << ", result: " << res << std::endl;
            for (int i=0; i<bestTopRow.size(); i++) {
                std::cout << bestTopRow[i] << " ";
            }
            std::cout << std::endl;
            tried.insert(toRemove);

            // If it's better, add the new things to the queue
            if (res < 0) {
                if (newTopRow.size() < bestTopRow.size()) {
                    bestTopRow = newTopRow;
                }
                for (int i=0; i<topRow.size(); i++) {
                    std::set<int> toAdd = toRemove;
                    toAdd.insert(i);
                    if (tried.find(toAdd) == tried.end()) {
                        queue.push_back(toAdd);
                    }
                }
            }

        }

        return 0;

    }

    // Brute force the minimal moment set TODO
    // 56 - removed from l3: 1 21 23 24 26 27 28 29 32 33 34 37 38 39 40 41 42 44 49 51 54 55 56 57 58 61 62 63 65 67 68 69 72 73 75 79 80 81 82 83 84 87 88 90 91 92 93 94 96 97 98 99 100 101 103 104 105 108 109 110 111 114 115 116 117 120 121 122 123 124 125 126 127 128 129 130 132 133 134 135 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252
    // 57 giving 5.0035: 1 <A2> <A3> <B1> <B2> <B3> <A1A2> <A1A3> <A1B1> <A1B2> <A1B3> <A2A1> <A2A3> <A2B1> <A2B2> <A2B3> <A3A1> <A3A2> <A3B1> <A3B2> <B1B2> <B1B3> <B2B1> <B2B3> <B3B1> <B3B2> <A1A2A1> <A1A2A3> <A1A2B1> <A1A2B2> <A1A2B3> <A1A3A2> <A1A3B1> <A1A3B2> <A1B1B2> <A1B1B3> <A1B2B1> <A1B2B3> <A1B3B1> <A1B3B2> <A2A1A2> <A2A1B1> <A2A1B2> <A2A1B3> <A2A3A1> <A2A3A2> <A2A3B2> <A2B1B2> <A2B2B3> <A2B3B1> <A2B3B2> <A3A1B1> <A3A1B2> <A3A2B1> <A3A2B2> <A3B1B2> <A3B2B1>
    if (testing == 9) {

        std::vector<Mon> monomialsInProblemVec;
        addSingleMonomials(monomialsInProblemVec, bellFunc);
        std::vector<Poly> possibleMonomials = generateMonomials(monomialsInProblemVec, level, 0);
        std::vector<Poly> topRow = {Poly(1)};
        for (int i=0; i<possibleMonomials.size(); i++) {
            topRow.push_back(possibleMonomials[i]);
        }
        //for (int i=0; i<momentMatrices[0].size(); i++) {
            //topRow.push_back(momentMatrices[0][i][0]);
        //}
        momentMatrices[0] = generateFromTopRow(topRow, use01);
        double resOG = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity);
        std::cout << std::setprecision(10) << "Initial: " << resOG << std::endl;

        std::vector<std::set<int>> queue;
        for (int i=1; i<topRow.size(); i++) {
            queue.push_back({i});
        }

        // For each thing in the queue, try removing those
        std::vector<Poly> bestTopRow = topRow;
        int bestSize = 100000;
        std::set<std::set<int>> tried;
        int iter = 0;
        while (queue.size() > 0) {
            iter++;

            // Get the next thing
            std::cout << "Queue size: " << queue.size() << std::endl;
            std::set<int> toRemove = queue.back();
            queue.pop_back();

            // Make sure we haven't tried this before
            if (tried.find(toRemove) != tried.end()) {
                continue;
            }

            // Generate the new top row
            std::vector<Poly> newTopRow;
            for (int i=0; i<topRow.size(); i++) {
                if (toRemove.find(i) == toRemove.end()) {
                    newTopRow.push_back(topRow[i]);
                }
            }

            // Generate the new moment matrix
            momentMatrices[0] = generateFromTopRow(newTopRow, use01);

            // Solve the problem
            std::cout << iter << " Removed ";
            for (auto& i : toRemove) {
                std::cout << i << " ";
            }
            double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, 0);
            std::cout << ", size: " << momentMatrices[0].size() << ", result: " << std::setprecision(10) << res << std::endl;
            tried.insert(toRemove);

            // If it's better, add the new things to the queue
            if (res <= resOG + 1e-10) {
                if (momentMatrices[0].size() < bestSize) {
                    bestTopRow = newTopRow;
                    bestSize = momentMatrices[0].size();
                }
                for (int i=1; i<topRow.size(); i++) {
                //for (int i=topRow.size()-1; i>0; i--) {
                    std::set<int> toAdd = toRemove;
                    toAdd.insert(i);
                    if (tried.find(toAdd) == tried.end()) {
                        queue.push_back(toAdd);
                    }
                }
            }

            // Output the best
            std::cout << "Best (" << bestSize << "): " << std::endl;
            for (int i=0; i<bestTopRow.size(); i++) {
                std::cout << bestTopRow[i] << " ";
            }
            std::cout << std::endl;


        }

        return 0;

    }

    // Annealing type approach TODO
    // 20 giving 5.00588: 1 <B2B3> <B1> <B2A3> <B2B1> <B1A2> <B3A3> <B2> <A3A1> <A2> <A3> <A1A3> <A2B2> <A1B1> <A1B3> <A2B3> <A3B1> <B2A1> <B3> <A1>
    // 5.5
    // 5.00376
    // 5.0035
    if (testing >= 10) {

        int matrixSizeLimit = testing;

        momentMatrices = generateAllMomentMatrices(bellFunc, {}, 1, verbosity);

        std::vector<Mon> monomialsInProblemVec;
        addSingleMonomials(monomialsInProblemVec, bellFunc);
        std::vector<Poly> possibleMonomials = generateMonomials(monomialsInProblemVec, level, 0);

        // Correlation matrix for each monomial
        // First val is average score when present
        // Second val is average score when not present
        std::vector<std::vector<std::pair<double,double>>> correlations(possibleMonomials.size(), std::vector<std::pair<double,double>>(possibleMonomials.size(), {0,0}));
        std::vector<std::vector<std::pair<double,double>>> correlationsAveraged(possibleMonomials.size(), std::vector<std::pair<double,double>>(possibleMonomials.size(), {0,0}));
        std::vector<std::vector<std::pair<int,int>>> correlationCounts(possibleMonomials.size(), std::vector<std::pair<int,int>>(possibleMonomials.size(), {0,0}));

        // For some number of iterations
        double bestRes = 1000000;
        std::vector<int> topRowLocs = {-1};
        std::vector<Poly> topRow = {Poly(1)};
        for (int i=0; i<matrixSizeLimit; i++) {
            topRow.push_back(*std::next(possibleMonomials.begin(), i));
            topRowLocs.push_back(i);
        }
        std::vector<Poly> topRowBest = topRow;
        for (int i=0; i<maxIters; i++) {

            std::vector<Poly> prevTopRow = topRow;

            // Average the correlations
            for (int i=0; i<correlations.size(); i++) {
                for (int j=i; j<correlations[i].size(); j++) {
                    if (correlationCounts[i][j].first > 0) {
                        correlationsAveraged[i][j].first = correlations[i][j].first / correlationCounts[i][j].first;
                    }
                    if (correlationCounts[i][j].second > 0) {
                        correlationsAveraged[i][j].second = correlations[i][j].second / correlationCounts[i][j].second;
                    }
                }
            }


            // If we're at the max matrix size
            //topRow = {Poly(1)};
            if (topRow.size() >= matrixSizeLimit) {

                // Remove some random moments
                if (i < 10000000) {
                    for (int j=0; j<2; j++) {
                        int toRemove = rand(1, topRow.size()-1);
                        topRow.erase(topRow.begin() + toRemove);
                        topRowLocs.erase(topRowLocs.begin() + toRemove);
                    }
                } else {

                    // Find the worse moment
                    int worseIndex = -1;
                    double worseScore = -10000;
                    for (int j=1; j<topRow.size(); j++) {
                        int mainLoc = topRowLocs[j];
                        double score = 0;
                        int count = 0;
                        for (int k=1; k<topRow.size(); k++) {
                            if (j == k) {
                                continue;
                            }
                            score += correlationsAveraged[mainLoc][topRowLocs[k]].first;
                            score -= correlationsAveraged[mainLoc][topRowLocs[k]].second;
                        }
                        if (score > worseScore) {
                            worseScore = score;
                            worseIndex = j;
                        }
                    }
                    if (worseIndex != -1) {
                        std::cout << "Removing " << worseIndex << " (" << topRow[worseIndex] << ")" << std::endl;
                        topRow.erase(topRow.begin() + worseIndex);
                        topRowLocs.erase(topRowLocs.begin() + worseIndex);
                    }

                }

            }

            // If we're less than the max matrix size, add a random moment
            while (topRow.size() < matrixSizeLimit) {
                int newLoc = rand(0, possibleMonomials.size()-1);
                Poly newMoment = *std::next(possibleMonomials.begin(), newLoc);
                if (std::find(topRowLocs.begin(), topRowLocs.end(), newLoc) == topRowLocs.end()) {
                    topRow.push_back(newMoment);
                    topRowLocs.push_back(newLoc);
                }
            }

            // Regenerate the moment matrix
            momentMatrices[0] = generateFromTopRow(topRow, use01);

            // Run the solver
            for (int j=0; j<topRowBest.size(); j++) {
                std::cout << topRowBest[j] << " ";
            }
            std::cout << std::endl;
            std::map<Mon, std::complex<double>> varVals;
            double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, 0, {-1, 1}, 0, &varVals);
            std::cout << i << " " << topRow.size() << " " << std::setprecision(10) << res << " " << bestRes << std::endl;

            // For each pair, update the correlation matrix
            for (int j=0; j<possibleMonomials.size(); j++) {
                for (int k=j; k<possibleMonomials.size(); k++) {
                    if (std::find(topRowLocs.begin(), topRowLocs.end(), j) != topRowLocs.end() && std::find(topRowLocs.begin(), topRowLocs.end(), k) != topRowLocs.end()) {
                        correlations[j][k].first += std::exp(res);
                        correlationCounts[j][k].first++;
                    } else {
                        correlations[j][k].second += std::exp(res);
                        correlationCounts[j][k].second++;
                    }
                }
            }

            // If it's worse, revert
            if (res > bestRes) {
                topRow = prevTopRow;
            } else {
                bestRes = res;
                topRowBest = topRow;
            }

            if (res <= 5.003501928) {
                break;
            }

        }

        // Average the correlations
        for (int i=0; i<correlations.size(); i++) {
            for (int j=i; j<correlations[i].size(); j++) {
                if (correlationCounts[i][j].first > 0) {
                    correlationsAveraged[i][j].first = correlations[i][j].first / correlationCounts[i][j].first;
                }
                if (correlationCounts[i][j].second > 0) {
                    correlationsAveraged[i][j].second = correlations[i][j].second / correlationCounts[i][j].second;
                }
            }
        }

        // Output the monomials
        std::cout << "Monomials: " << std::endl;
        for (int i=0; i<possibleMonomials.size(); i++) {
            std::cout << possibleMonomials[i] << " ";
        }
        std::cout << std::endl;
        
        // Pick the best
        std::vector<std::tuple<Poly,Poly,double,double>> bestPairs;
        for (int i=0; i<correlations.size(); i++) {
            for (int j=i; j<correlations[i].size(); j++) {
                if (correlationCounts[i][j].first > 0 && correlationCounts[i][j].second > 0) {
                    bestPairs.push_back(std::make_tuple(possibleMonomials[i], possibleMonomials[j], correlations[i][j].first, correlationsAveraged[i][j].second));
                }
            }
        }
        std::sort(bestPairs.begin(), bestPairs.end(), [](auto& left, auto& right) {
            return (std::get<2>(left) - std::get<3>(left)) < (std::get<2>(right) - std::get<3>(right));
        });

        // Output the best pairs
        std::cout << "50 Best pairs: " << std::endl;
        for (int i=0; i<std::min(50, int(bestPairs.size())); i++) {
            std::cout << std::get<0>(bestPairs[i]) << " " << std::get<1>(bestPairs[i]) << " " << std::get<2>(bestPairs[i]) << " " << std::get<3>(bestPairs[i]) << std::endl;
        }

        // Form a top row based on the first monomials 
        topRow = {Poly(1)};
        int k = 0;
        while (topRow.size() < matrixSizeLimit) {
            Poly mon1 = std::get<0>(bestPairs[k]);
            Poly mon2 = std::get<1>(bestPairs[k]);
            if (std::find(topRow.begin(), topRow.end(), mon1) == topRow.end()) {
                topRow.push_back(mon1);
            }
            if (std::find(topRow.begin(), topRow.end(), mon2) == topRow.end()) {
                topRow.push_back(mon2);
            }
            k++;
        }

        // Check the performance of this
        momentMatrices[0] = generateFromTopRow(topRow, use01);
        std::map<Mon, std::complex<double>> varVals;
        for (int j=0; j<topRow.size(); j++) {
            std::cout << topRow[j] << " ";
        }
        std::cout << std::endl;
        double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, 0, {-1, 1}, 0, &varVals);
        std::cout << "Final result: " << res << std::endl;

        return 0;

    }

    // Solve the problem
    std::pair<int,int> varBounds = {-1, 1};
    if (useDual) {
        varBounds = {-100000, 100000};
    }
    double res = 0;
    std::map<Mon, std::complex<double>> varVals;
    if (solver == "MOSEK") {
        if (verbosity >= 1) {
            std::cout << "Solving with MOSEK..." << std::endl;
        }
        res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varBounds, 0, &varVals);
    } else if (solver == "SCS") {
        if (verbosity >= 1) {
            std::cout << "Solving with SCS..." << std::endl;
        }
        res = solveSCS(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varBounds);
    } else if (solver == "Optim") {
        if (verbosity >= 1) {
            std::cout << "Solving with Optim..." << std::endl;
        }
        res = solveOptim(objective, constraintsZero, momentMatrices, varVals, verbosity, maxIters, numExtra, stepSize, tolerance);
    }
    if (useDual) {
        res = -res;
    }
    std::cout << "Result: " << std::setprecision(10) << res << std::endl;

    // If I3322, convert to the 0/1 version too
    if (problemName == "I3322" && !use01) {
        std::cout << "Result in 0/1: " << (res/4.0)-1.0 << std::endl;
    }

    // If Ising, round and provide a lower bound too
    if (problemName == "ISING") {

        std::map<Mon, std::complex<double>> roundedVals;
        for (auto& pair : varVals) {
            if (pair.first.size() == 1 && pair.first[0].first != 'z' && pair.first[0].first != 'y') {
                roundedVals[pair.first] = roundTo1(std::real(pair.second));
            }
        }
        if (roundedVals.size() == 0) {
            for (auto& pair : bellFunc) {
                roundedVals[pair.first] = roundTo1(std::real(pair.second));
            }
        }
        std::map<Mon, std::complex<double>> roundedValsAll = roundedVals;
        for (auto& pair : roundedVals) {
            for (auto& pair2 : roundedVals) {
                if (pair.first < pair2.first) {
                    roundedValsAll[pair.first*pair2.first] = roundedVals[pair.first]*roundedVals[pair2.first];
                }
            }
        }
        double lowerBound = bellFunc.eval(roundedValsAll).real();
        std::cout << "Lower bound: " << lowerBound << std::endl;

        // Anneal to try to get a better solution
        std::map<Mon, std::complex<double>> bestVals = roundedVals;
        std::map<Mon, std::complex<double>> newVals = roundedVals;
        std::map<Mon, std::complex<double>> prevVals = roundedVals;
        double bestVal = lowerBound;
        double temp = 1;
        int numIters = 1000;
        for (int i=0; i<numIters; i++) {
            prevVals = newVals;
            for (auto& pair : newVals) {
                if (rand(0, 1) < std::exp(-1.0/temp)) {
                    newVals[pair.first] = -newVals[pair.first];
                }
            }
            std::map<Mon, std::complex<double>> newValsAll = newVals;
            for (auto& pair : newVals) {
                for (auto& pair2 : newVals) {
                    if (pair.first < pair2.first) {
                        newValsAll[pair.first*pair2.first] = newVals[pair.first]*newVals[pair2.first];
                    }
                }
            }
            double newVal = bellFunc.eval(newValsAll).real();
            if (newVal > bestVal) {
                bestVal = newVal;
                bestVals = newVals;
            } else if (rand(0, 1) > std::exp((newVal-bestVal)/temp)) {
                newVals = prevVals;
            }
            temp -= 1.0/numIters;
        }
        std::cout << "After annealing: " << bestVal << std::endl;

    }

    // Start a timer
    std::chrono::steady_clock::time_point timeFinishedSolving = std::chrono::steady_clock::now();

    // Output timings
    if (verbosity >= 1) {
        int timeToGen = std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedGenerating - timeStart).count();
        int timeToSolve = std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedSolving - timeFinishedGenerating).count();
        std::cout << "Time to generate: " << timeToGen << "ms" << std::endl;
        std::cout << "Time to solve: " << timeToSolve << "ms" << std::endl;
    }

    // Exit without errors
    return 0;

}
