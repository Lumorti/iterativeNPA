// Standard includes
#include <iostream>
#include <vector>

// Import MOSEK
#include "fusion.h"

// Import Eigen
#include <Eigen/Dense>

// Define the type for the polynomials
typedef std::vector<std::pair<char,int>> monomial;
typedef std::vector<std::pair<double, monomial>> polynomial;
typedef std::vector<polynomial> polynomialVector;
typedef std::vector<std::vector<polynomial>> polynomialMatrix;

// Convert a constant polynomial to a double
double polynomialToDouble(polynomial p) {
    
    // Check if it's zero
    if (p.size() == 0 || (p.size() == 1 && p[0].first == 0)) {
        return 0;
    }

    // Check if it's a constant
    if (p.size() == 1 && p[0].second.size() == 0) {
        return p[0].first;
    }

    // Otherwise return 0
    return 0;

}

// Convert a constant polynomial matrix to a double matrix
std::vector<std::vector<double>> polynomialToDouble(polynomialMatrix m) {

    // Create the matrix
    std::vector<std::vector<double>> toReturn(m.size(), std::vector<double>(m[0].size(), 0));

    // Iterate through the matrix
    for (unsigned long int i=0; i<m.size(); i++) {
        for (unsigned long int j=0; j<m[i].size(); j++) {
            toReturn[i][j] = polynomialToDouble(m[i][j]);
        }
    }

    // Return the matrix
    return toReturn;

}

// Convert a constant polynomial vector to a double vector
std::vector<double> polynomialToDouble(polynomialVector r) {

    // Create the vector
    std::vector<double> toReturn(r.size(), 0);

    // Iterate through the vector
    for (unsigned long int i=0; i<r.size(); i++) {
        toReturn[i] = polynomialToDouble(r[i]);
    }

    // Return the vector
    return toReturn;

}

// Generate a monomial from a given string (e.g. "<A1B2>")
monomial stringToMonomial(std::string asString) {

    // If the string is empty, return an empty monomial
    if (asString == "" || asString.find("<") == std::string::npos) {
        return monomial();
    }

    // Iterate through the string
    monomial toReturn;
    char currentVariable = ' ';
    std::string currentIndex;
    for (unsigned long int i=0; i<asString.size(); i++) {

        // Check if this char is an integer
        if (asString[i] >= '0' && asString[i] <= '9') {
            currentIndex += asString[i];

        // Otherwise it's a variable
        } else if (asString[i] >= 'A' && asString[i] <= 'Z') {
            currentIndex += asString[i];
            if (currentVariable != ' ') {
                toReturn.push_back(std::make_pair(currentVariable, std::stoi(currentIndex)));
            }
            currentVariable = asString[i];
            currentIndex = "";
        }

    }

    // Add the last variable
    toReturn.push_back(std::make_pair(currentVariable, std::stoi(currentIndex)));

    // Return the monomial
    return toReturn;

}
        
// When printing a monomial, print it as a string
std::ostream& operator<<(std::ostream& os, const monomial& m) {

    // Iterate through the monomial
    if (m.size() > 0) {
        os << "<";
        for (unsigned long int i=0; i<m.size(); i++) {
            os << m[i].first << m[i].second;
        }
        os << ">";
    } else {
        os << "1";
    }

    // Return the stream
    return os;

}

// When printing a polynomial, print it as a string
std::ostream& operator<<(std::ostream& os, const polynomial& p) {

    // Check if it's zero
    if (p.size() == 0 || (p.size() == 1 && p[0].first == 0)) {
        os << "0";
        return os;
    }

    // Iterate through the polynomial
    for (unsigned long int i=0; i<p.size(); i++) {
        if (p[i].first == -1) {
            os << "-" << p[i].second;
        } else if (p[i].first == 1) {
            if (i == 0) {
                os << p[i].second;
            } else {
                os << "+" << p[i].second;
            }
        } else if (p[i].second.size() == 0 && p[i].first > 0) {
            if (i == 0) {
                os << p[i].first;
            } else {
                os << "+" << p[i].first;
            }
        } else if (p[i].second.size() == 0 && p[i].first < 0) {
            os << p[i].first;
        } else if (p[i].first > 0) {
            if (i == 0) {
                os << p[i].first << p[i].second;
            } else {
                os << "+" << p[i].first << p[i].second;
            }
        } else if (p[i].first != 0) {
            os << p[i].first << p[i].second;
        }
    }

    // Return the stream
    return os;

}

// When printing a polynomial matrix, print it as a string
std::ostream& operator<<(std::ostream& os, const polynomialMatrix& m) {

    // Determine the maximum width of each column
    std::vector<int> columnWidths(m[0].size(), 0);
    for (unsigned long int i=0; i<m.size(); i++) {
        for (unsigned long int j=0; j<m[i].size(); j++) {
            std::stringstream ss;
            ss << m[i][j];
            int sizeWhenWritten = ss.str().size();
            columnWidths[j] = std::max(columnWidths[j], sizeWhenWritten);
        }
    }

    // Iterate through the matrix
    for (unsigned long int i=0; i<m.size(); i++) {
        for (unsigned long int j=0; j<m[i].size(); j++) {
            std::stringstream ss;
            ss << m[i][j];
            int sizeWhenWritten = ss.str().size();
            for (unsigned long int k=0; k<columnWidths[j]-sizeWhenWritten; k++) {
                ss << " ";
            }
            ss << " ";
            os << ss.str();
        }
        os << std::endl;
    }

    // Return the stream
    return os;

}

// When printing a vector of polynomials, print it as a string
std::ostream& operator<<(std::ostream& os, const std::vector<polynomial>& v) {

    // Iterate through the vector
    for (unsigned long int i=0; i<v.size(); i++) {
        os << v[i] << std::endl;
    }

    // Return the stream
    return os;

}

// Generate a polynomial from a given string (e.g. "2<A1B2>-<A3A1>")
polynomial stringToPolynomial(std::string asString) {

    // Iterate through the string
    polynomial toReturn;
    std::string currentCoefficient;
    std::string currentMonomial;
    for (unsigned long int i=0; i<asString.size(); i++) {

        // If finished a monomial
        if (i > 0 && (asString[i] == '>' || ((asString[i] == '+' || asString[i] == '-') && asString[i-1] != '>'))) {

            // Ensure the coefficient is convertable
            if (currentCoefficient == "" || currentCoefficient == "+") {
                currentCoefficient = "1";
            } else if (currentCoefficient == "-") {
                currentCoefficient = "-1";
            }

            // Add the term
            if (asString[i] == '>') {
                currentMonomial += asString[i];
            }
            toReturn.push_back(std::make_pair(std::stod(currentCoefficient), stringToMonomial(currentMonomial)));
            currentMonomial = "";
            currentCoefficient = "";

            // The plus and the minus are for the next term
            if (asString[i] == '+' || asString[i] == '-') {
                currentCoefficient += asString[i];
            }

        // If starting or continuing a monomial
        } else if (asString[i] == '<' || currentMonomial.size() > 0) {
            currentMonomial += asString[i];

        // Otherwise it's for the coefficient
        } else if (asString[i] != ' ') {
            currentCoefficient += asString[i];

        }

    }

    // If there's a coefficient left over, add it
    if (currentCoefficient != "") {
        toReturn.push_back(std::make_pair(std::stod(currentCoefficient), stringToMonomial(currentMonomial)));
    }

    // Return the polynomial
    return toReturn;

}

// Given a polynomial and a variable list, evaluate the polynomial
double evaluatePolynomial(polynomial& p, std::vector<monomial>& variables, std::vector<double>& varVals) {

    // Replace each variable with its value
    double toReturn = 0;
    for (int i=0; i<p.size(); i++) {

        // Find this variable
        for (int j=0; j<variables.size(); j++) {
            if (p[i].second == variables[j]) {
                toReturn += p[i].first*varVals[j];
                break;
            }
        }

    }

    // Return the value
    return toReturn;

}

// Given a matrix and a variable list, return the matrix with the variables replaced
Eigen::MatrixXd replaceVariables(polynomialMatrix& momentMatrix, std::vector<monomial>& variables, std::vector<double>& varVals) {

    // Replace each variable with its value
    Eigen::MatrixXd momentMatrixEigen = Eigen::MatrixXd::Zero(momentMatrix.size(), momentMatrix.size());
    for (int i=0; i<momentMatrix.size(); i++) {
        for (int j=0; j<momentMatrix[i].size(); j++) {
            for (int k=0; k<momentMatrix[i][j].size(); k++) {

                // Find this variable
                bool found = false;
                for (int l=0; l<variables.size(); l++) {
                    if (momentMatrix[i][j][k].second == variables[l]) {
                        momentMatrixEigen(i, j) += momentMatrix[i][j][k].first*varVals[l];
                        found = true;
                        break;
                    }
                }

            }

        }
    }
    return momentMatrixEigen;

}

// Get the eigenvalues and vectors of a matrix after replacement
void getEigens(polynomialMatrix& momentMatrix, std::vector<monomial>& variables, std::vector<double>& varVals, std::vector<std::vector<double>>& eigenvectors, std::vector<double>& eigenvalues) {

    // Replace each variable with its value
    Eigen::MatrixXd momentMatrixEigen = replaceVariables(momentMatrix, variables, varVals);

    // Get the eigenvalues and vectors of this
    Eigen::EigenSolver<Eigen::MatrixXd> es(momentMatrixEigen);
    Eigen::MatrixXd eigenVectorsEigen = es.eigenvectors().real();
    Eigen::VectorXd eigenValuesEigen = es.eigenvalues().real();

    // Copy into the output vectors
    for (int i=0; i<eigenVectorsEigen.cols(); i++) {
        std::vector<double> eigenVector;
        for (int j=0; j<eigenVectorsEigen.rows(); j++) {
            eigenVector.push_back(eigenVectorsEigen(j, i));
        }
        eigenvectors.push_back(eigenVector);
    }
    for (int i=0; i<eigenValuesEigen.size(); i++) {
        eigenvalues.push_back(eigenValuesEigen(i));
    }

}

// Given a monomila, reduce it to its simplest form
monomial reduceMonomial(monomial mon_) {

    monomial mon = mon_;

    // Sort the monomial as much as we can
    for (int i=0; i<mon.size(); i++) {
        for (int j=0; j<int(mon.size())-1; j++) {
            if (mon[j].first != mon[j+1].first && mon[j] > mon[j+1]) {
                std::swap(mon[j], mon[j+1]);
            }
        }
    }

    // <A1A1> = <1>
    int i = 0;
    while (i < int(mon.size())-1) {
        if (mon[i] == mon[i+1]) {
            mon.erase(mon.begin()+i+1);
            mon.erase(mon.begin()+i);
            i = -1;
        }
        i++;
    }

    return mon;

}

// Combinations of a vector
void combo(std::vector<int> &alphabet, int n, std::vector<std::vector<int>> &result, std::vector<int> curr) {
    if (n == 0) {
        result.push_back(curr);
        return;
    }
    int largestSoFar = 0;
    if (curr.size() > 0) {
        largestSoFar = curr[curr.size()-1];
    }
    for (int i=largestSoFar; i<alphabet.size(); i++) {
        if (alphabet[i] > largestSoFar) {
            std::vector<int> newCurr = curr;
            newCurr.push_back(alphabet[i]);
            combo(alphabet, n - 1, result, newCurr);
        }
    }
    return;
}

// Sampled combinations of a vector
void combo(std::vector<int> &alphabet, int n, std::vector<std::vector<int>> &result, int numToSample) {
    int numFails = 0;
    while (result.size() < numToSample) {
        std::vector<int> curr;
        std::vector<int> alphabetCopy = alphabet;
        for (int j=0; j<n; j++) {
            int i = rand() % alphabetCopy.size();
            curr.push_back(alphabetCopy[i]);
            alphabetCopy.erase(alphabetCopy.begin()+i);
        }
        std::sort(curr.begin(), curr.end());
        if (std::find(result.begin(), result.end(), curr) == result.end()) {
            result.push_back(curr);
        } else {
            numFails++;
            if (numFails > numToSample*2) {
                break;
            }
        }
    }
    return;
}

// Add all single order monomials from a functional to a list of variables
void addSingleMonomials(std::vector<monomial>& variables, polynomial functional) {

    // Iterate through the functional
    for (long unsigned int i=0; i<functional.size(); i++) {

        // Iterate through the monomial
        for (long unsigned int j=0; j<functional[i].second.size(); j++) {
            monomial currentMonomial = {functional[i].second[j]};

            // Check if this moment is already in the list
            bool found = false;
            for (long unsigned int k=0; k<variables.size(); k++) {
                if (variables[k] == currentMonomial) {
                    found = true;
                    break;
                }
            }

            // If not, add it in the correct place
            if (!found) {
                bool hasBeenAdded = false;
                for (long unsigned int k=0; k<variables.size(); k++) {
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

// Generate a moment matrix given the top row
polynomialMatrix generateFromTopRow(std::vector<monomial> monomsInTopRow) {

    // Generate all combinations of the top row
    polynomialMatrix matrixToReturn = std::vector<std::vector<polynomial>>(monomsInTopRow.size(), std::vector<polynomial>(monomsInTopRow.size()));
    for (long unsigned int i=0; i<monomsInTopRow.size(); i++) {
        for (long unsigned int j=i; j<monomsInTopRow.size(); j++) {

            // Form the new moment
            monomial newMonomial;
            for (long unsigned int k=0; k<monomsInTopRow[i].size(); k++) {
                newMonomial.push_back(monomsInTopRow[i][k]);
            }
            for (int k=int(monomsInTopRow[j].size())-1; k>=0; k--) {
                newMonomial.push_back(monomsInTopRow[j][k]);
            }

            // Reduce the monomial
            newMonomial = reduceMonomial(newMonomial);

            // Set the matrix elements
            matrixToReturn[i][j] = polynomial(1, std::make_pair(1, newMonomial));
            matrixToReturn[j][i] = polynomial(1, std::make_pair(1, newMonomial));

        }
    }

    // Return the matrix
    return matrixToReturn;

}

// Generate a moment matrix given the top row as polynomials
polynomialMatrix generateFromTopRow(std::vector<polynomial> monomsInTopRow) {

    // Generate all combinations of the top row
    polynomialMatrix matrixToReturn = std::vector<std::vector<polynomial>>(monomsInTopRow.size(), std::vector<polynomial>(monomsInTopRow.size()));
    for (long unsigned int i=0; i<monomsInTopRow.size(); i++) {
        for (long unsigned int j=i; j<monomsInTopRow.size(); j++) {

            // Form the new polynomial
            polynomial newPolynomial;
            for (long unsigned int k=0; k<monomsInTopRow[i].size(); k++) {
                for (long unsigned int l=0; l<monomsInTopRow[j].size(); l++) {

                    // Form the new moment
                    monomial newMonomial;
                    for (long unsigned int m=0; m<monomsInTopRow[i][k].second.size(); m++) {
                        newMonomial.push_back(monomsInTopRow[i][k].second[m]);
                    }
                    for (int m=int(monomsInTopRow[j][l].second.size())-1; m>=0; m--) {
                        newMonomial.push_back(monomsInTopRow[j][l].second[m]);
                    }

                    // Reduce the monomial
                    newMonomial = reduceMonomial(newMonomial);

                    // Add to the polynomial
                    bool found = false;
                    double newCoefficient = monomsInTopRow[i][k].first * monomsInTopRow[j][l].first;
                    for (long unsigned int m=0; m<newPolynomial.size(); m++) {
                        if (newPolynomial[m].second == newMonomial) {
                            newPolynomial[m].first += newCoefficient;
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        newPolynomial.push_back(std::make_pair(newCoefficient, newMonomial));
                    }

                }
            }

            // Set the matrix elements
            matrixToReturn[i][j] = newPolynomial;
            matrixToReturn[j][i] = newPolynomial;

        }
    }

    // Return the matrix
    return matrixToReturn;

}
// Generate all moment matrices for a given level given a polynomial
std::vector<polynomialMatrix> generateAllMomentMatrices(polynomial functional, int level, int subLevel=-1, int numToSample=-1) {

    // First get the list of all moments used by iterating through the polynomial
    std::vector<monomial> variables;
    addSingleMonomials(variables, functional);

    // Generate all moments up to the given level
    std::vector<monomial> monomsInTopRow = {};
    if (level >= 1) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            monomial currentMonomial = {variables[i][0]};
            currentMonomial = reduceMonomial(currentMonomial);
            if (currentMonomial.size() > 0 && std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentMonomial) == monomsInTopRow.end()) {
                monomsInTopRow.push_back(currentMonomial);
            }
        }
    }
    if (level >= 2) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                monomial currentMonomial = {variables[i][0], variables[j][0]};
                currentMonomial = reduceMonomial(currentMonomial);
                if (currentMonomial.size() > 0 && std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentMonomial) == monomsInTopRow.end()) {
                    monomsInTopRow.push_back(currentMonomial);
                }
            }
        }
    }
    if (level >= 3) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                for (long unsigned int k=0; k<variables.size(); k++) {
                    monomial currentMonomial = {variables[i][0], variables[j][0], variables[k][0]};
                    currentMonomial = reduceMonomial(currentMonomial);
                    if (currentMonomial.size() > 0 && std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentMonomial) == monomsInTopRow.end()) {
                        monomsInTopRow.push_back(currentMonomial);
                    }
                }
            }
        }
    }
    if (level >= 4) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                for (long unsigned int k=0; k<variables.size(); k++) {
                    for (long unsigned int l=0; l<variables.size(); l++) {
                        monomial currentMonomial = {variables[i][0], variables[j][0], variables[k][0], variables[l][0]};
                        currentMonomial = reduceMonomial(currentMonomial);
                        if (currentMonomial.size() > 0 && std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentMonomial) == monomsInTopRow.end()) {
                            monomsInTopRow.push_back(currentMonomial);
                        }
                    }
                }
            }
        }
    }
    if (level >= 5) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                for (long unsigned int k=0; k<variables.size(); k++) {
                    for (long unsigned int l=0; l<variables.size(); l++) {
                        for (long unsigned int m=0; m<variables.size(); m++) {
                            monomial currentMonomial = {variables[i][0], variables[j][0], variables[k][0], variables[l][0], variables[m][0]};
                            currentMonomial = reduceMonomial(currentMonomial);
                            if (currentMonomial.size() > 0 && std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentMonomial) == monomsInTopRow.end()) {
                                monomsInTopRow.push_back(currentMonomial);
                            }
                        }
                    }
                }
            }
        }
    }
    if (level >= 6) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                for (long unsigned int k=0; k<variables.size(); k++) {
                    for (long unsigned int l=0; l<variables.size(); l++) {
                        for (long unsigned int m=0; m<variables.size(); m++) {
                            for (long unsigned int n=0; n<variables.size(); n++) {
                                monomial currentMonomial = {variables[i][0], variables[j][0], variables[k][0], variables[l][0], variables[m][0], variables[n][0]};
                                currentMonomial = reduceMonomial(currentMonomial);
                                if (currentMonomial.size() > 0 && std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentMonomial) == monomsInTopRow.end()) {
                                    monomsInTopRow.push_back(currentMonomial);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if  (level >= 7) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                for (long unsigned int k=0; k<variables.size(); k++) {
                    for (long unsigned int l=0; l<variables.size(); l++) {
                        for (long unsigned int m=0; m<variables.size(); m++) {
                            for (long unsigned int n=0; n<variables.size(); n++) {
                                for (long unsigned int o=0; o<variables.size(); o++) {
                                    monomial currentMonomial = {variables[i][0], variables[j][0], variables[k][0], variables[l][0], variables[m][0], variables[n][0], variables[o][0]};
                                    currentMonomial = reduceMonomial(currentMonomial);
                                    if (currentMonomial.size() > 0 && std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentMonomial) == monomsInTopRow.end()) {
                                        monomsInTopRow.push_back(currentMonomial);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Sublevel is capped at the size
    int maxSubLevel = 1+int(monomsInTopRow.size());
    if (subLevel < 0 || subLevel > maxSubLevel) {
        subLevel = maxSubLevel;
    }

    // Create the index array to be used for combinations
    std::vector<int> indices;
    for (int i=0; i<monomsInTopRow.size(); i++) {
        indices.push_back(i+1);
    }
    std::vector<std::vector<int>> combinations;

    // If we just want one big matrix, don't bother with combs
    if (subLevel == maxSubLevel) {
        combinations = {indices};

    // If we want all the combs 
    } else if (numToSample == -1) {
        combo(indices, subLevel-1, combinations, {});

    // If we just want to sample some of the combs
    } else {
        combo(indices, subLevel-1, combinations, numToSample);
    }

    // Each moment mat should start with 1
    monomsInTopRow.insert(monomsInTopRow.begin(), monomial());
    for (long unsigned int i=0; i<combinations.size(); i++) {
        combinations[i].insert(combinations[i].begin(), 0);
    }

    // Form the moment matrices
    std::vector<polynomialMatrix> matricesToReturn;
    for (int k=0; k<combinations.size(); k++) {

        // Get this top row combination
        std::vector<monomial> monomsInTopRowComb;
        for (long unsigned int i=0; i<combinations[k].size(); i++) {
            monomsInTopRowComb.push_back(monomsInTopRow[combinations[k][i]]);
        }

        // Form that matrix
        polynomialMatrix newMatrix = generateFromTopRow(monomsInTopRowComb);
        matricesToReturn.push_back(newMatrix);

    }

    // Return the moment matrix
    return matricesToReturn;

}

// Add variables from a moment matrix
void addVariables(std::vector<monomial>& variables, polynomialMatrix toAdd) {

    // Iterative through the matrix
    for (long unsigned int i=0; i<toAdd.size(); i++) {
        for (long unsigned int j=0; j<toAdd[i].size(); j++) {

            // Iterate through the polynomial
            for (long unsigned int k=0; k<toAdd[i][j].size(); k++) {
                monomial currentMonomial = toAdd[i][j][k].second;

                // Check if this monomial is already in the list
                bool found = false;
                for (long unsigned int l=0; l<variables.size(); l++) {
                    if (variables[l] == currentMonomial) {
                        found = true;
                        break;
                    }
                }

                // If not, add it
                if (!found) {
                    variables.push_back(currentMonomial);
                }

            }

        }
    }

}

// Add variables from a polynomial
void addVariables(std::vector<monomial>& variables, polynomial toAdd) {

    // Iterate through the monomial
    for (long unsigned int i=0; i<toAdd.size(); i++) {
        monomial currentMonomial = {toAdd[i].second};

        // Check if this monomial is already in the list
        bool found = false;
        for (long unsigned int j=0; j<variables.size(); j++) {
            if (variables[j] == currentMonomial) {
                found = true;
                break;
            }
        }

        // If not, add it
        if (!found) {
            variables.push_back(currentMonomial);
        }

    }

}

// Convert to MOSEK form and solve
double solveMOSEK(polynomial obj, std::vector<polynomialMatrix>& psd, std::vector<polynomial> constraintsZero, std::vector<polynomial> constraintsPositive, int verbosity, std::vector<monomial>& variables, std::vector<double>& variableValues, bool isDual=false) {

    // Get the list of variables
    int oneIndex = 0;
    variables = {monomial()};
    for (int i=0; i<psd.size(); i++) {
        addVariables(variables, psd[i]);
    }
    for (int i=0; i<constraintsZero.size(); i++) {
        addVariables(variables, constraintsZero[i]);
    }
    for (int i=0; i<constraintsPositive.size(); i++) {
        addVariables(variables, constraintsPositive[i]);
    }
    addVariables(variables, obj);

    // The c vector defining the objective
    std::vector<double> c(variables.size());
    for (int i=0; i<obj.size(); i++) {

        // Find the location of this variable
        monomial toFind = {obj[i].second};
        for (int j=0; j<variables.size(); j++) {
            if (variables[j] == toFind) {
                c[j] += obj[i].first;
                break;
            }
        }

    }
    auto cM = monty::new_array_ptr<double>(c);

    // The A matrix defining the equality constraints
    std::vector<int> ARows;
    std::vector<int> ACols;
    std::vector<double> AVals;
    for (int i=0; i<constraintsZero.size(); i++) {
        for (int j=0; j<constraintsZero[i].size(); j++) {

            // The row is the constraint number
            ARows.push_back(i);

            // Find the location of this variable
            monomial toFind = {constraintsZero[i][j].second};
            for (int k=0; k<variables.size(); k++) {
                if (variables[k] == toFind) {
                    ACols.push_back(k);
                    break;
                }
            }

            // The coefficient is the value
            AVals.push_back(constraintsZero[i][j].first);

        }
    }
    auto AM = mosek::fusion::Matrix::sparse(constraintsZero.size(), variables.size(), monty::new_array_ptr<int>(ARows), monty::new_array_ptr<int>(ACols), monty::new_array_ptr<double>(AVals));

    // The B vector defining the positivity constraints
    std::vector<int> BRows;
    std::vector<int> BCols;
    std::vector<double> BVals;
    for (int i=0; i<constraintsPositive.size(); i++) {
        for (int j=0; j<constraintsPositive[i].size(); j++) {

            // The row is the constraint number
            BRows.push_back(i);

            // Find the location of this variable
            monomial toFind = {constraintsPositive[i][j].second};
            for (int k=0; k<variables.size(); k++) {
                if (variables[k] == toFind) {
                    BCols.push_back(k);
                    break;
                }
            }

            // The coefficient is the value
            BVals.push_back(constraintsPositive[i][j].first);

        }
    }
    auto BM = mosek::fusion::Matrix::sparse(constraintsPositive.size(), variables.size(), monty::new_array_ptr<int>(BRows), monty::new_array_ptr<int>(BCols), monty::new_array_ptr<double>(BVals));

    // Output the variable list
    if (verbosity >= 2) {
        std::cout << "Variables:" << std::endl;
        for (int i=0; i<variables.size(); i++) {
            std::cout << i << " " << variables[i] << std::endl;
        }
    }

    // The vectors defining the PSD constraints
    std::vector<std::shared_ptr<monty::ndarray<int,1>>> indicesPSDPerMat;
    std::vector<std::shared_ptr<monty::ndarray<double,1>>> coeffsPSDPerMat;
    std::vector<std::pair<int,int>> matDims;
    for (int k=0; k<psd.size(); k++) {

        // Determine how many mats we need to sum to make this matrix
        int numMatsNeeded = 0;
        for (int i=0; i<psd[k].size(); i++) {
            for (int j=i; j<psd[k][i].size(); j++) {
                if (psd[k][i][j].size() > numMatsNeeded) {
                    numMatsNeeded = psd[k][i][j].size();
                }
            }
        }

        // For each part of the sum
        int sVecSize = psd[k].size() * (psd[k].size() + 1) / 2;
        std::vector<int> indicesPSD;
        std::vector<double> coeffsPSD;
        for (int l=0; l<numMatsNeeded; l++) {

            // The indices and coefficients for the svec
            for (int i=0; i<psd[k].size(); i++) {
                for (int j=i; j<psd[k][i].size(); j++) {

                    // If there are no more, this is zero
                    if (l >= psd[k][i][j].size()) {
                        indicesPSD.push_back(0);
                        coeffsPSD.push_back(0.0);
                        continue;
                    }

                    // Find this in the variable list
                    monomial toFind = {psd[k][i][j][l].second};
                    for (int k=0; k<variables.size(); k++) {
                        if (variables[k] == toFind) {
                            indicesPSD.push_back(k);
                            break;
                        }
                    }

                    // The coeff for mosek's svec
                    if (i != j) {
                        coeffsPSD.push_back(psd[k][i][j][l].first*std::sqrt(2.0));
                    } else {
                        coeffsPSD.push_back(psd[k][i][j][l].first);
                    }

                }

            }

        }

        // Convert to MOSEK form
        auto indicesPSDM = monty::new_array_ptr<int>(indicesPSD);
        auto coeffsPSDM = monty::new_array_ptr<double>(coeffsPSD);

        // Add to the list
        indicesPSDPerMat.push_back(indicesPSDM);
        coeffsPSDPerMat.push_back(coeffsPSDM);
        matDims.push_back({numMatsNeeded, sVecSize});

    }

    // Create a model
    mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});
    if (verbosity >= 2) {
        M->setLogHandler([=](const std::string & msg) {std::cout << msg << std::flush;});
    }

    // Create the main variable vector
    mosek::fusion::Variable::t xM = M->variable(variables.size());
    //mosek::fusion::Variable::t xM = M->variable(variables.size(), mosek::fusion::Domain::inRange(-1, 1));
    if (!isDual) {
        M->constraint(xM, mosek::fusion::Domain::inRange(-1.0, 1.0));
    }

    // The objective function
    if (isDual) {
        M->objective(mosek::fusion::ObjectiveSense::Minimize, mosek::fusion::Expr::dot(cM, xM));
    } else {
        M->objective(mosek::fusion::ObjectiveSense::Maximize, mosek::fusion::Expr::dot(cM, xM));
    }

    // The one variable should be fixed
    M->constraint(xM->index(oneIndex), mosek::fusion::Domain::equalsTo(1.0));

    // The matrix of this should be PSD
    for (int k=0; k<psd.size(); k++) {
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
    if (verbosity >= 1) {
        std::cout << "Objective value: " << objPrimal << std::endl;
    }

    // Get all of the variable values
    auto xMLevel = *(xM->level());
    variableValues.resize(variables.size());
    for (int i=0; i<variables.size(); i++) {
        variableValues[i] = xMLevel[i];
    }

    // Check the eigenvalues of each moment matrix
    if (verbosity >= 1) {
        for (int i=0; i<psd.size(); i++) {
            std::vector<std::vector<double>> eigenvectors;
            std::vector<double> eigenvalues;
            getEigens(psd[i], variables, variableValues, eigenvectors, eigenvalues);
            std::sort(eigenvalues.begin(), eigenvalues.end());
            std::cout << "Min eigenvalue: " << eigenvalues[0] << std::endl;
        }
    }

    // If superverbose, output all monomials
    if (verbosity >= 3) {
        std::cout << "Solution: " << std::endl;
        for (int i=0; i<variables.size(); i++) {
            std::cout << variables[i] << ": " << variableValues[i] << std::endl;
        }
        std::cout << "There are " << variables.size() << " variables." << std::endl;

        // Count the number of unique values
        std::vector<double> uniqueValues;
        for (int i=0; i<variableValues.size(); i++) {
            if (std::abs(variableValues[i]) > 1e-5) {
                bool found = false;
                for (int j=0; j<uniqueValues.size(); j++) {
                    if (std::abs((variableValues[i] - uniqueValues[j]) / variableValues[i]) < 1e-3) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    uniqueValues.push_back(variableValues[i]);
                }
            }
        }
        std::cout << "There are " << uniqueValues.size() << " unique values." << std::endl;

        Eigen::MatrixXd A = replaceVariables(psd[0], variables, variableValues);
        std::cout << "Moment matrix:" << std::endl;
        std::cout << psd[0] << std::endl;
        std::cout << "Moment matrix with vars replaced:" << std::endl;
        std::cout << A << std::endl;

    // If verbose, just output monomials that are in the objective
    } else if (verbosity >= 2) {
        std::cout << "Solution: " << std::endl;
        std::vector<monomial> variablesInObj = {};
        addVariables(variablesInObj, obj);
        for (int i=0; i<variablesInObj.size(); i++) {
            for (int j=0; j<variables.size(); j++) {
                if (variablesInObj[i] == variables[j]) {
                    std::cout << variablesInObj[i] << ": " << variableValues[j] << std::endl;
                    break;
                }
            }
        }
    } 

    // Return the objective
    return objPrimal;

}

// Wrapper in case vars aren't needed
double solveMOSEK(polynomial obj, std::vector<polynomialMatrix>& psd, std::vector<polynomial> constraintsZero, std::vector<polynomial> constraintsPositive, int verbosity, bool isDual=false) {
    std::vector<monomial> variables;
    std::vector<double> variableValues;
    return solveMOSEK(obj, psd, constraintsZero, constraintsPositive, verbosity, variables, variableValues, isDual);
}

// Encode a monomial into a format better suited for a neural net 
std::vector<double> encodeMonomial(monomial m, polynomialMatrix& mat) {

    // The output
    std::vector<double> output(2, 0.0);

    // Find the monomial in the matrix
    for (int i=0; i<mat.size(); i++) {
        for (int j=0; j<mat[i].size(); j++) {
            if (mat[i][j][0].second == m) {
                output[0] = i;
                output[1] = j;
                return output;
            }
        }

    }

    // Return the output
    return output;

}

// Get the adjugate of a matrix
Eigen::MatrixXd adjugate(const Eigen::MatrixXd& A) {

    // Quite expensive, there might be a better way
    int n = A.rows();
    Eigen::MatrixXd adjugateMatrix(n, n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {

            // Calculate the (n-1)x(n-1) submatrix Mij
            // Delete the empty rows so the determinant can be non zero
            Eigen::MatrixXd submatrix = Eigen::MatrixXd::Zero(n-1, n-1);
            int submatrixRow = 0;
            for (int k=0; k<n; k++) {
                if (k != i) {
                    int submatrixCol = 0;
                    for (int l=0; l<n; l++) {
                        if (l != j) {
                            submatrix(submatrixRow, submatrixCol) = A(k, l);
                            submatrixCol++;
                        }
                    }
                    submatrixRow++;
                }
            }

            // Calculate the cofactor Cij
            adjugateMatrix(j, i) = std::pow(-1, i + j) * submatrix.determinant();

        }
    }

    return adjugateMatrix;

}

// Generic entry function
int main(int argc, char* argv[]) {

    // Define the Bell scenario
    int numOutcomesA = 2;
    int numOutcomesB = 2;
    int level = 1;
    int subLevel = -1;
    polynomial bellFunc;
    int testing = 0;
    int numToSample = 1000;
    int verbosity = 1;
    int maxIters = 1000;
    std::string seed = "";

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
        } else if (argAsString == "--chsh") {
            bellFunc = stringToPolynomial("<A1B1>+<A1B2>+<A2B1>-<A2B2>");

        // I3322
        // however for now it's l1 = 5.5, l2 = 5.00631, l3 (88= 5.0035, l4 = 
        // l1 (7x7) = 5.5
        // l2 (28x28) = 5.00612
        // l3 (88x88) = 5.0035
        // l4 (244x244) = 5.0035
        } else if (argAsString == "--I3322") {
            bellFunc = stringToPolynomial("<A1>-<A2>+<B1>-<B2>-<A1B1>+<A1B2>+<A2B1>-<A2B2>+<A1B3>+<A2B3>+<A3B1>+<A3B2>");

        // Set the level
        } else if (argAsString == "-l") {
            level = std::stoi(argv[i+1]);
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
        } else if (argAsString == "-h") {
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
            std::cout << "  -v <num>        Verbosity level" << std::endl;
            std::cout << "  -t              Test the moment matrix" << std::endl;
            return 0;

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
    std::vector<polynomialMatrix> momentMatrices = generateAllMomentMatrices(bellFunc, level, subLevel, numToSample);
    if (verbosity >= 1) {
        std::cout << "Generated " << momentMatrices.size() << " moment matrices" << std::endl;
        int largestMomentMatrix = 0;
        for (int i=0; i<momentMatrices.size(); i++) {
            if (momentMatrices[i].size() > largestMomentMatrix) {
                largestMomentMatrix = momentMatrices[i].size();
            }
        }
        std::cout << "Largest moment matrix has size " << largestMomentMatrix << std::endl;
    }

    // Define other constraints
    polynomial objective = bellFunc;
    std::vector<polynomial> constraintsZero;
    std::vector<polynomial> constraintsPositive;

    // If testing TODO
    if (testing == 1) {

        // Put into into standard SDP form
        // s.t. X >= 0
        
        // C.X
        int matSize = momentMatrices[0].size();
        polynomialMatrix C = std::vector<polynomialVector>(matSize, polynomialVector(matSize, polynomial()));
        for (int i=0; i<objective.size(); i++) {
            monomial monToFind = objective[i].second;
            int locX = -1;
            int locY = -1;
            for (int j=0; j<momentMatrices[0].size(); j++) {
                for (int k=j+1; k<momentMatrices[0][j].size(); k++) {
                    if (momentMatrices[0][j][k].size() == 1 && momentMatrices[0][j][k][0].second == monToFind) {
                        locX = j;
                        locY = k;
                        break;
                    }
                }
            }
            C[locX][locY] = {std::make_pair(objective[i].first/2.0, monomial())};
            C[locY][locX] = {std::make_pair(objective[i].first/2.0, monomial())};
        }

        // A_i.X = b_i
        std::vector<polynomialMatrix> As;
        std::vector<polynomial> b;
        for (int i=0; i<momentMatrices[0].size(); i++) {
            polynomialMatrix A = std::vector<polynomialVector>(matSize, polynomialVector(matSize, polynomial()));
            A[i][i] = stringToPolynomial("1");
            As.push_back(A);
            b.push_back(stringToPolynomial("1"));
        }
        std::vector<monomial> alreadyDone;
        for (int i=0; i<momentMatrices[0].size(); i++) {
            for (int j=i+1; j<momentMatrices[0][i].size(); j++) {

                // Check if we've already found all versions of this monom
                monomial monToFind = momentMatrices[0][i][j][0].second;
                if (std::find(alreadyDone.begin(), alreadyDone.end(), monToFind) != alreadyDone.end()) {
                    continue;
                } else {
                    alreadyDone.push_back(monToFind);
                }

                // Find all monomials that are the same
                std::vector<std::pair<int,int>> equalMonomLocs;
                //equalMonomLocs.push_back(std::make_pair(j,i));
                for (int k=i; k<momentMatrices[0].size(); k++) {
                    for (int l=k+1; l<momentMatrices[0][k].size(); l++) {
                        if (k == i && l == j) {
                            continue;
                        }
                        if (momentMatrices[0][k][l].size() == 1 && momentMatrices[0][k][l][0].second == monToFind) {
                            equalMonomLocs.push_back(std::make_pair(k,l));
                            //equalMonomLocs.push_back(std::make_pair(l,k));
                        }
                    }
                }

                // For each equal monom, add a matrix
                std::cout << "For monomial " << monToFind << " at (" << i << "," << j << ") found " << equalMonomLocs.size() << " equal monomials" << std::endl;
                for (int k=0; k<equalMonomLocs.size(); k++) {
                    std::cout << "Adding constraint for monomial " << monToFind << " at (" << equalMonomLocs[k].first << "," << equalMonomLocs[k].second << ")" << std::endl;
                    polynomialMatrix A = std::vector<polynomialVector>(matSize, polynomialVector(matSize, polynomial()));
                    A[i][j] = stringToPolynomial("1");
                    A[equalMonomLocs[k].first][equalMonomLocs[k].second] = stringToPolynomial("-1");
                    As.push_back(A);
                    b.push_back(stringToPolynomial("0"));
                }

            }
        }

        // Output the primal in this form
        std::cout << "Moment mat:" << std::endl;
        std::cout << momentMatrices[0] << std::endl << std::endl;
        std::cout << "C: " << std::endl;
        std::cout << C << std::endl << std::endl;
        std::cout << "As: " << std::endl;
        for (int i=0; i<As.size(); i++) {
            std::cout << As[i] << std::endl << std::endl;
        }
        std::cout << "b: " << std::endl;
        for (int i=0; i<b.size(); i++) {
            std::cout << b[i] << std::endl;
        }

        // Convert everything to MOSEK form and test it
        std::vector<std::shared_ptr<monty::ndarray<double,2>>> AsM;
        std::vector<std::vector<std::vector<double>>> AsMat;
        for (int i=0; i<As.size(); i++) {
            AsMat.push_back(polynomialToDouble(As[i]));
            AsM.push_back(monty::new_array_ptr<double>(AsMat[i]));
        }
        std::vector<std::vector<double>> CMat = polynomialToDouble(C);
        std::shared_ptr<monty::ndarray<double,2>> CM = monty::new_array_ptr<double>(CMat);
        std::vector<double> bM = polynomialToDouble(b);

        // Create a model
        mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});
        if (verbosity >= 2) {
            M->setLogHandler([=](const std::string & msg) {std::cout << msg << std::flush;});
        }
        mosek::fusion::Variable::t XM = M->variable(matSize, matSize, mosek::fusion::Domain::inPSDCone());
        M->objective(mosek::fusion::ObjectiveSense::Maximize, mosek::fusion::Expr::dot(CM, XM));
        for (int i=0; i<As.size(); i++) {
            M->constraint(mosek::fusion::Expr::dot(AsM[i], XM), mosek::fusion::Domain::equalsTo(bM[i]));
        }

        // Output the primal objective value
        M->solve();
        double objPrimal = M->primalObjValue();
        if (verbosity >= 1) {
            std::cout << "Objective value: " << objPrimal << std::endl;
        }
        
        // Convert to the dual
        // b.y
        polynomial newObjective;
        for (int i=0; i<bM.size(); i++) {
            if (std::abs(bM[i]) > 1e-10) {
                newObjective.push_back(std::make_pair(bM[i], stringToMonomial("<C" + std::to_string(i) + ">")));
            }
        }
        std::cout << "New Objective: " << newObjective << std::endl;
        // s.t. C-\sum_i y_i A_i >= 0
        polynomialMatrix newMomentMat = polynomialMatrix(matSize, polynomialVector(matSize, polynomial()));
        for (int i=0; i<matSize; i++) {
            for (int j=0; j<matSize; j++) {
                newMomentMat[i][j].push_back(std::make_pair(CMat[i][j], monomial()));
                for (int k=0; k<AsMat.size(); k++) {
                    if (std::abs(AsMat[k][i][j]) > 1e-10) {
                        newMomentMat[i][j].push_back(std::make_pair(-AsMat[k][i][j], stringToMonomial("<C" + std::to_string(k) + ">")));
                    }
                }
            }
        }
        std::cout << "New moment matrix: " << std::endl;
        std::cout << newMomentMat << std::endl << std::endl;

        // Solve the dual
        objective = newObjective;
        momentMatrices = {newMomentMat};
        constraintsZero = {};
        constraintsPositive = {};
        solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, true);

        return 0;

	// Try the iterative approach TODO
    } else if (testing == 2) {

		// Run the main problem
		std::vector<double> dualSol;
        solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, false, dualSol);
		
		// Run the subproblem

		// Add the constraint

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

    // Convert to MOSEK form and solve
    solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity);

    // Exit without errors
    return 0;

}
