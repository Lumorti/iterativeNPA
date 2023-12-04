// Standard includes
#include <iostream>
#include <vector>

// Import MOSEK
#include "fusion.h"

// Optim
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"

// Import Eigen
#include <Eigen/Dense>

// Autodiff
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

// Define the type for the polynomials
typedef std::vector<std::pair<char,int>> monomial;
typedef std::vector<std::pair<double, monomial>> polynomial;
typedef std::vector<polynomial> polynomialVector;
typedef std::vector<std::vector<polynomial>> polynomialMatrix;
typedef std::vector<double> doubleVector;
typedef std::vector<std::vector<double>> doubleMatrix;

// Check if a polynomial is constant
bool polynomialIsConstant(polynomial p) {

    // Check if it's zero
    if (p.size() == 0 || (p.size() == 1 && p[0].first == 0)) {
        return true;
    }

    // Check if it's a constant
    if (p.size() == 1 && p[0].second.size() == 0) {
        return true;
    }

    // Otherwise return false
    return false;

}

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

// Given a polynomial, combine all of the terms of the same monomial
polynomial simplify(polynomial p) {

    // Create the return polynomial
    polynomial toReturn;

    // Iterate through the polynomial
    for (unsigned long int i=0; i<p.size(); i++) {

        // Check if this monomial is already in the return polynomial
        bool found = false;
        for (unsigned long int j=0; j<toReturn.size(); j++) {
            if (toReturn[j].second == p[i].second) {
                toReturn[j].first += p[i].first;
                found = true;
                break;
            }
        }

        // If not, add it
        if (!found) {
            toReturn.push_back(p[i]);
        }

    }

    // Return the simplified polynomial
    return toReturn;

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

// When printing a matrix of doubles, print it as a string
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double>>& m) {

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
    eigenvalues.clear();
    eigenvectors.clear();
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

// Given a monomial, reduce it to its simplest form
monomial reduceMonomial(monomial mon_, bool use01=false, bool trySwap=true) {

    // Sort the monomial as much as we can
    monomial mon = mon_;
    for (int i=0; i<mon.size(); i++) {
        for (int j=0; j<int(mon.size())-1; j++) {
            if (mon[j].first != mon[j+1].first && mon[j] > mon[j+1]) {
                std::swap(mon[j], mon[j+1]);
            }
        }
    }

    // <A1A1> = <A1>
    if (use01) {
        int i = 0;
        while (i < int(mon.size())-1) {
            if (mon[i] == mon[i+1]) {
                mon.erase(mon.begin()+i);
                i = -1;
            }
            i++;
        }

    // <A1A1> = <1>
    } else {
        int i = 0;
        while (i < int(mon.size())-1) {
            if (mon[i] == mon[i+1]) {
                mon.erase(mon.begin()+i+1);
                mon.erase(mon.begin()+i);
                i = -1;
            }
            i++;
        }
    }

    // Flip it to see if it's smaller
    if (trySwap) {
        monomial monFlipped = mon;
        std::reverse(monFlipped.begin(), monFlipped.end());
        if (monFlipped < mon) {
            mon = monFlipped;
        }
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
polynomialMatrix generateFromTopRow(std::vector<monomial> monomsInTopRow, bool use01=false) {

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
            newMonomial = reduceMonomial(newMonomial, use01);

            // Set the matrix elements
            matrixToReturn[i][j] = polynomial(1, std::make_pair(1, newMonomial));
            matrixToReturn[j][i] = polynomial(1, std::make_pair(1, newMonomial));

        }
    }

    // Return the matrix
    return matrixToReturn;

}

// Generate a moment matrix given the top row as polynomials
polynomialMatrix generateFromTopRow(std::vector<polynomial> monomsInTopRow, bool use01=false) {

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
                    newMonomial = reduceMonomial(newMonomial, use01);

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
std::vector<polynomialMatrix> generateAllMomentMatrices(polynomial functional, int level, int subLevel=-1, int numToSample=-1, bool use01=false) {

    // First get the list of all moments used by iterating through the polynomial
    std::vector<monomial> variables;
    addSingleMonomials(variables, functional);

    // Generate all moments up to the given level
    std::vector<monomial> monomsInTopRow = {};
    if (level >= 1) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            monomial currentMonomial = {variables[i][0]};
            currentMonomial = reduceMonomial(currentMonomial, use01, false);
            if (currentMonomial.size() > 0 && std::find(monomsInTopRow.begin(), monomsInTopRow.end(), currentMonomial) == monomsInTopRow.end()) {
                monomsInTopRow.push_back(currentMonomial);
            }
        }
    }
    if (level >= 2) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=0; j<variables.size(); j++) {
                monomial currentMonomial = {variables[i][0], variables[j][0]};
                currentMonomial = reduceMonomial(currentMonomial, use01, false);
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
                    currentMonomial = reduceMonomial(currentMonomial, use01, false);
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
                        currentMonomial = reduceMonomial(currentMonomial, use01, false);
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
                            currentMonomial = reduceMonomial(currentMonomial, use01, false);
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
                                currentMonomial = reduceMonomial(currentMonomial, use01, false);
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
                                    currentMonomial = reduceMonomial(currentMonomial, use01, false);
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
        polynomialMatrix newMatrix = generateFromTopRow(monomsInTopRowComb, use01);
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
double solveMOSEK(polynomial obj, std::vector<polynomialMatrix>& psd, std::vector<polynomial> constraintsZero, std::vector<polynomial> constraintsPositive, int verbosity, std::vector<monomial>& variables, std::vector<double>& variableValues, bool use01=false) {

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

    // If we use 0/1 rather than +1/-1, constrain to be in [0,1]
    if (use01) {
        M->constraint(xM, mosek::fusion::Domain::inRange(0.0, 1.0));
    }

    // The objective function
    M->objective(mosek::fusion::ObjectiveSense::Maximize, mosek::fusion::Expr::dot(cM, xM));

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

// The structure used to hold the data for the optimisation
struct optimData {
    int numD;
    double penalty;
    polynomial objective;
    std::vector<polynomial> equalityCons;
};

// Given a monomial and an autodiff vector, evaluate the monomial
autodiff::real evalCommutingMonom(const double& coeff, const monomial& monom, const autodiff::ArrayXreal& x, const int& numD) {
    autodiff::real term = coeff;
    if (monom.size() > 0) {
        for (int j=0; j<monom.size(); j++) {
            if (monom[j].first == 'D') {
                term *= x(monom[j].second);
            } else {
                term *= x(monom[j].second + numD);
            }
        }
    }
    return term;
}

// Given a monomial and an Eigen vector, evaluate the monomial
double evalCommutingMonom(const double& coeff, const monomial& monom, const Eigen::VectorXd& x, const int& numD) {
    double term = coeff;
    if (monom.size() > 0) {
        for (int j=0; j<monom.size(); j++) {
            if (monom[j].first == 'D') {
                term *= x(monom[j].second);
            } else {
                term *= x(monom[j].second + numD);
            }
        }
    }
    return term;
}

// Given a monomial and an Eigen vector, evaluate the monomial without one of the vars
double evalCommutingMonomWithout(const double& coeff, const monomial& monom, const Eigen::VectorXd& x, const int& numD, const int& var) {
    double term = coeff;
    if (monom.size() > 0) {
        for (int j=0; j<monom.size(); j++) {
            if (j != var) {
                if (monom[j].first == 'D') {
                    term *= x(monom[j].second);
                } else {
                    term *= x(monom[j].second + numD);
                }
            }
        }
    }
    return term;
}


autodiff::real optFunctionAD(const autodiff::ArrayXreal& x, const optimData* optDataRecast) {

    // The cost is a combination of the objective and the equality constraints
    autodiff::real cost = 0.0;
    for (int i=0; i<optDataRecast->objective.size(); i++) {
        cost -= evalCommutingMonom(optDataRecast->objective[i].first, optDataRecast->objective[i].second, x, optDataRecast->numD);
    }
    std::vector<double> equalityConVals(optDataRecast->equalityCons.size(), 0.0);
    for (int i=0; i<optDataRecast->equalityCons.size(); i++) {
        autodiff::real equalityConVal = 0.0;
        for (int j=0; j<optDataRecast->equalityCons[i].size(); j++) {
            equalityConVal += evalCommutingMonom(optDataRecast->equalityCons[i][j].first, optDataRecast->equalityCons[i][j].second, x, optDataRecast->numD);
        }
        //DEBUG
        //std::cout << "Equality constraint " << i << ": " << equalityConVals[i] << std::endl;
        cost += equalityConVal * equalityConVal * optDataRecast->penalty;
    }

    return cost;

}

// Optim function with auto diff
double optFunction(const Eigen::VectorXd& x, Eigen::VectorXd* grad_out, void* optData) {

    // Recast this generic pointer into the correct format
    optimData* optDataRecast = reinterpret_cast<optimData*>(optData);

    // Using the types from autodiff
    autodiff::real u;
    autodiff::ArrayXreal xd = x.eval();

    // If requested, compute the gradient alongside the function using auto diff
    if (grad_out) {
        Eigen::VectorXd grad_tmp = autodiff::gradient(optFunctionAD, autodiff::wrt(xd), autodiff::at(xd, optDataRecast), u);
        *grad_out = grad_tmp;

    // Otherwise, just evaluate the function
    } else {
        u = optFunctionAD(xd, optDataRecast);
    }

    // Convert to double and return
    return u.val();

}

// Cost/gradient function for optim
static double gradFunction(const Eigen::VectorXd& x, Eigen::VectorXd* gradOut, void* optData) {
    
    // Recast this generic pointer into the correct format
    optimData* optDataRecast = reinterpret_cast<optimData*>(optData);

    // Output each variables DEBUG
    //std::cout << "Variables: " << std::endl;
    //for (int i=0; i<x.size(); i++) {
        //if (i < optDataRecast->numD) {
            //std::cout << "D" << i << ": " << x(i) << std::endl;
        //} else {
            //std::cout << "C" << i - optDataRecast->numD << ": " << x(i) << std::endl;
        //}
    //}

    // The cost is a combination of the objective and the equality constraints
    double cost = 0.0;
    for (int i=0; i<optDataRecast->objective.size(); i++) {
        cost -= evalCommutingMonom(optDataRecast->objective[i].first, optDataRecast->objective[i].second, x, optDataRecast->numD);
    }
    std::vector<double> equalityConVals(optDataRecast->equalityCons.size(), 0.0);
    for (int i=0; i<optDataRecast->equalityCons.size(); i++) {
        equalityConVals[i] = 0.0;
        for (int j=0; j<optDataRecast->equalityCons[i].size(); j++) {
            equalityConVals[i] += evalCommutingMonom(optDataRecast->equalityCons[i][j].first, optDataRecast->equalityCons[i][j].second, x, optDataRecast->numD);
        }
        //DEBUG
        //std::cout << "Equality constraint " << i << ": " << equalityConVals[i] << std::endl;
        cost += equalityConVals[i] * equalityConVals[i] * optDataRecast->penalty;
    }

    // Calculate the gradient
    if (gradOut) {

        // Reset the gradient
        *gradOut = Eigen::VectorXd::Zero(x.size());

        // Check the contribution from the objective
        for (int j=0; j<optDataRecast->objective.size(); j++) {
            if (optDataRecast->objective[j].second.size() > 0) {
                int ind = optDataRecast->objective[j].second[0].second;
                if (optDataRecast->objective[j].second[0].first == 'C') {
                    ind += optDataRecast->numD;
                }
                (*gradOut)(ind) -= optDataRecast->objective[j].first;
            }
        }

        // Check the contribution from the equality constraints
        for (int i=0; i<optDataRecast->equalityCons.size(); i++) {
            for (int j=0; j<optDataRecast->equalityCons[i].size(); j++) {

                // Check if it's non-constant
                if (optDataRecast->equalityCons[i][j].second.size() >= 2) {

                    // If it's a square of a variable
                    if (optDataRecast->equalityCons[i][j].second[0] == optDataRecast->equalityCons[i][j].second[1]) {

                        // Get the index of the variable
                        int ind = optDataRecast->equalityCons[i][j].second[0].second;
                        if (optDataRecast->equalityCons[i][j].second[0].first == 'C') {
                            ind += optDataRecast->numD;
                        }

                        // Add the contribution
                        (*gradOut)(ind) += 4.0 * x(ind) * equalityConVals[i] * optDataRecast->penalty;

                    // Otherwise we need to evaluate the rest of the monomial
                    } else {

                        // Get the index of the variable
                        int ind1 = optDataRecast->equalityCons[i][j].second[0].second;
                        if (optDataRecast->equalityCons[i][j].second[0].first == 'C') {
                            ind1 += optDataRecast->numD;
                        }
                        int ind2 = optDataRecast->equalityCons[i][j].second[1].second;
                        if (optDataRecast->equalityCons[i][j].second[1].first == 'C') {
                            ind2 += optDataRecast->numD;
                        }

                        // Evaluate that monomial without each variable
                        double without1 = evalCommutingMonomWithout(optDataRecast->objective[j].first, optDataRecast->objective[j].second, x, optDataRecast->numD, 0);
                        double without2 = evalCommutingMonomWithout(optDataRecast->objective[j].first, optDataRecast->objective[j].second, x, optDataRecast->numD, 1);
                        double coeff = optDataRecast->equalityCons[i][j].first;

                        // Add the contribution
                        (*gradOut)(ind1) += 2.0 * coeff * without1 * equalityConVals[i] * optDataRecast->penalty;
                        (*gradOut)(ind2) += 2.0 * coeff * without2 * equalityConVals[i] * optDataRecast->penalty;

                    }

                } else if (optDataRecast->equalityCons[i][j].second.size() == 1) {

                    // Get the index of the variable
                    int ind = optDataRecast->equalityCons[i][j].second[0].second;
                    if (optDataRecast->equalityCons[i][j].second[0].first == 'C') {
                        ind += optDataRecast->numD;
                    }

                    double coeff = optDataRecast->equalityCons[i][j].first;

                    // Add the contribution
                    (*gradOut)(ind) += 2.0 * coeff * equalityConVals[i] * optDataRecast->penalty;

                }

            }
        }

        // Output the gradient DEBUG
        //std::cout << "Gradient: " << std::endl;
        //for (int i=0; i<gradOut->size(); i++) {
            //if (i < optDataRecast->numD) {
                //std::cout << "D" << i << ": " << (*gradOut)(i) << std::endl;
            //} else {
                //std::cout << "C" << i - optDataRecast->numD << ": " << (*gradOut)(i) << std::endl;
            //}
        //}

    }

    // Return just the cost
    return cost;

}

// Given a location in the Cholesky matrix, return the index in the vector
int cholLocToInd(int i, int j, int width) {
    if (i <= j) {
        return i*width + j - (i*(i+1))/2;
    } else {
        return cholLocToInd(j, i, width);
    }
}

// Generic entry function
int main(int argc, char* argv[]) {

    // Define the Bell scenario
    int numOutcomesA = 2;
    int numOutcomesB = 2;
    int level = 1;
    int subLevel = -1;
    polynomial bellFunc = stringToPolynomial("<A1B1>+<A1B2>+<A2B1>-<A2B2>");
    int testing = 0;
    int numToSample = 1000;
    int verbosity = 1;
    int maxIters = 1000;
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
            bellFunc = stringToPolynomial("<A1B1>+<A1B2>+<A2B1>-<A2B2>");
            problemName = "CHSH";

        // If told to use a 0/1 output
        } else if (argAsString == "-01") {
            use01 = true;

        // I3322 TODO
        // l1 (7x7) = 5.5
        // l2 (28x28) = 5.00612
        // l3 (88x88) = 5.0035
        // l4 (244x244) = 5.0035
        } else if (argAsString == "--I3322" || argAsString == "--i3322") {
            if (use01) {
                bellFunc = stringToPolynomial("-<A2>-<B1>-2<B2>+<A1B1>+<A1B2>+<A2B1>+<A2B2>-<A1B3>+<A2B3>-<A3B1>+<A3B2>");
            } else {
                bellFunc = stringToPolynomial("<A1>-<A2>+<B1>-<B2>-<A1B1>+<A1B2>+<A2B1>-<A2B2>+<A1B3>+<A2B3>+<A3B1>+<A3B2>");
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
            std::cout << "  -01             Use 0/1 output instead of -1/1" << std::endl;
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
    std::vector<polynomialMatrix> momentMatrices = generateAllMomentMatrices(bellFunc, level, subLevel, numToSample, use01);
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

    // If told to add extra to the top row
    if (extraMoments.size() > 0) {
        std::vector<polynomial> topRow = momentMatrices[0][0];
        for (int i=0; i<extraMoments.size(); i++) {
            topRow.push_back(stringToPolynomial(extraMoments[i]));
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
    polynomial objective = bellFunc;
    std::vector<polynomial> constraintsZero;
    std::vector<polynomial> constraintsPositive;

    // If testing
    if (testing == 1) {

        // Put into into standard SDP form
        // min C.X s.t. A.X = b, X >= 0
        
        // Define C
        int matSize = momentMatrices[0].size();
        std::vector<std::vector<double>> C = std::vector<std::vector<double>>(matSize, std::vector<double>(matSize, 0));
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
            C[locX][locY] = -objective[i].first/2.0;
            C[locY][locX] = -objective[i].first/2.0;
        }

        // Define the A matrices and b vector
        std::vector<doubleMatrix> As;
        std::vector<double> b;
        for (int i=0; i<momentMatrices[0].size(); i++) {
            for (int j=i; j<momentMatrices[0][i].size(); j++) {

                // Trying to find the matrix relating this to other elements
                doubleMatrix A(matSize, doubleVector(matSize, 0));
                if (i == j) { 
                    A[i][j] = 1;
                } else {
                    A[i][j] = 0.5;
                }
                double bVal = 0;
                bool somethingFound = false;

                // For each term in this poly, try to find an early term 
                // that would reduce this to a constant
                for (int k=0; k<momentMatrices[0][i][j].size(); k++) {

                    // If it's a constant
                    if (momentMatrices[0][i][j][k].second.size() == 0) {
                        bVal += momentMatrices[0][i][j][k].first;
                        somethingFound = true;
                        continue;
                    }

                    // Otherwise find an earlier variable
                    monomial monomToFind = momentMatrices[0][i][j][k].second;
                    bool found = false;
                    for (int l=0; l<momentMatrices[0].size(); l++) {
                        for (int m=l; m<momentMatrices[0][l].size(); m++) {
                            if (l*momentMatrices[0].size() + m >= i*momentMatrices[0].size() + j) {
                                break;
                            }
                            if (momentMatrices[0][l][m].size() == 1 && momentMatrices[0][l][m][0].second == monomToFind) {
                                if (l == m) { 
                                    A[l][m] = -momentMatrices[0][i][j][k].first;
                                } else {
                                    A[l][m] = -momentMatrices[0][i][j][k].first / 2.0;
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

                // Add this matrix
                if (somethingFound) {
                    As.push_back(A);
                    b.push_back(bVal);
                }

            }
        }

        // Output the A matrices
        //std::cout << "A matrices: " << std::endl;
        //for (int i=0; i<As.size(); i++) {
            //std::cout << As[i] << std::endl;
        //}

        // Convert the constraints to the dual objective b.y
        polynomial newObjective;
        for (int i=0; i<b.size(); i++) {
            if (std::abs(b[i]) > 1e-10) {
                newObjective.push_back(std::make_pair(b[i], stringToMonomial("<D" + std::to_string(i) + ">")));
            }
        }
        std::cout << "Dual objective: " << newObjective << std::endl;

        // Convert the objective to the dual constraints C-\sum_i y_i A_i >= 0
        polynomialMatrix newMomentMat = polynomialMatrix(matSize, polynomialVector(matSize, polynomial()));
        for (int i=0; i<matSize; i++) {
            for (int j=i; j<matSize; j++) {
                newMomentMat[i][j].push_back(std::make_pair(C[i][j], monomial()));
                for (int k=0; k<As.size(); k++) {
                    if (std::abs(As[k][i][j]) > 1e-10) {
                        newMomentMat[i][j].push_back(std::make_pair(-As[k][i][j], stringToMonomial("<D" + std::to_string(k) + ">")));
                    }
                }
            }
        }

        // Symmetrize the matrix
        for (int i=0; i<matSize; i++) {
            for (int j=i+1; j<matSize; j++) {
                newMomentMat[j][i] = newMomentMat[i][j];
            }
        }
        std::cout << "Dual moment matrix: " << std::endl;
        std::cout << newMomentMat << std::endl << std::endl;

        // Solve the dual
        polynomial objectiveDual = newObjective;
        std::vector<polynomialMatrix> momentMatricesDual = {newMomentMat};
        std::vector<polynomial> constraintsZeroDual = {};
        std::vector<polynomial> constraintsPositiveDual = {};
        std::vector<double> varVals;
        std::vector<monomial> varNames;
        solveMOSEK(objectiveDual, momentMatricesDual, constraintsZeroDual, constraintsPositiveDual, verbosity, varNames, varVals, use01);

        // Sort the varNames array by the names, based on the x.second for each x
        std::sort(varNames.begin(), varNames.end(), [](const monomial& a, const monomial& b) {
            if (a.size() > 0 && b.size() > 0) {
                return a[0].second < b[0].second;
            } else {
                return a.size() < b.size();
            }
        });

        // Get the eigenvalues if we start with all zeros
        std::vector<std::vector<double>> eigenvectors;
        std::vector<double> eigenvalues;
        varVals[0] = 1;
        for (int i=1; i<varVals.size(); i++) {
            varVals[i] = 0.0;
        }
        getEigens(newMomentMat, varNames, varVals, eigenvectors, eigenvalues);
        double minEigenvalue = 100000000;
        for (int i=0; i<eigenvalues.size(); i++) {
            minEigenvalue = std::min(minEigenvalue, eigenvalues[i]);
        }
        std::cout << "Minimum eigenvalue: " << minEigenvalue << std::endl;
        for (int i=0; i<newObjective.size(); i++) {
            for (int j=0; j<varNames.size(); j++) {
                if (varNames[j] == newObjective[i].second) {
                    varVals[j] += minEigenvalue - 1e-8;
                    break;
                }
            }
        }

        // Now recalculate the min eigenvalue
        getEigens(newMomentMat, varNames, varVals, eigenvectors, eigenvalues);
        minEigenvalue = 100000000;
        for (int i=0; i<eigenvalues.size(); i++) {
            minEigenvalue = std::min(minEigenvalue, eigenvalues[i]);
        }
        std::cout << "Minimum eigenvalue: " << minEigenvalue << std::endl;

        // Now calculate the Cholesky decomposition
        Eigen::MatrixXd A = replaceVariables(newMomentMat, varNames, varVals);
        Eigen::MatrixXd L = A.llt().matrixL();
        std::cout << "Cholesky decomposition: " << std::endl;
        std::cout << L.transpose() << std::endl;

        // Recostruct the moment matrix
        Eigen::MatrixXd testMat = L * L.transpose();
        std::cout << "Reconstructed moment matrix: " << std::endl;
        std::cout << testMat << std::endl;

        // Get the generic LL^T matrix
        std::vector<polynomial> equals;
        polynomialMatrix LLT = polynomialMatrix(matSize, polynomialVector(matSize, polynomial()));
        for (int i=0; i<matSize; i++) {
            for (int j=i; j<matSize; j++) {

                // Generate this element of the generic LLT matrix
                polynomial thisLLT;
                for (int k=0; k<=i; k++) {
                    thisLLT.push_back(std::make_pair(1.0, stringToMonomial("<C" + std::to_string(cholLocToInd(k, i, matSize)) + "> " + "<C" + std::to_string(cholLocToInd(k, j, matSize)) + ">")));
                }

                // Now subtract the polynomial from the dual moment matrix
                for (int k=0; k<momentMatricesDual[0][i][j].size(); k++) {
                    thisLLT.push_back(std::make_pair(-momentMatricesDual[0][i][j][k].first, momentMatricesDual[0][i][j][k].second));
                }

                // Combine any like terms
                thisLLT = simplify(thisLLT);

                // Add this polynomial to the list
                equals.push_back(thisLLT);

                // Output this polynomial
                std::cout << thisLLT << std::endl;

            }
        }

        // Count the vars needed, ignoring 1
        int numD = varNames.size();
        if (varNames[0] == monomial()) {
            numD--;
        }
        int numVars = numD + matSize*(matSize+1)/2;
        std::cout << "Num D: " << numD << std::endl;
        std::cout << "Num vars: " << numVars << std::endl;

        // Set up the optimisation
        Eigen::VectorXd x = Eigen::VectorXd::Zero(numVars);
        optimData optData;
        optData.objective = objectiveDual;
        optData.equalityCons = equals;
        optData.penalty = 1e8;
        optData.numD = numD;

        // Output all the variable values
        std::cout << "Variable values: " << std::endl;
        for (int i=0; i<varNames.size(); i++) {
            std::cout << varNames[i] << ": " << varVals[i] << std::endl;
        }

        // Set the initial x
        for (int i=0; i<numD; i++) {
            x[i] = varVals[i+1];
        }
        for (int i=0; i<matSize; i++) {
            for (int j=i; j<matSize; j++) {
                x[numD + cholLocToInd(i, j, matSize)] = L(j, i);
            }
        }

        // Check the initial cost
        Eigen::VectorXd gradData = Eigen::VectorXd::Zero(numVars);
        double startCost = gradFunction(x, &gradData, &optData);
        std::cout << "Starting cost: " << startCost << std::endl;

        // Settings for the optimiser
        optim::algo_settings_t settings;
        settings.print_level = verbosity;
        settings.iter_max = maxIters;
        settings.rel_objfn_change_tol = 1e-10;
        settings.gd_settings.method = 6;
        settings.lbfgs_settings.par_M = 10;
        settings.lbfgs_settings.wolfe_cons_1 = 1e-3;
        settings.lbfgs_settings.wolfe_cons_2 = 0.8;

        // Level 3 uses 500MB and takes about 2 seconds
        // Level 4 uses 10GB

        // Optimise TODO
        if (maxIters > 0) {
            optData.penalty = 1e3;
            for (int i=0; i<10; i++) {
                std::cout << "Optimising with " << numVars << " variables and penalty " << optData.penalty << std::endl;
                bool success = optim::lbfgs(x, optFunction, &optData, settings);
                std::cout << "Result = " << gradFunction(x, &gradData, &optData) << std::endl;
                optData.penalty *= 2;
            }
        }

		return 0;

    } else if (testing == 2) {

        std::cout << "Primal moment matrix: " << std::endl;
        std::cout << momentMatrices[0] << std::endl << std::endl;
        std::cout << "Primal objective: " << std::endl;
        std::cout << objective << std::endl << std::endl;

        // Keep iterating, for now just a fixed number of times
        for (int iter=0; iter<maxIters; iter++) {

            // Put into into standard SDP form
            // min C.X s.t. A.X = b, X >= 0

            // Define C
            int matSize = momentMatrices[0].size();
            std::vector<std::vector<double>> C = std::vector<std::vector<double>>(matSize, std::vector<double>(matSize, 0));
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
                C[locX][locY] = -objective[i].first/2.0;
                C[locY][locX] = -objective[i].first/2.0;
            }

            // Define the A matrices and b vector
            std::vector<doubleMatrix> As;
            std::vector<double> b;
            for (int i=0; i<momentMatrices[0].size(); i++) {
                for (int j=i; j<momentMatrices[0][i].size(); j++) {

                    // Trying to find the matrix relating this to other elements
                    doubleMatrix A(matSize, doubleVector(matSize, 0));
                    if (i == j) { 
                        A[i][j] = 1;
                    } else {
                        A[i][j] = 0.5;
                    }
                    double bVal = 0;
                    bool somethingFound = false;

                    // For each term in this poly, try to find an early term 
                    // that would reduce this to a constant
                    for (int k=0; k<momentMatrices[0][i][j].size(); k++) {

                        // If it's a constant
                        if (momentMatrices[0][i][j][k].second.size() == 0) {
                            bVal += momentMatrices[0][i][j][k].first;
                            somethingFound = true;
                            continue;
                        }

                        // Otherwise find an earlier variable
                        monomial monomToFind = momentMatrices[0][i][j][k].second;
                        bool found = false;
                        for (int l=0; l<momentMatrices[0].size(); l++) {
                            for (int m=l; m<momentMatrices[0][l].size(); m++) {
                                if (l*momentMatrices[0].size() + m >= i*momentMatrices[0].size() + j) {
                                    break;
                                }
                                if (momentMatrices[0][l][m].size() == 1 && momentMatrices[0][l][m][0].second == monomToFind) {
                                    if (l == m) { 
                                        A[l][m] = -momentMatrices[0][i][j][k].first;
                                    } else {
                                        A[l][m] = -momentMatrices[0][i][j][k].first / 2.0;
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

                    // Add this matrix
                    if (somethingFound) {
                        As.push_back(A);
                        b.push_back(bVal);
                    }

                }
            }
                
            // Convert the constraints to the dual objective b.y
            polynomial newObjective;
            for (int i=0; i<b.size(); i++) {
                if (std::abs(b[i]) > 1e-10) {
                    newObjective.push_back(std::make_pair(b[i], stringToMonomial("<D" + std::to_string(i) + ">")));
                }
            }
            std::cout << "Dual objective: " << newObjective << std::endl;

            // Convert the objective to the dual constraints C-\sum_i y_i A_i >= 0
            polynomialMatrix newMomentMat = polynomialMatrix(matSize, polynomialVector(matSize, polynomial()));
            for (int i=0; i<matSize; i++) {
                for (int j=i; j<matSize; j++) {
                    newMomentMat[i][j].push_back(std::make_pair(C[i][j], monomial()));
                    for (int k=0; k<As.size(); k++) {
                        if (std::abs(As[k][i][j]) > 1e-10) {
                            newMomentMat[i][j].push_back(std::make_pair(-As[k][i][j], stringToMonomial("<D" + std::to_string(k) + ">")));
                        }
                    }
                }
            }

            // Symmetrize the matrix
            for (int i=0; i<matSize; i++) {
                for (int j=i+1; j<matSize; j++) {
                    newMomentMat[j][i] = newMomentMat[i][j];
                }
            }
            std::cout << "Dual moment matrix: " << std::endl;
            std::cout << newMomentMat << std::endl << std::endl;

            // Solve the dual
            polynomial objectiveDual = newObjective;
            std::vector<polynomialMatrix> momentMatricesDual = {newMomentMat};
            std::vector<polynomial> constraintsZeroDual = {};
            std::vector<polynomial> constraintsPositiveDual = {};
            std::vector<double> varVals;
            std::vector<monomial> varNames;
            solveMOSEK(objectiveDual, momentMatricesDual, constraintsZeroDual, constraintsPositiveDual, verbosity, varNames, varVals);

            // Using this dual info, generate the sub-problem
            std::vector<polynomialMatrix> momentMatricesSub = generateAllMomentMatrices(bellFunc, level+1, subLevel, numToSample, use01);
            Eigen::MatrixXd Y = replaceVariables(momentMatricesDual[0], varNames, varVals);
            std::cout << "Y: " << std::endl;
            std::cout << Y << std::endl << std::endl;
            polynomial objectiveSub;
            for (int i=0; i<Y.rows()-iter; i++) {
                for (int j=i; j<Y.cols()-iter; j++) {
                    for (int k=0; k<momentMatricesSub[0][i][j].size(); k++) {
                        if (i == j) {
                            objectiveSub.push_back(std::make_pair(-Y(i,j), momentMatricesSub[0][i][j][k].second));
                        } else {
                            objectiveSub.push_back(std::make_pair(-2.0*Y(i,j), momentMatricesSub[0][i][j][k].second));
                        }
                    }
                }
            }
            std::cout << "Sub problem objective: " << std::endl;
            objectiveSub = simplify(objectiveSub);
            std::cout << objectiveSub << std::endl << std::endl;
            std::vector<polynomial> constraintsZeroSub = {};
            std::vector<polynomial> constraintsPositiveSub = {};

            // Solve the sub problem
            std::vector<double> varValsSub;
            std::vector<monomial> varNamesSub;
            double valSub = solveMOSEK(objectiveSub, momentMatricesSub, constraintsZeroSub, constraintsPositiveSub, verbosity, varNamesSub, varValsSub);

            // Check for convergence
            if (std::abs(valSub) < 1e-5) {
                std::cout << "Converged due to sub-problem giving zero" << std::endl;
                break;
            }

            // Output the vars
            //std::cout << "Vars: " << std::endl;
            //for (int i=0; i<varNamesSub.size(); i++) {
                //std::cout << varNamesSub[i] << " = " << varValsSub[i] << std::endl;
            //}

            // Generate the new constraint 
            polynomial newPositivityCon = objectiveSub;
            newPositivityCon.push_back(std::make_pair(-valSub, monomial()));
            //for (int i=0; i<objectiveSub.size(); i++) {
                //for (int j=0; j<varNamesSub.size(); j++) {
                    //if (varNamesSub[j] == objectiveSub[i].second) {
                        //newPositivityCon.push_back(std::make_pair(objectiveSub[i].first, objectiveSub[i].second));
                        //break;
                    //}
                //}
            //}
            newPositivityCon = simplify(newPositivityCon);
            for (int i=0; i<newPositivityCon.size(); i++) {
                newPositivityCon[i].first *= -1;
            }
            std::cout << "New positivity constraint: " << std::endl;
            std::cout << newPositivityCon << std::endl << std::endl;

            // Add this to the original main moment matrix, expanding by one
            momentMatrices[0].push_back(polynomialVector(matSize+1, stringToPolynomial("0")));
            for (int i=0; i<matSize; i++) {
                momentMatrices[0][i].push_back(stringToPolynomial("0"));
            }
            momentMatrices[0][matSize][matSize] = newPositivityCon;

            // Output the new moment matrix
            std::cout << "New moment matrix: " << std::endl;
            std::cout << momentMatrices[0] << std::endl << std::endl;

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

    // Convert to MOSEK form and solve
    std::vector<monomial> varNames;
    std::vector<double> varVals;
    double res = solveMOSEK(objective, momentMatrices, constraintsZero, constraintsPositive, verbosity, varNames, varVals, use01);
    std::cout << "Result: " << res << std::endl;

    // If I3322, convert to the 0/1 version too
    if (problemName == "I3322" && !use01) {
        std::cout << "Result in 0/1: " << (res/4.0)-1.0 << std::endl;
    }

    // Exit without errors
    return 0;

}
