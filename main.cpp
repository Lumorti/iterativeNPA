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
typedef std::vector<std::vector<polynomial>> polynomialMatrix;

// Generate a monomial from a given string (e.g. "<A1B2>")
monomial stringToMonomial(std::string asString) {

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
        } else {
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
        if (asString[i] == '>') {

            // Ensure the coefficient is convertable
            if (currentCoefficient == "+") {
                currentCoefficient = "1";
            } else if (currentCoefficient == "-") {
                currentCoefficient = "-1";
            }

            // Add the term
            currentMonomial += asString[i];
            toReturn.push_back(std::make_pair(std::stod(currentCoefficient), stringToMonomial(currentMonomial)));
            currentMonomial = "";
            currentCoefficient = "";

        // If starting or continuing a monomial
        } else if (asString[i] == '<' || currentMonomial.size() > 0) {
            currentMonomial += asString[i];

        // Otherwise it's for the coefficient
        } else if (asString[i] != ' ') {
            currentCoefficient += asString[i];

        }

    }

    // Return the polynomial
    return toReturn;

}

// Generate the moment matrix for a given level given a polynomial
polynomialMatrix generateMomentMatrix(polynomial functional, int level) {

    // First get the list of all moments used by iterating through the polynomial
    std::vector<monomial> variables;
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

            // If not, add it
            if (!found) {
                variables.push_back(currentMonomial);
            }

        }

    }

    // Generate all moments up to the given level
    std::vector<monomial> monomsInTopRow = {monomial()};
    if (level >= 1) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            monomsInTopRow.push_back(variables[i]);
        }
    }
    if (level >= 2) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=i; j<variables.size(); j++) {
                monomial currentMonomial = {variables[i][0], variables[i][1], variables[j][0], variables[j][1]};
                monomsInTopRow.push_back(currentMonomial);
            }
        }
    }
    if (level >= 3) {
        for (long unsigned int i=0; i<variables.size(); i++) {
            for (long unsigned int j=i; j<variables.size(); j++) {
                for (long unsigned int k=j; k<variables.size(); k++) {
                    monomial currentMonomial = {variables[i][0], variables[i][1], variables[j][0], variables[j][1], variables[k][0], variables[k][1]};
                    monomsInTopRow.push_back(currentMonomial);
                }
            }
        }
    }

    // Form the moment matrix
    polynomialMatrix matrixToReturn = std::vector<std::vector<polynomial>>(monomsInTopRow.size(), std::vector<polynomial>(monomsInTopRow.size()));
    for (long unsigned int i=0; i<monomsInTopRow.size(); i++) {
        for (long unsigned int j=i; j<monomsInTopRow.size(); j++) {

            // If it's diagonal it's just 1
            if (i == j) {
                matrixToReturn[i][j] = polynomial(1, std::make_pair(1, monomial()));
                continue;
            }

            // Form the new moment
            monomial newMoment;
            for (long unsigned int k=0; k<monomsInTopRow[i].size(); k++) {
                newMoment.push_back(monomsInTopRow[i][k]);
            }
            for (long unsigned int k=0; k<monomsInTopRow[j].size(); k++) {
                newMoment.push_back(monomsInTopRow[j][k]);
            }

            // Set the matrix element
            matrixToReturn[i][j] = polynomial(1, std::make_pair(1, newMoment));
            matrixToReturn[j][i] = polynomial(1, std::make_pair(1, newMoment));

        }
    }

    // Return the moment matrix
    return matrixToReturn;

}

// Add variables from a moment matrix
void addVariables(std::vector<monomial>& variables, polynomialMatrix toAdd) {

    // Iterative through the matrix
    for (long unsigned int i=0; i<toAdd.size(); i++) {
        for (long unsigned int j=0; j<toAdd[i].size(); j++) {

            // Iterate through the polynomial
            for (long unsigned int k=0; k<toAdd[i][j].size(); k++) {
                monomial currentMonomial = toAdd[i][j][k].second;

                // Check if this moment is already in the list
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

        // Check if this moment is already in the list
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
void solveMOSEK(polynomial obj, polynomialMatrix psd, std::vector<polynomial> constraintsZero, std::vector<polynomial> constraintsPositive) {

    // Get the list of variables
    int oneIndex = 0;
    std::vector<monomial> variables = {monomial()};
    addVariables(variables, obj);
    addVariables(variables, psd);

    // The c vector defining the objective
    std::vector<double> c(variables.size());
    for (int i=0; i<obj.size(); i++) {

        // Find the location of this variable
        monomial toFind = {obj[i].second};
        for (int j=0; j<variables.size(); j++) {
            if (variables[j] == toFind) {
                c[j] = obj[i].first;
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

    // The vectors defining the PSD constraints
    std::vector<int> indicesPSD;
    std::vector<double> coeffsPSD;
    for (int i=0; i<psd.size(); i++) {
        for (int j=i; j<psd[i].size(); j++) {

            // Find this in the variable list
            monomial toFind = {psd[i][j][0].second};
            for (int k=0; k<variables.size(); k++) {
                if (variables[k] == toFind) {
                    indicesPSD.push_back(k);
                    break;
                }
            }

            // The coeff for mosek's svec
            if (i != j) {
                coeffsPSD.push_back(std::sqrt(2.0));
            } else {
                coeffsPSD.push_back(1.0);
            }

        }
    }

    // Create a model TODO get CHSH working
    mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});
    M->setLogHandler([=](const std::string & msg) {std::cout << msg << std::flush;});

    // Create the main variable vector
    mosek::fusion::Variable::t xM = M->variable(variables.size(), mosek::fusion::Domain::inRange(-1, 1));

    // The objective function
    M->objective(mosek::fusion::ObjectiveSense::Maximize, mosek::fusion::Expr::dot(cM, xM));

    // The one variable should be fixed
    M->constraint(xM->index(oneIndex), mosek::fusion::Domain::equalsTo(1.0));

    // The matrix of this should be PSD
    M->constraint(mosek::fusion::Expr::mulElm(monty::new_array_ptr<double>(coeffsPSD), xM->pick(monty::new_array_ptr<int>(indicesPSD))), mosek::fusion::Domain::inSVecPSDCone());

    // Linear equality constraints
    M->constraint(mosek::fusion::Expr::mul(AM, xM), mosek::fusion::Domain::equalsTo(0.0));

    // Linear positivity constraints
    M->constraint(mosek::fusion::Expr::mul(BM, xM), mosek::fusion::Domain::greaterThan(0.0));

    // Solve the problem
    M->solve();

    // Output the solution
    double objPrimal = M->primalObjValue();
    std::cout << "Objective value: " << objPrimal << std::endl;

}

// Generic entry function
int main(int argc, char* argv[]) {

    // Define the Bell scenario
    int numOutcomesA = 2;
    int numOutcomesB = 2;
    int level = 1;

    // CHSH (c = 2, q = 2sqrt(2) w/ level 1)
    polynomial bellFunc = stringToPolynomial("-<A1>-<B1>+<A1B1>+<A1B2>+<A2B1>-<A2B2>");

    // I3322
    //polynomial bellFunc = stringToPolynomial("-<A2>-<B1>-2<B2>+<A1B1>+<A1B2>+<A2B1>+<A2B2>-<A1B3>+<A2B3>-<A3B1>+<A3B2>");

    // Define the moment matrix
    polynomialMatrix momentMatrix = generateMomentMatrix(bellFunc, level);

    // Define other constraints
    std::vector<polynomial> constraintsZero;
    std::vector<polynomial> constraintsPositive;

    // Output the problem
    if (bellFunc.size() > 0) {
        std::cout << "Bell functional: " << std::endl;
        std::cout << bellFunc << std::endl << std::endl;
    }
    if (momentMatrix.size() > 0) {
        std::cout << "Moment matrix: " << std::endl;
        std::cout << momentMatrix << std::endl << std::endl;
    }
    if (constraintsZero.size() > 0) {
        std::cout << "Zero constraints: " << std::endl;
        std::cout << constraintsZero << std::endl << std::endl;
    }
    if (constraintsPositive.size() > 0) {
        std::cout << "Positive constraints: " << std::endl;
        std::cout << constraintsPositive << std::endl << std::endl;
    }

    // Convert to MOSEK form and solve
    solveMOSEK(bellFunc, momentMatrix, constraintsZero, constraintsPositive);

    // Exit without errors
    return 0;

}
