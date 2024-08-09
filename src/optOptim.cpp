#include "optOptim.h"
#include "utils.h"
#include <chrono>

// Optim
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"

// Import Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>

// The structure used to hold the data for the optimisation
struct optimData {
    Poly objective;
    Eigen::SparseMatrix<double> A;
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>* linSolver;
    Eigen::VectorXd b;
    std::vector<std::vector<Poly>> momentMatrix;
    std::vector<Mon> varList;
    double tol;
    int maxIters;
    Eigen::VectorXd xOg;
    bool perIterOutput;
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
    double minEig = eigVals.minCoeff();

    // The cost function to minimize
    double cost = std::pow(errorLin,2) + std::pow(minEig,2);

    // Calculate the gradient
    if (gradOut) {

        // SDP projection
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

        // Linear projection
        Eigen::VectorXd projX = optDataRecast->linSolver->solveWithGuess(optDataRecast->b, xPos);

        // Reset the gradient
        *gradOut = Eigen::VectorXd::Zero(x.size());

        // Set the gradient
        for (int i=0; i<gradOut->size(); i++) {
            (*gradOut)(i) = x(i) - projX(i);
        }

        // Early convergence
        if (errorLin < optDataRecast->tol && minEig > -optDataRecast->tol) {
            *gradOut = Eigen::VectorXd::Zero(x.size());
        }

        // Output the objective of the x + grad
        double newCost = gradFunction(x + *gradOut, nullptr, optData);
        std::cout << "new cost: " << newCost << std::endl;

        // Per iteration output
        if (optDataRecast->perIterOutput) {
            //std::cout << "obj=" << obj << "  lin=" << errorLin << "  eig=" << minEig << "       \r" << std::flush;
            std::cout << "obj=" << obj << "  lin=" << errorLin << "  eig=" << minEig << "       \n" << std::flush;
        }

    }

    // Return the cost
    return cost;

}

// Attempt to solve using Optim
double solveOptim(Poly& objective, std::vector<Poly>& constraintsZero, std::vector<std::vector<std::vector<Poly>>>& momentMatrices, std::map<Mon, std::complex<double>>& startVals, int verbosity, int maxIters, int numExtra, double distance, double tolerance) {

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
    if (verbosity >= 1) {
        std::cout << "Starting distance: " << distance << std::endl;
    }
    std::map<Mon, std::complex<double>> varVals;
    if (startVals.size() == 0) {
        for (auto& mon : vars) {
            varVals[mon] = 0;
        }
    } else {
        for (auto& mon : vars) {
            varVals[mon] = startVals[mon];
        }
    }
    for (auto& mon : vars) {
        if (objective.contains(mon)) {
            if (objective[mon].real() > 0) {
                varVals[mon] += distance;
            } else {
                varVals[mon] -= distance;
            }
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

    // Set up the linear solver
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> linSolver;
    linSolver.setTolerance(tolerance/100.0);
    linSolver.compute(A);

    // Convert to a vector
    Eigen::VectorXd x = Eigen::VectorXd::Zero(varList.size());
    for (int i=0; i<varList.size(); i++) {
        x(i) = varVals[varList[i]].real();
    }
    Eigen::VectorXd xOg = x;

    // Test With Optim
    optimData optData;
    optData.objective = objective;
    optData.A = A;
    optData.linSolver = &linSolver;
    optData.b = b;
    optData.momentMatrix = momentMatrices[0];
    optData.varList = varList;
    optData.tol = tolerance;
    optData.xOg = xOg;
    optim::algo_settings_t settings;
    settings.print_level = 0;
    settings.iter_max = maxIters;
    settings.grad_err_tol = 1e-12;
    settings.rel_sol_change_tol = 1e-12;
    settings.rel_objfn_change_tol = 1e-12;
    settings.lbfgs_settings.par_M = 10;
    settings.lbfgs_settings.wolfe_cons_1 = 1e-3;
    settings.lbfgs_settings.wolfe_cons_2 = 0.9;
    settings.gd_settings.method = 0;

    // Solve
    Eigen::VectorXd projX = linSolver.solveWithGuess(b, x);
    optData.perIterOutput = verbosity >= 1;
    bool success = optim::lbfgs(projX, gradFunction, &optData, settings);
    //bool success2 = optim::gd(projX, gradFunction, &optData, settings);
    //bool success3 = optim::lbfgs(projX, gradFunction, &optData, settings);
    std::cout << std::endl;
    for (int i=0; i<varList.size(); i++) {
        varVals[varList[i]] = projX(i);
    }
    x = projX;

    // Travel a bit in the objective direction
    double prevObj = 10000000;
    distance = x.size() / 10.0;
    std::vector<Eigen::VectorXd> xList;
    optData.perIterOutput = 0;
    tolerance *= 1000;
    double reductionFactor = std::pow(tolerance / (distance+1), 1.0/numExtra);
    if (numExtra > 0 && verbosity >= 1) {
        std::cout << "Reduction factor: " << reductionFactor << std::endl;
    }
    for (int extraIter=0; extraIter<numExtra; extraIter++) {

        // Move a bit in the direction
        Eigen::VectorXd delta = xOg - x;
        delta.normalize();
        x += delta*distance;

        // Do something based on the previous x TODO

        // The last iteration should be longer
        if (extraIter == numExtra-1) {
            settings.iter_max = 1000000;
        } else {
            settings.iter_max = maxIters;
        }

        // Solve
        std::chrono::steady_clock::time_point timeStartLin = std::chrono::steady_clock::now();
        Eigen::VectorXd projX = linSolver.solveWithGuess(b, x);
        std::chrono::steady_clock::time_point timeFinishedLin = std::chrono::steady_clock::now();
        bool success = optim::lbfgs(projX, gradFunction, &optData, settings);
        std::chrono::steady_clock::time_point timeFinishedOpt = std::chrono::steady_clock::now();
        std::map<Mon, std::complex<double>> varValsNew;
        for (int i=0; i<varList.size(); i++) {
            varValsNew[varList[i]] = projX(i);
        }

        // Calculate the objective
        double newObj = -objective.eval(varValsNew).real();
        std::cout << "iter=" << extraIter << "  dis=" << distance << "  obj=" << newObj << "  timeLin=" << std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedLin - timeStartLin).count() << "ms  timeOpt=" << std::chrono::duration_cast<std::chrono::milliseconds>(timeFinishedOpt - timeFinishedLin).count() << "ms" << std::endl;

        // Then reduce the distance
        distance *= reductionFactor;

        // Adjust the distance
        double objDiff = std::abs(newObj - prevObj);
        x = projX;
        prevObj = newObj;

    }

    // Update the varVals
    for (int i=0; i<varList.size(); i++) {
        varVals[varList[i]] = x(i);
    }
    varVals[Mon()] = 1;

    // Calculate the final objective
    double finalObjective = objective.eval(varVals).real();

    // Verbose output of the final solution
    if (verbosity >= 1) {
        std::cout << std::endl;
        int numIters = settings.opt_iter;
        std::cout << "Iterations needed: " << numIters << std::endl;
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

    // Return the values
    if (startVals.size() == 0) {
        startVals = varVals;
    }

    // Return the objective
    return finalObjective;

}
