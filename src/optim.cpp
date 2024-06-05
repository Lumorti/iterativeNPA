#include "optim.h"
#include "utils.h"

// Optim
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"

// Import Eigen
#include <Eigen/Dense>

// The structure used to hold the data for the optimisation
struct optimData {
    double penalty;
    Poly objective;
    std::vector<Poly> equalityCons;
    std::vector<std::vector<Poly>> momentMatrix;
    std::vector<Mon> varList;
};

// Cost/gradient function for optim
static double gradFunction(const Eigen::VectorXd& x, Eigen::VectorXd* gradOut, void* optData) {
    
    // Recast this generic pointer into the correct format
    optimData* optDataRecast = reinterpret_cast<optimData*>(optData);

    std::map<Mon, std::complex<double>> xMap;
    for (int i=0; i<x.size(); i++) {
        xMap[optDataRecast->varList[i]] = x(i);
    }

    // The cost is a combination of the objective and the equality constraints
    double cost = 0.0;
    cost -= optDataRecast->objective.eval(xMap).real();
    for (size_t i=0; i<optDataRecast->equalityCons.size(); i++) {
        cost += optDataRecast->penalty * std::pow(optDataRecast->equalityCons[i].eval(xMap).real(), 2);
    }
    Eigen::MatrixXd momentMatrix = replaceVariables(optDataRecast->momentMatrix, xMap).real();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(momentMatrix);
    Eigen::VectorXd eigenvalues = es.eigenvalues();
    double minEigen = eigenvalues.minCoeff();
    double logDet = 0;
    if (minEigen < 0) {
        logDet = -1e10;
    } else {
        logDet = std::log(std::sqrt(minEigen));
    }
    cost -= optDataRecast->penalty * logDet;

    //std::cout << std::endl;
    //std::cout << "Min eigen = " << minEigen << std::endl;
    //std::cout << "Cost = " << cost << std::endl;
    //std::cout << "x vector = " << std::endl;
    //std::cout << x << std::endl;
    //std::cout << "Moment mat = " << std::endl;
    //std::cout << momentMatrix << std::endl;

    // Calculate the gradient
    if (gradOut) {

        // Reset the gradient
        *gradOut = Eigen::VectorXd::Zero(x.size());

        // Gradient of the matrix
        Eigen::MatrixXd gradMat = Eigen::MatrixXd::Zero(momentMatrix.rows(), momentMatrix.cols());
        for (int i=0; i<momentMatrix.rows(); i++) {
            gradMat(i, i) = -1e-3;
        }

        // Plus the inverse
        Eigen::MatrixXd invMat = momentMatrix.inverse();
        gradMat -= optDataRecast->penalty * invMat;

        // For each element of the grad mat, split it between the terms in the moment matrix
        double gradScaling = 1e-8;
        for (int i=0; i<optDataRecast->momentMatrix.size(); i++) {
            for (int j=0; j<optDataRecast->momentMatrix[i].size(); j++) {
                for (auto term : optDataRecast->momentMatrix[i][j]) {
                    if (!term.first.isConstant()) {
                        int varLoc = std::find(optDataRecast->varList.begin(), optDataRecast->varList.end(), term.first) - optDataRecast->varList.begin();
                        (*gradOut)(varLoc) += gradScaling * std::real(term.second) * gradMat(i, j);
                    }
                }
            }
        }

        //std::cout << "Grad mat = " << std::endl;
        //std::cout << gradMat << std::endl;
        //std::cout << "Var list = " << std::endl;
        //std::cout << optDataRecast->varList << std::endl;
        //std::cout << "Grad = " << std::endl;
        //std::cout << *gradOut << std::endl;

    }

    // Return just the cost
    return cost;

}

double solveOptim(Poly objective, std::vector<Poly> cons, std::vector<std::vector<Poly>> momentMatrix, int verbosity) {

    // Get the list of variables
    std::set<Mon> varList;
    addVariables(varList, momentMatrix);
    addVariables(varList, objective);
    for (size_t i=0; i<cons.size(); i++) {
        addVariables(varList, cons[i]);
    }
    std::vector<Mon> varListVec(varList.begin(), varList.end());
    for (size_t i=0; i<varListVec.size(); i++) {
        if (varListVec[i].isConstant()) {
            varListVec.erase(varListVec.begin() + i);
            i--;
        }
    }
    int numVars = varListVec.size();

    // Set up the optimisation
    optimData optData;
    optData.objective = objective;
    optData.equalityCons = cons;
    optData.momentMatrix = momentMatrix;
    optData.varList = varListVec;
    optData.penalty = 1;

    std::map<Mon, std::complex<double>> allZeros;
    for (int i=0; i<varListVec.size(); i++) {
        allZeros[varListVec[i]] = 0;
    }
    Eigen::MatrixXcd allZerosMat = replaceVariables(momentMatrix, allZeros);
    std::cout << "All zeros moment mat = " << std::endl;
    std::cout << allZerosMat << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(allZerosMat);
    Eigen::VectorXd eigenvalues = es.eigenvalues();
    double minEigen = eigenvalues.minCoeff();

    // Starting point
    Eigen::VectorXd x = Eigen::VectorXd::Zero(numVars);
    for (int i=0; i<momentMatrix.size(); i++) {
        if (!momentMatrix[i][i].isConstant()) {
            int varLoc = std::find(varListVec.begin(), varListVec.end(), momentMatrix[i][i].getKey()) - varListVec.begin();
            x(varLoc) = 1.1*minEigen;
        }
    }

    // Check the initial cost
    Eigen::VectorXd gradData = Eigen::VectorXd::Zero(numVars);
    double startCost = gradFunction(x, &gradData, &optData);
    std::cout << "Starting cost: " << startCost << std::endl;
    std::map<Mon, std::complex<double>> xMapStart;
    for (int i=0; i<x.size(); i++) {
        xMapStart[optData.varList[i]] = x(i);
    }
    std::cout << "Starting obj = " << objective.eval(xMapStart) << std::endl;
    Eigen::MatrixXcd MStart = replaceVariables(momentMatrix, xMapStart);
    std::cout << "Starting moment mat = " << std::endl;
    std::cout << MStart << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esStart(MStart);
    std::cout << "Starting eigenvalues = " << std::endl;
    std::cout << esStart.eigenvalues() << std::endl;

    // Settings for the optimiser
    optim::algo_settings_t settings;
    settings.print_level = verbosity - 1;
    settings.iter_max = 10000;
    settings.rel_objfn_change_tol = 1e-5;
    settings.gd_settings.method = 0;
    settings.lbfgs_settings.par_M = 10;
    settings.lbfgs_settings.wolfe_cons_1 = 1e-3;
    settings.lbfgs_settings.wolfe_cons_2 = 0.8;
    //settings.vals_bound = true;
    //settings.lower_bounds = -minEigen*Eigen::VectorXd::Ones(numVars);
    //settings.upper_bounds = minEigen*Eigen::VectorXd::Ones(numVars);

    //optData.penalty = 1e2;
    //bool success = optim::nm(x, gradFunction, &optData, settings);

    // Optimise
    optData.penalty = 1;
    for (int iter=0; iter<10000; iter++) {
    //for (int iter=0; iter<100000; iter++) {
        std::cout << "Optimising with " << numVars << " variables and penalty " << optData.penalty << std::endl;
        bool success = optim::gd(x, gradFunction, &optData, settings);
        //bool success = optim::de(x, gradFunction, &optData, settings);
        //bool success = optim::pso(x, gradFunction, &optData, settings);
        //bool success = optim::lbfgs(x, gradFunction, &optData, settings);
        //bool success = optim::nm(x, gradFunction, &optData, settings);
        std::cout << "Result = " << gradFunction(x, &gradData, &optData) << std::endl;
        optData.penalty /= 5;

        std::map<Mon, std::complex<double>> xMap;
        for (int i=0; i<x.size(); i++) {
            xMap[optData.varList[i]] = x(i);
        }
        std::cout << "Obj = " << objective.eval(xMap) << std::endl;
        Eigen::MatrixXcd M = replaceVariables(momentMatrix, xMap);
        //std::cout << "Moment mat = " << std::endl;
        //std::cout << M << std::endl;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(M);
        std::cout << "Min eigenvalue = " << es.eigenvalues().minCoeff() << std::endl;

        if (optData.penalty < 1e-5) {
            break;
        }

    }

    optData.penalty = 0;
    double result = gradFunction(x, &gradData, &optData);
    std::cout << "Optim result = " << result << std::endl;

    return result;

}
