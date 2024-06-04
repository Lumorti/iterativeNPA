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
    std::set<Mon> varList;
};

// Cost/gradient function for optim
static double gradFunction(const Eigen::VectorXd& x, Eigen::VectorXd* gradOut, void* optData) {
    
    // Recast this generic pointer into the correct format
    optimData* optDataRecast = reinterpret_cast<optimData*>(optData);

    std::map<Mon, std::complex<double>> xMap;
    int i = 0;
    for (auto mon : optDataRecast->varList) {
        xMap[mon] = std::real(x(i));
        i++;
    }

    //std::cout << "X = " << x << std::endl;

    // The cost is a combination of the objective and the equality constraints
    double cost = 0.0;
    cost -= optDataRecast->objective.eval(xMap).real();
    for (size_t i=0; i<optDataRecast->equalityCons.size(); i++) {
        cost += optDataRecast->penalty * std::pow(optDataRecast->equalityCons[i].eval(xMap).real(), 2);
    }
    Eigen::MatrixXcd momentMatrix = replaceVariables(optDataRecast->momentMatrix, xMap);
    cost -= optDataRecast->penalty * std::log(momentMatrix.determinant().real());

    // Calculate the gradient
    if (gradOut) {

        // Reset the gradient
        *gradOut = Eigen::VectorXd::Zero(x.size());

    }

    // Return just the cost
    return cost;

}

double solveOptim(Poly objective, std::vector<Poly> cons, std::vector<std::vector<Poly>> momentMatrix) {

    // Get the list of variables
    std::set<Mon> varList;
    addVariables(varList, momentMatrix);
    addVariables(varList, objective);
    for (size_t i=0; i<cons.size(); i++) {
        addVariables(varList, cons[i]);
    }
    int numVars = varList.size();

    // Set up the optimisation
    Eigen::VectorXd x = Eigen::VectorXd::Zero(numVars);
    optimData optData;
    optData.objective = objective;
    optData.equalityCons = cons;
    optData.momentMatrix = momentMatrix;
    optData.varList = varList;
    optData.penalty = 1e3;

    // Check the initial cost
    Eigen::VectorXd gradData = Eigen::VectorXd::Zero(numVars);
    double startCost = gradFunction(x, &gradData, &optData);
    std::cout << "Starting cost: " << startCost << std::endl;

    // Settings for the optimiser
    optim::algo_settings_t settings;
    settings.print_level = 1;
    settings.iter_max = 1000;
    settings.rel_objfn_change_tol = 1e-10;
    settings.gd_settings.method = 6;
    settings.lbfgs_settings.par_M = 10;
    settings.lbfgs_settings.wolfe_cons_1 = 1e-3;
    settings.lbfgs_settings.wolfe_cons_2 = 0.8;
    settings.vals_bound = true;
    settings.lower_bounds = -1*Eigen::VectorXd::Ones(numVars);
    settings.upper_bounds = 1*Eigen::VectorXd::Ones(numVars);

    //optData.penalty = 1e2;
    //bool success = optim::nm(x, gradFunction, &optData, settings);

    // Optimise
    optData.penalty = 1e1;
    for (int iter=0; iter<10; iter++) {
        std::cout << "Optimising with " << numVars << " variables and penalty " << optData.penalty << std::endl;
        //bool success = optim::lbfgs(x, gradFunction, &optData, settings);
        bool success = optim::nm(x, gradFunction, &optData, settings);
        std::cout << "Result = " << gradFunction(x, &gradData, &optData) << std::endl;
        optData.penalty *= 2;

        std::map<Mon, std::complex<double>> xMap;
        int i = 0;
        for (auto mon : varList) {
            xMap[mon] = std::real(x(i));
            i++;
        }
        std::cout << "Obj = " << objective.eval(xMap) << std::endl;
        Eigen::MatrixXcd M = replaceVariables(momentMatrix, xMap);
        std::cout << "Moment mat = " << std::endl;
        std::cout << M << std::endl;
        std::cout << "Det = " << M.determinant() << std::endl;

    }

    optData.penalty = 0;
    double result = gradFunction(x, &gradData, &optData);
    std::cout << "Result = " << result << std::endl;

    return 0.0;

}
